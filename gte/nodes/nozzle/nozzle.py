from typing import Any, Dict, Tuple, Union

from numpy import isnan, nan
from scipy.optimize import root
from substance import Substance
from thermodynamics import adiabatic_index

try:
    from ...checks import check_mass_flow, check_pressure, check_temperature
    from ...config import EPSREL
    from ...config import parameters as gtep
    from ...errors import TYPE_ERROR
    from ...utils import integral_average
    from ..node import Node
except ImportError:
    import os
    import sys

    sys.path.insert(0, os.getcwd())

    from gte.checks import check_mass_flow, check_pressure, check_temperature
    from gte.config import EPSREL
    from gte.config import parameters as gtep
    from gte.errors import TYPE_ERROR
    from gte.nodes.node import Node
    from gte.utils import integral_average


class Nozzle(Node):
    """Выходное устройство"""

    variables: Tuple[str, str] = (gtep.eff_speed, gtep.pipi, gtep.force)
    n_vars: int = 2

    __slots__ = ()  # нет новых атрибутов

    def __init__(self, parameters: Dict[str, float], name: str = "Nozzle"):
        Node.__init__(self, parameters, name)

    @classmethod
    def _equations(cls, x: Tuple[float], args: Dict[str, Any]) -> Tuple[float, float]:
        """
        pi* = P*_outlet / P*_inlet
        c_outlet = eff_speed * (2 * hcp * T*_outlet * (1 - pi* ** ((k - 1) / k))) ** 0.5
        force = m * c_outlet
        """
        eff_speed, pipi, force = args.get(gtep.eff_speed), args.get(gtep.pipi), args.get(gtep.force)
        outlet_PP = x[0]
        if eff_speed is None:
            eff_speed = x[1]
        elif pipi is None:
            pipi = x[1]
        elif force is None:
            force = x[1]
        else:
            pass  # validate

        inlet, outlet = args["inlet"], args["outlet"]

        hcp, _ = integral_average(
            inlet.functions[gtep.hcp],
            **{
                gtep.TT: (inlet.parameters[gtep.TT], outlet.parameters[gtep.TT]),
                gtep.PP: (inlet.parameters[gtep.PP], outlet_PP),
                gtep.eo: (inlet.parameters.get(gtep.eo, 0), outlet.parameters.get(gtep.eo, 0)),
            },
        )
        k = adiabatic_index(inlet.parameters[gtep.gc], hcp)

        c = eff_speed * (2 * hcp * outlet.parameters[gtep.TT] * (1 - pipi ** ((k - 1) / k))) ** 0.5

        return (
            pipi - outlet_PP / inlet.parameters[gtep.PP],
            force - outlet.parameters[gtep.m] * c,
        )

    @classmethod
    def predict(cls, parameters: Dict[str, Union[float, int]], inlet: Substance) -> Tuple[Dict[str, float], Substance]:
        """Начальные приближения"""
        Node.validate_substance(inlet)

        if not isinstance(parameters, dict):
            raise TypeError(TYPE_ERROR.format(f"{type(parameters)=}", dict))
        if len(parameters) != cls.n_vars:
            raise ValueError(f"{len(parameters)=} must be {cls.n_vars}")
        for parameter, value in parameters.items():
            if parameter not in cls.variables:
                raise ValueError(f"{parameter=} not in {cls.variables}")
            if not isinstance(value, (float, int)):
                raise TypeError(TYPE_ERROR.format(f"{type(value)=}", float))

        gc_i = inlet.functions[gtep.gc](inlet.parameters)
        hcp_i = inlet.functions[gtep.hcp](inlet.parameters)
        k_i = adiabatic_index(gc_i, hcp_i)

        outlet = Substance(
            inlet.name,
            inlet.composition,
            parameters={
                gtep.m: inlet.parameters[gtep.m],
                gtep.TT: inlet.parameters[gtep.TT],
            },
            functions=inlet.functions,
        )
        if gtep.eo in inlet.parameters:
            outlet.parameters["oxidizer"] = inlet.parameters["oxidizer"]
            outlet.parameters[gtep.eo] = inlet.parameters[gtep.eo]

        vars = {}

        if gtep.eff_speed not in parameters:
            outlet.parameters[gtep.PP] = inlet.parameters[gtep.PP] * parameters[gtep.pipi]
            outlet.parameters[gtep.c] = parameters[gtep.force] / inlet.parameters[gtep.m]
            vars[gtep.eff_speed] = outlet.parameters[gtep.c] / (2 * hcp_i * outlet.parameters[gtep.TT] * (1 - parameters[gtep.pipi] ** ((k_i - 1) / k_i))) ** 0.5
        elif gtep.pipi not in parameters:
            outlet.parameters[gtep.c] = parameters[gtep.force] / inlet.parameters[gtep.m]
            vars[f"{gtep.pipi}"] = (1 - ((outlet.parameters[gtep.c] / parameters[gtep.eff_speed]) ** 2) / (2 * hcp_i * outlet.parameters[gtep.TT])) ** (k_i / (k_i - 1))
            outlet.parameters[gtep.PP] = inlet.parameters[gtep.PP] * vars[f"{gtep.pipi}"]
        elif gtep.force not in parameters:
            outlet.parameters[gtep.PP] = inlet.parameters[gtep.PP] * parameters[gtep.pipi]
            outlet.parameters[gtep.c] = parameters[gtep.eff_speed] * (2 * hcp_i * outlet.parameters[gtep.TT] * (1 - parameters[gtep.pipi] ** ((k_i - 1) / k_i))) ** 0.5
            vars[gtep.force] = inlet.parameters[gtep.m] * outlet.parameters[gtep.c]
        else:
            raise ArithmeticError(f"{parameters=}")

        return vars, outlet

    @classmethod
    def calculate(cls, parameters: Dict[str, Union[float, int]], inlet: Substance) -> Tuple[Dict[str, float], Substance]:
        parameters_ = parameters.copy()  # необходима копия т.к. словарь передается по ссылке и при добавлении искоромй vars становится нерасчетным
        prediction, outlet_ = cls.predict(parameters_, inlet)

        outlet = Substance(
            inlet.name,
            inlet.composition,
            parameters={
                gtep.m: inlet.parameters[gtep.m],
                gtep.TT: inlet.parameters[gtep.TT],
            },
            functions=inlet.functions,
        )
        if gtep.eo in inlet.parameters:
            outlet.parameters["oxidizer"] = inlet.parameters["oxidizer"]
            outlet.parameters[gtep.eo] = inlet.parameters[gtep.eo]

        args: Dict[str, Any] = {"inlet": inlet, "outlet": outlet, **parameters_}  # НУ
        x0 = [outlet_.parameters[gtep.PP]] + [prediction[v] for v in cls.variables if v not in parameters_]

        result = root(cls._equations, x0, args, method="lm")
        outlet.parameters[gtep.PP] = float(result.x[0])
        for v in cls.variables:
            if v not in parameters_:
                parameters_[v] = float(result.x[1])

        outlet = Node.calculate_substance(outlet)

        ranges = {
            gtep.TT: (inlet.parameters[gtep.TT], outlet.parameters[gtep.TT]),
            gtep.PP: (inlet.parameters[gtep.PP], outlet.parameters[gtep.PP]),
            gtep.eo: (inlet.parameters.get(gtep.eo, 0), outlet.parameters.get(gtep.eo, 0)),
        }
        hcp, _ = integral_average(inlet.functions[gtep.hcp], **ranges)
        k = adiabatic_index(inlet.parameters[gtep.gc], hcp)

        outlet.parameters[gtep.c] = parameters_[gtep.eff_speed] * (2 * hcp * outlet.parameters[gtep.TT] * (1 - parameters_[gtep.pipi] ** ((k - 1) / k))) ** 0.5
        outlet.parameters[gtep.T] = outlet.parameters[gtep.TT] - outlet.parameters[gtep.c] ** 2 / (2 * outlet.parameters[gtep.hcp])

        eff_speed = parameters_.get(gtep.eff_speed, cls.efficiency_speed(inlet, outlet))
        pipi = parameters_.get(gtep.pipi, cls.total_pressure_ratio(inlet, outlet))
        force = parameters_.get(gtep.force, cls.force(inlet, outlet))

        return {gtep.pipi: pipi, gtep.eff_speed: eff_speed, gtep.force: force}, outlet

    @classmethod
    def validate(cls, inlet: Substance, outlet: Substance, epsrel: float = EPSREL) -> Dict[int, float]:
        eff_speed = cls.efficiency_speed(inlet, outlet)
        pipi = cls.total_pressure_ratio(inlet, outlet)
        force = cls.force(inlet, outlet)

        x0 = (outlet.parameters[gtep.PP], eff_speed, pipi, force)
        args = {"inlet": inlet, "outlet": outlet, gtep.eff_speed: eff_speed, gtep.pipi: pipi, gtep.force: force}

        result: Dict[int, float] = {}
        for i, null in enumerate(cls._equations(x0, args)):
            if isnan(null) or abs(null) > epsrel:
                result[i] = null

        return result

    @classmethod
    def check_real(cls, inlet: Substance, outlet: Substance) -> str:
        if not check_mass_flow(inlet.parameters[gtep.m]):
            return f"inlet {gtep.m} {inlet.parameters[gtep.m]}"
        if not check_mass_flow(outlet.parameters[gtep.m]):
            return f"outlet {gtep.m} {outlet.parameters[gtep.m]}"

        if not check_temperature(inlet.parameters[gtep.TT]):
            return f"inlet {gtep.TT} {inlet.parameters[gtep.TT]}"
        if not check_temperature(outlet.parameters[gtep.TT]):
            return f"outlet {gtep.TT} {outlet.parameters[gtep.TT]}"

        if not check_pressure(inlet.parameters[gtep.PP]):
            return f"inlet {gtep.PP} {inlet.parameters[gtep.PP]}"
        if not check_pressure(outlet.parameters[gtep.PP]):
            return f"outlet {gtep.PP} {outlet.parameters[gtep.PP]}"

        return ""

    @classmethod
    def efficiency_speed(cls, inlet: Substance, outlet: Substance) -> float:
        """КПД сохранения скорости"""
        pipi = cls.total_pressure_ratio(inlet, outlet)  # + проверка
        ranges = {
            gtep.TT: (inlet.parameters[gtep.TT], outlet.parameters[gtep.TT]),
            gtep.PP: (inlet.parameters[gtep.PP], outlet.parameters[gtep.PP]),
            gtep.eo: (inlet.parameters.get(gtep.eo, 0), outlet.parameters.get(gtep.eo, 0)),
        }
        hcp, _ = integral_average(inlet.functions[gtep.hcp], **ranges)
        k = adiabatic_index(inlet.parameters[gtep.gc], hcp)
        return outlet.parameters.get(gtep.c, nan) / (2 * hcp * outlet.parameters[gtep.TT] * (1 - pipi ** ((k - 1) / k))) ** 0.5

    @classmethod
    def total_pressure_ratio(cls, inlet: Substance, outlet: Substance) -> float:
        """Степень повышения полного давления"""
        if not isinstance(inlet, Substance):
            raise TypeError(TYPE_ERROR.format(f"{type(inlet)=}", Substance))
        if not isinstance(outlet, Substance):
            raise TypeError(TYPE_ERROR.format(f"{type(outlet)=}", Substance))

        return outlet.parameters.get(gtep.PP, nan) / inlet.parameters.get(gtep.PP, nan)

    @classmethod
    def force(cls, inlet: Substance, outlet: Substance) -> float:
        """Реактивная сила"""
        if not isinstance(inlet, Substance):
            raise TypeError(TYPE_ERROR.format(f"{type(inlet)=}", Substance))
        if not isinstance(outlet, Substance):
            raise TypeError(TYPE_ERROR.format(f"{type(outlet)=}", Substance))

        return outlet.parameters.get(gtep.m, nan) * outlet.parameters.get(gtep.c, nan)


if __name__ == "__main__":
    from gte.fixtures import exhaust as inlet

    inlet.parameters[gtep.TT] = 1200
    inlet.parameters[gtep.PP] = 101325 * 2

    test_cases = (
        {"parameters": {gtep.pipi: 1 / 1.8, gtep.eff_speed: 0.99}},
        {"parameters": {gtep.pipi: 1 / 1.8, gtep.force: 31_000}},
        {"parameters": {gtep.eff_speed: 0.99, gtep.force: 31_000}},
    )
    for test_case in test_cases:
        n = Nozzle(test_case["parameters"], name="test")
        print(f"{n.is_solvable=}")

        vars, outlet = n.calculate(n.parameters, inlet)

        for k, v in outlet.parameters.items():
            print(f"{k:<25}: {v:.4f}")
        print(vars)

        print(f"{Nozzle.validate(inlet, outlet) = }")
        print(f"{Nozzle.check_real(inlet, outlet) = }")
        print()
