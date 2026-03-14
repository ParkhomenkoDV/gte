import os
from typing import Any, Dict, Tuple, Union

from numpy import isnan, nan
from scipy.optimize import root
from substance import Substance
from thermodynamics import adiabatic_index
from thermodynamics import parameters as tdp

try:
    from ...checks import check_mass_flow, check_temperature
    from ...config import EPSREL
    from ...config import parameters as gtep
    from ...errors import TYPE_ERROR
    from ...utils import call_with_kwargs, integral_average
    from ..node import GTENode
except ImportError:
    import os
    import sys

    sys.path.insert(0, os.getcwd())

    from gte.checks import check_mass_flow, check_temperature
    from gte.config import EPSREL
    from gte.config import parameters as gtep
    from gte.errors import TYPE_ERROR
    from gte.nodes.node import GTENode
    from gte.utils import call_with_kwargs, integral_average


class Nozzle(GTENode):
    """Выходное устройство"""

    variables: Tuple[str, str] = (gtep.eff_speed, gtep.pipi, gtep.force)
    n_vars: int = 2
    models: Dict[str, Any] = {}
    figure: Tuple[Tuple[float, ...], Tuple[float, ...]] = (
        (+0.4, -0.4, -0.4, +0.4),
        (+0.4, +0.4, -0.4, -0.4),
    )

    def __init__(self, parameters: Dict[str, float], name: str = "Nozzle"):
        """Инициализация объекта сопла"""
        GTENode.__init__(self, parameters, name)

    @classmethod
    def predict(cls, inlet: Substance, parameters: Dict[str, Union[float, int]]) -> Tuple[Dict[str, float], Substance]:
        """Начальные приближения"""
        GTENode.validate_substance(inlet)

        if not isinstance(parameters, dict):
            raise TypeError(TYPE_ERROR.format(f"{type(parameters)=}", dict))
        if len(parameters) != cls.n_vars:
            raise ValueError(f"{len(parameters)=} must be {cls.n_vars}")
        for parameter, value in parameters.items():
            if parameter not in cls.variables:
                raise ValueError(f"{parameter=} not in {cls.variables}")
            if not isinstance(value, (float, int)):
                raise TypeError(TYPE_ERROR.format(f"{type(value)=}", float))

        inlet_params = {tdp.t: inlet.parameters[gtep.TT], tdp.p: inlet.parameters[gtep.PP], tdp.eo: inlet.parameters.get(gtep.eo)}
        gc_i = call_with_kwargs(inlet.functions[gtep.gc], inlet_params)
        hcp_i = call_with_kwargs(inlet.functions[gtep.hcp], inlet_params)
        k_i = adiabatic_index(gc_i, hcp_i)

        outlet = Substance(
            inlet.name,
            inlet.composition,
            parameters={
                gtep.m: inlet.parameters[gtep.m],
                gtep.eo: inlet.parameters.get(gtep.eo),
                gtep.TT: inlet.parameters[gtep.TT],
            },
            functions=inlet.functions,
        )
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
                tdp.t: (inlet.parameters[gtep.TT], outlet.parameters[gtep.TT]),
                tdp.p: (inlet.parameters[gtep.PP], outlet_PP),
                tdp.eo: (inlet.parameters.get(gtep.eo), outlet.parameters.get(gtep.eo)),
            },
        )
        k = adiabatic_index(inlet.parameters[gtep.gc], hcp)

        c = eff_speed * (2 * hcp * outlet.parameters[gtep.TT] * (1 - pipi ** ((k - 1) / k))) ** 0.5

        return (
            pipi - outlet_PP / inlet.parameters[gtep.PP],
            force - outlet.parameters[gtep.m] * c,
        )

    @classmethod
    def calculate(cls, inlet: Substance, parameters: Dict[str, Union[float, int]]) -> Tuple[Dict[str, float], Substance]:
        prediction, outlet_ = cls.predict(inlet, parameters)

        outlet = Substance(
            inlet.name,
            inlet.composition,
            parameters={
                gtep.m: inlet.parameters[gtep.m],
                gtep.eo: inlet.parameters.get(gtep.eo),
                gtep.TT: inlet.parameters[gtep.TT],
            },
            functions=inlet.functions,
        )

        args: Dict[str, Any] = {"inlet": inlet, "outlet": outlet, **parameters}  # НУ
        x0 = [outlet_.parameters[gtep.PP]] + [prediction[v] for v in cls.variables if v not in parameters]

        result = root(cls._equations, x0, args, method="lm")
        outlet.parameters[gtep.PP] = float(result.x[0])
        for v in cls.variables:
            if v not in parameters:
                parameters[v] = float(result.x[1])

        outlet = GTENode.calculate_substance(outlet)

        ranges = {
            tdp.t: (inlet.parameters[gtep.TT], outlet.parameters[gtep.TT]),
            tdp.p: (inlet.parameters[gtep.PP], outlet.parameters[gtep.PP]),
            tdp.eo: (inlet.parameters.get(gtep.eo), outlet.parameters.get(gtep.eo)),
        }
        hcp, _ = integral_average(inlet.functions[gtep.hcp], **ranges)
        k = adiabatic_index(inlet.parameters[gtep.gc], hcp)

        outlet.parameters[gtep.c] = parameters[gtep.eff_speed] * (2 * hcp * outlet.parameters[gtep.TT] * (1 - parameters[gtep.pipi] ** ((k - 1) / k))) ** 0.5

        # outlet.parameters[gtep.T] = outlet.parameters[gtep.TT] - outlet.parameters[gtep.c] ** 2 / (2 * outlet.parameters[gtep.hcp])
        eff_speed = parameters.get(gtep.eff_speed, cls.efficiency_speed(inlet, outlet))
        pipi = parameters.get(gtep.pipi, cls.total_pressure_ratio(inlet, outlet))
        force = parameters.get(gtep.force, cls.force(inlet, outlet))

        return {gtep.pipi: pipi, gtep.eff_speed: eff_speed, gtep.force: force}, outlet

    @classmethod
    def validate(cls, inlet: Substance, outlet: Substance, epsrel: float = EPSREL) -> Dict[int, float]:
        pipi = cls.total_pressure_ratio(inlet, outlet)
        eff_speed = cls.speed_efficiency(inlet, outlet)

        x0 = (outlet.parameters[gtep.TT], outlet.parameters[gtep.PP])
        args = {"inlet": inlet, "outlet": outlet, gtep.p_eff: pipi, gtep.eff_speed: eff_speed}

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

        return ""

    @classmethod
    def efficiency_speed(cls, inlet: Substance, outlet: Substance) -> float:
        """КПД сохранения скорости"""
        pipi = cls.total_pressure_ratio(inlet, outlet)  # + проверка
        ranges = {
            tdp.t: (inlet.parameters[gtep.TT], outlet.parameters[gtep.TT]),
            tdp.p: (inlet.parameters[gtep.PP], outlet.parameters[gtep.PP]),
            tdp.eo: (inlet.parameters.get(gtep.eo), outlet.parameters.get(gtep.eo)),
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

        return inlet.parameters.get(gtep.PP, nan) / outlet.parameters.get(gtep.PP, nan)

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
        vars, outlet = Nozzle.calculate(inlet, parameters=test_case["parameters"])

        for k, v in outlet.parameters.items():
            print(f"{k:<40}: {v}")
        print(vars)

        # print(f"{Nozzle.validate(inlet, outlet) = }")
        # print(f"{Nozzle.check_real(inlet, outlet) = }")
        print()

    n = Nozzle({})
