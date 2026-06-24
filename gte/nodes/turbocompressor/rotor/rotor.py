from typing import Any, Dict, Tuple, Union

from numpy import isnan, nan
from scipy.optimize import root
from substance import Substance
from thermodynamics import adiabatic_index

try:
    from ....checks import check_mass_flow, check_temperature
    from ....config.config import EPSREL
    from ....config.config import parameters as gtep
    from ....errors import TYPE_ERROR
    from ....utils.utils import integral_average
    from ...node import Node
except ImportError:
    import os
    import sys

    sys.path.insert(0, os.getcwd())

    from gte.checks import check_mass_flow, check_temperature
    from gte.config.config import EPSREL
    from gte.config.config import parameters as gtep
    from gte.errors import TYPE_ERROR
    from gte.nodes.node import Node
    from gte.utils.utils import integral_average


class Rotor(Node):
    """Ротор"""

    variables: Tuple[str, str, str] = (gtep.effeff, gtep.titi, gtep.pipi)
    n_vars: int = 2

    __slots__ = ()  # нет новых атрибутов

    def __init__(self, parameters: Dict[str, float], name: str = "Rotor"):
        Node.__init__(self, parameters, name)

    @classmethod
    def _equations(cls, x: Tuple[float], args: Dict[str, Any]) -> Tuple[float, float, float]:
        """
        T*_outlet = T*_inlet * (1 + (pi* ** ((k-1) / k) - 1) / eff*)
        ti* = T*_outlet / T*_inlet
        pi* = P*_outlet / P*_inlet
        """
        outlet_TT, outlet_PP = x[0], x[1]
        effeff, titi, pipi = args.get(gtep.effeff, x[2]), args.get(gtep.titi, x[2]), args.get(gtep.pipi, x[2])

        inlet, outlet = args["inlet"], args["outlet"]

        ranges = {
            gtep.TT: (inlet.parameters[gtep.TT], outlet_TT),
            gtep.PP: (inlet.parameters[gtep.PP], outlet_PP),
            gtep.eo: (inlet.parameters.get(gtep.eo, 0), outlet.parameters.get(gtep.eo, 0)),
        }
        gc, _ = integral_average(inlet.functions[gtep.gc], **ranges)
        hcp, _ = integral_average(inlet.functions[gtep.hcp], **ranges)
        k = adiabatic_index(gc, hcp)

        return (
            outlet_TT - inlet.parameters[gtep.TT] * (1 + (pipi ** ((k - 1) / k) - 1) / effeff),
            titi - outlet_TT / inlet.parameters[gtep.TT],
            pipi - outlet_PP / inlet.parameters[gtep.PP],
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

        gc_i: float = inlet.functions[gtep.gc](inlet.parameters)
        hcp_i: float = inlet.functions[gtep.hcp](inlet.parameters)
        k_i: float = adiabatic_index(gc_i, hcp_i)

        outlet = Substance(
            inlet.name,
            inlet.composition,
            parameters={
                gtep.m: inlet.parameters[gtep.m],
            },
            functions=inlet.functions,
        )
        if gtep.eo in inlet.parameters:
            outlet.parameters["oxidizer"] = inlet.parameters["oxidizer"]
            outlet.parameters[gtep.eo] = inlet.parameters[gtep.eo]

        if gtep.titi in parameters and gtep.pipi in parameters:
            outlet.parameters[gtep.TT] = inlet.parameters[gtep.TT] * parameters[gtep.titi]
            outlet.parameters[gtep.PP] = inlet.parameters[gtep.PP] * parameters[gtep.pipi]
        elif gtep.effeff in parameters and gtep.titi in parameters:
            outlet.parameters[gtep.TT] = inlet.parameters[gtep.TT] * parameters[gtep.titi]
            outlet.parameters[gtep.PP] = inlet.parameters[gtep.PP] * ((outlet.parameters[gtep.TT] / inlet.parameters[gtep.TT] - 1) * parameters[gtep.effeff] + 1) ** (k_i / (k_i - 1))
        elif gtep.effeff in parameters and gtep.pipi in parameters:
            outlet.parameters[gtep.TT] = inlet.parameters[gtep.TT] * (1 + (parameters[gtep.pipi] ** ((k_i - 1) / k_i) - 1) / parameters[gtep.effeff])
            outlet.parameters[gtep.PP] = inlet.parameters[gtep.PP] * parameters[gtep.pipi]
        else:
            raise ArithmeticError(f"{parameters=}")

        outlet = Node.calculate_substance(outlet)

        return {gtep.effeff: cls.total_efficiency(inlet, outlet), gtep.titi: cls.total_temperature_ratio(inlet, outlet), gtep.pipi: cls.total_pressure_ratio(inlet, outlet)}, outlet

    @classmethod
    def calculate(cls, parameters: Dict[str, Union[float, int]], inlet: Substance) -> Tuple[Dict[str, float], Substance]:
        prediction, outlet_ = cls.predict(parameters, inlet)

        outlet = Substance(
            inlet.name,
            inlet.composition,
            parameters={
                gtep.m: inlet.parameters[gtep.m],
            },
            functions=inlet.functions,
        )

        if gtep.eo in inlet.parameters:
            outlet.parameters["oxidizer"] = inlet.parameters["oxidizer"]
            outlet.parameters[gtep.eo] = inlet.parameters[gtep.eo]

        args: Dict[str, Any] = {"inlet": inlet, "outlet": outlet, **parameters}  # НУ
        x0 = [outlet_.parameters[gtep.TT], outlet_.parameters[gtep.PP]] + [prediction[v] for v in cls.variables if v not in parameters]

        result = root(cls._equations, x0, args, method="hybr")

        outlet.parameters[gtep.TT], outlet.parameters[gtep.PP] = float(result.x[0]), float(result.x[1])
        outlet = Node.calculate_substance(outlet)

        effeff = parameters.get(gtep.effeff, cls.total_efficiency(inlet, outlet))
        titi = parameters.get(gtep.titi, cls.total_temperature_ratio(inlet, outlet))
        pipi = parameters.get(gtep.pipi, cls.total_pressure_ratio(inlet, outlet))

        return {gtep.effeff: effeff, gtep.titi: titi, gtep.pipi: pipi, gtep.power: cls.power(inlet, outlet)}, outlet

    @classmethod
    def validate(cls, inlet: Substance, outlet: Substance, epsrel: float = EPSREL) -> Dict[int, float]:
        effeff = cls.total_efficiency(inlet, outlet)
        titi = cls.total_temperature_ratio(inlet, outlet)
        pipi = cls.total_pressure_ratio(inlet, outlet)

        x0 = (outlet.parameters[gtep.TT], outlet.parameters[gtep.PP], effeff, titi, pipi)
        args = {"inlet": inlet, "outlet": outlet, gtep.effeff: effeff, gtep.titi: titi, gtep.pipi: pipi}

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
    def total_efficiency(cls, inlet: Substance, outlet: Substance) -> float:
        """Адиабатический КПД"""
        if not isinstance(inlet, Substance):
            raise TypeError(TYPE_ERROR.format(f"{type(inlet)=}", Substance))
        if not isinstance(outlet, Substance):
            raise TypeError(TYPE_ERROR.format(f"{type(outlet)=}", Substance))

        titi = outlet.parameters[gtep.TT] / inlet.parameters[gtep.TT]
        pipi = cls.total_pressure_ratio(inlet, outlet)

        ranges = {
            gtep.TT: (inlet.parameters[gtep.TT], outlet.parameters[gtep.TT]),
            gtep.PP: (inlet.parameters[gtep.PP], outlet.parameters[gtep.PP]),
            gtep.eo: (inlet.parameters.get(gtep.eo, 0), outlet.parameters.get(gtep.eo, 0)),
        }
        gc, _ = integral_average(inlet.functions[gtep.gc], **ranges)
        hcp, _ = integral_average(inlet.functions[gtep.hcp], **ranges)
        k = adiabatic_index(gc, hcp)

        return (pipi ** ((k - 1) / k) - 1) / (titi - 1)

    @classmethod
    def total_temperature_ratio(cls, inlet: Substance, outlet: Substance) -> float:
        """Степень повышения полной температуры"""
        if not isinstance(inlet, Substance):
            raise TypeError(TYPE_ERROR.format(f"{type(inlet)=}", Substance))
        if not isinstance(outlet, Substance):
            raise TypeError(TYPE_ERROR.format(f"{type(outlet)=}", Substance))

        return outlet.parameters.get(gtep.TT, nan) / inlet.parameters.get(gtep.TT, nan)

    @classmethod
    def total_pressure_ratio(cls, inlet: Substance, outlet: Substance) -> float:
        """Степень повышения полного давления"""
        if not isinstance(inlet, Substance):
            raise TypeError(TYPE_ERROR.format(f"{type(inlet)=}", Substance))
        if not isinstance(outlet, Substance):
            raise TypeError(TYPE_ERROR.format(f"{type(outlet)=}", Substance))

        return outlet.parameters.get(gtep.PP, nan) / inlet.parameters.get(gtep.PP, nan)

    @classmethod
    def power(cls, inlet: Substance, outlet: Substance) -> float:
        """Мощность"""
        if not isinstance(inlet, Substance):
            raise TypeError(TYPE_ERROR.format(f"{type(inlet)=}", Substance))
        if not isinstance(outlet, Substance):
            raise TypeError(TYPE_ERROR.format(f"{type(outlet)=}", Substance))

        ranges = {
            gtep.TT: (inlet.parameters[gtep.TT], outlet.parameters[gtep.TT]),
            gtep.PP: (inlet.parameters[gtep.PP], outlet.parameters[gtep.PP]),
            gtep.eo: (inlet.parameters.get(gtep.eo, 0), outlet.parameters.get(gtep.eo, 0)),
        }
        hcp, _ = integral_average(inlet.functions[gtep.hcp], **ranges)

        return inlet.parameters[gtep.m] * hcp * (outlet.parameters[gtep.TT] - inlet.parameters[gtep.TT])


if __name__ == "__main__":
    from gte.fixtures import air, exhaust

    test_cases = (
        # compressor
        {"parameters": {gtep.pipi: 6, gtep.effeff: 0.85}, "inlet": air},
        {"parameters": {gtep.pipi: 6, gtep.titi: 1.775}, "inlet": air},
        {"parameters": {gtep.effeff: 0.85, gtep.titi: 1.775}, "inlet": air},
        # turbine
        {"parameters": {gtep.effeff: 1 / 0.9, gtep.titi: 1 / 1.775}, "inlet": exhaust},
        {"parameters": {gtep.pipi: 1 / 4.0, gtep.titi: 1 / 1.775}, "inlet": exhaust},
        {"parameters": {gtep.effeff: 1 / 0.9, gtep.pipi: 1 / 4.0}, "inlet": exhaust},
    )
    for test_case in test_cases:
        r = Rotor(test_case["parameters"], "test")
        print(f"{r.is_solvable=}")

        vars, outlet = r.calculate(r.parameters, test_case["inlet"])

        for k, v in outlet.parameters.items():
            print(f"{k:25}: {v:.4f}")
        print(vars)

        print(f"{Rotor.validate(test_case["inlet"], outlet) = }")
        print(f"{Rotor.check_real(test_case["inlet"], outlet) = }")
        print()
