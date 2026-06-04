from typing import Any, Dict, Tuple, Union

from numpy import isnan, nan
from substance import Substance

try:
    from ...checks import check_mass_flow, check_pressure, check_temperature
    from ...config import EPSREL
    from ...config import parameters as gtep
    from ...errors import TYPE_ERROR
    from ..node import GTENode
except ImportError:
    import os
    import sys

    sys.path.insert(0, os.getcwd())

    from gte.checks import check_mass_flow, check_pressure, check_temperature
    from gte.config import EPSREL
    from gte.config import parameters as gtep
    from gte.errors import TYPE_ERROR
    from gte.nodes.node import GTENode


class Channel(GTENode):
    """Канал"""

    variables: Tuple[str, str] = (gtep.titi, gtep.pipi)
    n_vars: int = 2

    __slots__ = ()  # нет новых атрибутов

    def __init__(self, parameters: Dict[str, float], name: str = "Channel"):
        GTENode.__init__(self, parameters, name)

    @classmethod
    def predict(cls, parameters: Dict[str, Union[float, int]], inlet: Substance) -> Tuple[Dict[str, float], Substance]:
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

        outlet = Substance(
            inlet.name,
            inlet.composition,
            parameters={
                gtep.m: inlet.parameters[gtep.m],
                gtep.TT: inlet.parameters[gtep.TT] * parameters[gtep.titi],
                gtep.PP: inlet.parameters[gtep.PP] * parameters[gtep.pipi],
            },
            functions=inlet.functions,
        )
        if gtep.eo in inlet.parameters:
            outlet.parameters["oxidizer"] = inlet.parameters["oxidizer"]
            outlet.parameters[gtep.eo] = inlet.parameters[gtep.eo]

        vars = {}

        return vars, outlet

    @classmethod
    def _equations(cls, x: Tuple[float], args: Dict[str, Any]) -> Tuple[float, float]:
        """
        ti* = T*_outlet / T*_inlet
        pi* = P*_outlet / P*_inlet
        """
        titi, pipi = args.get(gtep.titi), args.get(gtep.pipi)
        outlet_TT, outlet_PP = x[0], x[1]

        inlet, _ = args["inlet"], args["outlet"]

        return (
            titi - outlet_TT / inlet.parameters[gtep.TT],
            pipi - outlet_PP / inlet.parameters[gtep.PP],
        )

    @classmethod
    def calculate(cls, parameters: Dict[str, Union[float, int]], inlet: Substance) -> Tuple[Dict[str, float], Substance]:
        _, outlet = cls.predict(parameters, inlet)

        outlet = GTENode.calculate_substance(outlet)

        titi = parameters.get(gtep.titi, cls.total_temperature_ratio(inlet, outlet))
        pipi = parameters.get(gtep.pipi, cls.total_pressure_ratio(inlet, outlet))

        return {gtep.titi: titi, gtep.pipi: pipi}, outlet

    @classmethod
    def validate(cls, inlet: Substance, outlet: Substance, epsrel: float = EPSREL) -> Dict[int, float]:
        titi = cls.total_temperature_ratio(inlet, outlet)
        pipi = cls.total_pressure_ratio(inlet, outlet)

        x0 = (titi, pipi)
        args = {"inlet": inlet, "outlet": outlet, gtep.titi: titi, gtep.pipi: pipi}

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


if __name__ == "__main__":
    from gte.fixtures import air as inlet

    test_cases = ({"parameters": {gtep.pipi: 0.95, gtep.titi: 1.05}},)

    for test_case in test_cases:
        ch = Channel(test_case["parameters"], name="test")
        print(f"{ch.is_solvable=}")

        vars, outlet = ch.calculate(ch.parameters, inlet)

        for k, v in outlet.parameters.items():
            print(f"{k:<25}: {v:.4f}")
        print(vars)

        print(f"{Channel.validate(inlet, outlet) = }")
        print(f"{Channel.check_real(inlet, outlet) = }")
        print()
