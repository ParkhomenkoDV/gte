import os
from typing import Any, Dict, Tuple, Union

from numpy import isnan, nan
from scipy.optimize import root
from substance import Substance
from thermodynamics import adiabatic_index
from thermodynamics import parameters as tdp

try:
    from ...checks import check_efficiency, check_mass_flow, check_pressure_ratio, check_temperature
    from ...config import EPSREL
    from ...config import parameters as gtep
    from ...errors import TYPE_ERROR
    from ...utils import call_with_kwargs, integral_average
    from ..node import GTENode
except ImportError:
    import os
    import sys

    sys.path.insert(0, os.getcwd())

    from gte.checks import check_efficiency, check_mass_flow, check_pressure_ratio, check_temperature
    from gte.config import EPSREL
    from gte.config import parameters as gtep
    from gte.errors import TYPE_ERROR
    from gte.nodes.node import GTENode
    from gte.utils import call_with_kwargs, integral_average


class Compressor(GTENode):
    """Компрессор"""

    variables: Tuple[str, str, str] = (gtep.effeff, gtep.pipi, gtep.power)
    n_vars: int = 2
    models: Dict[str, Any] = {}
    figure: Tuple[Tuple[float, ...], Tuple[float, ...]] = (
        (-0.4, +0.4, +0.4, -0.4, -0.4),
        (+0.4, +0.2, -0.2, -0.4, +0.4),
    )

    __slots__ = ()  # нет новых атрибутов

    def __init__(self, parameters: Dict[str, float], name: str = "Compressor"):
        """Инициализация объекта компрессора"""
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

        inlet_params: Dict[str, float] = {tdp.t: inlet.parameters[gtep.TT], tdp.p: inlet.parameters[gtep.PP], tdp.eo: inlet.parameters.get(gtep.eo)}
        gc_i: float = call_with_kwargs(inlet.functions[gtep.gc], inlet_params)
        hc_i: float = call_with_kwargs(inlet.functions[gtep.hcp], inlet_params)
        k_i: float = adiabatic_index(gc_i, hc_i)

        outlet = Substance(
            inlet.name,
            inlet.composition,
            parameters={
                gtep.m: inlet.parameters[gtep.m],
            },
            functions=inlet.functions,
        )
        vars: Dict[str, float] = {}

        if gtep.pipi not in parameters:
            outlet.parameters[gtep.TT] = inlet.parameters[gtep.TT] + parameters[gtep.power] / inlet.parameters[gtep.m] / hc_i
            outlet.parameters[gtep.PP] = inlet.parameters[gtep.PP] * ((outlet.parameters[gtep.TT] / inlet.parameters[gtep.TT] - 1) * parameters[gtep.effeff] + 1) ** (k_i / (k_i - 1))
            vars[gtep.pipi] = outlet.parameters[gtep.PP] / inlet.parameters[gtep.PP]
        elif gtep.effeff not in parameters:
            outlet.parameters[gtep.TT] = inlet.parameters[gtep.TT] + parameters[gtep.power] / inlet.parameters[gtep.m] / hc_i
            outlet.parameters[gtep.PP] = inlet.parameters[gtep.PP] * parameters[gtep.pipi]
            vars[gtep.effeff] = (parameters[gtep.pipi] ** ((k_i - 1) / k_i) - 1) / (outlet.parameters[gtep.TT] / inlet.parameters[gtep.TT] - 1)
        elif gtep.power not in parameters:
            outlet.parameters[gtep.TT] = inlet.parameters[gtep.TT] * (1 + parameters[gtep.pipi] ** ((k_i - 1) / k_i)) / parameters[gtep.effeff]
            outlet.parameters[gtep.PP] = inlet.parameters[gtep.PP] * parameters[gtep.pipi]
            vars[gtep.power] = inlet.parameters[gtep.m] * hc_i * (outlet.parameters[gtep.TT] - inlet.parameters[gtep.TT])
        else:
            raise ArithmeticError(f"{parameters=}")

        outlet = GTENode.calculate_substance(outlet)

        return vars, outlet

    @classmethod
    def _equations(cls, x: Tuple[float], args: Dict[str, Any]) -> Tuple[float, float, float]:
        """
        power = m * hcp * (T*_outlet - T*_inlet)
        T*_outlet = T*_inlet * (1 + (pi* ** ((k-1) / k) - 1) / eff*)
        pi* = P*_outlet / P*_inlet
        """
        effeff, pipi, power = args.get(gtep.effeff), args.get(gtep.pipi), args.get(gtep.power)
        outlet_TT, outlet_PP = x[0], x[1]
        if effeff is None:
            effeff = x[2]
        elif pipi is None:
            pipi = x[2]
        elif power is None:
            power = x[2]
        else:
            pass  # validate

        inlet, outlet = args["inlet"], args["outlet"]

        ranges = {
            tdp.t: (inlet.parameters[gtep.TT], outlet_TT),
            tdp.p: (inlet.parameters[gtep.PP], outlet_PP),
            tdp.eo: (inlet.parameters.get(gtep.eo), outlet.parameters.get(gtep.eo)),
        }
        gc, _ = integral_average(inlet.functions[gtep.gc], **ranges)
        hcp, _ = integral_average(inlet.functions[gtep.hcp], **ranges)
        k = adiabatic_index(gc, hcp)

        return (
            power - inlet.parameters[gtep.m] * hcp * (outlet_TT - inlet.parameters[gtep.TT]),
            outlet_TT - inlet.parameters[gtep.TT] * (1 + (pipi ** ((k - 1) / k) - 1) / effeff),
            pipi - outlet_PP / inlet.parameters[gtep.PP],
        )

    @classmethod
    def calculate(cls, inlet: Substance, parameters: Dict[str, Union[float, int]]) -> Tuple[Dict[str, float], Substance]:
        prediction, outlet_ = cls.predict(inlet, parameters)

        outlet = Substance(
            inlet.name,
            inlet.composition,
            parameters={
                gtep.m: inlet.parameters[gtep.m],
            },
            functions=inlet.functions,
        )

        args: Dict[str, Any] = {"inlet": inlet, "outlet": outlet, **parameters}  # НУ
        x0 = [outlet_.parameters[gtep.TT], outlet_.parameters[gtep.PP]] + [prediction[v] for v in cls.variables if v not in parameters]

        result = root(cls._equations, x0, args, method="lm")

        outlet.parameters[gtep.TT], outlet.parameters[gtep.PP] = float(result.x[0]), float(result.x[1])
        outlet = GTENode.calculate_substance(outlet)

        effeff = parameters.get(gtep.effeff, cls.total_efficiency(inlet, outlet))
        pipi = parameters.get(gtep.pipi, cls.total_pressure_ratio(inlet, outlet))
        power = parameters.get(gtep.power, cls.power(inlet, outlet))

        return {gtep.effeff: effeff, gtep.pipi: pipi, gtep.power: power}, outlet

    @classmethod
    def validate(cls, inlet: Substance, outlet: Substance, epsrel: float = EPSREL) -> Dict[int, float]:
        effeff = cls.total_efficiency(inlet, outlet)
        pipi = cls.total_pressure_ratio(inlet, outlet)
        power = cls.power(inlet, outlet)

        x0 = (outlet.parameters[gtep.TT], outlet.parameters[gtep.PP])
        args = {"inlet": inlet, "outlet": outlet, gtep.effeff: effeff, gtep.pipi: pipi, gtep.power: power}

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

        effeff = cls.total_efficiency(inlet, outlet)
        if not check_efficiency(effeff):
            return f"{effeff}"

        pipi = cls.total_pressure_ratio(inlet, outlet)
        if not check_pressure_ratio(pipi):
            return f"{pipi}"

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
            tdp.t: (inlet.parameters[gtep.TT], outlet.parameters[gtep.TT]),
            tdp.p: (inlet.parameters[gtep.PP], outlet.parameters[gtep.PP]),
            tdp.eo: (inlet.parameters.get(gtep.eo), outlet.parameters.get(gtep.eo)),
        }
        gc, _ = integral_average(inlet.functions[gtep.gc], **ranges)
        hcp, _ = integral_average(inlet.functions[gtep.hcp], **ranges)
        k = adiabatic_index(gc, hcp)

        return (pipi ** ((k - 1) / k) - 1) / (titi - 1)

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
            tdp.t: (inlet.parameters[gtep.TT], outlet.parameters[gtep.TT]),
            tdp.p: (inlet.parameters[gtep.PP], outlet.parameters[gtep.PP]),
            tdp.eo: (inlet.parameters.get(gtep.eo), outlet.parameters.get(gtep.eo)),
        }
        hcp, _ = integral_average(inlet.functions[gtep.hcp], **ranges)

        return inlet.parameters[gtep.m] * hcp * (outlet.parameters[gtep.TT] - inlet.parameters[gtep.TT])


if __name__ == "__main__":
    from gte.fixtures import air as inlet

    test_cases = (
        {"parameters": {gtep.pipi: 6, gtep.effeff: 0.85}},
        {"parameters": {gtep.pipi: 6, gtep.power: 12 * 10**6}},
        {"parameters": {gtep.effeff: 0.85, gtep.power: 12 * 10**6}},
    )
    for test_case in test_cases:
        c = Compressor(test_case["parameters"], "test")
        print(f"{c.is_solvable=}")

        vars, outlet = c.calculate(inlet, parameters=c.parameters)

        for k, v in outlet.parameters.items():
            print(f"{k:25}: {v:.4f}")
        print(vars)

        print(f"{Compressor.validate(inlet, outlet) = }")
        print(f"{Compressor.check_real(inlet, outlet) = }")
        print()
