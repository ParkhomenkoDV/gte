import os
from typing import Any, Callable, Dict, Tuple, Union

import matplotlib.pyplot as plt
from numpy import isnan, nan
from scipy.optimize import root
from substance import Substance
from thermodynamics import adiabatic_index
from thermodynamics import parameters as tdp

try:
    from ...checks import check_characteristic, check_efficiency, check_mass_flow, check_pressure_ratio, check_temperature
    from ...config import EPSREL
    from ...config import parameters as gtep
    from ...errors import TYPE_ERROR
    from ...utils import call_with_kwargs, integral_average
    from ..node import GTENode
except ImportError:
    import os
    import sys

    sys.path.insert(0, os.getcwd())

    from gte.checks import check_characteristic, check_efficiency, check_mass_flow, check_pressure_ratio, check_temperature
    from gte.config import EPSREL
    from gte.config import parameters as gtep
    from gte.errors import TYPE_ERROR
    from gte.nodes.node import GTENode
    from gte.utils import call_with_kwargs, integral_average


class Compressor(GTENode):
    """Компрессор"""

    variables: Tuple[str, str, str] = (gtep.effeff, gtep.pipi, gtep.power)
    models: Dict[str, Any] = {}
    figure: Tuple[Tuple[float, ...], Tuple[float, ...]] = (
        (-0.4, +0.4, +0.4, -0.4, -0.4),
        (+0.4, +0.2, -0.2, -0.4, +0.4),
    )

    def __init__(self, total_efficiency: Callable, total_pressure_ratio: Callable, name: str = "Compressor"):
        """Инициализация объекта компрессора"""
        for function in (total_efficiency, total_pressure_ratio):
            check_characteristic(function, {gtep.rf, gtep.m, gtep.TT, gtep.PP})

        GTENode.__init__(self, {gtep.effeff: total_efficiency, gtep.pipi: total_pressure_ratio}, name)

    def plot_characteristic(
        self,
        rotation_frequency=(0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10),
        mass_flow=(0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1),
    ) -> None:
        fg = plt.figure(figsize=(8, 10))
        gs = fg.add_gridspec(2, 1)  # строки, столбцы

        for i, func in enumerate(self.characteristic.values()):
            fg.add_subplot(gs[i, 0])
            plt.axis("equal")
            plt.xlabel(gtep.m, fontsize=12)
            plt.ylabel(func.target, fontsize=12)

            for rf in rotation_frequency:
                x, y = [], []
                for m in mass_flow:
                    x.append(m)
                    y.append(func(**{gtep.rf: rf, gtep.m: m}))
                plt.plot(x, y, label=f"{rf=:.4f}")

            plt.grid()
            plt.legend()

        plt.tight_layout()
        plt.show()

    @classmethod
    def predict(cls, inlet: Substance, parameters: Dict[str, Union[float, int]], use_ml: bool = True) -> Dict[str, float]:
        """Начальные приближения"""
        GTENode.validate_substance(inlet)

        if not isinstance(parameters, dict):
            raise TypeError(TYPE_ERROR.format(f"{type(parameters)=}", dict))
        if len(parameters) != 2:
            raise ValueError(f"{len(parameters)=} must be 2")
        for parameter, value in parameters.items():
            if parameter not in (gtep.pipi, gtep.effeff, gtep.power):
                raise ValueError(f"{parameter=} not in {(gtep.pipi, gtep.effeff, gtep.power)}")
            if not isinstance(value, (float, int)):
                raise TypeError(TYPE_ERROR.format(f"{type(value)=}", float))

        if not isinstance(use_ml, bool):
            raise TypeError(TYPE_ERROR.format(f"{type(use_ml)=}", bool))

        inlet_params = {tdp.t: inlet.parameters[gtep.TT], tdp.p: inlet.parameters[gtep.PP], tdp.eo: inlet.parameters.get(gtep.eo)}
        gc_i = call_with_kwargs(inlet.functions[gtep.gc], inlet_params)
        hc_i = call_with_kwargs(inlet.functions[gtep.hcp], inlet_params)
        k_i = adiabatic_index(gc_i, hc_i)

        prediction: Dict[str, float] = {}
        if gtep.pipi not in parameters:
            prediction[f"outlet_{gtep.TT}"] = inlet.parameters[gtep.TT] + parameters[gtep.power] / inlet.parameters[gtep.m] / hc_i
            prediction[f"outlet_{gtep.PP}"] = inlet.parameters[gtep.PP] * ((prediction[f"outlet_{gtep.TT}"] / inlet.parameters[gtep.TT] - 1) * parameters[gtep.effeff] + 1) ** (k_i / (k_i - 1))
            if use_ml:  # TODO
                prediction[gtep.pipi] = prediction[f"outlet_{gtep.PP}"] / inlet.parameters[gtep.PP]
            else:
                prediction[gtep.pipi] = prediction[f"outlet_{gtep.PP}"] / inlet.parameters[gtep.PP]
        elif gtep.effeff not in parameters:
            prediction[f"outlet_{gtep.TT}"] = inlet.parameters[gtep.TT] + parameters[gtep.power] / inlet.parameters[gtep.m] / hc_i
            prediction[f"outlet_{gtep.PP}"] = inlet.parameters[gtep.PP] * parameters[gtep.pipi]
            if use_ml:  # TODO
                prediction[gtep.effeff] = (parameters[gtep.pipi] ** ((k_i - 1) / k_i) - 1) / (prediction[f"outlet_{gtep.TT}"] / inlet.parameters[gtep.TT] - 1)
            else:
                prediction[gtep.effeff] = (parameters[gtep.pipi] ** ((k_i - 1) / k_i) - 1) / (prediction[f"outlet_{gtep.TT}"] / inlet.parameters[gtep.TT] - 1)
        elif gtep.power not in parameters:
            prediction[f"outlet_{gtep.TT}"] = inlet.parameters[gtep.TT] * (1 + parameters[gtep.pipi] ** ((k_i - 1) / k_i)) / parameters[gtep.effeff]
            prediction[f"outlet_{gtep.PP}"] = inlet.parameters[gtep.PP] * parameters[gtep.pipi]
            if use_ml:  # TODO
                prediction[gtep.power] = inlet.parameters[gtep.m] * hc_i * (prediction[f"outlet_{gtep.TT}"] - inlet.parameters[gtep.TT])
            else:
                prediction[gtep.power] = inlet.parameters[gtep.m] * hc_i * (prediction[f"outlet_{gtep.TT}"] - inlet.parameters[gtep.TT])
        else:
            raise ArithmeticError(f"{parameters=}")

        order = (f"outlet_{gtep.TT}", f"outlet_{gtep.PP}", gtep.effeff, gtep.pipi, gtep.power)  # порядок выдачи и получения признаков

        return {k: prediction[k] for k in order if k in prediction}

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
    def calculate(cls, inlet: Substance, parameters: Dict[str, Union[float, int]]) -> Substance:
        prediction: Dict[str, float] = cls.predict(inlet, parameters, use_ml=False)

        outlet = Substance(
            inlet.name,
            inlet.composition,
            parameters={
                gtep.m: inlet.parameters[gtep.m],
            },
            functions=inlet.functions,
        )

        args: Dict[str, Any] = {"inlet": inlet, "outlet": outlet, **parameters}  # НУ
        x0 = tuple(prediction.values())

        result = root(cls._equations, x0, args, method="lm")

        outlet.parameters[gtep.TT], outlet.parameters[gtep.PP] = float(result.x[0]), float(result.x[1])
        outlet = GTENode.calculate_substance(outlet)

        return outlet

    def solve(self, inlet: Substance, rotation_frequency: float) -> Dict[str, Any]:
        GTENode.validate_substance(inlet)

        effeff = self.characteristic[gtep.effeff](rotation_frequency, inlet.parameters[gtep.m])
        pipi = self.characteristic[gtep.pipi](rotation_frequency, inlet.parameters[gtep.m])

        outlet = self.calculate(inlet, parameters={gtep.effeff: effeff, gtep.pipi: pipi})

        hcp, _ = integral_average(
            inlet.functions[gtep.hcp],
            **{
                tdp.t: (inlet.parameters[gtep.TT], outlet.parameters[gtep.TT]),
                tdp.p: (inlet.parameters[gtep.PP], outlet.parameters[gtep.PP]),
                tdp.eo: (inlet.parameters.get(gtep.eo), outlet.parameters.get(gtep.eo)),
            },
        )

        power = inlet.parameters[gtep.m] * hcp * (outlet.parameters[gtep.TT] - inlet.parameters[gtep.TT])

        return {"outlet": outlet, gtep.power: power}

    @classmethod
    def validate(cls, inlet: Substance, outlet: Substance, epsrel: float = EPSREL) -> Dict[int, float]:
        effeff = cls.total_efficiency(inlet, outlet)
        pipi = cls.pressure_ratio(inlet, outlet)
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

        pipi = cls.pressure_ratio(inlet, outlet)
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
        pipi = cls.pressure_ratio(inlet, outlet)

        ranges = {
            tdp.t: (inlet.parameters[gtep.TT], outlet.parameters[gtep.TT]),
            tdp.p: (inlet.parameters[gtep.PP], outlet.parameters[gtep.PP]),
            tdp.eo: (inlet.parameters.get(gtep.eo), outlet.parameters.get(gtep.eo)),
        }
        gc, _ = integral_average(inlet.functions[gtep.gc], **ranges)
        hc, _ = integral_average(inlet.functions[gtep.hcp], **ranges)
        k = adiabatic_index(gc, hc)

        return (pipi ** ((k - 1) / k) - 1) / (titi - 1)

    @classmethod
    def pressure_ratio(cls, inlet: Substance, outlet: Substance) -> float:
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
        hc, _ = integral_average(inlet.functions[gtep.hcp], **ranges)

        return inlet.parameters[gtep.m] * hc * (outlet.parameters[gtep.TT] - inlet.parameters[gtep.TT])


if __name__ == "__main__":
    from gte.fixtures import air as inlet

    test_cases = (
        {"parameters": {gtep.pipi: 6, gtep.effeff: 0.85}},
        {"parameters": {gtep.pipi: 6, gtep.power: 12 * 10**6}},
        {"parameters": {gtep.effeff: 0.85, gtep.power: 12 * 10**6}},
    )
    for test_case in test_cases:
        outlet = Compressor.calculate(inlet, parameters=test_case["parameters"])

        for k, v in outlet.parameters.items():
            print(f"{k:<40}: {v}")

        print(f"{Compressor.validate(inlet, outlet) = }")
        print(f"{Compressor.check_real(inlet, outlet) = }")
        print()

    c = Compressor(
        total_efficiency=lambda rotation_frequency, mass: 0.85,
        total_pressure_ratio=lambda rotation_frequency, mass: 6,
    )
    result = c.solve(inlet, 0)
    print(result)
