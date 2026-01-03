import inspect
import os
import pickle
from typing import Any, Callable, Dict, Tuple, Union

from numpy import isnan, nan
from scipy.optimize import root
from substance import Substance
from thermodynamics import adiabatic_index
from thermodynamics import parameters as tdp

try:
    from .checks import check_efficiency, check_mass_flow, check_pressure_ratio, check_temperature
    from .config import EPSREL
    from .config import parameters as gtep
    from .node import GTENode
    from .utils import call_with_kwargs, integral_average
except ImportError:
    from checks import check_efficiency, check_mass_flow, check_pressure_ratio, check_temperature
    from config import EPSREL
    from config import parameters as gtep
    from node import GTENode
    from utils import call_with_kwargs, integral_average


models = {}
for model in (gtep.pipi, gtep.effeff, gtep.power):
    path = f"gte/models/compressor_{model}.pkl.xz"
    if os.path.exists(path):
        with open(path, "rb") as file:
            models[model] = pickle.load(file)
    else:
        print(f"'{path}' not found!")


class Turbine(GTENode):
    """Турбина"""

    variables = (gtep.effeff, gtep.pipi, gtep.power)
    models: Dict[str, Any] = models

    def __init__(self, name: str = "Turbine", characteristic: Dict[str, Callable] = None):
        GTENode.__init__(self, name=name)

        assert isinstance(characteristic, dict), TypeError(f"{type(characteristic)=} must be dict")
        assert gtep.effeff in characteristic, KeyError(f"{gtep.effeff} not in {characteristic=}")
        assert gtep.pipi in characteristic, KeyError(f"{gtep.pipi} not in {characteristic=}")
        effeff, pipi = characteristic[gtep.effeff], characteristic[gtep.pipi]
        signature = inspect.signature(effeff)
        assert gtep.rf in signature.parameters, KeyError()
        assert gtep.mf in signature.parameters, KeyError()
        signature = inspect.signature(pipi)
        assert gtep.rf in signature.parameters, KeyError()
        assert gtep.mf in signature.parameters, KeyError()
        self.characteristic: Dict[str, Callable] = {gtep.effeff: effeff, gtep.pipi: pipi}

    def solve(self, inlet: Substance, rotation_frequency: float) -> Dict[str, Any]:
        GTENode.validate_substance(inlet)

        effeff = self.characteristic[gtep.effeff](rotation_frequency, inlet.parameters[gtep.mf])
        pipi = self.characteristic[gtep.pipi](rotation_frequency, inlet.parameters[gtep.mf])

        outlet = self.calculate(inlet, parameters={gtep.effeff: effeff, gtep.pipi: pipi})

        hc, _ = integral_average(
            inlet.functions[gtep.hcp],
            {
                tdp.t: (inlet.parameters[gtep.TT], outlet.parameters[gtep.TT]),
                tdp.p: (inlet.parameters[gtep.PP], outlet.parameters[gtep.PP]),
                tdp.eo: (inlet.parameters[gtep.eo], outlet.parameters[gtep.eo]),
            },
        )

        power = inlet.parameters[gtep.mf] * hc * (outlet.parameters[gtep.TT] - inlet.parameters[gtep.TT])  # TODO

        return {"outlet": outlet, gtep.power: power}

    @classmethod
    def predict(cls, inlet: Substance, parameters: Dict[str, float | int], use_ml: bool = True) -> Dict[str, float]:
        """Начальные приближения"""
        GTENode.validate_substance(inlet)

        inlet_params = {tdp.t: inlet.parameters.get(gtep.TT), tdp.p: inlet.parameters.get(gtep.PP), tdp.eo: inlet.parameters.get(gtep.eo)}
        gc_i = call_with_kwargs(inlet.functions[gtep.gc], inlet_params)
        hc_i = call_with_kwargs(inlet.functions[gtep.hcp], inlet_params)
        k_i = adiabatic_index(gc_i, hc_i)

        prediction: Dict[str, float] = {}
        if gtep.pipi not in parameters:
            prediction[f"outlet_{gtep.TT}"] = inlet.parameters[gtep.TT] - parameters[gtep.power] / inlet.parameters[gtep.mf] / hc_i
            prediction[f"outlet_{gtep.PP}"] = inlet.parameters[gtep.PP] * (1 - (1 - prediction[f"outlet_{gtep.TT}"] / inlet.parameters[gtep.TT]) / parameters[gtep.effeff]) ** (k_i / (k_i - 1))
            if use_ml:
                prediction[gtep.pipi] = (1 - (1 - prediction[f"outlet_{gtep.TT}"] / inlet.parameters[gtep.TT]) / parameters[gtep.effeff]) ** (k_i / (1 - k_i))
            else:
                prediction[gtep.pipi] = (1 - (1 - prediction[f"outlet_{gtep.TT}"] / inlet.parameters[gtep.TT]) / parameters[gtep.effeff]) ** (k_i / (1 - k_i))
        elif gtep.effeff not in parameters:
            prediction[f"outlet_{gtep.TT}"] = inlet.parameters[gtep.TT] - parameters[gtep.power] / inlet.parameters[gtep.mf] / hc_i
            prediction[f"outlet_{gtep.PP}"] = inlet.parameters[gtep.PP] / parameters[gtep.pipi]
            if use_ml:
                prediction[gtep.effeff] = (1 - prediction[f"outlet_{gtep.TT}"] / inlet.parameters[gtep.TT]) / (1 - parameters[gtep.pipi] ** ((1 - k_i) / k_i))
            else:
                prediction[gtep.effeff] = (1 - prediction[f"outlet_{gtep.TT}"] / inlet.parameters[gtep.TT]) / (1 - parameters[gtep.pipi] ** ((1 - k_i) / k_i))
        elif gtep.power not in parameters:
            prediction[f"outlet_{gtep.TT}"] = inlet.parameters[gtep.TT] * (1 - parameters[gtep.effeff] * (1 - parameters[gtep.pipi] ** ((1 - k_i) / k_i)))
            prediction[f"outlet_{gtep.PP}"] = inlet.parameters[gtep.PP] / parameters[gtep.pipi]
            if use_ml:
                prediction[gtep.power] = inlet.parameters[gtep.mf] * hc_i * (inlet.parameters[gtep.TT] - prediction[f"outlet_{gtep.TT}"])
            else:
                prediction[gtep.power] = inlet.parameters[gtep.mf] * hc_i * (inlet.parameters[gtep.TT] - prediction[f"outlet_{gtep.TT}"])
        else:
            raise ArithmeticError(f"{parameters=}")

        order = (f"outlet_{gtep.TT}", f"outlet_{gtep.PP}", gtep.effeff, gtep.pipi, gtep.power)  # порядок выдачи и получения признаков

        return {k: prediction[k] for k in order if k in prediction}

    @classmethod
    def _equations(cls, x: Tuple[float], args: Dict[str, Any]) -> Tuple[float, float, float]:
        """
        power = mf * hcp * (T*_inlet - T*_outlet)
        T*_outlet = T*_inlet * (1 - eff* * (1 - pi* ** ((1-k) / k)))
        pi* = P*_inlet / P*_outlet
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
        hc, _ = integral_average(inlet.functions[gtep.hcp], **ranges)
        k = adiabatic_index(gc, hc)

        return (
            power - inlet.parameters[gtep.mf] * hc * (inlet.parameters[gtep.TT] - outlet_TT),
            outlet_TT - inlet.parameters[gtep.TT] * (1 - effeff * (1 - pipi ** ((1 - k) / k))),
            pipi - inlet.parameters[gtep.PP] / outlet_PP,
        )

    @classmethod
    def calculate(cls, inlet: Substance, parameters: Dict[str, Union[float, int]]) -> Substance:  # TODO cooling
        prediction: Dict[str, float] = cls.predict(inlet, parameters, use_ml=False)

        outlet = Substance(
            inlet.name,
            inlet.composition,
            parameters={
                gtep.mf: inlet.parameters[gtep.mf],
                gtep.eo: inlet.parameters[gtep.eo],  # TODO посчитать через массу!
            },
            functions=inlet.functions,
        )

        args: Dict[str, Any] = {"inlet": inlet, "outlet": outlet, **parameters}  # НУ
        x0 = tuple(prediction.values())

        result = root(cls._equations, x0, args, method="lm")

        outlet.parameters[gtep.TT], outlet.parameters[gtep.PP] = float(result.x[0]), float(result.x[1])
        outlet = GTENode.calculate_substance(outlet)

        return outlet

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
        if not check_mass_flow(inlet.parameters[gtep.mf]):
            return f"inlet {gtep.mf} {inlet.parameters[gtep.mf]}"
        if not check_mass_flow(outlet.parameters[gtep.mf]):
            return f"outlet {gtep.mf} {outlet.parameters[gtep.mf]}"

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
        assert isinstance(inlet, Substance), TypeError(f"{type(inlet)=} must be {Substance}")
        assert isinstance(outlet, Substance), TypeError(f"{type(outlet)=} must be {Substance}")

        titi = outlet.parameters[gtep.TT] / inlet.parameters[gtep.TT]
        pipi = cls.pressure_ratio(inlet, outlet)

        ranges = {
            tdp.t: (inlet.parameters[gtep.TT], outlet.parameters[gtep.TT]),
            tdp.p: (inlet.parameters[gtep.PP], outlet.parameters[gtep.PP]),
            tdp.eo: (inlet.parameters.get(gtep.eo), outlet.parameters.get(gtep.eo)),
        }
        gc, _ = integral_average(inlet.functions[gtep.gc], **ranges)
        hcp, _ = integral_average(inlet.functions[gtep.hcp], **ranges)
        k = adiabatic_index(gc, hcp)

        return (1 - titi) / (1 - pipi ** ((1 - k) / k))

    @classmethod
    def pressure_ratio(cls, inlet: Substance, outlet: Substance) -> float:
        """Степень повышения полного давления"""
        assert isinstance(inlet, Substance), TypeError(f"{type(inlet)=} must be {Substance}")
        assert isinstance(outlet, Substance), TypeError(f"{type(outlet)=} must be {Substance}")

        return inlet.parameters.get(gtep.PP, nan) / outlet.parameters.get(gtep.PP, nan)

    @classmethod
    def power(cls, inlet: Substance, outlet: Substance) -> float:
        """Мощность"""
        assert isinstance(inlet, Substance), TypeError(f"{type(inlet)=} must be {Substance}")
        assert isinstance(outlet, Substance), TypeError(f"{type(outlet)=} must be {Substance}")

        ranges = {
            tdp.t: (inlet.parameters[gtep.TT], outlet.parameters[gtep.TT]),
            tdp.p: (inlet.parameters[gtep.PP], outlet.parameters[gtep.PP]),
            tdp.eo: (inlet.parameters.get(gtep.eo), outlet.parameters.get(gtep.eo)),
        }
        hc, _ = integral_average(inlet.functions[gtep.hcp], **ranges)

        return inlet.parameters[gtep.mf] * hc * (inlet.parameters[gtep.TT] - outlet.parameters[gtep.TT])


if __name__ == "__main__":
    from fixtures import exhaust

    test_cases = (
        {"parameters": {gtep.pipi: 3, gtep.effeff: 0.9}},
        {"parameters": {gtep.pipi: 3, gtep.power: 12 * 10**6}},
        {"parameters": {gtep.effeff: 0.9, gtep.power: 12 * 10**6}},
    )

    for test_case in test_cases:
        outlet = Turbine.calculate(exhaust, parameters=test_case["parameters"])

        for k, v in outlet.parameters.items():
            print(f"{k:<40}: {v}")

        print(f"{Turbine.validate(exhaust, outlet) = }")
        print(f"{Turbine.check_real(exhaust, outlet) = }")
        print()
