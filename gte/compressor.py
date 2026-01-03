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
    if os.path.isfile(path):
        with open(path, "rb") as file:
            models[model] = pickle.load(file)
    else:
        print(f"'{path}' not found!")


class Compressor(GTENode):
    """Компрессор"""

    variables = (gtep.effeff, gtep.pipi, gtep.power)
    models: Dict[str, Any] = models

    def __init__(self, name: str = "Compressor", characteristic: Dict[str, Callable] = None):
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

        hcp, _ = integral_average(
            inlet.functions[gtep.hcp],
            {
                tdp.t: (inlet.parameters[gtep.TT], outlet.parameters[gtep.TT]),
                tdp.p: (inlet.parameters[gtep.PP], outlet.parameters[gtep.PP]),
                tdp.eo: (inlet.parameters[gtep.eo], outlet.parameters[gtep.eo]),
            },
        )

        power = inlet.parameters[gtep.mf] * hcp * (outlet.parameters[gtep.TT] - inlet.parameters[gtep.TT])

        return {"outlet": outlet, gtep.power: power}

    @classmethod
    def predict(cls, inlet: Substance, parameters: Dict[str, float | int], use_ml: bool = True) -> Dict[str, float]:
        """Начальные приближения"""
        GTENode.validate_substance(inlet)

        assert isinstance(parameters, dict), TypeError(f"{type(parameters)=} must be dict")
        assert len(parameters) == 2, ArithmeticError(f"{len(parameters)=} must be 2")
        for key, value in parameters.items():
            assert key in (gtep.pipi, gtep.effeff, gtep.power)
            assert isinstance(value, (float, int)), TypeError(f"{type(value)=} must be numeric")

        assert isinstance(use_ml, bool), TypeError(f"{type(use_ml)=} must be bool")

        inlet_params = {tdp.t: inlet.parameters.get(gtep.TT), tdp.p: inlet.parameters.get(gtep.PP), tdp.eo: inlet.parameters.get(gtep.eo)}
        gc_i = call_with_kwargs(inlet.functions[gtep.gc], inlet_params)
        hc_i = call_with_kwargs(inlet.functions[gtep.hcp], inlet_params)
        k_i = adiabatic_index(gc_i, hc_i)

        prediction: Dict[str, float] = {}
        if gtep.pipi not in parameters:
            prediction[f"outlet_{gtep.TT}"] = inlet.parameters[gtep.TT] + parameters[gtep.power] / inlet.parameters[gtep.mf] / hc_i
            prediction[f"outlet_{gtep.PP}"] = inlet.parameters[gtep.PP] * ((prediction[f"outlet_{gtep.TT}"] / inlet.parameters[gtep.TT] - 1) * parameters[gtep.effeff] + 1) ** (k_i / (k_i - 1))
            if use_ml:
                prediction[gtep.pipi] = prediction[f"outlet_{gtep.PP}"] / inlet.parameters[gtep.PP]
            else:
                prediction[gtep.pipi] = prediction[f"outlet_{gtep.PP}"] / inlet.parameters[gtep.PP]
        elif gtep.effeff not in parameters:
            prediction[f"outlet_{gtep.TT}"] = inlet.parameters[gtep.TT] + parameters[gtep.power] / inlet.parameters[gtep.mf] / hc_i
            prediction[f"outlet_{gtep.PP}"] = inlet.parameters[gtep.PP] * parameters[gtep.pipi]
            if use_ml:
                prediction[gtep.effeff] = (parameters[gtep.pipi] ** ((k_i - 1) / k_i) - 1) / (prediction[f"outlet_{gtep.TT}"] / inlet.parameters[gtep.TT] - 1)
            else:
                prediction[gtep.effeff] = (parameters[gtep.pipi] ** ((k_i - 1) / k_i) - 1) / (prediction[f"outlet_{gtep.TT}"] / inlet.parameters[gtep.TT] - 1)
        elif gtep.power not in parameters:
            prediction[f"outlet_{gtep.TT}"] = inlet.parameters[gtep.TT] * (1 + parameters[gtep.pipi] ** ((k_i - 1) / k_i)) / parameters[gtep.effeff]
            prediction[f"outlet_{gtep.PP}"] = inlet.parameters[gtep.PP] * parameters[gtep.pipi]
            if use_ml:
                prediction[gtep.power] = inlet.parameters[gtep.mf] * hc_i * (prediction[f"outlet_{gtep.TT}"] - inlet.parameters[gtep.TT])
            else:
                prediction[gtep.power] = inlet.parameters[gtep.mf] * hc_i * (prediction[f"outlet_{gtep.TT}"] - inlet.parameters[gtep.TT])
        else:
            raise ArithmeticError(f"{parameters=}")

        order = (f"outlet_{gtep.TT}", f"outlet_{gtep.PP}", gtep.effeff, gtep.pipi, gtep.power)  # порядок выдачи и получения признаков

        return {k: prediction[k] for k in order if k in prediction}

    @classmethod
    def _equations(cls, x: Tuple[float], args: Dict[str, Any]) -> Tuple[float, float, float]:
        """
        power = mf * hcp * (T*_outlet - T*_inlet)
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

        inlet = args["inlet"]

        ranges = {
            tdp.t: (inlet.parameters[gtep.TT], outlet_TT),
            tdp.p: (inlet.parameters[gtep.PP], outlet_PP),
        }
        gc, _ = integral_average(inlet.functions[gtep.gc], **ranges)
        hc, _ = integral_average(inlet.functions[gtep.hcp], **ranges)
        k = adiabatic_index(gc, hc)

        return (
            power - inlet.parameters[gtep.mf] * hc * (outlet_TT - inlet.parameters[gtep.TT]),
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
                gtep.mf: inlet.parameters[gtep.mf],
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
        hc, _ = integral_average(inlet.functions[gtep.hcp], **ranges)
        k = adiabatic_index(gc, hc)

        return (pipi ** ((k - 1) / k) - 1) / (titi - 1)

    @classmethod
    def pressure_ratio(cls, inlet: Substance, outlet: Substance) -> float:
        """Степень повышения полного давления"""
        assert isinstance(inlet, Substance), TypeError(f"{type(inlet)=} must be {Substance}")
        assert isinstance(outlet, Substance), TypeError(f"{type(outlet)=} must be {Substance}")

        return outlet.parameters.get(gtep.PP, nan) / inlet.parameters.get(gtep.PP, nan)

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

        return inlet.parameters[gtep.mf] * hc * (outlet.parameters[gtep.TT] - inlet.parameters[gtep.TT])


if __name__ == "__main__":
    from fixtures import air

    test_cases = (
        {"parameters": {gtep.pipi: 6, gtep.effeff: 0.85}},
        {"parameters": {gtep.pipi: 6, gtep.power: 12 * 10**6}},
        {"parameters": {gtep.effeff: 0.85, gtep.power: 12 * 10**6}},
    )
    for test_case in test_cases:
        outlet = Compressor.calculate(inlet=air, parameters=test_case["parameters"])

        for k, v in outlet.parameters.items():
            print(f"{k:<40}: {v}")

        print(f"{Compressor.validate(air, outlet) = }")
        print(f"{Compressor.check_real(air, outlet) = }")
        print()
