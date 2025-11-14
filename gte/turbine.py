import os
import pickle
from copy import deepcopy
from typing import Any, Dict, Tuple

from numpy import isnan, nan
from scipy.optimize import root
from substance import Substance
from thermodynamics import adiabatic_index
from thermodynamics import parameters as tdp

try:
    from .checks import check_efficiency, check_pressure_ratio, check_temperature
    from .config import EPSREL
    from .config import parameters as gtep
    from .node import GTENode
    from .utils import call_with_kwargs, integral_average
except ImportError:
    from checks import check_efficiency, check_pressure_ratio, check_temperature
    from config import EPSREL
    from config import parameters as gtep
    from node import GTENode
    from utils import call_with_kwargs, integral_average


models = {}
for model in (gtep.TT, gtep.PP, gtep.pipi, gtep.effeff, gtep.power):
    path = f"gte/models/compressor_{model}.pkl"
    if os.path.exists(path):
        with open(path, "rb") as file:
            models[model] = pickle.load(file)
    else:
        print(f"'{path}' not found!")


class Turbine(GTENode):
    """Турбина"""

    models: Dict[str, Any] = models

    __slots__ = GTENode.__slots__ + [gtep.pipi, gtep.effeff, gtep.power]

    def __init__(self, name: str = "Turbine"):
        GTENode.__init__(self, name=name)

        setattr(self, gtep.pipi, nan)
        setattr(self, gtep.effeff, nan)
        setattr(self, gtep.power, nan)

    @property
    def variables(self) -> Dict[str, float]:
        return {
            gtep.pipi: getattr(self, gtep.pipi),
            gtep.effeff: getattr(self, gtep.effeff),
            gtep.power: getattr(self, gtep.power),
        }

    def predict(self) -> Dict[str, float]:
        """Начальные приближения"""
        prediction = {
            f"outlet_{gtep.TT}": self.inlet.parameters[gtep.TT] * 0.7,  # TODO model
            f"outlet_{gtep.PP}": self.inlet.parameters[gtep.PP] / 2,  # TODO model
        }
        for k, v in self.variables.items():
            if not isnan(v):
                continue
            if k == gtep.pipi:
                prediction[k] = 3  # TODO: model or formula
            elif k == gtep.effeff:
                prediction[k] = 1.0
            elif k == gtep.power:
                prediction[k] = 24 * 10**6  # TODO: model or formula
        return prediction

    def _equations(self, x: Tuple[float], args: Dict[str, Any]) -> Tuple:
        """Уравнения"""
        if not len(x):
            return (nan, nan, nan)

        self.outlet.parameters[gtep.TT] = x[0]
        self.outlet.parameters[gtep.PP] = x[1]
        idx = 2
        for k in self.variables:
            if k not in args:
                setattr(self, k, x[idx])
                idx += 1

        mf = (self.inlet.parameters[gtep.mf] + self.outlet.parameters[gtep.mf]) / 2
        ranges = {
            tdp.t: (self.inlet.parameters[gtep.TT], self.outlet.parameters[gtep.TT]),
            tdp.p: (self.inlet.parameters[gtep.PP], self.outlet.parameters[gtep.PP]),
            tdp.eo: (self.inlet.parameters.get(gtep.eo), self.outlet.parameters.get(gtep.eo)),
        }
        gc, _ = integral_average(self.inlet.functions[gtep.gc], **ranges)
        Cp, _ = integral_average(self.inlet.functions[gtep.Cp], **ranges)
        k = adiabatic_index(gc, Cp)

        return (
            getattr(self, gtep.power) - mf * Cp * (self.inlet.parameters[gtep.TT] - self.outlet.parameters[gtep.TT]),
            self.outlet.parameters[gtep.TT] - self.inlet.parameters[gtep.TT] * (1 - getattr(self, gtep.effeff) * (1 - getattr(self, gtep.pipi) ** ((1 - k) / k))),
            getattr(self, gtep.pipi) - self.inlet.parameters[gtep.PP] / self.outlet.parameters[gtep.PP],
        )

    def solve(self, substance_inlet: Substance, x0: Dict = None) -> Substance:  # TODO *inlet_substances
        GTENode.validate_substance(self, substance_inlet)
        self.inlet = deepcopy(substance_inlet)
        self.outlet = Substance(
            self.inlet.name,
            self.inlet.composition,
            parameters={gtep.mf: self.inlet.parameters[gtep.mf] - self.leak},
            functions=self.inlet.functions,
        )

        self.outlet.parameters[gtep.eo] = self.inlet.parameters[gtep.eo]  # TODO посчитать через массу!

        if x0 is None:
            x0 = tuple(self.predict().values())
        else:
            assert isinstance(x0, dict), TypeError(f"type x0 must be dict, but has {type(x0)}")
            for k, v in self.predict().items():
                if k not in x0:
                    x0[k] = v
        args = {k: v for k, v in self.variables.items() if not isnan(v)}
        count_variables = sum(1 if k not in args else 0 for k in self.variables)
        count_equations = len(self._equations(x0, args)) - 2  # outlet_TT, outlet_PP

        if count_variables < count_equations:
            raise ArithmeticError(f"{count_variables=} < {count_equations=}")
        elif count_variables > count_equations:
            raise ArithmeticError(f"{count_variables=} > {count_equations=}")

        root(self._equations, x0, args, method="lm")

        outlet_parameters = {tdp.t: self.outlet.parameters.get(gtep.TT), tdp.p: self.outlet.parameters.get(gtep.PP), tdp.eo: self.outlet.parameters.get(gtep.eo)}
        self.outlet.parameters[gtep.gc] = call_with_kwargs(self.outlet.functions[gtep.gc], outlet_parameters)
        self.outlet.parameters[gtep.Cp] = call_with_kwargs(self.outlet.functions[gtep.Cp], outlet_parameters)
        self.outlet.parameters[gtep.DD] = self.outlet.parameters[gtep.PP] / (self.outlet.parameters[gtep.gc] * self.outlet.parameters[gtep.TT])
        self.outlet.parameters[gtep.k] = adiabatic_index(self.outlet.parameters[gtep.gc], self.outlet.parameters[gtep.Cp])

        return self.outlet

    def validate(self, epsrel: float = EPSREL) -> bool:
        """Проверка найденного решения"""
        x0 = (self.outlet.parameters[gtep.TT], self.outlet.parameters[gtep.PP])
        args = {k: v for k, v in self.variables.items() if not isnan(v)}

        result = True
        for i, null in enumerate(self._equations(x0, args)):
            if isnan(null) or abs(null) > epsrel:
                result = False
                print(f"{i}: {null:.6f}")

        return result

    @property
    def is_real(self):
        checks = (
            check_efficiency(getattr(self, gtep.effeff)),
            check_temperature(self.outlet.parameters[gtep.TT]),
            check_pressure_ratio(getattr(self, gtep.pipi)),
        )
        return all(checks)


if __name__ == "__main__":
    from colorama import Fore
    from fixtures import exhaust

    test_cases = (
        {"name": "1", "turbine": {gtep.pipi: 3, gtep.effeff: 0.9, "leak": 0.03}},
        {"name": "2", "turbine": {gtep.pipi: 3, gtep.power: 12 * 10**6, "leak": 0.03}},
        {"name": "3", "turbine": {gtep.effeff: 0.9, gtep.power: 12 * 10**6, "leak": 0.03}},
    )

    for test_case in test_cases:
        t = Turbine(test_case["name"])

        for k, v in test_case["turbine"].items():
            setattr(t, k, v)

        t.solve(exhaust)

        for k, v in t.summary.items():
            print(f"{k:<40}: {v}")

        print(Fore.GREEN + f"{t.validate() = }" + Fore.RESET)
        print(Fore.GREEN + f"{t.is_real = }" + Fore.RESET)
        print()
