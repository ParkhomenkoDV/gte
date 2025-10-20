from copy import deepcopy

from numpy import isnan, nan
from scipy.optimize import root
from substance import Substance
from thermodynamics import adiabatic_index, сritical_sonic_velocity

try:
    from .checks import check_efficiency, check_temperature
    from .config import EPSREL
    from .config import parameters as gtep
    from .node import GTENode
    from .utils import call_with_kwargs, integral_average
except ImportError:
    from checks import check_efficiency, check_temperature
    from config import EPSREL
    from config import parameters as gtep
    from node import GTENode
    from utils import call_with_kwargs, integral_average


class Compressor(GTENode):
    """Компрессор"""

    __slots__ = (gtep.pipi, gtep.effeff, gtep.power)

    def __init__(self, name="Compressor"):
        GTENode.__init__(self, name=name)

        setattr(self, gtep.pipi, nan)
        setattr(self, gtep.effeff, nan)
        setattr(self, gtep.power, nan)

    @property
    def variables(self) -> dict[str:float]:
        return {
            gtep.pipi: getattr(self, gtep.pipi),
            gtep.effeff: getattr(self, gtep.effeff),
            gtep.power: getattr(self, gtep.power),
        }

    @property
    def _x0(self) -> dict[str:float]:
        """Начальные приближения"""
        x0 = {
            f"outlet_{gtep.TT}": self.inlet.parameters[gtep.TT],
            f"outlet_{gtep.PP}": self.inlet.parameters[gtep.PP],
        }
        for k, v in self.variables.items():
            if not isnan(v):
                continue
            if k == gtep.pipi:
                x0[k] = 6  # TODO: model or formula
            elif k == gtep.effeff:
                x0[k] = 1.0
            elif k == gtep.power:
                x0[k] = 20 * 10**6  # TODO: model or formula
        return x0

    def equations(self, x: tuple, args: dict) -> tuple:
        """Уравнения"""
        self.outlet.parameters[gtep.TT] = x[0]
        self.outlet.parameters[gtep.PP] = x[1]
        idx = 2
        for k in self.variables:
            if k not in args:
                setattr(self, k, x[idx])
                idx += 1

        mf = (self.inlet.parameters[gtep.mf] + self.outlet.parameters[gtep.mf]) / 2
        gc, _ = integral_average(
            self.inlet.functions[gtep.gc],
            **{
                gtep.TT: (self.inlet.parameters[gtep.TT], self.outlet.parameters[gtep.TT]),
                gtep.PP: (self.inlet.parameters[gtep.PP], self.outlet.parameters[gtep.PP]),
            },
        )
        Cp, _ = integral_average(
            self.inlet.functions[gtep.Cp],
            **{
                gtep.TT: (self.inlet.parameters[gtep.TT], self.outlet.parameters[gtep.TT]),
                gtep.PP: (self.inlet.parameters[gtep.PP], self.outlet.parameters[gtep.PP]),
            },
        )
        k = adiabatic_index(gc, Cp)

        return (
            getattr(self, gtep.power) - mf * Cp * (self.outlet.parameters[gtep.TT] - self.inlet.parameters[gtep.TT]),
            self.outlet.parameters[gtep.TT] - self.inlet.parameters[gtep.TT] * (1 + (getattr(self, gtep.pipi) ** ((k - 1) / k) - 1) / getattr(self, gtep.effeff)),
            getattr(self, gtep.pipi) - self.outlet.parameters[gtep.PP] / self.inlet.parameters[gtep.PP],
        )

    def calculate(self, substance_inlet: Substance, x0: dict = None) -> Substance:
        GTENode.validate_substance(self, substance_inlet)
        self.inlet = deepcopy(substance_inlet)
        self.outlet.name = self.inlet.name
        self.outlet.functions = self.inlet.functions

        self.outlet.parameters[gtep.mf] = self.inlet.parameters[gtep.mf] - self.mass_flow_leak

        if x0 is None:
            x0 = tuple(self._x0.values())
        else:
            assert isinstance(x0, dict), TypeError(f"type x0 must be dict, but has {type(x0)}")
            for k, v in self._x0.items():
                if k not in x0:
                    x0[k] = v
        args = {k: v for k, v in self.variables.items() if not isnan(v)}
        count_variables = sum(1 if k not in args else 0 for k in self.variables)
        count_equations = len(self.equations(x0, args)) - 2  # outlet_TT, outlet_PP

        if count_variables < count_equations:
            raise ArithmeticError(f"{count_variables=} < {count_equations=}")
        elif count_variables > count_equations:
            raise ArithmeticError(f"{count_variables=} > {count_equations=}")

        root(self.equations, x0, args, method="lm")

        self.outlet.parameters[gtep.gc] = call_with_kwargs(self.outlet.functions[gtep.gc], self.outlet.parameters)
        self.outlet.parameters[gtep.Cp] = call_with_kwargs(self.outlet.functions[gtep.Cp], self.outlet.parameters)
        self.outlet.parameters[gtep.DD] = self.outlet.parameters[gtep.PP] / (self.outlet.parameters[gtep.gc] * self.outlet.parameters[gtep.TT])
        self.outlet.parameters[gtep.k] = adiabatic_index(self.outlet.parameters[gtep.gc], self.outlet.parameters[gtep.Cp])
        self.outlet.parameters[gtep.a_critical] = сritical_sonic_velocity(self.outlet.parameters[gtep.k], self.outlet.parameters[gtep.gc], self.outlet.parameters[gtep.TT])

        return self.outlet

    def validate(self, epsrel: float = EPSREL) -> bool:
        """Проверка найденного решения"""
        x0 = (self.outlet.parameters[gtep.TT], self.outlet.parameters[gtep.PP])
        args = {k: v for k, v in self.variables.items() if not isnan(v)}

        result = True
        for i, null in enumerate(self.equations(x0, args)):
            if abs(null) > epsrel:
                result = False
                print(f"{i}: {null:.6f}")

        return result

    @property
    def is_real(self):
        checks = (
            check_efficiency(getattr(self, gtep.effeff)),
            check_temperature(self.outlet.parameters[gtep.TT]),
        )
        return all(checks)


if __name__ == "__main__":
    from colorama import Fore
    from fixtures import air

    for k, v in gtep.items():
        print(f"{k:<10}: {v}")

    test_cases = (
        {"name": "1", "compressor": {gtep.pipi: 6, gtep.effeff: 0.85, "mass_flow_leak": 0.03}},
        {"name": "2", "compressor": {gtep.pipi: 6, gtep.power: 12 * 10**6, "mass_flow_leak": 0.03}},
        {"name": "3", "compressor": {gtep.effeff: 0.85, gtep.power: 12 * 10**6, "mass_flow_leak": 0.03}},
    )

    for test_case in test_cases:
        c = Compressor(test_case["name"])
        c.summary

        for k, v in test_case["compressor"].items():
            setattr(c, k, v)

        c.calculate(air)

        c.summary

        print(Fore.GREEN + f"{c.validate() = }" + Fore.RESET)
        print(Fore.GREEN + f"{c.is_real = }" + Fore.RESET)
        print()
