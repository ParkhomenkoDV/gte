from copy import deepcopy

from mathematics import eps
from numpy import isinf, isnan, nan
from scipy.optimize import root
from substance import Substance
from thermodynamics import T0, adiabatic_index, сritical_sonic_velocity

try:
    from .checks import check_efficiency, check_temperature
    from .config import EPSREL
    from .config import parameters as gtep
    from .node import GTENode
    from .utils import call_with_kwargs, enthalpy
except ImportError:
    from checks import check_efficiency, check_temperature
    from config import EPSREL
    from config import parameters as gtep
    from node import GTENode
    from utils import call_with_kwargs, enthalpy


class CombustionChamber(GTENode):
    """Камера сгорания"""

    __slots__ = (gtep.effburn, gtep.peff, "fuel")

    def __init__(self, name="CombustionChamber"):
        GTENode.__init__(self, name=name)

        setattr(self, gtep.effburn, nan)
        setattr(self, gtep.peff, nan)

        self.fuel = Substance("fuel")

    @property
    def variables(self) -> dict:
        return {
            gtep.effburn: getattr(self, gtep.effburn),
            gtep.peff: getattr(self, gtep.peff),
        }

    @property
    def _x0(self) -> dict[str:float]:
        """Начальные приближения"""
        x0 = {
            f"outlet_{gtep.TT}": self.inlet.parameters[gtep.TT],
            f"outlet_{gtep.PP}": self.inlet.parameters[gtep.PP],
            f"outlet_{gtep.eo}": 3,  # 1 == clean_exhaust => prohibited
        }
        for k, v in self.variables.items():
            if not isnan(v):
                continue
            if k == "efficiency_burn":
                x0[k] = 1.0
            elif k == gtep.peff:
                x0[k] = 0.95  # TODO: model or formula
        return x0

    def equations(self, x: tuple, args: dict) -> tuple:
        """Уравнения"""
        self.outlet.parameters[gtep.TT] = x[0]
        self.outlet.parameters[gtep.PP] = x[1]
        self.outlet.parameters[gtep.eo] = x[2]
        idx = 3
        for k in self.variables:
            if k not in args:
                setattr(self, k, x[idx])
                idx += 1

        mf_i = self.inlet.parameters[gtep.mf]
        PP_i = self.inlet.parameters[gtep.PP]
        e_i = enthalpy(self.inlet.functions[gtep.Cp], **{gtep.TT: (T0 + 15, self.inlet.parameters[gtep.TT]), gtep.PP: (101325, self.inlet.parameters[gtep.PP])})

        mf_f = self.fuel.parameters[gtep.mf]
        e_f = enthalpy(self.fuel.functions[gtep.C], **{gtep.TT: (T0 + 15, self.fuel.parameters[gtep.TT])})

        mf_o = self.outlet.parameters[gtep.mf]

        return (
            (
                mf_o
                * enthalpy(
                    self.outlet.functions[gtep.Cp],
                    **{
                        gtep.TT: (T0 + 15, self.outlet.parameters[gtep.TT]),
                        gtep.PP: (101325, self.outlet.parameters[gtep.PP]),
                        gtep.eo: (1, self.outlet.parameters[gtep.eo]),
                    },
                )
            )
            - (mf_i * e_i)
            - (mf_f * (e_f + self.efficiency_burn * self.fuel.parameters.get("lower_heating_value", nan))),
            self.outlet.parameters[gtep.PP] - PP_i * getattr(self, gtep.peff),
            1 - (mf_f / mf_o) * self.fuel.parameters.get("stoichiometry", nan) * self.outlet.parameters[gtep.eo],
        )

    def calculate(self, substance_inlet: Substance, fuel: Substance, x0=None) -> Substance:
        GTENode.validate_substance(self, substance_inlet)
        GTENode.validate_substance(self, fuel)
        assert fuel.parameters.get("stoichiometry") is not None, KeyError("fuel has not parameter 'stoichiometry'")
        assert isinstance(fuel.parameters["stoichiometry"], (int, float)), TypeError(f"type fuel stoichiometry must be numeric, but has {type(fuel.parameters['stoichiometry'])}")
        assert fuel.parameters.get("lower_heating_value") is not None, KeyError("fuel has not parameter 'lower_heating_value'")
        assert isinstance(fuel.parameters["lower_heating_value"], (int, float)), TypeError(f"type fuel lower_heating_value must be numeric, but has {type(fuel.parameters['lower_heating_value'])}")
        assert fuel.functions.get(gtep.gc) is not None, KeyError(f"fuel has not function '{gtep.gc}'")

        self.inlet = deepcopy(substance_inlet)
        self.fuel = deepcopy(fuel)
        self.outlet = Substance("exhaust") + fuel

        self.outlet.functions[gtep.gc], self.outlet.functions[gtep.Cp] = self.fuel.functions[gtep.gc], self.fuel.functions[gtep.Cp]

        self.outlet.parameters[gtep.mf] = self.inlet.parameters[gtep.mf] + self.fuel.parameters[gtep.mf] - self.mass_flow_leak

        if x0 is None:
            x0 = tuple(self._x0.values())
        args = {k: v for k, v in self.variables.items() if not isnan(v)}
        count_variables = sum(1 if k not in args else 0 for k in self.variables)
        count_equations = len(self.equations(x0, args)) - 3  # outlet_TT, outlet_PP, outlet_eo

        if count_variables < count_equations:
            raise ArithmeticError(f"{count_variables=} < {count_equations=}")
        elif count_variables > count_equations:
            raise ArithmeticError(f"{count_variables=} > {count_equations=}")

        root(self.equations, x0, args, method="lm")

        self.outlet.parameters[gtep.gc] = call_with_kwargs(self.outlet.functions[gtep.gc], self.outlet.parameters)
        self.outlet.parameters[gtep.Cp] = self.outlet.functions[gtep.Cp](self.outlet.parameters[gtep.TT])
        self.outlet.parameters[gtep.DD] = self.outlet.parameters[gtep.PP] / (self.outlet.parameters[gtep.gc] * self.outlet.parameters[gtep.TT])
        self.outlet.parameters[gtep.k] = adiabatic_index(self.outlet.parameters[gtep.gc], self.outlet.parameters[gtep.Cp])
        self.outlet.parameters[gtep.a_critical] = сritical_sonic_velocity(self.outlet.parameters[gtep.k], self.outlet.parameters[gtep.gc], self.outlet.parameters[gtep.TT])

        return self.outlet

    def validate(self, epsrel: float = EPSREL) -> bool:
        """Проверка найденного решения"""
        x0 = (self.outlet.parameters[gtep.TT], self.outlet.parameters[gtep.PP], self.outlet.parameters[gtep.eo])
        args = {k: v for k, v in self.variables.items() if not isnan(v)}

        result = True
        for i, null in enumerate(self.equations(x0, args)):
            epsilon = eps("rel", null, 0)
            if epsilon > epsrel and not isinf(epsilon):
                result = False
                print(f"{i}: {null:.6f}")

        return result

    @property
    def is_real(self):
        checks = (
            check_efficiency(getattr(self, "efficiency_burn")),
            check_temperature(self.outlet.parameters[gtep.TT]),
            self.inlet.parameters[gtep.TT] <= self.outlet.parameters[gtep.TT],
            0 <= self.outlet.parameters[gtep.eo],
        )
        return all(checks)


if __name__ == "__main__":
    from colorama import Fore
    from fixtures import air, kerosene

    air.parameters[gtep.TT] = 600
    air.parameters[gtep.PP] = 101325 * 6

    for k, v in gtep.items():
        print(f"{k:<10}: {v}")

    cc = CombustionChamber()
    cc.summary

    cc.efficiency_burn = 0.99
    setattr(cc, gtep.peff, 0.95)
    cc.mass_flow_leak = 0

    cc.calculate(air, kerosene)

    cc.summary

    print(Fore.GREEN + f"{cc.validate() = }" + Fore.RESET)
    print(Fore.GREEN + f"{cc.is_real = }" + Fore.RESET)
    print()
