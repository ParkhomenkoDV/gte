import os
import pickle
from copy import deepcopy
from typing import Any, Dict, Tuple

from numpy import isnan, nan
from scipy.optimize import root
from substance import Substance
from thermodynamics import T0, adiabatic_index, heat_capacity_p, heat_capacity_p_exhaust, сritical_sonic_velocity
from thermodynamics import parameters as tdp

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

models = {}
for model in (gtep.TT, gtep.PP, gtep.pipi, gtep.effeff, gtep.power):
    path = f"models/compressor_{model}.pickle"
    if os.path.exists(path):
        models[model] = pickle.load(path)
    else:
        print(f"'{path}' not found!")


class CombustionChamber(GTENode):
    """Камера сгорания"""

    models: Dict[str, Any] = models

    __slots__ = (gtep.effburn, gtep.peff, "fuel")

    def __init__(self, name="CombustionChamber"):
        GTENode.__init__(self, name=name)

        setattr(self, gtep.effburn, nan)
        setattr(self, gtep.peff, nan)

        self.fuel = Substance("fuel")

    @property
    def variables(self) -> Dict[str, float]:
        return {
            gtep.effburn: getattr(self, gtep.effburn),
            gtep.peff: getattr(self, gtep.peff),
        }

    @property
    def x0(self) -> Dict[str, float]:
        """Начальные приближения"""
        x0 = {
            f"outlet_{gtep.TT}": self.inlet.parameters[gtep.TT],
            f"outlet_{gtep.PP}": self.inlet.parameters[gtep.PP],
        }
        for k, v in self.variables.items():
            if not isnan(v):
                continue
            if k == "efficiency_burn":
                x0[k] = 1.0
            elif k == gtep.peff:
                x0[k] = 0.95  # TODO: model or formula
        return x0

    def equations(self, x: Tuple, args: Dict) -> Tuple:
        """Уравнения"""
        self.outlet.parameters[gtep.TT] = x[0]
        self.outlet.parameters[gtep.PP] = x[1]
        idx = 2
        for k in self.variables:
            if k not in args:
                setattr(self, k, x[idx])
                idx += 1

        return (
            # (mf_i * enthalpy_i) + mf_f * (Q*eff + enthalpy_f) = (mf_i + mf_f) * enthalpy_o
            (
                self.outlet.parameters[gtep.mf]
                * enthalpy(
                    self.outlet.functions[gtep.Cp],
                    **{
                        tdp.t: (T0 + 15, self.outlet.parameters[gtep.TT]),
                        tdp.p: (101325, self.outlet.parameters[gtep.PP]),
                        tdp.eo: (self.outlet.parameters[gtep.eo], self.outlet.parameters[gtep.eo]),
                    },
                )
            )
            - (self.inlet.parameters[gtep.mf] * self.inlet.parameters.get("enthalpy", nan))
            - self.fuel.parameters[gtep.mf] * (self.fuel.parameters.get("enthalpy", nan) + self.efficiency_burn * self.fuel.parameters.get("lower_heat", nan)),
            self.outlet.parameters[gtep.PP] - self.inlet.parameters[gtep.PP] * getattr(self, gtep.peff),
        )

    def __validate_fuel(self, fuel: Substance) -> None:
        """Проверка параметров горючего"""
        assert len(fuel.composition) > 0, ValueError(f"{fuel.composition = }")

        stoichiometry = fuel.parameters.get("stoichiometry")
        assert stoichiometry is not None, KeyError("fuel has not parameter 'stoichiometry'")
        assert isinstance(stoichiometry, (int, float)), TypeError(f"type fuel stoichiometry must be numeric, but has {type(stoichiometry)}")

        lower_heat = fuel.parameters.get("lower_heat")
        assert lower_heat is not None, KeyError("fuel has not parameter 'lower_heat'")
        assert isinstance(lower_heat, (int, float)), TypeError(f"type fuel lower_heat must be numeric, but has {type(lower_heat)}")

        assert fuel.functions.get(gtep.gc) is not None, KeyError(f"fuel has not function '{gtep.gc}'")

    def calculate(self, substance_inlet: Substance, fuel: Substance, x0=None) -> Substance:
        GTENode.validate_substance(self, substance_inlet)
        GTENode.validate_substance(self, fuel)
        self.__validate_fuel(fuel)

        self.inlet = deepcopy(substance_inlet)
        self.fuel = deepcopy(fuel)
        self.outlet = Substance("exhaust")

        self.outlet.functions[gtep.gc] = self.fuel.functions[gtep.gc]

        stoichiometry = fuel.parameters["stoichiometry"]
        H2O = self.fuel.composition.get("H2O", 0)  # массовая доля волы в смеси
        self.outlet.functions[gtep.Cp] = lambda temperature, excess_oxidizing: (
            (1 - H2O) * (heat_capacity_p_exhaust(temperature, fuel.composition) + self.inlet.functions[gtep.Cp](temperature) * excess_oxidizing * stoichiometry) + H2O * excess_oxidizing * stoichiometry * heat_capacity_p("H2O", temperature)
        ) / (1 - H2O + excess_oxidizing * stoichiometry)

        self.outlet.parameters[gtep.mf] = self.inlet.parameters[gtep.mf] + self.fuel.parameters[gtep.mf] - self.leak
        self.outlet.parameters[gtep.eo] = self.inlet.parameters[gtep.mf] / self.fuel.parameters[gtep.mf] / self.fuel.parameters["stoichiometry"]

        if x0 is None:
            x0 = tuple(self.x0.values())
        args = {k: v for k, v in self.variables.items() if not isnan(v)}
        count_variables = sum(1 if k not in args else 0 for k in self.variables)
        count_equations = len(self.equations(x0, args)) - 2  # outlet_TT, outlet_PP

        if count_variables < count_equations:
            raise ArithmeticError(f"{count_variables=} < {count_equations=}")
        elif count_variables > count_equations:
            raise ArithmeticError(f"{count_variables=} > {count_equations=}")

        self.inlet.parameters["enthalpy"] = enthalpy(self.inlet.functions[gtep.Cp], **{tdp.t: (T0 + 15, self.inlet.parameters[gtep.TT]), tdp.p: (101325, self.inlet.parameters[gtep.PP])})
        self.fuel.parameters["enthalpy"] = enthalpy(self.fuel.functions[gtep.hc], **{tdp.t: (T0 + 15, self.fuel.parameters[gtep.TT])})

        root(self.equations, x0, args, method="lm")

        outlet_parameters = {tdp.t: self.outlet.parameters.get(gtep.TT), tdp.p: self.outlet.parameters.get(gtep.PP), tdp.eo: self.outlet.parameters.get(gtep.eo)}
        self.outlet.parameters[gtep.gc] = call_with_kwargs(self.outlet.functions[gtep.gc], outlet_parameters)
        self.outlet.parameters[gtep.Cp] = call_with_kwargs(self.outlet.functions[gtep.Cp], outlet_parameters)
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
            if isnan(null) or abs(null) > epsrel:
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
    cc.leak = 0

    cc.calculate(air, kerosene)

    cc.summary

    print(Fore.GREEN + f"{cc.validate() = }" + Fore.RESET)
    print(Fore.GREEN + f"{cc.is_real = }" + Fore.RESET)
    print()
