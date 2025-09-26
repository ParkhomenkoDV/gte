from copy import deepcopy

from mathematics import eps
from numpy import isinf, isnan, nan
from scipy.optimize import fsolve
from substance import Substance
from thermodynamics import (
    T0,
    adiabatic_index,
    gas_const,
    heat_capacity_at_constant_pressure,
    lower_heating_value,
    stoichiometry,
)

from src.checks import check_efficiency, check_temperature
from src.config import EPSREL
from src.config import parameters as gtep
from src.nodes.node import GTENode
from src.utils import call_with_kwargs, enthalpy


class CombustionChamber(GTENode):
    """Камера сгорания"""

    def __init__(self, name="CombustionChamber"):
        GTENode.__init__(self, name=name)

        self.efficiency_burn = nan
        setattr(self, gtep.peff, nan)

    @property
    def variables(self) -> dict:
        return {
            "efficiency_burn": self.efficiency_burn,
            gtep.peff: getattr(self, gtep.peff),
        }

    @property
    def _x0(self) -> dict[str:float]:
        """Начальные приближения"""
        x0 = {
            f"outlet_{gtep.TT}": self.inlet.parameters[gtep.TT],
            f"outlet_{gtep.PP}": self.inlet.parameters[gtep.PP],
            f"outlet_{gtep.eo}": 1,
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
        idx = 2
        for k in self.variables:
            if k not in args:
                setattr(self, k, x[idx])
                idx += 1

        mf_i = self.inlet.parameters[gtep.mf]
        TT_i = self.inlet.parameters[gtep.TT]
        PP_i = self.inlet.parameters[gtep.PP]
        e_i = enthalpy(self.inlet.functions[gtep.Cp], **{gtep.TT: (0, TT_i), gtep.PP: (0, PP_i)})

        mf_f = self.fuel.parameters[gtep.mf]
        TT_f = self.fuel.parameters[gtep.TT]
        e_f = enthalpy(self.fuel.functions[gtep.Cp], **{gtep.TT: (0, TT_f)})

        mf_o = self.outlet.parameters[gtep.mf]

        return (
            (
                mf_o
                * enthalpy(
                    self.outlet.functions[gtep.Cp],
                    **{
                        gtep.TT: (0, self.outlet.parameters[gtep.TT]),
                        gtep.eo: (0, self.outlet.parameters[gtep.eo]),
                    },
                )
            )
            - (mf_i * e_i)
            - (mf_f * (e_f + self.efficiency_burn * lower_heating_value(self.fuel.name))),
            self.outlet.parameters[gtep.PP] - PP_i * getattr(self, gtep.peff),
            1 - (mf_f / mf_o) * stoichiometry(self.fuel.name) * self.outlet.parameters[gtep.eo],
        )

    def calculate(self, substance_inlet: Substance, fuel: Substance, x0=None) -> Substance:
        GTENode.validate_substance(self, substance_inlet)
        GTENode.validate_substance(self, fuel)
        self.inlet = deepcopy(substance_inlet)
        self.fuel = deepcopy(fuel)
        self.outlet = Substance("exhaust") + fuel

        self.outlet.functions[gtep.gc] = lambda excess_oxidizing: gas_const("EXHAUST", excess_oxidizing, fuel=fuel.name)
        self.outlet.functions[gtep.Cp] = lambda total_temperature, excess_oxidizing: heat_capacity_at_constant_pressure(
            "EXHAUST",
            total_temperature,
            excess_oxidizing,
            fuel=fuel.name,
        )

        self.outlet.parameters[gtep.mf] = self.inlet.parameters[gtep.mf] + self.fuel.parameters[gtep.mf] - self.mass_flow_leak
        self.outlet.parameters[gtep.eo] = 1

        x0 = tuple(self._x0.values())
        args = {k: v for k, v in self.variables.items() if not isnan(v)}
        count_variables = sum(1 if k not in args else 0 for k in self.variables)
        count_equations = len(self.equations(x0, args)) - 3  # outlet_TT, outlet_PP, outlet_eo

        if count_variables < count_equations:
            raise ArithmeticError(f"{count_variables=} < {count_equations=}")
        elif count_variables > count_equations:
            raise ArithmeticError(f"{count_variables=} > {count_equations=}")

        fsolve(self.equations, x0, args)

        assert self.inlet.parameters[gtep.TT] <= self.outlet.parameters[gtep.TT], Exception("inlet temperature > outlet temperature")
        assert 0 < self.outlet.parameters[gtep.eo], Exception(f"outlet {gtep.eo} < 0")

        self.outlet.parameters[gtep.gc] = call_with_kwargs(self.outlet.functions[gtep.gc], self.outlet.parameters)
        self.outlet.parameters[gtep.Cp] = self.outlet.functions[gtep.Cp](self.outlet.parameters[gtep.TT], self.outlet.parameters[gtep.eo])
        self.outlet.parameters[gtep.DD] = self.outlet.parameters[gtep.PP] / (self.outlet.parameters[gtep.gc] * self.outlet.parameters[gtep.TT])
        self.outlet.parameters[gtep.k] = adiabatic_index(self.outlet.parameters[gtep.gc], self.outlet.parameters[gtep.Cp])

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
        )
        return all(checks)


if __name__ == "__main__":
    for k, v in gtep.items():
        print(f"{k:<10}: {v}")

    substance_inlet = Substance(
        "air",
        parameters={
            gtep.gc: 287.14,
            gtep.TT: 300 * 2,
            gtep.PP: 101_325 * 6,
            gtep.mf: 100,
            gtep.Cp: 1006,
            gtep.k: 1.4,
            gtep.c: 0,
        },
        functions={
            gtep.gc: lambda total_temperature: gas_const("AIR"),
            gtep.Cp: lambda total_temperature: heat_capacity_at_constant_pressure("AIR", total_temperature),
        },
    )

    fuel = Substance(
        "kerosene",
        parameters={
            gtep.TT: 40 + T0,
            gtep.PP: 101_325,
            gtep.mf: 1,
        },
        functions={
            gtep.Cp: lambda total_temperature: 200,
        },
    )

    cc = CombustionChamber()
    cc.summary

    cc.efficiency_burn = 0.99
    setattr(cc, gtep.peff, 0.95)
    cc.mass_flow_leak = 0

    cc.calculate(substance_inlet, fuel)

    cc.summary

    print(f"{cc.validate() = }")
    print(f"{cc.is_real = }")
