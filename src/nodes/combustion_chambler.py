from copy import deepcopy

import numpy as np
from numpy import isnan, nan
from scipy import integrate
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

from src.config import parameters as gtep
from src.nodes import GTENode


def average_integral(f, *borders) -> float:
    """Среднеинтегральное значение"""
    for border in borders:
        assert isinstance(border, (tuple, list))
        for b in border:
            assert isinstance(b, (int, float, np.number)), f"{type(b)}"

    if borders[0][0] == borders[0][1]:
        return f(borders[0][0])
    else:
        return integrate.quad(f, borders[0][0], borders[0][1])[0] / (
            borders[0][1] - borders[0][0]
        )


class CombustionChamber(GTENode):
    """Камера сгорания"""

    def __init__(self, name="CombustionChamber"):
        GTENode.__init__(self, name=name)

        self.fuel = Substance("fuel")

        self.efficiency_burn = nan
        self.total_pressure_loss = 0

    @property
    def __x0(self) -> dict[str:float]:
        """Начальные приближения"""
        x0 = {
            "outlet_" + gtep.TT: self.inlet.parameters[gtep.TT],
            "outlet_" + gtep.PP: self.inlet.parameters[gtep.PP],
        }

        if isnan(self.outlet.parameters[gtep.eo]):
            x0["outlet_" + gtep.eo] = 2  # TODO: model
        else:
            raise "недоопределено"

        return x0

    def equations(self, x0: tuple, args: dict) -> tuple:
        """Уравнения"""
        self.outlet.parameters[gtep.TT] = x0[0]
        self.outlet.parameters[gtep.PP] = x0[1]
        if gtep.eo not in args:
            self.outlet.parameters[gtep.eo] = x0[2]
        elif gtep.eo in args:
            pass

        mf_i = self.inlet.parameters[gtep.mf]
        TT_i = self.inlet.parameters[gtep.TT]
        PP_i = self.inlet.parameters[gtep.PP]
        f_Cp_i = self.inlet.functions[gtep.Cp]
        e_i = average_integral(f_Cp_i, (T0, TT_i))

        name_f = self.fuel.name
        mf_f = self.fuel.parameters[gtep.mf]
        TT_f = self.fuel.parameters[gtep.TT]
        f_Cp_f = self.fuel.functions[gtep.Cp]
        e_f = average_integral(f_Cp_f, (T0, TT_f))

        mf_o = self.outlet.parameters[gtep.mf]
        f_Cp_o = self.outlet.functions[gtep.Cp]

        return (
            (mf_o * average_integral(f_Cp_o, (T0, self.outlet.parameters[gtep.TT])))
            - (mf_i * e_i)
            - (mf_f * (e_f + self.efficiency_burn * lower_heating_value(name_f))),
            self.outlet.parameters[gtep.PP] - PP_i * (1 - self.total_pressure_loss),
            1 - (mf_f / mf_o) * stoichiometry(name_f) * self.outlet.parameters[gtep.eo],
        )

    def calculate(
        self, substance_inlet: Substance, fuel: Substance, x0=None
    ) -> Substance:
        GTENode.validate_substance(self, substance_inlet)
        GTENode.validate_substance(self, fuel)
        self.inlet = deepcopy(substance_inlet)
        self.fuel = deepcopy(fuel)
        self.outlet = Substance("exhaust") + fuel

        self.outlet.functions[gtep.gc] = lambda excess_oxidizing: gas_const(
            "EXHAUST", excess_oxidizing, fuel.name
        )
        self.outlet.functions[gtep.Cp] = (
            lambda TT, excess_oxidizing: heat_capacity_at_constant_pressure(
                "EXHAUST",
                TT,
                excess_oxidizing=excess_oxidizing,
                fuel=fuel.name,
            )
        )
        self.outlet.parameters[gtep.eo] = nan

        self.outlet.parameters[gtep.mf] = (
            self.inlet.parameters[gtep.mf]
            + self.fuel.parameters[gtep.mf]
            - self.mass_flow_leak
        )

        args = {}
        if not isnan(self.efficiency_burn) and not isnan(self.total_pressure_loss):
            args.update(
                {
                    "efficiency_burn": self.efficiency_burn,
                    "total_pressure_loss": self.total_pressure_loss,
                }
            )
        else:
            raise

        fsolve(self.equations, tuple(self.__x0.values()), args)

        assert self.inlet.parameters[gtep.TT] <= self.outlet.parameters[gtep.TT], (
            Exception("inlet temperature > outlet temperature")
        )
        assert 0 < self.outlet.parameters[gtep.eo], Exception(f"outlet {gtep.eo} < 0")

        self.outlet.parameters[gtep.gc] = self.outlet.functions[gtep.gc](
            self.outlet.parameters[gtep.eo]
        )
        self.outlet.parameters[gtep.PP] = (
            self.inlet.parameters[gtep.PP] * self.total_pressure_loss
        )
        self.outlet.parameters[gtep.DD] = self.outlet.parameters[gtep.PP] / (
            self.outlet.parameters[gtep.gc] * self.outlet.parameters[gtep.TT]
        )
        self.outlet.parameters[gtep.Cp] = self.outlet.functions[gtep.Cp](
            self.outlet.parameters[gtep.TT], self.outlet.parameters[gtep.eo]
        )
        self.outlet.parameters[gtep.k] = adiabatic_index(
            self.outlet.parameters[gtep.gc], self.outlet.parameters[gtep.Cp]
        )

        return self.outlet


if __name__ == "__main__":
    for k, v in gtep.items():
        print(f"{k:<10}: {v}")

    substance_inlet = Substance(
        "air",
        parameters={
            gtep.gc: 287,
            gtep.TT: 300 * 2,
            gtep.PP: 101_325 * 6,
            gtep.mf: 50,
            gtep.Cp: 1006,
            gtep.k: 1.4,
        },
        functions={
            gtep.Cp: lambda T: 1006,
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
            gtep.Cp: lambda T: 200,
        },
    )

    cc = CombustionChamber()
    cc.summary

    cc.total_pressure_loss = 0.04
    cc.efficiency_burn = 0.99
    cc.mass_flow_leak = 0

    cc.calculate(substance_inlet, fuel)

    cc.summary

    for k, v in cc.summary.items():
        print(f"{k:<50}: {v}")
