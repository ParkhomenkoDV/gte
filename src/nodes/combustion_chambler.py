from copy import deepcopy

from node import GTENode
from numpy import nan
from substance import Substance
from thermodynamics import (
    T0,
    # Cp,
    # Qa1,
    gas_const,
    # atmosphere_standard,
    stoichiometry,
)

from src.config import parameters as gtep


class CombustionChamber(GTENode):
    """Камера сгорания"""

    def __init__(self, name="CombustionChamber"):
        GTENode.__init__(self, name=name)

        self.total_pressure_loss = nan
        self.efficiency_burn = nan
        self.total_pressure_loss = 0

    def get_outlet_parameters(self) -> None:
        """if hasattr(self, "TT3"):  # для КС
            g_fuel3 = (
                Cp(
                    self.substance,
                    T=self.TT3,
                    P=self.PP3,
                    a_ox=getattr(self, "a_ox1", None),
                    fuel=fuel,
                )
                * self.TT3
            )
            g_fuel3 -= self.Cp1 * self.TT1
            g_fuel3 /= (
                Qa1(fuel) * self.η_burn
                - (1 + l_stoichiometry(fuel))
                * (
                    Cp("EXHAUST", fuel=fuel, a_ox=1, T=self.TT3) * self.TT3
                    - Cp("EXHAUST", fuel=fuel, a_ox=1, T=(T0 + 15)) * (T0 + 15)
                )
                + l_stoichiometry(fuel) * (Cp(self.substance, self.TT3)) * self.TT3
                - Cp(self.substance, (T0 + 15)) * (T0 + 15)
                + (
                    Cp(fuel, T=self.T_fuel) * self.T_fuel
                    - Cp(fuel, T=T0 + 15) * (T0 + 15)
                )
            )
            self.a_ox3 = 1 / g_fuel3 / l_stoichiometry(fuel)
            if self.a_ox3 < 1:
                self.warnings[1].add("a_ox < 1! Топливо не сгорает полностью!")
        elif hasattr(self, "a_ox3"):  # для ФК
            self.g_fuel = 1 / l_stoichiometry(fuel) / self.a_ox3
        elif hasattr(self, "g_fuel"):
            self.a_ox3 = 1 / self.g_fuel / l_stoichiometry(fuel)

        self.Cp3 = Cp("EXHAUST", T=self.TT3, P=self.PP3, a_ox=self.a_ox3, fuel=fuel)
        self.k3 = self.Cp3 / (self.Cp3 - self.R_gas3)"""

    def calculate(
        self,
        substance_inlet: Substance,
        fuel: Substance,
        epsrel: float = 0.01,
        niter: int = 10,
        **kwargs,
    ) -> Substance:
        GTENode.validate_substance(self, substance_inlet)
        GTENode.validate_substance(self, fuel)
        self.inlet = deepcopy(substance_inlet)
        self.outlet = Substance("exhaust") + fuel

        # define functions
        self.outlet.functions[gtep.gc] = lambda excess_oxidizing: gas_const(
            "EXHAUST", excess_oxidizing, fuel.name
        )
        self.outlet.functions[gtep.Cp] = lambda TT: nan

        self.outlet.parameters[gtep.mf] = (
            self.inlet.parameters[gtep.mf] + fuel.parameters[gtep.mf]
        )

        if "TT_max" in kwargs:
            self.outlet.parameters[gtep.TT] = kwargs["TT_max"]
            assert self.inlet.parameters[gtep.TT] <= self.outlet.parameters[gtep.TT], (
                ArithmeticError("inlet temperature > outlet temperature")
            )
            self.outlet.parameters[gtep.eo] = 1 / (
                fuel.parameters[gtep.mf]
                / self.outlet.parameters[gtep.mf]
                * stoichiometry(fuel.name)
            )
        elif gtep.eo in kwargs:
            self.outlet.parameters[gtep.eo] = kwargs[gtep.eo]

        else:
            raise AssertionError("CombustionChamber must have TT_max")

        self.outlet.parameters[gtep.gc] = self.outlet.functions[gtep.gc](
            self.outlet.parameters[gtep.eo]
        )
        self.outlet.parameters[gtep.PP] = (
            self.inlet.parameters[gtep.PP] * self.total_pressure_loss
        )
        self.outlet.parameters[gtep.DD] = self.outlet.parameters[gtep.PP] / (
            self.outlet.parameters[gtep.gc] * self.outlet.parameters[gtep.TT]
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
        functions={},
    )

    cc = CombustionChamber()
    cc.total_pressure_loss = 0.96
    cc.efficiency_burn = 0.99
    cc.mass_flow_leak = 0

    cc.calculate(substance_inlet, fuel, TT_max=1800)

    for k, v in cc.summary.items():
        print(f"{k:<40}: {v}")
