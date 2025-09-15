from copy import deepcopy

from mathematics import eps
from node import GTENode
from numpy import inf, isnan, nan
from scipy import integrate
from substance import Substance
from thermodynamics import (
    Cp,
    R_gas,
    adiabatic_index,
    mixing_param,
    η_polytropic,
)

from config import EPSREL, NITER
from config import parameters as gtep


class Turbine:
    """Турбина"""

    def __init__(self, name="Turbine"):
        GTENode.__init__(self, name=name)

        self.pipi = nan
        self.eff = nan
        self.power = nan

    def equations(self, x, *args: tuple) -> tuple:
        """Уравнения"""
        TT_o, PP_o, power, pipi, gc, Cp = x
        mf_i, TT_i, PP_i, f_gas_const, f_Cp = args
        effeff = 0.90
        return (
            power - mf_i * Cp * (TT_o - TT_i),
            TT_o - TT_i * (1 + (pipi ** ((k - 1) / k) - 1) / effeff),
            PP_i - PP_i * pipi,
            gc - integrate.quad(f_gas_const, TT_i, TT_o)[0],
            Cp - integrate.quad(f_Cp, TT_i, TT_o)[0],
            k - adiabatic_index(gc, Cp),
        )

    def get_outlet_parameters(self, error=0.01, Niter=100, **kwargs):
        """Расчет параметров после"""
        if hasattr(self, "a_ox1"):
            self.a_ox3 = self.a_ox1 * (1 + self.g_cool)
        self.R_gas3 = R_gas(
            self.substance, a_ox=getattr(self, "a_ox1", None), fuel=fuel
        )

        assert hasattr(self, "ηη"), f"{type(self).__name__} object has no attribute ηη!"
        self.k3 = self.k1  # нулевое приближение
        if self.L:
            for iteration in range(Niter):
                k2 = 0.5 * (self.k1 + self.k3)
                Cp2 = k2 / (k2 - 1) * R_gas(self.substance, a_ox=self.a_ox3, fuel=fuel)
                assert hasattr(self, "ηη"), (
                    f"{type(self).__name__} object has no attribute ηη!"
                )
                self.ππ = (1 - self.L / (Cp2 * self.TT1 * self.ηη)) ** (k2 / (1 - k2))
                self.TT3 = self.TT1 - self.L / Cp2
                # if self.g_cool: self.TT3 = mixing_param([self.TT3, *TT_cools], [], [])
                self.Cp3 = Cp(self.substance, T=self.TT3, a_ox=self.a_ox3, fuel=fuel)
                if isnan(k2):
                    self.warnings[3].add(
                        f"T*3 < 0 [К]! Не хватает теплоперепада в турбине_{self}!"
                    )
                    return
                if (
                    abs(eps("rel", self.Cp3 / (self.Cp3 - self.R_gas3), self.k3))
                    <= error
                ):
                    break
                self.k3 = self.Cp3 / (self.Cp3 - self.R_gas3)
            else:
                print(
                    f"Iteration limit in module {Turbine.__name__} in function {self.solve.__name__} L!"
                )
            self.PP3 = self.PP1 / self.ππ
        elif hasattr(self, "PP3"):
            self.ππ = self.PP1 / self.PP3
            for iteration in range(Niter):
                k2 = 0.5 * (self.k1 + self.k3)
                self.TT3 = self.TT1 * (1 - (1 - self.ππ ** ((1 - k2) / k2)) * self.ηη)
                self.Cp3 = Cp(self.substance, T=self.TT3, a_ox=self.a_ox3, fuel=fuel)
                if (
                    abs(eps("rel", self.Cp3 / (self.Cp3 - self.R_gas3), self.k3))
                    <= error
                ):
                    break
                self.k3 = self.Cp3 / (self.Cp3 - self.R_gas3)
            else:
                print(
                    f"Iteration limit in module {Turbine.__name__} in function {self.solve.__name__} PP3!"
                )
            self.L = self.Cp1 * self.TT1 - self.Cp3 * self.TT3
        k2 = 0.5 * (self.k1 + self.k3)
        self.ηn = η_polytropic(what="E", ππ=self.ππ, ηη=self.ηη, k=k2)
        self.ρρ3 = self.PP3 / (self.R_gas3 * self.TT3)
        assert hasattr(self, "g_leak"), (
            f"{type(self).__name__} object has no attribute g_leak!"
        )
        self.g3 = self.g1 - self.g_leak

    def calculate(
        self,
        substance_inlet: Substance,
        epsrel: float = EPSREL,
        niter: int = NITER,
        **kwargs,
    ) -> Substance:
        GTENode.validate_substance(self, substance_inlet)
        self.inlet = deepcopy(substance_inlet)
        self.outlet = deepcopy(self.inlet)

        self.outlet.parameters[gtep.PP] = self.inlet.parameters[gtep.PP] / self.pipi


if __name__ == "__main__":
    for k, v in gtep.items():
        print(f"{k:<10}: {v}")

    substance_inlet = Substance(
        "exhaust",
        parameters={
            gtep.gc: 287,
            gtep.TT: 1700,
            gtep.PP: 101_325 * 5,
            gtep.mf: 51,
            gtep.Cp: 1206,
            gtep.k: 1.33,
        },
        functions={
            gtep.Cp: lambda T: 1006,
        },
    )

    t = Turbine()

    t.power = 16_000_000
    t.eff = 0.91

    t.calculate(substance_inlet)
    for k, v in t.__dict__.items():
        print(k, "=", v)
