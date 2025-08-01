from copy import deepcopy

from mathematics import eps
from node import GTENode
from numpy import inf, isnan, nan
from substance import Substance
from thermodynamics import Cp, R_gas, g_cool_BMSTU, mixing_param, η_polytropic

from config import EPSREL, NITER
from config import parameters as gtep


class Turbine:
    """Турбина"""

    def __init__(self, name="Turbine"):
        GTENode.__init__(self, name=name)

        self.pipi = nan
        self.eff = nan
        self.power = nan

    def get_inlet_parameters(self) -> None:
        """Расчет параметров перед"""

        if not hasattr(self, "g1"):
            self.g1 = 1
        self.R_gas1 = R_gas(
            self.substance, a_ox=getattr(self, "a_ox1", None), fuel=fuel
        )
        self.ρρ1 = self.PP1 / (self.R_gas1 * self.TT1)
        self.Cp1 = Cp(
            self.substance,
            T=self.TT1,
            P=self.PP1,
            a_ox=getattr(self, "a_ox1", None),
            fuel=fuel,
        )
        self.k1 = self.Cp1 / (self.Cp1 - self.R_gas1)

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
        epsrel: float = EPSREL,
        niter: int = NITER,
        **kwargs,
    ) -> Substance:
        self.get_inlet_parameters()

        if self._shafts:
            m = kwargs.get("m", {})
            assert m, f"{type(self).__name__} object has no attribute m!"
            L_C = N = 0
            η = 1
            for shaft in self._shafts:
                for tp, els in shaft.items():
                    if tp == "-":
                        for el in els:
                            if hasattr(el, "L"):
                                L_C += el.L * m[find_node_in_scheme(scheme, el)[0]]
                            elif hasattr(el, "N"):
                                N += shaft[-1].N
                            else:
                                raise "!"
                    elif tp == "0":
                        for el in els:
                            if hasattr(el, "η"):
                                η *= el.η
                    else:
                        raise "!"
            self.L = (
                L_C / η / self.η_mechanical
            )  # удельная работа от компрессоров TODO сделать КПД мех атрибутом вала
            G = kwargs.get("G", inf)
            if isnan(G) or G != 0:
                G = inf
            L_N = (
                (N / (G / m[find_node_in_scheme(scheme, self)[0]]))
                / η
                / self.η_mechanical
            )  # удельная работа от нагрузок
            self.L += L_N  # удельная работа от нагрузок
            assert hasattr(self, "L") or L_N, (
                f"{type(self).__name__} object has no attributes L and/or N!"
            )
            self.L /= m[find_node_in_scheme(scheme, self)[0]]
        assert hasattr(self, "L") or (hasattr(self, "N") and hasattr(self, "G")), (
            f"{type(self).__name__} object has no attributes L and/or (N and G)"
        )

        g_leaks = self.g_leak
        g_fuels = 0
        if scheme:
            for node in range(n):
                g_leaks += scheme[c][node].g_leak  # суммарные утечки до турбины
                g_fuels += (
                    scheme[c][node].g_fuel if hasattr(scheme[c][node], "g_fuel") else 0
                )
        assert hasattr(self, "T_lim"), (
            f"{type(self).__name__} object has no attribute T_lim!"
        )
        self.g_cool = g_cool_BMSTU(self.TT1, T_lim=self.T_lim)
        self.g_cool = self.g_cool * (1 - g_leaks) / (1 + self.g_cool - g_fuels)

        self.get_outlet_parameters(how=how, error=epsrel, Niter=niter)


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

    t.L = 250000
    t.ηη = 0.91

    t.calculate(substance_inlet)
    for k, v in t.__dict__.items():
        print(k, "=", v)
