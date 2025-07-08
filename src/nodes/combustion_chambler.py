from thermodynamics import (
    T0,
    atmosphere_standard,
    Cp,
    R_gas,
    Qa1,
    l_stoichiometry,
    g_cool_BMSTU,
)
from colorama import Fore

from node import Node

from tools import export2, isnum, COOR, Axis, angle, rounding, eps, dist, dist2line
from decorators import timeit


def find_node_in_scheme(scheme, node2find) -> tuple:
    """Поиск положения узла в схеме"""
    for contour in scheme:
        for i, node in enumerate(scheme[contour]):
            if scheme[contour][i] is node2find:
                return contour, i


class CombustionChamber(Node):
    """Камера сгорания"""

    def __init__(self, name="CombustionChamber"):
        Node.init(self)
        self.name = name
        self.warnings = {0: set(), 1: set(), 2: set(), 3: set()}  # предупреждения

    def __str__(self) -> str:
        return self.name

    def set_combination(self, combination, combustionchamber_main) -> None:
        """Установка комбинации"""
        varible_params = [
            key
            for key, value in combustionchamber_main.__dict__.items()
            if type(value) is list and len(value) and not key.startswith("_")
        ]
        positions = [0] * len(varible_params)

        for i in range(combination):
            for j, param in enumerate(varible_params):
                if (
                    positions[j]
                    == len(getattr(combustionchamber_main, varible_params[j])) - 1
                ):
                    positions[j] = 0
                else:
                    positions[j] += 1
                    continue
        for j, param in enumerate(varible_params):
            setattr(
                self,
                varible_params[j],
                getattr(combustionchamber_main, varible_params[j])[positions[j]],
            )

    def input_parameters(self):
        # TODO сделать ввод либо ТТ3 либо a_ox3
        # TT3 <-> g_fuel <-> a_ox
        correct_input = False
        while not correct_input:
            self.TT3_var = [TT3 for TT3 in input("T*3 [К] = ").split()]
            for TT3 in self.TT3_var:
                if not isnum(TT3) or float(TT3) < 0:
                    print(Fore.RED + "T*3 must be a number > 0")
                    correct_input = False
                    break
                correct_input = True
        self.TT3_var = [float(TT3) for TT3 in self.TT3_var]

        correct_input = False
        while not correct_input:
            self.η_burn_var = [burn for burn in input("η горения [] = ").split()]
            for burn in self.η_burn_var:
                if not isnum(burn) or float(burn) < 0 or float(burn) > 1:
                    print(Fore.RED + "η_burn must be a number in [0..1]!")
                    correct_input = False
                    break
                correct_input = True
        self.η_burn_var = [float(burn) for burn in self.η_burn_var]

        correct_input = False
        while not correct_input:
            self.T_fuel_var = [T_fuel for T_fuel in input("T горючего [К] = ").split()]
            for T_fuel in self.T_fuel_var:
                if not isnum(T_fuel) or float(T_fuel) < 0:
                    print(Fore.RED + "T_fuel must be a number > 0")
                    correct_input = False
                    break
                correct_input = True
        self.T_fuel_var = [float(Tf) for Tf in self.T_fuel_var]

        correct_input = False
        while not correct_input:
            self.σ_var = [sigma for sigma in input("σ [] = ").split()]
            for sigma in self.σ_var:
                if not isnum(sigma) or float(sigma) < 0 or float(sigma) > 1:
                    print(Fore.RED + "σ must be a number in [0..1]!")
                    correct_input = False
                    break
                correct_input = True
        self.σ_var = [float(sigma) for sigma in self.σ_var]

        correct_input = False
        while not correct_input:
            self.g_leak_var = [gleak for gleak in input("g утечки [] = ").split()]
            for gleak in self.g_leak_var:
                if not isnum(gleak) or float(gleak) < 0 or float(gleak) > 1:
                    print(Fore.RED + "g_leak must be a number in [0..1]!")
                    correct_input = False
                    break
                correct_input = True
        self.g_leak_var = [float(sigma) for sigma in self.g_leak_var]

    def get_inlet_parameters(self) -> None:
        """Расчет параметров перед"""
        if scheme:
            if hasattr(scheme[c][n - 1], "a_ox3"):
                self.a_ox1 = scheme[c][n - 1].a_ox3
            self.TT1 = scheme[c][n - 1].TT3
            self.PP1 = scheme[c][n - 1].PP3
            self.g1 = scheme[c][n - 1].g3
        assert hasattr(self, "TT1"), (
            f"{type(self).__name__} object has no attribute TT1!"
        )
        assert hasattr(self, "PP1"), (
            f"{type(self).__name__} object has no attribute PP1!"
        )
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

    def get_outlet_parameters(self, how="all") -> None:
        """Расчет параметров после"""
        assert hasattr(self, "σ"), f"{type(self).__name__} object has no attribute σ!"
        self.PP3 = self.PP1 * self.σ

        if hasattr(self, "TT3"):  # для КС
            if self.TT1 > self.TT3:
                self.warnings[3].add("T*1_КС > T*3_КС!")
                return
            assert hasattr(self, "η_burn"), (
                f"{type(self).__name__} object has no attribute η_burn!"
            )
            assert hasattr(self, "T_fuel"), (
                f"{type(self).__name__} object has no attribute T_fuel!"
            )
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
        else:
            raise f"{type(self).__name__} object has no attributes TT3 or a_ox3 or g_fuel!"

        assert hasattr(self, "g_leak"), (
            f"{type(self).__name__} object has no attribute g_leak!"
        )
        g_leaks = self.g_leak  # суммарные утечки до КС (включительно)
        if scheme:
            for node in range(n):
                g_leaks += scheme[c][node].g_leak
        assert hasattr(self, "T_lim"), (
            f"{type(self).__name__} object has no attribute T_lim!"
        )
        g_cools = g_cool_BMSTU(
            self.TT3, T_lim=self.T_lim
        )  # суммарные массовый расход на охлаждение после КС
        self.g_fuel = g_fuel3 * (1 - g_cools - g_leaks) / (1 - g_fuel3)

        self.g3 = self.g1 - self.g_leak - g_cools + self.g_fuel
        self.R_gas3 = R_gas("EXHAUST", a_ox=self.a_ox3, fuel=fuel)
        self.ρρ3 = self.PP3 / (self.R_gas3 * self.TT3)
        self.Cp3 = Cp("EXHAUST", T=self.TT3, P=self.PP3, a_ox=self.a_ox3, fuel=fuel)
        self.k3 = self.Cp3 / (self.Cp3 - self.R_gas3)
        self.substance = "EXHAUST"

    def __calculate(self, how="all", **kwargs) -> None:
        global c, n, scheme, fuel
        c = n = None
        scheme = kwargs.get("scheme", {})

        if scheme:
            c, n = find_node_in_scheme(scheme, self)
            self.substance = scheme[c][n - 1].substance
        else:
            self.substance = kwargs.get("substance", "")  # рабочее тело
        assert self.substance, (
            f"{type(self).__name__} object has no attribute substance!"
        )
        fuel = kwargs.get("fuel", "")  # горючее
        assert fuel, (
            f"{CombustionChamber.__name__}.{CombustionChamber.solve.__name__}() method has no arguments fuel!"
        )

        self.get_inlet_parameters()
        self.get_outlet_parameters(how=how)

    def solve(self, how="all", error=0.01, Niter=100, **kwargs) -> None:
        self.__calculate(how=how, error=error, Niter=Niter, **kwargs)


if __name__ == "__main__":
    cc = CombustionChamber()

    cc.TT1 = 600
    cc.PP1 = 101325 * 6

    cc.T_fuel = 40 + T0
    cc.η_burn = 0.99
    cc.TT3 = 1600
    cc.T_lim = 1800
    cc.σ = 0.94
    cc.g_leak = 0

    cc.solve(substance="AIR", fuel="КЕРОСИН")
    for k, v in cc.__dict__.items():
        print(k, "=", v)
