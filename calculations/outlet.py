from numpy import nan, inf
from thermodynamics import T0, Cp, R_gas
from tools import eps, isnum
from colorama import Fore

from combustion_chambler import CombustionChamber  # камера сгорания


def find_node_in_scheme(scheme, node) -> tuple:
    """Поиск положения узла в схеме"""
    for contour in scheme:
        for i, _ in enumerate(scheme[contour]):
            if scheme[contour][i] is node:
                return contour, i


class Outlet:
    """Выходное устройство"""

    def __init__(self):
        self.TT1 = nan  # полная температура перед
        self.TT3 = nan  # полная температура после
        self.PP1 = nan  # полное давление перед
        self.PP3 = nan  # полное давление после
        self.ρρ1 = nan  # полная плотность перед
        self.ρρ3 = nan  # полная плотность после

        self.T1 = nan  # статическая температура перед
        self.T3 = nan  # статическая температура после
        self.P1 = nan  # статическое давление перед
        self.P3 = nan  # статическое давление после
        self.ρ1 = nan  # статическая плотность перед
        self.ρ3 = nan  # статическая плотность после
        self.c1 = nan  # абсолютная скорость перед
        self.c3 = nan  # абсолютная скорость после
        self.g1 = nan  # относительный массовый расход перед
        self.g3 = nan  # относительный массовый расход после

        self.Cp1 = nan  # теплоемкость при постоянном давлении перед
        self.Cp3 = nan  # теплоемкость при постоянном давлении после
        self.R_gas1 = nan  # газовая постоянная перед
        self.R_gas3 = nan  # газовая постоянная после
        self.k1 = nan  # показатель адиабаты перед
        self.k3 = nan  # показатель адиабаты после
        self.a_ox1 = inf  # коэффициент избытка окислителя перед
        self.a_ox3 = inf  # коэффициент избытка окислителя после

        self.warnings = set()  # предупреждения

    def get_variability(self):
        result = 1
        for attr in (self.σ, self.g_leak):
            if type(attr) is list: result *= len(attr)
        return result

    def set_combination(self, combination, outlet_main):
        positions = [0] * 2
        for i in range(combination):
            if positions[0] == len(outlet_main.σ) - 1:
                positions[0] = 0
            else:
                positions[0] += 1
                continue

            if positions[1] == len(outlet_main.g_leak) - 1:
                positions[1] = 0
            else:
                positions[1] += 1
                continue

        self.σ = outlet_main.σ[positions[0]]
        self.g_leak = outlet_main.g_leak[positions[1]]

    def input_parameters(self):
        correct_input = False
        while not correct_input:
            self.σ_var = [sigma for sigma in input('σ [] = ').split()]
            for sigma in self.σ_var:
                if not isnum(sigma) or float(sigma) < 0 or float(sigma) > 1:
                    print(Fore.RED + 'σ must be a number in [0..1]!')
                    correct_input = False
                    break
                correct_input = True
        self.σ_var = [float(sigma) for sigma in self.σ_var]

        correct_input = False
        while not correct_input:
            self.g_leak_var = [gleak for gleak in input('g утечки [] = ').split()]
            for gleak in self.g_leak_var:
                if not isnum(gleak) or float(gleak) < 0 or float(gleak) > 1:
                    print(Fore.RED + 'g_leak must be a number in [0..1]!')
                    correct_input = False
                    break
                correct_input = True
        self.g_leak_var = [float(gleak) for gleak in self.g_leak_var]

    def solve(self, error=0.01, Niter=100, **kwargs):

        substance = kwargs.get('substance', '')  # рабочее тело
        fuel = kwargs.get('fuel', '')  # горючее

        self.g1 = 1
        self.g3 = 1 - self.g_leak

        scheme = kwargs.get('scheme', {})
        if scheme:
            c, n = find_node_in_scheme(scheme, self)
            self.a_ox1 = scheme[c][n - 1].a_ox3
            self.TT1 = scheme[c][n - 1].TT3
            self.PP1 = scheme[c][n - 1].PP3
            for node in scheme[c][:n]:
                self.g1 -= node.g_leak
                self.g1 += node.g_fuel if hasattr(node, 'g_fuel') else 0
                self.g3 -= node.g_leak
                self.g3 += node.g_fuel if hasattr(node, 'g_fuel') else 0

            substance = 'AIR'
            for node in scheme[c][:n]:
                if type(node) is CombustionChamber:
                    substance = 'EXHAUST'
                    break

        self.a_ox1 = kwargs.get('a_ox1', nan) if self.a_ox1 is nan else self.a_ox1  # коэффициент избытка окислителя
        # перед
        self.R_gas1 = R_gas(substance, a_ox=self.a_ox1, fuel=fuel)
        self.TT1 = kwargs.get('TT3', nan) if self.TT1 is nan else self.TT1
        self.PP1 = kwargs.get('PP3', nan) if self.PP1 is nan else self.PP1
        self.ρρ1 = self.PP1 / (self.R_gas1 * self.TT1)
        self.Cp1 = Cp(substance, T=self.TT1, a_ox=self.a_ox1, fuel=fuel)
        self.k1 = self.Cp1 / (self.Cp1 - self.R_gas1)
        self.c1 = kwargs.get('c3', nan)
        # после
        self.a_ox3 = self.a_ox1
        self.R_gas3 = R_gas(substance, a_ox=self.a_ox3, fuel=fuel)
        self.TT3 = self.TT1
        self.PP3 = self.PP1 * self.σ
        self.ρρ3 = self.PP1 / (self.R_gas3 * self.TT1)
        self.Cp3 = Cp(substance, T=self.TT3, a_ox=self.a_ox3, fuel=fuel)
        self.k3 = self.Cp3 / (self.Cp3 - self.R_gas3)
        self.c3 = self.c1 * nan

        return {'T*1': self.TT1, 'T*3': self.TT3,
                'P*1': self.PP1, 'P*3': self.PP3,
                'ρ*1': self.ρρ1, 'ρ*3': self.ρρ3,
                'g1': self.g1, 'g3': self.g3,
                'Cp1': self.Cp1, 'Cp3': self.Cp3,
                'R_gas1': self.R_gas1, 'R_gas3': self.R_gas3,
                'k1': self.k1, 'k3': self.k3,
                'a_ox1': self.a_ox1, 'a_ox3': self.a_ox3,
                'P1': self.P1, 'P3': self.P3,
                'T1': self.T1, 'T3': self.T3,
                'ρ1': self.ρ1, 'ρ3': self.ρ3,
                'c1': self.c1, 'c3': self.c3,
                'σ': self.σ, 'g_leak': self.g_leak}


if __name__ == '__main__':
    outlet = Outlet()
    outlet.σ = 0.98
    outlet.g_leak = 0.0001
    n.TT1 = 1500
    n.PP1 = 600_000
    n.a_ox = 3
    print(outlet.solve(substance='EXHAUST', fuel='КЕРОСИН'))
