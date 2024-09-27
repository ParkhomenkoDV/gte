from numpy import nan
from thermodynamics import T0, Cp, R_gas, η_polytropic
from tools import isnum, eps

from combustion_chambler import CombustionChamber  # камера сгорания


def find_node_in_scheme(scheme, node2find) -> tuple:
    """Поиск положения узла в схеме"""
    for contour in scheme:
        for i, node in enumerate(scheme[contour]):
            if scheme[contour][i] is node2find: return contour, i


class HeatExchanger:
    """Теплообменный аппарат"""

    def __init__(self, name='HeatExchanger'):
        self.name = name
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
        self.a_ox1 = nan  # коэффициент избытка окислителя перед
        self.a_ox3 = nan  # коэффициент избытка окислителя после

        self.τ = nan  #
        self.σ = nan  # коэффициент потери полного давления
        self.g_leak = nan  # относительные массовые утечки

        self.warnings = set()  # предупреждения

    def get_variability(self):
        result = 1
        for attr in (self.τ, self.σ, self.g_leak):
            if type(attr) is list: result *= len(attr)
        return result

    def set_combination(self, combination, heatexchanger_main):
        varible_params = ('τ', 'σ', 'g_leak')
        positions = [0] * len(varible_params)

        for i in range(combination):
            for j, param in enumerate(varible_params):
                if positions[j] == len(getattr(heatexchanger_main, varible_params[j])) - 1:
                    positions[j] = 0
                else:
                    positions[j] += 1
                    continue
        for j, param in enumerate(varible_params):
            setattr(self, varible_params[j], getattr(heatexchanger_main, varible_params[j])[positions[j]])

    def input_parameters(self):
        pass

    def solve(self, **kwargs):
        substance = kwargs.get('substance', '')  # рабочее тело
        fuel = kwargs.get('fuel', '')  # горючее

        scheme = kwargs.get('scheme', {})
        if scheme:
            c, n = find_node_in_scheme(scheme, self)
            self.a_ox = scheme[c][n - 1].a_ox3
            self.TT1 = scheme[c][n - 1].TT3
            self.PP1 = scheme[c][n - 1].PP3
            self.g1 = scheme[c][n - 1].g3 - self.g_leak
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

        self.a_ox = kwargs.get('a_ox', nan) if self.a_ox is nan else self.a_ox  # коэффициент избытка окислителя
        self.R_gas1 = R_gas(substance, a_ox=self.a_ox, fuel=fuel)
        self.TT1 = kwargs.get('TT3', nan) if self.TT1 is nan else self.TT1
        self.PP1 = kwargs.get('PP3', nan) if self.PP1 is nan else self.PP1
        self.ρρ1 = self.PP1 / (self.R_gas1 * self.TT1)
        self.Cp1 = Cp(substance, T=self.TT1, P=self.PP1, a_ox=self.a_ox, fuel=fuel)
        self.k1 = self.Cp1 / (self.Cp1 - self.R_gas1)

        self.R_gas3 = R_gas(substance, a_ox=self.a_ox, fuel=fuel)
        self.TT3 = self.TT1 * self.τ
        self.PP3 = self.PP1 * self.σ
        self.ρρ3 = self.PP3 / (self.R_gas3 * self.TT3)
        self.Cp3 = Cp(substance, T=self.TT3, P=self.PP3, a_ox=self.a_ox, fuel=fuel)
        self.k3 = self.Cp3 / (self.Cp3 - self.R_gas3)


if __name__ == '__main__':
    he = HeatExchanger()
    he.τ = 0.67
    he.σ = 0.98
    he.g_leak = 0
    print(he.solve(substance='EXHAUST', fuel='КЕРОСИН', TT3=1550, PP3=600_000, a_ox=3))
