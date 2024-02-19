import sys
from numpy import nan, inf, linspace
from colorama import Fore
from scipy.optimize import fsolve

from thermodynamics import Substance, atmosphere_standard, R_gas

sys.path.append('D:/Programming/Python')

from decorators import timeit


def Cp(T, g):
    return sum()


def entalphy_t(T0, T1, g):
    return


def temperature_i(i, T0, g):
    return


def pressure_t(T0, T1, q):
    return


def temperature_p(p, T0, g):
    return


def find_node_in_scheme(scheme, node2find) -> tuple:
    """Поиск положения узла в схеме"""
    for contour in scheme:
        for i, node in enumerate(scheme[contour]):
            if scheme[contour][i] is node2find: return contour, i


def get_combination(combination: int, max_list: list[int]) -> list[int]:
    n = max_list.copy()
    for i, max_element in enumerate(max_list):
        n[i] = combination % max_element
        combination //= max_element
    return n


class Inlet:
    """Входное устройство"""

    def get_inlet_parameters(self, **kwargs) -> None:
        """Расчет параметров перед"""
        # Невозмущенные параметры
        self.R_gas_i = self.substance.R()
        if not hasattr(self, 'T1') and self.substance.strip().upper() == 'AIR':
            self.T1 = atmosphere_standard(kwargs.get('H', None))['T']
        if not hasattr(self, 'P1') and self.substance.strip().upper() == 'AIR':
            self.P1 = atmosphere_standard(kwargs.get('H', None))['P']
        assert self.T1 and self.P1, f'{type(self).__name__} object has no attributes T1 and P1!'
        self.ρ1 = self.P1 / (self.R_gas1 * self.T1)
        self.g1 = 1
        self.Cp1 = Cp(self.substance, T=self.T1, P=self.P1, a_ox=getattr(self, 'a_ox1', None), fuel=fuel)
        self.k1 = self.Cp1 / (self.Cp1 - self.R_gas1)
        self.Mc1 = 0
        self.c1 = 0
        self.TT1 = self.T1
        self.PP1 = self.P1
        self.ρρ1 = self.ρ1

    def calculate(self, G, *args, **kwargs):
        """Расчет параметров перед"""

        T1 = 288
        self.T_o = T1
        P1 = 100_000
        P3 = P1
        return {'T_o': self.T_o, 'P3': P3}


class Compressor:

    def calculate(self, G, *args, **kwargs):
        Cp = 1004
        T1 = scheme[self.place['contour']][self.place['node']-1].T_o
        #T1 = 288
        T3 = 800
        Lc = Cp * (T3 - T1)
        return {'N': Lc * G, 'G': G * 0.999}


class CombustionChamber:

    def calculate(self, G, *args, **kwargs):
        a_ox3 = 1.2
        l0 = 14
        g_fuel = 1 / l0 / a_ox3


class Turbine:
    def calculate(self, G, ππ, *args, **kwargs):
        k = 1.33
        T1 = 1650
        T2 = 1300
        eff = 0.9
        Lt = -k / (k - 1) * 288 * T1 * (1 - 1 / (ππ ** ((k - 1) / k))) * eff
        return {'N': Lt * G, 'G': G}


class Nozzle:
    def calculate(self, G, *args, **kwargs):
        Cp2 = 1200
        TT3 = 600
        ππ = 1.4
        k2 = 1.33
        v_ = 0.99
        c3 = v_ * (2 * Cp2 * TT3 * (1 - ππ ** ((1 - k2) / k2))) ** 0.5
        self.R = c3 * G
        return {'R': self.R}


class Load:
    def calculate(self, *args, **kwargs):
        return {'N': 10_000}


class GTE:
    """ГТД"""

    def __init__(self, name='GTE') -> None:
        self.name = name
        self.scheme = {}  # схема ГТД
        self.shafts = []  # валы ГТД
        self.contouring = {1: 1}  # степень контурности []

        self.M = nan  # Мах полета []
        self.v = nan  # скорость полета [м/с]

    def __str__(self) -> str:
        return self.name

    def eq(self, points, *args, **kwargs):
        res = []

        # уравнения неразрывности
        '''for i in range(1):
            res.append(Turbine(points[0], points[1])['G'] - Compressor.calculate(points[0])['G'])'''

        # Баланс мощностей
        for contour, shaft in self.shafts.items():
            res.append(sum([node.calculate(*points, *args, **kwargs)['N'] for i, node in enumerate(shaft)]))

        res.append(self.scheme[1][-1].calculate(points[0])['R'] - self.R)

        return res

    def get_varibles(self):
        # Массовй расход
        vars = {'G': 300}

        # Степени понижения полного давления в турбинах
        for contour in self.scheme:
            for i, node in enumerate(self.scheme[contour]):
                if isinstance(node, Turbine):
                    vars['pipi' + str(i)] = 3
        print(f'points0: {vars}')
        return vars

    @timeit()
    def calculate(self, *args, **kwargs):
        vars0 = self.get_varibles()
        vars_list = fsolve(self.eq, list(vars0.values()))
        vars = dict()
        for i, key in enumerate(vars0.keys()):
            vars[key] = vars_list[i]
        print(Fore.GREEN + f'points: {vars}')
        print(Fore.CYAN + f'eq: {gte.eq(list(vars.values()))}')
        return vars


if __name__ == '__main__':
    gte = GTE('Jumo 004b')
    if 1:
        print(Fore.CYAN + f'{gte}' + Fore.RESET)
        gte.scheme = {1: [Inlet(), Compressor(), CombustionChamber(), Turbine(), Nozzle()]}
        gte.shafts = {1: [gte.scheme[1][1], gte.scheme[1][3], Load()]}

        gte.m = {1: 1}

        gte.H = 0
        gte.M = 0
        gte.R = 10_000

        gte.substance = Substance({'AIR': 1.0})
        gte.fuel = Substance({'KEROSENE': 1.0})

        gte.scheme[1][0].place = {'contour': 1, 'node': 0}
        gte.scheme[1][0].σ = 0.98
        gte.scheme[1][0].g_leak = 0.005

        gte.scheme[1][1].place = {'contour': 1, 'node': 1}
        gte.scheme[1][1].ππ = list(linspace(3, 43, 40 + 1))
        gte.scheme[1][1].ηη = [0.86]
        gte.scheme[1][1].g_leak = [0.05]

        '''
        gte.scheme[1][2].T_fuel = [40 + 273.15]
        gte.scheme[1][2].η_burn = [0.99]
        gte.scheme[1][2].TT3 = list(linspace(800, 1200, 8 + 1))
        gte.scheme[1][2].T_lim = [1000]
        gte.scheme[1][2].σ = [0.94]
        gte.scheme[1][2].g_leak = [0]

        gte.scheme[1][3]._shafts = [{'-': [gte.scheme[1][1]]}]
        gte.scheme[1][3].ηη = [0.92]
        gte.scheme[1][3].η_mechanical = [0.99]
        gte.scheme[1][3].T_lim = [1000]
        gte.scheme[1][3].g_leak = [0.05]

        gte.scheme[1][4].PP3 = [101325]
        gte.scheme[1][4].ηη = [0.96]
        gte.scheme[1][4].v_ = [0.98]
        gte.scheme[1][4].g_leak = [0.001]

        gte.validate_scheme()
        gte.export_gte_main()
        gte.solve(how=how, error=error, Niter=Niter, file_type='xlsx')'''
    scheme = gte.scheme
    gte.calculate(H=gte.H, M=gte.M, substance=gte.substance, fuel=gte.fuel, scheme=gte.scheme)
