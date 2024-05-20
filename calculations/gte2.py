import sys
from numpy import nan, inf, linspace, sqrt
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
        self.gas_const_i = self.substance.gas_const[0]
        assert hasattr(self, "T_i"), 'hasattr(self, "T_i")'
        assert hasattr(self, "P_i"), 'hasattr(self, "P_i")'
        self.ρ_i = self.P_i / (self.gas_const_i * self.T_i)

        self.Cp_i = 1004  # Cp(self.substance, T=self.T1, P=self.P1) #TODO
        self.k_i = self.Cp_i / (self.Cp_i - self.gas_const_i)
        self.Mc_i = 0
        self.c_i = 0
        self.F_i = inf

        # полные параметры
        self.TT_i = self.T_i
        self.PP_i = self.P_i
        self.ρρ_i = self.ρ_i

    def get_outlet_parameters(self, **kwargs):
        """Расчет параметров после"""
        G = kwargs.get('G', nan)

        self.gas_const_o = self.substance.gas_const[0]
        self.TT_o = self.TT_i
        self.PP_o = self.PP_i * self.σ
        self.ρρ_o = self.PP_o / (self.gas_const_o * self.TT_o)

        self.Mc_o = kwargs.get('M', nan)
        self.c_o = self.Mc_o * sqrt(self.gas_const_o)  # TODO
        self.F_o = G / (self.ρρ_o * self.c_o)

    def calculate(self, *args, **kwargs) -> dict[str: int | float]:
        """Расчет параметров"""
        self.substance = kwargs.get('substance', None)

        self.get_inlet_parameters(**kwargs)
        self.get_outlet_parameters(**kwargs)


class Compressor:
    """Компрессор"""

    def get_inlet_parameters(self, **kwargs) -> None:
        """Расчет параметров перед"""
        self.gas_const_i = self.substance.gas_const[0]
        self.TT_i = scheme[self.place['contour']][self.place['pos'] - 1].TT_o
        self.PP_i = scheme[self.place['contour']][self.place['pos'] - 1].PP_o

    def get_outlet_parameters(self, **kwargs):
        """Расчет параметров после"""
        self.TT_o = 800

    def calculate(self, *args, **kwargs) -> dict[str:int | float]:
        self.substance = kwargs.get('substance', None)
        scheme = kwargs.get('scheme', dict())
        G = kwargs.get('G', nan)

        self.get_inlet_parameters(**kwargs)
        self.get_outlet_parameters(**kwargs)

        Cp = 1004
        Lc = Cp * (self.TT_o - self.TT_i)
        return {'N': Lc * G, 'G': G * 0.999}


class CombustionChamber:
    """Камера сгорания"""

    def calculate(self, *args, **kwargs):
        a_ox3 = 1.2
        l0 = 14
        g_fuel = 1 / l0 / a_ox3


class Turbine:
    """Турбина"""

    def calculate(self, *args, **kwargs):
        G = kwargs.get('G', nan)
        ππ = kwargs.get('pipi_1_3', nan)

        k = 1.33
        T1 = 1650
        T2 = 1300
        eff = 0.9
        Lt = -k / (k - 1) * 288 * T1 * (1 - 1 / (ππ ** ((k - 1) / k))) * eff
        return {'N': Lt * G, 'G': G}


class Nozzle:

    def calculate(self, *args, **kwargs):
        G = kwargs.get('G', nan)

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


class GTE_mode:
    def __init__(self):
        self.H = nan  # высота полета [м]
        self.M = nan  # Мах полета []

    def __call__(self):
        return {'H': self.H, 'M': self.M}


class GTE_scheme:
    GTE_NODES = (Inlet, Compressor, CombustionChamber, Turbine, Nozzle)

    def __init__(self, scheme: dict):
        assert type(scheme) is dict, 'type(scheme) is dict'
        scheme = dict(sorted(scheme.items(), key=lambda item: item[1]))
        keys, values = map(list, (scheme.keys(), scheme.values()))
        assert all(map(lambda key: type(key) is int, keys)), 'all(map(lambda key: type(key) is int, keys))'
        """assert all(map(lambda value: type(value) in self.GTE_NODES, values)), \
            'all(map(lambda value: type(value) in self.GTE_NODES, values))'"""
        assert 1 in keys, '1 in keys'
        assert all(map(lambda key: abs(keys[i] - keys[i + 1]) == 1, keys[:-1]))

        self.__scheme = scheme


class GTE:
    """ГТД"""

    def __init__(self, name='GTE') -> None:
        self.name = name
        self.scheme = GTE_scheme({1: []})  # схема ГТД
        self.shafts = []  # валы ГТД
        self.contouring = {1: 1}  # степень контурности []

        self.mode = GTE_mode()  # режим работы

        self.v = nan  # скорость полета [м/с]

    def __str__(self) -> str:
        return self.name

    def equations(self, points: tuple | list, *args, **kwargs) -> list:
        """СНЛАУ"""

        p = {key: points[i] for i, key in enumerate(self.__vars.keys())}  # преобразование списка параметров в словарь

        eq = list()  # список НЛАУ

        # уравнения неразрывности
        '''for i in range(1):
            res.append(Turbine(points[0], points[1])['G'] - Compressor.calculate(points[0])['G'])'''

        # Баланс мощностей
        for contour, shaft in self.shafts.items():
            eq.append(sum([node.calculate(**p, scheme=self.scheme, substance=self.substance)['N']
                           for i, node in enumerate(shaft)]))

        # требования
        eq.append(self.scheme[1][-1].calculate(**p, **kwargs)['R'] - self.R)

        return eq

    # TODO: обучить модель
    def get_varibles(self) -> dict[str:int | float]:
        """Начальные приближения"""

        # Массовый расход
        vars0 = {'G': 300}

        # Степени понижения полного давления в турбинах
        for contour in self.scheme:
            for i, node in enumerate(self.scheme[contour]):
                if isinstance(node, Turbine):
                    vars0[f'pipi_{contour}_{i}'] = 3

        print(f'points0: {vars0}')
        self.__vars = vars0
        return vars0

    @timeit(6)
    def calculate(self, Niter: int = 10, epsilon=0.01, *args, **kwargs):
        """Решение СНЛАУ"""

        assert type(Niter) is int, 'type(Niter) is int'
        assert 1 <= Niter, '1 <= Niter'

        self.placement()

        vars0 = self.get_varibles()
        for i in range(Niter):

            # расчет ГТД в строчку
            for contour in self.scheme:
                for node in self.scheme[contour]:
                    node.calculate(vars0, **kwargs)

            vars_list = fsolve(self.equations, tuple(vars0.values()))
            vars = {key: vars_list[i] for i, key in enumerate(vars0.keys())}

            print(Fore.GREEN + f'points: {vars}' + Fore.RESET)
            print(Fore.CYAN + f'eq: {gte.equations(list(vars.values()))}' + Fore.RESET)
            if all(map(lambda x0, x: abs(x - x0) / x0 <= epsilon, vars0.values(), vars.values())): break
            vars0 = vars  # обновление параметров
        else:
            print(Fore.RED + 'Решение не найдено!' + Fore.RESET)
        return vars

    def placement(self):
        for contour in self.scheme:
            for i, node in enumerate(self.scheme[contour]):
                node.place = {'contour': contour, 'pos': i}


if __name__ == '__main__':
    gte = GTE('Jumo 004b')
    if 1:
        print(Fore.CYAN + f'{gte}' + Fore.RESET)
        gte.scheme = {1: [Inlet(), Compressor(), CombustionChamber(), Turbine(), Nozzle()]}
        gte.shafts = {1: [gte.scheme[1][1], gte.scheme[1][3], Load()]}

        gte.m = {1: 1}

        gte.mode.H = 0
        gte.mode.M = 0

        gte.R = 10_000

        gte.substance = Substance({'N2': 0.755, 'O2': 0.2315, 'Ar': 0.01292, 'Ne': 0.000014, 'H': 0.000008})
        gte.fuel = Substance({'KEROSENE': 1.0})

        gte.scheme[1][0].σ = 0.98
        gte.scheme[1][0].g_leak = 0.005
        gte.scheme[1][0].T_i = 288
        gte.scheme[1][0].P_i = 100000

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
        gte.solve(how=how, error=error, Niter=Niter, file_type='xlsx')
        '''
    scheme = gte.scheme
    gte.calculate(scheme=gte.scheme, **gte.mode(), substance=gte.substance, fuel=gte.fuel)
