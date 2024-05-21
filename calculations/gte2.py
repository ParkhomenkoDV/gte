import sys
from numpy import nan, inf, linspace, sqrt
from colorama import Fore
from scipy.optimize import fsolve

from thermodynamics import Substance, atmosphere_standard, GDF, R_gas

sys.path.append('D:/Programming/Python')

import decorators


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


"""
Порядок расчета ТД параметров:
R -> T -> P -> ro -> Cp -> k
"""


class Inlet:
    """Входное устройство"""

    def get_inlet_parameters(self, **kwargs) -> dict[str:int | float]:
        """Расчет параметров перед"""

        # Невозмущенные параметры
        self.gas_const_i = self.substance.gas_const[0]
        assert hasattr(self, "T_i"), 'hasattr(self, "T_i")'
        assert hasattr(self, "P_i"), 'hasattr(self, "P_i")'
        self.ρ_i = self.P_i / (self.gas_const_i * self.T_i)

        # полные параметры
        self.TT_i = self.T_i
        self.PP_i = self.P_i
        self.ρρ_i = self.ρ_i
        self.Cp_i = self.substance.Cp(T=self.TT_i, P=self.PP_i)[0]
        self.k_i = self.Cp_i / (self.Cp_i - self.gas_const_i)

        self.Mc_i = 0
        self.c_i = 0
        self.F_i = inf
        return self.__dict__

    def get_outlet_parameters(self, **kwargs) -> dict[str:int | float]:
        """Расчет параметров после"""
        G = kwargs.get('G', nan)
        self.Mc_o = kwargs.get('M', nan)

        self.gas_const_o = self.substance.gas_const[0]
        self.TT_o = self.TT_i
        self.PP_o = self.PP_i * self.σ
        self.ρρ_o = self.PP_o / (self.gas_const_o * self.TT_o)
        self.Cp_o = self.substance.Cp(T=self.TT_o, P=self.PP_o)[0]
        self.k_o = self.Cp_o / (self.Cp_o - self.gas_const_o)

        self.lc = sqrt(((self.k_o + 1) / 2 * self.Mc_o ** 2) / (1 + (self.k_o - 1) / 2 * self.Mc_o ** 2))
        self.T_o = self.TT_o * GDF('T', self.lc, self.k_o)
        self.P_o = self.PP_o * GDF('P', self.lc, self.k_o)
        self.ρ_o = self.P_o / (self.gas_const_o * self.T_o)
        self.c_o = self.Mc_o * sqrt(self.k_o * self.gas_const_o * self.T_o)
        self.F_o = G / (self.ρρ_o * self.c_o) if self.c_o != 0 else inf
        return self.__dict__

    def calculate(self, *args, **kwargs) -> dict[str: int | float]:
        """Расчет параметров"""
        self.substance = kwargs.get('substance', None)

        self.get_inlet_parameters(**kwargs)
        self.get_outlet_parameters(**kwargs)


class Compressor:
    """Компрессор"""

    def get_inlet_parameters(self, **kwargs) -> dict[str:int | float]:
        """Расчет параметров перед"""
        scheme = kwargs.get('scheme', dict())

        self.gas_const_i = self.substance.gas_const[0]
        self.TT_i = scheme[self.place['contour']][self.place['pos'] - 1].TT_o
        self.PP_i = scheme[self.place['contour']][self.place['pos'] - 1].PP_o
        self.ρρ_i = self.PP_i / (self.gas_const_i * self.TT_i)
        self.Cp_i = self.substance.Cp(T=self.TT_i, P=self.PP_i)[0]
        self.k_i = self.Cp_i / (self.Cp_i - self.gas_const_i)
        return self.__dict__

    def get_outlet_parameters(self, **kwargs) -> dict[str:int | float]:
        """Расчет параметров после"""
        self.gas_const_o = self.substance.gas_const[0]
        self.TT_o = self.TT_i * (1 + (self.ππ ** ((self.k_i - 1) / self.k_i) - 1) / self.eff)
        self.PP_o = self.PP_i * self.ππ
        self.ρρ_o = self.PP_o / (self.gas_const_o * self.TT_o)
        self.Cp_o = self.substance.Cp(T=self.TT_o, P=self.PP_o)[0]  # TODO
        self.k_o = self.Cp_o / (self.Cp_o - self.gas_const_o)
        return self.__dict__

    def calculate(self, *args, **kwargs) -> dict[str:int | float]:
        self.substance = kwargs.get('substance', None)

        G = kwargs.get('G', nan)

        self.get_inlet_parameters(**kwargs)
        self.get_outlet_parameters(**kwargs)

        Cp = 1004
        self.L = Cp * (self.TT_i - self.TT_o)  # потребитель
        self.N = self.L * G
        return self.__dict__


class CombustionChamber:
    """Камера сгорания"""

    def get_inlet_parameters(self, **kwargs) -> dict[str:int | float]:
        """Расчет параметров перед"""
        scheme = kwargs.get('scheme', dict())

        self.gas_const_i = self.substance.gas_const[0]
        self.TT_i = scheme[self.place['contour']][self.place['pos'] - 1].TT_o
        self.PP_i = scheme[self.place['contour']][self.place['pos'] - 1].PP_o
        self.ρρ_i = self.PP_i / (self.gas_const_i * self.TT_i)
        self.Cp_i = self.substance.Cp(T=self.TT_i, P=self.PP_i)[0]
        self.k_i = self.Cp_i / (self.Cp_i - self.gas_const_i)
        return self.__dict__

    def get_outlet_parameters(self, **kwargs) -> dict[str:int | float]:
        """Расчет параметров после"""
        self.gas_const_o = self.substance.gas_const[0]

        if hasattr(self, 'TT_o'):
            self.a_ox = 1
        elif hasattr(self, 'G_fuel'):
            self.TT_o = 1800
        elif hasattr(self, 'a_ox'):
            self.TT_o = 1800
        elif hasattr(self, 'g_fuel'):
            self.TT_o = 1800
        else:
            raise Exception('2222222222222')

        self.PP_o = self.PP_i * self.σ
        self.ρρ_o = self.PP_o / (self.gas_const_o * self.TT_o)
        self.Cp_o = self.substance.Cp(T=self.TT_o, P=self.PP_o)[0]  # TODO
        self.k_o = self.Cp_o / (self.Cp_o - self.gas_const_o)
        return self.__dict__

    def calculate(self, *args, **kwargs):
        self.substance = kwargs.get('substance', None)

        self.get_inlet_parameters(**kwargs)
        self.get_outlet_parameters(**kwargs)

        a_ox3 = 1.2
        l0 = 14
        g_fuel = 1 / l0 / a_ox3


class Turbine:
    """Турбина"""

    def get_inlet_parameters(self, **kwargs) -> dict[str:int | float]:
        """Расчет параметров перед"""
        scheme = kwargs.get('scheme', dict())

        self.gas_const_i = self.substance.gas_const[0]
        self.TT_i = scheme[self.place['contour']][self.place['pos'] - 1].TT_o
        self.PP_i = scheme[self.place['contour']][self.place['pos'] - 1].PP_o
        self.ρρ_i = self.PP_i / (self.gas_const_i * self.TT_i)
        self.Cp_i = self.substance.Cp(T=self.TT_i, P=self.PP_i)[0]
        self.k_i = self.Cp_i / (self.Cp_i - self.gas_const_i)
        return self.__dict__

    def get_outlet_parameters(self, **kwargs) -> dict[str:int | float]:
        """Расчет параметров после"""
        self.gas_const_o = self.substance.gas_const[0]
        self.TT_o = self.TT_i * (1 - (1 - self.ππ ** ((1 - self.k_i) / self.k_i)) * self.eff)
        self.PP_o = self.PP_i / self.ππ
        self.ρρ_o = self.PP_o / (self.gas_const_o * self.TT_o)
        self.Cp_o = self.substance.Cp(T=self.TT_o, P=self.PP_o)[0]  # TODO
        self.k_o = self.Cp_o / (self.Cp_o - self.gas_const_o)

        return self.__dict__

    def calculate(self, *args, **kwargs):
        G = kwargs.get('G', nan)
        self.ππ = kwargs.get('pipi_1_3', nan)

        self.substance = kwargs.get('substance', None)

        self.get_inlet_parameters(**kwargs)
        self.get_outlet_parameters(**kwargs)

        self.L = self.Cp_i * (self.TT_i - self.TT_o)
        self.N = self.L * G
        return self.__dict__


class Outlet:
    """Выходное устройство"""

    def get_inlet_parameters(self, **kwargs) -> dict[str:int | float]:
        """Расчет параметров перед"""
        scheme = kwargs.get('scheme', dict())

        self.gas_const_i = self.substance.gas_const[0]
        self.TT_i = scheme[self.place['contour']][self.place['pos'] - 1].TT_o
        self.PP_i = scheme[self.place['contour']][self.place['pos'] - 1].PP_o
        self.ρρ_i = self.PP_i / (self.gas_const_i * self.TT_i)
        self.Cp_i = self.substance.Cp(T=self.TT_i, P=self.PP_i)[0]
        self.k_i = self.Cp_i / (self.Cp_i - self.gas_const_i)
        return self.__dict__

    def get_outlet_parameters(self, **kwargs) -> dict[str:int | float]:
        """Расчет параметров после"""
        self.gas_const_o = self.substance.gas_const[0]
        self.TT_o = self.TT_i

        if hasattr(self, 'ππ'):
            self.PP_o = self.PP_i / self.ππ
        elif hasattr(self, 'PP_o'):
            self.ππ = self.PP_i / self.PP_o
        else:
            raise Exception('111111111')

        self.ρρ_o = self.PP_o / (self.gas_const_o * self.TT_o)
        self.Cp_o = self.substance.Cp(T=self.TT_o, P=self.PP_o)[0]  # TODO
        self.k_o = self.Cp_o / (self.Cp_o - self.gas_const_o)

        self.c_o = self.v_ * (2 * self.Cp_o * self.TT_o * (1 - self.ππ ** ((1 - self.k_o) / self.k_o))) ** 0.5

        return self.__dict__

    def calculate(self, *args, **kwargs) -> dict[str:int | float]:
        G = kwargs.get('G', nan)

        self.substance = kwargs.get('substance', None)
        self.get_inlet_parameters(**kwargs)
        self.get_outlet_parameters(**kwargs)

        self.R = self.c_o * G
        return self.__dict__


class Load:
    def calculate(self, *args, **kwargs) -> dict[str:int | float]:
        return self.__dict__


class GTE_mode:
    def __init__(self):
        self.H = nan  # высота полета [м]
        self.M = nan  # Мах полета []

    def __call__(self):
        return {'H': self.H, 'M': self.M}


class GTE_scheme:
    """Схема ГТД"""

    GTE_NODES = (Inlet, Compressor, CombustionChamber, Turbine, Outlet)

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
            eq.append(sum([node.calculate(**p, scheme=self.scheme, substance=self.substance)['N'] for node in shaft]))

        # требования
        eq.append(self.scheme[1][-1].calculate(**p, scheme=self.scheme, substance=self.substance)['R'] - self.R)

        return eq

    # TODO: обучить модель
    def get_varibles(self) -> dict[str:int | float]:
        """Начальные приближения"""

        # Массовый расход
        vars0 = {'G': 30}

        # Степени понижения полного давления в турбинах
        for contour in self.scheme:
            for i, node in enumerate(self.scheme[contour]):
                if isinstance(node, Turbine):
                    vars0[f'pipi_{contour}_{i}'] = 3

        print(f'points0: {vars0}')
        self.__vars = vars0
        return vars0

    @decorators.timeit(6)
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
                    node.calculate(**vars0, **kwargs)

            vars_list = fsolve(self.equations, tuple(vars0.values()),
                               xtol=epsilon, maxfev=100 * (len(vars0) + 1))
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

    # TODO:
    def solve(self):
        pass


if __name__ == '__main__':
    gte = GTE('Jumo 004b')
    if 1:
        print(Fore.CYAN + f'{gte}' + Fore.RESET)
        gte.scheme = {1: [Inlet(), Compressor(), CombustionChamber(), Turbine(), Outlet()]}
        gte.shafts = {1: [gte.scheme[1][1], gte.scheme[1][3]]}

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

        gte.scheme[1][1].ππ = 6  # list(linspace(3, 43, 40 + 1))
        gte.scheme[1][1].eff = 0.86
        gte.scheme[1][1].g_leak = 0.05

        gte.scheme[1][2].T_fuel = 40 + 273.15
        gte.scheme[1][2].η_burn = 0.99
        gte.scheme[1][2].TT_o = 1000
        gte.scheme[1][2].T_lim = 1000
        gte.scheme[1][2].σ = 0.94
        gte.scheme[1][2].g_leak = 0

        gte.scheme[1][3].eff = 0.92
        gte.scheme[1][3].η_mechanical = 0.99
        gte.scheme[1][3].T_lim = 1000
        gte.scheme[1][3].g_leak = 0.05

        gte.scheme[1][4].PP_o = 101325
        gte.scheme[1][4].eff = 0.96
        gte.scheme[1][4].v_ = 0.98
        gte.scheme[1][4].g_leak = 0.001

        '''
        gte.validate_scheme()
        gte.export_gte_main()
        gte.solve(how=how, error=error, Niter=Niter, file_type='xlsx')
        '''

    gte.calculate(scheme=gte.scheme, **gte.mode(), substance=gte.substance, fuel=gte.fuel)
    for contour in gte.scheme:
        for node in gte.scheme[contour]:
            print(f'{node.__class__.__name__}: {node.__dict__}')
