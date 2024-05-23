import sys
from copy import deepcopy
from tqdm import tqdm
from colorama import Fore

import pandas as pd
from numpy import nan, inf, linspace, sqrt, cos, sin, radians, prod
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

from thermodynamics import Substance, atmosphere_standard, GDF, R_gas, l_stoichiometry

sys.path.append('D:/Programming/Python')

import decorators


def It(T0, T1, g):
    return


def Ti(i, T0, g):
    return


def Pt(T0, T1, q):
    return


def Tp(p, T0, g):
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

        mode = kwargs.pop('mode', None)
        assert hasattr(mode, "T") and hasattr(mode, "P")

        self.gas_const_i = self.substance.gas_const[0]
        self.TT_i = mode.T  # TODO
        self.PP_i = mode.P  # TODO
        self.ρρ_i = self.PP_i / (self.gas_const_i * self.TT_i)
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
        fuel = kwargs.pop('fuel', '')

        self.gas_const_o = self.substance.gas_const[0]

        if hasattr(self, 'TT_o'):
            self.a_ox = 1
        elif hasattr(self, 'G_fuel'):
            self.TT_o = 1800
        elif hasattr(self, 'a_ox'):
            g_fuel = 1 / (self.a_ox * l_stoichiometry(fuel))  # приведена ко входу в КС
            self.TT_o = 1800
        elif hasattr(self, 'g_fuel'):
            self.a_ox = 1 / (self.g_fuel * l_stoichiometry(fuel))  # приведена ко входу в ГТД
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


class HeatExchanger:
    """Теплообменный аппарат"""
    pass


class GTE_mode:

    def __setattr__(self, key, value):
        # атмосферные условия
        if key == 'T':  # статическая окружающая температура [К]
            assert type(value) in (int, float)
            assert 0 < value
        elif key == 'P':  # статическое окружающее давление [Па]
            assert type(value) in (int, float)
            assert 0 <= value

        # высотно-скоростные характеристики
        elif key == 'H':  # высота полета [м]
            assert type(value) in (int, float)
        elif key == 'M':  # Мах полета []
            assert type(value) in (int, float)
            assert 0 <= value
        else:
            raise AttributeError('"T", "P", "H", "M"')

        object.__setattr__(self, key, value)


GTE_NODES = (Inlet, Compressor, CombustionChamber, Turbine, Outlet, Load)


class GTE_scheme(dict):
    """Схема ГТД"""

    def __init__(self, scheme: dict):
        assert type(scheme) is dict

        scheme = dict(sorted(scheme.items(), key=lambda item: item[0]))
        contours, contour_nodes = map(tuple, (scheme.keys(), scheme.values()))

        assert all(map(lambda contour: type(contour) is int, contours))
        for nodes in contour_nodes:
            assert all(map(lambda node: type(node) in GTE_NODES, nodes))

        assert 1 in contours
        for i in range(len(contours) - 1):
            assert contours[i + 1] - contours[i] == 1

        super(GTE_scheme, self).__init__(scheme)

    # никаких append/pop/insert! Только перезапись

    @staticmethod
    def Figures(node, **kwargs) -> tuple:
        x0 = kwargs.get('x0', 0)
        y0 = kwargs.get('y0', 0)
        x, y = [], []

        if type(node) is Inlet:
            x = [x0 - 0.4, x0 + 0.4, x0 + 0.4, x0 - 0.4]
            y = [y0 + 0.4, y0 + 0.4, y0 - 0.4, y0 - 0.4]
        elif type(node) == Compressor:
            x = [x0 - 0.4, x0 + 0.4, x0 + 0.4, x0 - 0.4, x0 - 0.4]
            y = [y0 + 0.4, y0 + 0.2, y0 - 0.2, y0 - 0.4, y0 + 0.4]
        elif type(node) == CombustionChamber:
            x = [0.4 * cos(alpha) + x0 for alpha in linspace(0, radians(360), 360)]
            y = [0.4 * sin(alpha) + y0 for alpha in linspace(0, radians(360), 360)]
        elif type(node) == Turbine:
            x = [x0 - 0.4, x0 + 0.4, x0 + 0.4, x0 - 0.4, x0 - 0.4]
            y = [y0 + 0.2, y0 + 0.4, y0 - 0.4, y0 - 0.2, y0 + 0.2]
        elif type(node) == Outlet:
            x = [x0 + 0.4, x0 - 0.4, x0 - 0.4, x0 + 0.4]
            y = [y0 + 0.4, y0 + 0.4, y0 - 0.4, y0 - 0.4]
        elif type(node) == HeatExchanger:
            x = [x0 - 0.4, x0 + 0.4, x0 + 0.4, x0 - 0.4, x0 - 0.4]
            y = [y0 + 0.4, y0 + 0.4, y0 - 0.4, y0 - 0.4, y0 + 0.4]
        elif type(node) == Load:
            x = [x0 - 0.4, x0, x0 + 0.4, x0 - 0.4]
            y = [y0 - 0.4, y0 + 0.4, y0 - 0.4, y0 - 0.4]
        return x, y

    def show(self, **kwargs):
        fg = plt.figure(figsize=kwargs.get('figsize', (max(map(len, self.values())) * 2, (len(self) + 1 + 2) * 2)))
        fg.suptitle('GTE scheme', fontsize=14, fontweight='bold')
        gs = fg.add_gridspec(len(self) + 1, 1)  # строки, столбцы

        for contour in self:
            fg.add_subplot(gs[len(self) - contour, 0])
            plt.grid(True)
            plt.axis('square')
            # plt.title('contour ' + to_roman(contour) + ' | ' + 'контур ' + to_roman(contour), fontsize=14)
            plt.xlim(0, len(self[contour]))
            plt.ylim(0, 1)
            plt.xticks(linspace(0, len(self[contour]), len(self[contour]) + 1))
            plt.yticks(linspace(0, 1, 1 + 1))

            x0 = y0 = 0.5

            for i, node in enumerate(self[contour]):
                plt.plot(*self.Figures(node, x0=x0, y0=y0), color='black', linewidth=3,
                         label=f'{contour}.{i + 1}: {node.__class__.__name__}')
                plt.text(x0, y0, f'{contour}.{i + 1}', fontsize=12, fontweight='bold',
                         ha='center', va='center')
                x0 += 1

        fg.add_subplot(gs[len(self), 0])
        plt.axis('off')
        plt.grid(False)
        plt.xlim(0, max(map(len, self.values())))
        plt.ylim(0, 1)
        plt.plot([0, max(map(len, self.values()))], [0.5, 0.5], color='black', linewidth=1.5, linestyle='dashdot')

        fg.legend(title='Specification', title_fontsize=14, alignment='center',
                  loc='lower center', fontsize=12, ncols=len(self),
                  frameon=True, framealpha=1.0, facecolor='white', edgecolor='black',
                  draggable=True)

        plt.show()


class GTE_shaft(list):
    """Вал ГТД"""

    def __init__(self, shaft: tuple | list):
        assert type(shaft) in (list, tuple)
        assert all(map(lambda node: type(node) in GTE_NODES))
        super(GTE_shaft, self).__init__(shaft)


class GTE:
    """ГТД"""

    @classmethod
    def version(cls) -> str:
        version = 7.0
        next_version = ('камера смешения',
                        'переход к массивам в местах постоянства диапазона значений'
                        'теплак плак-плак',
                        'переходный канал',
                        'type(node) is class -> isinstance(node, class)'
                        'ТД параметры через PTM 1677 as {"O2": 98, "H2": 2}',
                        'соотношение соответствующих относительных расходов к своим контурам',
                        'охлаждение турбины',
                        'get_inlet_parameters() for all nodes',
                        'get_outlet_parameters() for all nodes',
                        'продолжение расчета',
                        'multiprocessing',
                        'ускорение расчета до 6000 [ГТД/с]')
        print(f'{cls.__name__} version: {Fore.GREEN}{version}')
        for i, v in enumerate(next_version): print(cls.__name__ + ' version:', int(version) + i + 1, v)
        return str(version)

    def __init__(self, name='GTE', scheme=None) -> None:

        assert type(name) is str

        self.name = name  # название ГТД
        self.scheme = scheme  # GTE_scheme(scheme) if scheme is not None else GTE_scheme({1: tuple()})  # схема ГТД
        self.shafts = []  # валы ГТД
        self.contouring = {1: 1}  # степень контурности []

        self.mode = GTE_mode()  # режим работы

        self.v = nan  # скорость полета [м/с]

    def __setattr__(self, key, value):
        if key == 'name':
            assert type(value) is str
            object.__setattr__(self, key, value)
        elif key == 'scheme':
            if type(value) is dict:
                object.__setattr__(self, key, GTE_scheme(value))
            elif value is None:
                object.__setattr__(self, key, GTE_scheme({1: tuple()}))
            else:
                raise AttributeError('type(scheme) is dict')
        else:
            object.__setattr__(self, key, value)

    def describe(self) -> None:
        """Описание ГТД"""
        print(f'name: {self.name}')
        print()
        print('scheme:')
        for contour in self.scheme:
            print('\t' + f'contour: {contour}')
            for node in self.scheme[contour]:
                print('\t\t' + f'node: {node}')
                print('\t\t\t' + f'parameters: {dict(sorted(node.__dict__.items(), key=lambda item: item[0]))}')
        print()

    def summary(self):
        pass

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
        eq.append(sum([self.scheme[contour][-1].calculate(**p, scheme=self.scheme, substance=self.substance)['R']
                       for contour in self.scheme]) - self.R)

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
        """Расстановка мест положений в ГТД"""
        for contour in self.scheme:
            for i, node in enumerate(self.scheme[contour]):
                node.place = {'contour': contour, 'pos': i}

    def get_variability(self) -> int:
        """Максимальное количество комбинаций варьируемых параметров"""
        variability = prod([len(value) for key, value in self.mode.items()
                            if type(value) in (tuple, list)])
        return variability
        return prod([len(value) for key, value in self.__dict__.items()
                     if type(value) is list and len(value) and not key.startswith('_')])

    def gte_generator(self):
        """Генератор объектов ГТД с заданными варьируемыми параметрами"""
        max_combinations = [self.get_variability()]  # для ГТД
        for node in self.__path: max_combinations.append(node.get_variability())  # для узлов ГТД
        combinations = [1] * (1 + len(self.__path))

        for comb in tqdm(range(prod(max_combinations)), desc='Calculation', ncols=70):
            gte_var = deepcopy(self)
            gte_var.__set_combination(combinations[0], self)
            for i, node in enumerate(gte_var.__path): node.set_combination(combinations[i + 1], self.__path[i])
            gte_var.__update_combination(self, combinations, max_combinations)
            yield gte_var

    # TODO:
    def solve(self):
        for gte_var in self.gte_generator():
            gte_var.calculate(scheme=self.scheme, mode=self.mode, substance=self.substance, fuel=self.fuel)


if __name__ == '__main__':

    if 1:
        gte = GTE('Jumo 004b')
        print(Fore.CYAN + f'{gte.name}' + Fore.RESET)
        gte.scheme = {1: [Inlet(), Compressor(), CombustionChamber(), Turbine(), Outlet()]}
        gte.scheme.show()
        gte.shafts = {1: [gte.scheme[1][1], gte.scheme[1][3]]}

        gte.m = {1: 1}

        gte.mode.T = 288
        gte.mode.P = 101325
        gte.mode.H = 0
        gte.mode.M = 0

        gte.R = 10_000

        gte.substance = Substance({'N2': 0.755, 'O2': 0.2315, 'Ar': 0.01292, 'Ne': 0.000014, 'H2': 0.000008})
        gte.fuel = Substance({'KEROSENE': 1.0})

        gte.scheme[1][0].σ = 0.98
        gte.scheme[1][0].g_leak = 0.005

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

    if 0:
        gte = GTE('CFM-56')
        print(Fore.CYAN + f'{gte.name}' + Fore.RESET)
        gte.scheme = {1: [Inlet(), Compressor(), CombustionChamber(), Turbine(), Outlet()],
                      2: [Inlet(), Compressor(), Outlet()]}
        gte.scheme.show()
        gte.shafts = {1: [gte.scheme[1][1], gte.scheme[2][1], gte.scheme[1][3]]}

        gte.m = {1: 1}

        gte.mode.T = 288
        gte.mode.P = 101325
        gte.mode.H = 0
        gte.mode.M = 0

        gte.R = 10_000

        gte.substance = Substance({'N2': 0.755, 'O2': 0.2315, 'Ar': 0.01292, 'Ne': 0.000014, 'H2': 0.000008})
        gte.fuel = Substance({'KEROSENE': 1.0})

        gte.scheme[1][0].σ = 0.98
        gte.scheme[1][0].g_leak = 0.005

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

        gte.scheme[2][0].σ = 0.98
        gte.scheme[2][0].g_leak = 0.005

        gte.scheme[2][1].ππ = 6  # list(linspace(3, 43, 40 + 1))
        gte.scheme[2][1].eff = 0.86
        gte.scheme[2][1].g_leak = 0.05

        gte.scheme[2][2].PP_o = 101325
        gte.scheme[2][2].eff = 0.96
        gte.scheme[2][2].v_ = 0.98
        gte.scheme[2][2].g_leak = 0.001

    gte.calculate(scheme=gte.scheme, mode=gte.mode, substance=gte.substance, fuel=gte.fuel)

    gte.describe()
