from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from colorama import Fore
from numpy import cos, inf, linspace, nan, prod, radians, sin, sqrt
from scipy.optimize import fsolve
from thermodynamics import (
    Substance,
    atmosphere_standard,
    gas_const,
    gdf,
    stoichiometry,
)
from tqdm import tqdm

from src.config import parameters as gtep


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
            if scheme[contour][i] is node2find:
                return contour, i


# TODO: обучить модели регрессии по предсказанию НУ расчета
# TODO: gte.describe() # неизвестные параметры и необходимые уравнения
# TODO: __iter__ вместо gte.generator


class Variability:
    """Варьируемость"""

    @staticmethod
    def varible_parameters(obj) -> dict[str : tuple | list]:
        """Словарь с варьируемыми параметрами и их итераторами значений"""
        return {key: value for key, value in obj.__dict__.items() if type(value) in (tuple, list) and len(value) and not key.startswith("_")}

    def variability(self) -> int:
        """Максимальное количество комбинаций варьируемых параметров"""
        return int(prod([len(value) for value in self.varible_parameters(self).values()]))

    @staticmethod
    def get_combination(combination: int, max_list: list[int]) -> list[int]:
        n = max_list.copy()
        for i, max_element in enumerate(max_list):
            n[i] = combination % max_element
            combination //= max_element
        return n

    def _set_combination(self, combination: int, main_obj: object) -> None:
        """Установка комбинации"""
        varible_params = list(self.varible_parameters(main_obj).keys())
        positions = [0] * len(varible_params)

        for _ in range(combination):
            for j, param in enumerate(varible_params):
                if positions[j] == len(getattr(main_obj, varible_params[j])) - 1:
                    positions[j] = 0
                else:
                    positions[j] += 1
                    continue
        for j, param in enumerate(varible_params):
            setattr(self, param, getattr(main_obj, param)[positions[j]])

    def _update_combination(self, main_obj: object, combination: int, max_combination: int) -> None:
        if max_combination > 1:  # если есть варьируемые параметры
            if combination < max_combination:  # если не конец варьирования параметров
                self._set_combination(combination, main_obj)  # установка текущего параметра варьирования
            else:
                self._set_combination(0, main_obj)


class GTE_node:
    def get_inlet_parameters(self, **kwargs) -> dict[str : int | float]:
        """Расчет параметров перед"""
        scheme = kwargs.pop("scheme", dict())

        self.gas_const_i = self.substance.gas_const[0]
        if not hasattr(self, "TT_i"):
            self.TT_i = scheme[self._place["contour"]][self._place["pos"] - 1].TT_o
        if not hasattr(self, "PP_i"):
            self.PP_i = scheme[self._place["contour"]][self._place["pos"] - 1].PP_o
        self.ρρ_i = self.PP_i / (self.gas_const_i * self.TT_i)
        self.CpCp_i = self.substance.Cp(T=self.TT_i, P=self.PP_i)[0]
        self.kk_i = self.CpCp_i / (self.CpCp_i - self.gas_const_i)
        return self.__dict__


class Inlet(Variability):
    """Входное устройство"""

    def get_inlet_parameters(self, **kwargs) -> dict[str : int | float]:
        """Расчет параметров перед"""
        G = kwargs.get("G", nan)
        mode = kwargs.pop("mode", None)
        assert hasattr(mode, "T") and hasattr(mode, "P") and hasattr(mode, "M")

        self.gas_const_i = self.substance.gas_const[0]

        self.T_i = mode.T
        self.P_i = mode.P
        self.ρ_i = self.P_i / (self.gas_const_i * self.T_i)
        self.Cp_i = self.substance.Cp(T=self.T_i, P=self.P_i)[0]
        self.k_i = self.Cp_i / (self.Cp_i - self.gas_const_i)
        self.M_c_i = mode.M
        self.a = sqrt(self.k_i * self.gas_const_i * self.T_i)
        self.c_i = self.M_c_i * self.a
        self.F_i = G / (self.ρ_i * self.c_i) if self.c_i != 0 else inf

        self.TT_i = self.T_i * (1 - (self.k_i - 1) / 2 * self.M_c_i**2)
        self.PP_i = self.P_i * (1 - (self.k_i - 1) / 2 * self.M_c_i**2) ** (self.k_i / (self.k_i - 1))
        self.ρρ_i = self.PP_i / (self.gas_const_i * self.TT_i)
        self.CpCp_i = self.substance.Cp(T=self.TT_i, P=self.PP_i)[0]
        self.kk_i = self.CpCp_i / (self.CpCp_i - self.gas_const_i)

        return self.__dict__

    def get_outlet_parameters(self, **kwargs) -> dict[str : int | float]:
        """Расчет параметров после"""
        self.gas_const_o = self.substance.gas_const[0]
        self.TT_o = self.TT_i
        self.PP_o = self.PP_i * self.σ
        self.ρρ_o = self.PP_o / (self.gas_const_o * self.TT_o)
        self.CpCp_o = self.substance.Cp(T=self.TT_o, P=self.PP_o)[0]
        self.kk_o = self.CpCp_o / (self.CpCp_o - self.gas_const_o)
        return self.__dict__

    def calculate(self, *args, **kwargs) -> dict[str : int | float]:
        """Расчет параметров"""
        self.substance = kwargs.get("substance", None)

        self.get_inlet_parameters(**kwargs)
        self.get_outlet_parameters(**kwargs)


class Compressor(Variability, GTE_node):
    """Компрессор"""

    def get_outlet_parameters(self, **kwargs) -> dict[str : int | float]:
        """Расчет параметров после"""
        self.gas_const_o = self.substance.gas_const[0]
        self.TT_o = self.TT_i * (1 + (self.ππ ** ((self.kk_i - 1) / self.kk_i) - 1) / self.effeff)
        self.PP_o = self.PP_i * self.ππ
        self.ρρ_o = self.PP_o / (self.gas_const_o * self.TT_o)
        self.CpCp_o = self.substance.Cp(T=self.TT_o, P=self.PP_o)[0]
        self.kk_o = self.CpCp_o / (self.CpCp_o - self.gas_const_o)
        self.G_o = self.G_i * (1 - self.g_leak)
        return self.__dict__

    def calculate(self, *args, **kwargs) -> dict[str : int | float]:
        self.substance = kwargs.pop("substance", None)
        self.G_i = kwargs.pop("G", nan)

        self.get_inlet_parameters(**kwargs)
        self.get_outlet_parameters(**kwargs)

        self.L = self.substance.Cp(T=[self.TT_i, self.TT_o])[0] * (self.TT_i - self.TT_o)
        self.N = self.L * self.G_i
        return self.__dict__


class CombustionChamber(Variability, GTE_node):
    """Камера сгорания"""

    def get_outlet_parameters(self, **kwargs) -> dict[str : int | float]:
        """Расчет параметров после"""
        fuel = kwargs.pop("fuel", "")

        self.gas_const_o = self.substance.gas_const[0]

        if hasattr(self, "TT_o"):
            self.a_ox = 1
        elif hasattr(self, "G_fuel"):
            self.TT_o = 1800
        elif hasattr(self, "a_ox"):
            g_fuel = 1 / (self.a_ox * l_stoichiometry(fuel))  # приведена ко входу в КС
            self.TT_o = 1800
        elif hasattr(self, "g_fuel"):
            self.a_ox = 1 / (self.g_fuel * l_stoichiometry(fuel))  # приведена ко входу в ГТД
            self.TT_o = 1800
        else:
            raise Exception("2222222222222")

        self.PP_o = self.PP_i * self.σ
        self.ρρ_o = self.PP_o / (self.gas_const_o * self.TT_o)
        self.CpCp_o = self.substance.Cp(T=self.TT_o, P=self.PP_o)[0]  # TODO
        self.kk_o = self.CpCp_o / (self.CpCp_o - self.gas_const_o)
        return self.__dict__

    def calculate(self, *args, **kwargs):
        self.substance = kwargs.get("substance", None)

        self.get_inlet_parameters(**kwargs)
        self.get_outlet_parameters(**kwargs)

        a_ox3 = 1.2
        l0 = 14
        g_fuel = 1 / l0 / a_ox3


class Turbine(Variability, GTE_node):
    """Турбина"""

    def get_outlet_parameters(self, **kwargs) -> dict[str : int | float]:
        """Расчет параметров после"""
        self.gas_const_o = self.substance.gas_const[0]
        self.TT_o = self.TT_i * (1 - (1 - self.ππ ** ((1 - self.kk_i) / self.kk_i)) * self.effeff)
        self.PP_o = self.PP_i / self.ππ
        self.ρρ_o = self.PP_o / (self.gas_const_o * self.TT_o)
        self.CpCp_o = self.substance.Cp(T=self.TT_o, P=self.PP_o)[0]
        self.kk_o = self.CpCp_o / (self.CpCp_o - self.gas_const_o)
        self.G_o = self.G_i * (1 - self.g_leak)
        return self.__dict__

    def calculate(self, *args, **kwargs):
        self.G_i = kwargs.pop("G", nan)
        self.ππ = kwargs.pop("pipi_1_3", nan)

        self.substance = kwargs.get("substance", None)

        self.get_inlet_parameters(**kwargs)
        self.get_outlet_parameters(**kwargs)

        self.L = self.substance.Cp(T=[self.TT_i, self.TT_o])[0] * (self.TT_i - self.TT_o)
        self.N = self.L * self.G_i
        return self.__dict__


class Outlet(Variability, GTE_node):
    """Выходное устройство"""

    def get_outlet_parameters(self, **kwargs) -> dict[str : int | float]:
        """Расчет параметров после"""
        self.gas_const_o = self.substance.gas_const[0]
        self.TT_o = self.TT_i

        if hasattr(self, "ππ"):
            self.PP_o = self.PP_i / self.ππ
        elif hasattr(self, "PP_o"):
            self.ππ = self.PP_i / self.PP_o
        else:
            raise Exception("111111111")

        self.ρρ_o = self.PP_o / (self.gas_const_o * self.TT_o)
        self.Cp_o = self.substance.Cp(T=self.TT_o, P=self.PP_o)[0]
        self.k_o = self.Cp_o / (self.Cp_o - self.gas_const_o)
        self.G_o = self.G_i * (1 - self.g_leak)

        self.c_o = self.v_ * (2 * self.Cp_o * self.TT_o * (1 - self.ππ ** ((1 - self.k_o) / self.k_o))) ** 0.5

        return self.__dict__

    def calculate(self, *args, **kwargs) -> dict[str : int | float]:
        self.G_i = kwargs.get("G", nan)

        self.substance = kwargs.get("substance", None)
        self.get_inlet_parameters(**kwargs)
        self.get_outlet_parameters(**kwargs)

        self.R = self.c_o * self.G_o
        return self.__dict__


class Load(Variability):
    def calculate(self, *args, **kwargs) -> dict[str : int | float]:
        return self.__dict__


class GTE_mode(Variability):
    def __setattr__(self, key, value):
        """# атмосферные условия
        if key == 'T':  # статическая окружающая температура [К]
            assert type(value) in (int, float, tuple, list)
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

        elif key == 'R':
            assert type(value) in (int, float)

        else:
            raise AttributeError('"T", "P", "H", "M"')"""

        object.__setattr__(self, key, value)


GTE_NODES = (Inlet, Compressor, CombustionChamber, Turbine, Outlet, Load)


class GTE_scheme(dict):
    """Схема ГТД"""

    def __init__(self, scheme: dict):
        assert type(scheme) is dict

        scheme = dict(sorted(scheme.items(), key=lambda item: item[0]))  # сортировка по контурам по возрастанию
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
        x0 = kwargs.get("x0", 0)
        y0 = kwargs.get("y0", 0)
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
        fg = plt.figure(figsize=kwargs.get("figsize", (max(map(len, self.values())) * 2, (len(self) + 1 + 2) * 2)))
        fg.suptitle("GTE scheme", fontsize=14, fontweight="bold")
        gs = fg.add_gridspec(len(self) + 1, 1)  # строки, столбцы

        for contour in self:
            fg.add_subplot(gs[len(self) - contour, 0])
            plt.grid(True)
            plt.axis("square")
            # plt.title('contour ' + to_roman(contour) + ' | ' + 'контур ' + to_roman(contour), fontsize=14)
            plt.xlim(0, len(self[contour]))
            plt.ylim(0, 1)
            plt.xticks(linspace(0, len(self[contour]), len(self[contour]) + 1))
            plt.yticks(linspace(0, 1, 1 + 1))

            x0 = y0 = 0.5

            for i, node in enumerate(self[contour]):
                plt.plot(
                    *self.Figures(node, x0=x0, y0=y0),
                    color="black",
                    linewidth=3,
                    label=f"{contour}.{i + 1}: {node.__class__.__name__}",
                )
                plt.text(
                    x0,
                    y0,
                    f"{contour}.{i + 1}",
                    fontsize=12,
                    fontweight="bold",
                    ha="center",
                    va="center",
                )
                x0 += 1

        fg.add_subplot(gs[len(self), 0])
        plt.axis("off")
        plt.grid(False)
        plt.xlim(0, max(map(len, self.values())))
        plt.ylim(0, 1)
        plt.plot(
            [0, max(map(len, self.values()))],
            [0.5, 0.5],
            color="black",
            linewidth=1.5,
            linestyle="dashdot",
        )

        fg.legend(
            title="Specification",
            title_fontsize=14,
            alignment="center",
            loc="lower center",
            fontsize=12,
            ncols=len(self),
            frameon=True,
            framealpha=1.0,
            facecolor="white",
            edgecolor="black",
            draggable=True,
        )

        plt.show()


class GTE_shaft(list):
    """Вал ГТД"""

    def __init__(self, shaft: tuple | list):
        assert type(shaft) in (list, tuple)
        assert all(map(lambda node: type(node) in GTE_NODES))
        super(GTE_shaft, self).__init__(shaft)


class GTE(Variability):
    """ГТД"""

    @classmethod
    @property
    def version(cls) -> str:
        version = 8.0
        next_version = (
            "камера смешения",
            "переход к массивам в местах постоянства диапазона значенийтеплак плак-плак",
            "переходный канал",
            "type(node) is class -> isinstance(node, class)соотношение соответствующих относительных расходов к своим контурам",
            "охлаждение турбины",
            "продолжение расчета",
            "multiprocessing",
            "ускорение расчета до 6000 [ГТД/с]",
        )
        print(f"{cls.__name__} version: {Fore.GREEN}{version}")
        for i, v in enumerate(next_version):
            print(cls.__name__ + " version:", int(version) + i + 1, v)
        return str(version)

    @classmethod
    @property
    def author(cls) -> str:
        return "Daniil Vitalievich Andryushin"

    def __init__(self, name="GTE", scheme=None) -> None:
        assert type(name) is str

        self.name = name  # название ГТД
        self.scheme = scheme  # GTE_scheme(scheme) if scheme is not None else GTE_scheme({1: tuple()})  # схема ГТД
        self.shafts = []  # валы ГТД
        self.contouring = {1: 1}  # степень контурности []

        self.mode = GTE_mode()  # режим работы

    def __setattr__(self, key, value):
        if key == "name":
            assert type(value) is str
            object.__setattr__(self, key, value)
        elif key == "scheme":
            if type(value) is dict:
                object.__setattr__(self, key, GTE_scheme(value))
            elif value is None:
                object.__setattr__(self, key, GTE_scheme({1: tuple()}))
            else:
                raise AttributeError("type(scheme) is dict")
        else:
            object.__setattr__(self, key, value)

    def describe(self) -> None:
        """Выявление неизвестных данных и необходимых уравнений"""
        pass

    def summary(self) -> None:
        """Описание ГТД"""
        print(f"name: {self.name}")
        print()
        print("scheme:")
        for contour in self.scheme:
            print(2 * " " + f"contour: {contour}")
            for node in self.scheme[contour]:
                print(4 * " " + f"node: {node.__class__.__name__}")
                for key, value in dict(sorted(node.__dict__.items(), key=lambda item: item[0])).items():
                    if not key.startswith("_"):
                        print(6 * " " + f"{key}: {value}")
        print()
        print(f"substance: {self.substance}")
        print(f"fuel: {self.fuel}")

    def equations(self, points: tuple | list, *args, **kwargs) -> list:
        """СНЛАУ"""

        p = {key: points[i] for i, key in enumerate(self.__vars.keys())}  # преобразование списка параметров в словарь

        eq = list()  # список НЛАУ

        # уравнения неразрывности
        """for i in range(1):
            res.append(Turbine(points[0], points[1])['G'] - Compressor.calculate(points[0])['G'])"""

        # Баланс мощностей
        for contour, shaft in self.shafts.items():
            eq.append(sum([node.calculate(**p, scheme=self.scheme, substance=self.substance)["N"] for node in shaft]))

        # требования
        eq.append(sum([self.scheme[contour][-1].calculate(**p, scheme=self.scheme, substance=self.substance)["R"] for contour in self.scheme]) - self.mode.R)

        return eq

    # TODO: обучить модель
    def get_varibles(self, log=True) -> dict[str : int | float]:
        """Начальные приближения"""

        # Массовый расход
        vars0 = {"G": 30}

        # Степени понижения полного давления в турбинах
        for contour in self.scheme:
            for i, node in enumerate(self.scheme[contour]):
                if isinstance(node, Turbine):
                    vars0[f"pipi_{contour}_{i}"] = 3

        if log:
            print(f"points0: {vars0}")
        self.__vars = vars0
        return vars0

    # @decorators.warns("ignore")  # при отсутствии решения
    def __calculate(self, Niter: int = 10, xtol: float = 0.01, **kwargs):
        """Решение СНЛАУ"""

        log = kwargs.pop("log", False)

        vars0 = self.get_varibles(log=log)
        for i in range(Niter):
            # расчет ГТД в строчку
            for contour in self.scheme:
                for node in self.scheme[contour]:
                    node.calculate(**vars0, **kwargs)

            vars_list = fsolve(
                self.equations,
                tuple(vars0.values()),
                xtol=xtol,
                maxfev=100 * (len(vars0) + 1),
            )
            vars = {key: vars_list[i] for i, key in enumerate(vars0.keys())}

            if log:
                print(Fore.GREEN + f"variables: {vars}" + Fore.RESET)
                print(Fore.CYAN + f"zeros: {self.equations(list(vars.values()))}" + Fore.RESET)

            if all(
                map(
                    lambda x0, x: abs(x - x0) / x0 <= xtol,
                    vars0.values(),
                    vars.values(),
                )
            ):
                break
            vars0 = vars  # обновление параметров
        else:
            print(Fore.RED + "Решение не найдено!" + Fore.RESET)
        return vars

    def placement(self):
        """Расстановка мест положений в ГТД"""
        for contour in self.scheme:
            for i, node in enumerate(self.scheme[contour]):
                node._place = {"contour": contour, "pos": i}

    def gte_generator(self):  # TODO __iter__
        """Генератор объектов ГТД с заданными варьируемыми параметрами"""
        list_count_combinations = list()  # список из количеств комбинаций
        list_count_combinations.append(self.mode.variability())  # для режима работы ГТД
        for contour in self.scheme:
            for node in self.scheme[contour]:
                list_count_combinations.append(node.variability())  # для узлов ГТД

        for comb in tqdm(range(prod(list_count_combinations)), desc="Calculation", ncols=70):
            gte_var = deepcopy(self)  # обнуление параметров
            combinations = self.get_combination(comb, list_count_combinations)

            # установка параметров
            gte_var.mode._set_combination(combinations[0], self.mode)  # для режима работы ГТД
            k = 1
            for contour in gte_var.scheme:
                for node, main_node in zip(gte_var.scheme[contour], self.scheme[contour]):
                    node._set_combination(combinations[k], main_node)  # для узлов ГТД
                    k += 1

            yield gte_var

    def solve(self, Niter: int = 10, xtol: float = 0.01, log=False) -> list[object]:
        """Расчет ГТД"""

        assert type(Niter) is int
        assert 1 <= Niter

        assert type(xtol) is float
        assert 0 < xtol < 1

        self.placement()  # расстановка мест положений в ГТД

        result = list()  # TODO: multiprocessing
        for gte_var in self.gte_generator():
            gte_var.__calculate(
                scheme=gte_var.scheme,
                mode=gte_var.mode,
                substance=gte_var.substance,
                fuel=gte_var.fuel,
            )
            result.append(deepcopy(gte_var))
            if log:
                gte_var.describe()

        return result

    def dataframe(self) -> pd.DataFrame:
        result = dict()
        for key, value in self.__dict__.items():
            if type(value) in (int, float, str):
                result[key] = value
            elif type(value) is object:
                for k, v in value.__dict__.items():
                    if type(v) in (int, float, np.float64, str):
                        result[k] = v
            else:
                pass

        return pd.DataFrame([result])


if __name__ == "__main__":
    from src.config import parameters as gtep

    print(gtep)

    if 0:
        gte = GTE("Jumo 004b")
        gte.scheme = {1: [Inlet(), Compressor(), CombustionChamber(), Turbine(), Outlet()]}
        gte.shafts = {1: [gte.scheme[1][1], gte.scheme[1][3]]}

        gte.contouring = {1: 1}

        gte.mode.T = [288]
        gte.mode.P = 101325
        gte.mode.H = 0
        gte.mode.M = 0

        gte.mode.R = 10_000

        gte.substance = Substance({"N2": 0.755, "O2": 0.2315, "Ar": 0.01292, "Ne": 0.000014, "H2": 0.000008})
        gte.fuel = Substance({"C2H8N2": 1.0})

        gte.scheme[1][0].σ = [0.98]
        gte.scheme[1][0].g_leak = 0.005

        gte.scheme[1][1].ππ = list(linspace(3, 9, 6 + 1))
        gte.scheme[1][1].effeff = 0.86
        gte.scheme[1][1].g_leak = 0.05

        gte.scheme[1][2].T_fuel = 40 + 273.15
        gte.scheme[1][2].η_burn = 0.99
        gte.scheme[1][2].TT_o = 1200  # list(linspace(800, 1200, 4 + 1))
        gte.scheme[1][2].T_lim = 1000
        gte.scheme[1][2].σ = 0.94
        gte.scheme[1][2].g_leak = 0

        gte.scheme[1][3].effeff = 0.92
        gte.scheme[1][3].η_mechanical = 0.99
        gte.scheme[1][3].T_lim = 1000
        gte.scheme[1][3].g_leak = 0.05

        gte.scheme[1][4].PP_o = 101325
        gte.scheme[1][4].eff = 0.96
        gte.scheme[1][4].v_ = 0.98
        gte.scheme[1][4].g_leak = 0.001

        """
        gte.validate_scheme()
        gte.export_gte_main()
        """

    if 0:
        gte = GTE("CFM-56")
        gte.scheme = {
            1: [Inlet(), Compressor(), CombustionChamber(), Turbine(), Outlet()],
            2: [Inlet(), Compressor(), Outlet()],
        }
        gte.scheme.show()
        gte.shafts = {1: [gte.scheme[1][1], gte.scheme[2][1], gte.scheme[1][3]]}

        gte.m = {1: 1}

        gte.mode.T = 288
        gte.mode.P = 101325
        gte.mode.H = 0
        gte.mode.M = 0

        gte.R = 10_000

        gte.substance = Substance({"N2": 0.755, "O2": 0.2315, "Ar": 0.01292, "Ne": 0.000014, "H2": 0.000008})
        gte.fuel = Substance({"KEROSENE": 1.0})

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

    """gte.describe()
    # gte.scheme.show()

    for e in gte.solve():
        e.summary()
        print(e.dataframe())
    """
