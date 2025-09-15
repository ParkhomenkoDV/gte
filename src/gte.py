import sys
from copy import deepcopy

import pandas as pd
from colorama import Fore
from numpy import arange, cos, inf, isnan, linspace, nan, prod, radians, sin
from tqdm import tqdm

sys.path.append("D:/Programming/Python/")

# from thermodynamics import FUELS, OXIDIZERS, Cp, R_gas
# from tools import eps, export2, isiter, isnum, to_roman

"""from nodes.combustion_chambler import CombustionChamber  # камера сгорания
from nodes.compressor import Compressor  # компрессор
from nodes.gear import Gear  # коробка приводов
from nodes.heat_exchanger import HeatExchanger  # теплообменный аппарат

# узлы ГТД
from nodes.inlet import Inlet  # вход
from nodes.load import Load  # нагрузка
from nodes.mixing_chamber import MixingChamber  # камера смешения
from nodes.nozzle import Nozzle  # сопло
from nodes.outlet import Outlet  # выход
from nodes.propeller import Propeller  # пропеллер
from nodes.turbine import Turbine  # турбина"""

# 1 = перед, 2 = посередине, 3 = после


def find_node_in_scheme(scheme, node2find) -> tuple:
    """Поиск положения узла в схеме"""
    for contour in scheme:
        for i, node in enumerate(scheme[contour]):
            if scheme[contour][i] is node2find:
                return contour, i


def get_combination(combination: int, max_list: list[int]) -> list[int]:
    n = max_list.copy()
    for i, max_element in enumerate(max_list):
        n[i] = combination % max_element
        combination //= max_element
    return n


class GTE_scheme:
    """Схема ГТД"""

    def __init__(self, scheme=None):
        if scheme is None:
            self.__scheme = {1: []}

    def __call__(self) -> dict:
        return self.__scheme

    def validate(self, scheme) -> bool:
        assert isinstance(scheme, dict), "type(scheme) is dict!"
        return True

    def show(self):
        pass


class GTE:
    """ГТД"""

    @classmethod
    def version(cls) -> str:
        version = 6.0
        next_version = (
            "камера смешения",
            "переход к массивам в местах постоянства диапазона значений"
            "теплак плак-плак",
            "упразднение outlet -> nozzle and outlet = переходный канал",
            "type(node) is class -> isinstance(node, class)"
            'ТД параметры через PTM 1677 as {"O2": 98, "H2": 2}',
            "соотношение соответствующих относительных расходов к своим контурам",
            "охлаждение турбины",
            "get_inlet_parameters() for all nodes",
            "get_outlet_parameters() for all nodes",
            "продолжение расчета",
            "multiprocessing",
            "ускорение расчета до 6000 [ГТД/с]",
        )
        print(f"{cls.__name__} version: {Fore.GREEN}{version}")
        for i, v in enumerate(next_version):
            print(cls.__name__ + " version:", int(version) + i + 1, v)
        return str(version)

    def __init__(self, name="GTE") -> None:
        self.name = name
        self.scheme = GTE_scheme()
        self.loads = []  # нагрузки/отборы мощности ГТД
        self.__path = list()  # путь обхода расчета ГТД
        self.m = {1: 1}  # степень контурности []
        self.G = {1: nan}  # абсолютный массовый расход [кг/с]
        self.N = nan  # мощность ГТД [Вт]

        self.M = nan  # Мах полета []
        self.v = nan  # скорость полета [м/с]

        self.warnings = {0: set(), 1: set(), 2: set(), 3: set()}  # предупреждения

    def __len__(self) -> int:
        """Количество узлов в ГТД"""
        i = len(self.loads)
        for contour in self.scheme:
            i += len(self.scheme[contour])
        return i

    def __str__(self) -> str:
        return self.name

    def get_variability(self) -> int:
        """Максимальное количество комбинаций варьируемых параметров"""
        return prod(
            [
                len(value)
                for key, value in self.__dict__.items()
                if type(value) is list and len(value) and not key.startswith("_")
            ]
        )

    def set_comb(self, gte_main, combination, combinations, max_combinations) -> None:
        pass

    def __set_combination(self, combination, gte_main) -> None:
        """Установка комбинации"""
        varible_params = [
            key
            for key, value in gte_main.__dict__.items()
            if type(value) is list and len(value) and not key.startswith("_")
        ]
        positions = [0] * len(varible_params)

        for i in range(combination):
            for j, param in enumerate(varible_params):
                if positions[j] == len(getattr(gte_main, varible_params[j])) - 1:
                    positions[j] = 0
                else:
                    positions[j] += 1
                    continue
        for j, param in enumerate(varible_params):
            setattr(
                self,
                varible_params[j],
                getattr(gte_main, varible_params[j])[positions[j]],
            )

    def gte_generator(self):
        """Генератор объектов ГТД с заданными варьируемыми параметрами"""
        max_combinations = [self.get_variability()]  # для ГТД
        for node in self.__path:
            max_combinations.append(node.get_variability())  # для узлов ГТД
        combinations = [1] * (1 + len(self.__path))

        for comb in tqdm(range(prod(max_combinations)), desc="Calculation", ncols=70):
            gte_var = deepcopy(self)
            gte_var.__set_combination(combinations[0], self)
            for i, node in enumerate(gte_var.__path):
                node.set_combination(combinations[i + 1], self.__path[i])
            gte_var.__update_combination(self, combinations, max_combinations)
            yield gte_var

        """for contour in self.scheme:
            for node in self.scheme[contour]:
                for attr in node.attrs:
                    if isiter(attr):
                        yield"""

    def __update_combination(self, gte_main, combinations, max_combinations) -> None:
        self.warnings = {
            0: set(),
            1: set(),
            2: set(),
            3: set(),
        }  # обнуление предупреждений
        if max_combinations[0] > 1:  # если есть варьируемые параметры ГТД
            if (
                combinations[0] < max_combinations[0]
            ):  # если не конец варьирования параметров ГТД
                self.__set_combination(
                    combinations[0], gte_main
                )  # установка текущих параметров варьирования
                combinations[0] += 1  # задание следующего номера комбинации
                return
            else:
                self.__set_combination(0, gte_main)
                combinations[0] = 1

        for i, node in enumerate(self.__path):
            node.warnings = {
                0: set(),
                1: set(),
                2: set(),
                3: set(),
            }  # обнуление предупреждений
            if (
                max_combinations[i + 1] > 1
            ):  # проверка наличия варьируемых параметров узлов ГТД
                if (
                    combinations[i + 1] < max_combinations[i + 1]
                ):  # проверка конца варьирования параметров узлов ГТД
                    node.set_combination(combinations[i + 1], gte_main.__path[i])
                    combinations[i + 1] += 1
                    return
                else:
                    node.set_combination(0, gte_main.__path[i])
                    combinations[i + 1] = 1

    @classmethod
    def input_action(cls) -> str:
        """Ввод действия"""
        while True:
            print(f"{Fore.YELLOW}Input action:", end=" ")
            print(""""a": add; "d": del; "e": end.""")
            action = input("action: ").strip().lower()
            if action not in ("a", "add", "d", "del", "delete", "e", "end"):
                print(f"{Fore.RED}No such action")
                continue
            return action

    def find_node_in_GTE(self, find_node) -> int:
        i = 0
        for contour in self.scheme:
            for node in self.scheme[contour]:
                i += 1
                if node is find_node:
                    return i
        for node in self.loads:
            i += 1
            if node is find_node:
                return i

    def input_contour(self) -> int:
        """Ввод существующего контура ГТД"""
        contour = 1
        if len(self.scheme) > 1:
            while True:  # ввод контура
                contour = input(
                    "Input the contour in scheme GTE from 1 to "
                    + str(len(self.scheme))
                    + ": "
                )
                if isnum(contour, type_num="int") and 1 <= int(contour) <= len(
                    self.scheme
                ):
                    contour = int(contour)
                    break
                else:
                    print(
                        Fore.RED
                        + "contour is int num in [1..{}]!".format(len(self.scheme))
                    )
        return contour

    def validate_scheme(self) -> bool:
        """Проверка на правильность ввода схемы"""
        print()
        if len(self) == 0:
            print(Fore.RED + "Empty GTE scheme!")
            return False

        try:
            self.set_path()
        except:
            print(Fore.RED + "Incorrect GTE solve path!")
            return False

        if len(self.__path) != len(self):
            print(Fore.RED + "Incorrect count items of GTE!")
            return False

        return True  # выполнены все проверки

    # TODO добавить камеру смешения
    def set_path(self) -> None:
        """Путь обхода расчета ГТД"""
        self.__path = []  # обнуление пути обхода
        for contour in self.scheme:  # для каждого контура
            for node in self.scheme[contour]:  # для каждого узла в данном контуре
                if type(node) is Turbine:  # если узел - это турбина
                    if node._shafts:  # если у турбины есть валы
                        for shaft in node._shafts:  # для каждого вала данной турбины
                            for (
                                tp,
                                els,
                            ) in shaft.items():  # для каждой связи данного вала
                                if tp in ("-", "+"):
                                    for el in els:
                                        if type(el) is Compressor:
                                            c, n = find_node_in_scheme(self.scheme, el)
                                            for i in range(n + 1):
                                                if self.scheme[c][i] not in self.__path:
                                                    self.__path.append(
                                                        self.scheme[c][i]
                                                    )
                                        elif type(el) is Load or type(el) is Propeller:
                                            if el not in self.__path:
                                                self.__path.append(el)
                                        else:
                                            raise "Invalid turbine shafts"
                                if tp == "0":
                                    for el in els:
                                        if el not in self.__path:
                                            self.__path.append(el)
                    else:
                        raise "Turbine object has no attribute shafts!"
                    self.__path.append(node)  # добавить турбину в путь обхода расчета
                else:  # иначе если узел не турбина и данного узла нет в пути обхода расчета
                    if node not in self.__path:
                        self.__path.append(node)

        for i, node in enumerate(self.__path):
            node.name = i
        self.output_scheme(numerate=True)
        print(
            f"\n{Fore.CYAN}Путь обхода расчета ГТД:",
            list(f"{node.name}: " + type(node).__name__ for node in self.__path),
        )

    def input_scheme(self) -> None:
        """Ввод схемы ГТД"""
        while True:
            print(f"\n{Fore.YELLOW}INPUT GTE SCHEME")
            self.scheme = self.scheme if self.scheme else {1: []}
            while True:
                print()
                self.output_scheme(numerate=True)  # вывод схемы
                action = self.input_action()  # ввод действия
                if action in ("a", "add"):
                    node = self.input_node()
                    if node in ("CONTOUR", "КОНТУР"):
                        self.scheme[max(self.scheme.keys()) + 1] = []
                    elif type(node) in (
                        Inlet,
                        Compressor,
                        CombustionChamber,
                        Turbine,
                        Nozzle,
                        Outlet,
                        HeatExchanger,
                    ):
                        contour = self.input_contour()
                        place = 0
                        if self.scheme[contour]:
                            while True:  # ввод положения в контуре
                                place = input(
                                    "Input the place in scheme GTE from 1 to "
                                    + str(len(self.scheme[contour]) + 1)
                                    + ": "
                                )
                                if (
                                    isnum(place, type_num="int")
                                    and 1 <= int(place) <= len(self.scheme[contour]) + 1
                                ):
                                    place = int(place)
                                    break
                                else:
                                    print(
                                        Fore.RED
                                        + "place is int num in [1..{}]!".format(
                                            len(self.scheme[contour]) + 1
                                        )
                                    )
                        if place == 0 or place == len(self.scheme[contour]) + 1:
                            self.scheme[contour].insert(place - 1, node)
                        else:
                            self.scheme[contour][place - 1] = node
                    elif type(node) in (Gear, Load, Propeller):
                        place = 0
                        if self.loads:
                            while True:
                                place = input(
                                    "Input the place in GTE loads from 1 to "
                                    + str(len(self.loads) + 1)
                                    + ": "
                                )
                                if (
                                    isnum(place, type_num="int")
                                    and 1 <= int(place) <= len(self.loads) + 1
                                ):
                                    place = int(place)
                                    break
                                else:
                                    print(
                                        Fore.RED
                                        + "place is int num in [1..{}]!".format(
                                            len(self.loads) + 1
                                        )
                                    )
                        if place == 0 or place == len(self.loads) + 1:
                            self.loads.insert(place - 1, node)
                        else:
                            self.loads[place - 1] = node
                elif action in ("d", "del", "delete"):
                    if self.count_items() == 0:
                        print(f"{Fore.RED}No place to delete!")
                        continue
                    node = self.input_node()
                    if node in ("CONTOUR", "КОНТУР"):
                        contour = self.input_contour()
                        try:
                            del self.scheme[contour]
                        except:
                            print(Fore.RED + "No contour to delete!")
                    else:
                        pos = self.count_items()
                        if pos > 1:
                            self.output_scheme(numerate=True)  # вывод схемы
                            while True:
                                pos = input(
                                    "Input item position from [1..{}]: ".format(
                                        self.count_items()
                                    )
                                )
                                if isnum(pos) and 1 <= int(pos) <= self.count_items():
                                    pos = int(pos)
                                    break
                                else:
                                    print(
                                        Fore.RED
                                        + "Position is int num in [1..{}]!".format(
                                            self.count_items()
                                        )
                                    )
                        i = 0
                        for c in self.scheme:
                            for _, n in enumerate(self.scheme[c]):
                                i += 1
                                if i == pos:
                                    self.scheme[c].pop(_)
                        for _, n in enumerate(self.loads):
                            i += 1
                            if i == pos:
                                self.loads.pop(_)
                elif action in ("e", "end"):
                    break

            self.input_shafts()  # ввод валов (связей)

            if self.validate_scheme:
                break

    def input_mixing_nodes(self) -> None:
        """Ввод узлов смешения"""
        print(Fore.YELLOW + "\nINPUT MIXING NODES")

    def input_shafts(self) -> None:
        """Ввод валов (связей)"""
        print(Fore.YELLOW + "\nINPUT GTE SHAFTS")

        d = []

        for contour in self.scheme:
            for i_node, node in enumerate(self.scheme[contour]):
                if type(node) is Turbine:
                    s = []
                    while True:
                        print()
                        self.output_scheme(numerate=True)
                        print(
                            f"Contour {to_roman(contour)} Turbine position {i_node + 1}"
                        )

                        while True:
                            pos = input(f"Input item position from [1..{len(self)}]: ")
                            if isnum(pos) and 1 <= int(pos) <= len(self):
                                pos = int(pos)
                                break
                            else:
                                print(
                                    Fore.RED
                                    + f"Position is int num in [1..{len(self)}]!"
                                )

                        i = 0
                        for c in self.scheme:
                            for n in self.scheme[c]:
                                i += 1
                                if i == pos:
                                    item = n
                        for n in self.loads:
                            i += 1
                            if i == pos:
                                item = n

                        if item not in d:
                            s = s if s else []
                            if type(item) is Gear:
                                s.append(item)
                                d.append(item)
                                continue
                            elif type(item) is Compressor:
                                s.append(item)
                                d.append(item)
                            elif type(item) is Load or type(item) is Propeller:
                                s.append(item)
                                d.append(item)
                            else:
                                print(Fore.RED + "Such item with Turbine!")
                                continue
                        else:
                            print(Fore.RED + "Such node already!")
                            continue

                        node.shafts.append(s)
                        s = []
                        break

    def input_contouring(self) -> None:
        """Ввод степеней контурности контуров"""
        print()
        for contour in self.scheme:
            if contour == 1:
                self.m.append({1: 1})
                continue
            else:
                correct_input = False
                while not correct_input:
                    m_list = [
                        m
                        for m in input(
                            "contour " + to_roman(contour) + ": m = "
                        ).split()
                    ]
                    for m in m_list:
                        if not isnum(m) or float(m) < 0:
                            print(Fore.RED + "m must be a number >= 0!")
                            correct_input = False
                            break
                        correct_input = True
                        self.m.append({contour: float(m)})

    def output_scheme(self, numerate=False) -> None:
        """Вывод схемы ГТД"""
        i = 0
        for contour in self.scheme:
            print(Fore.MAGENTA + "contour " + to_roman(contour), ":", sep="", end=" ")
            if numerate is True:
                temp_list = []
                for node in self.scheme[contour]:
                    i += 1
                    temp_list.append("{}: ".format(i) + type(node).__name__)
            else:
                temp_list = list(type(node).__name__ for node in self.scheme[contour])
            print(temp_list)
        if numerate is True:
            temp_list = []
            for node in self.loads:
                i += 1
                temp_list.append("{}: ".format(i) + type(node).__name__)
        else:
            temp_list = list(type(node).__name__ for node in self.loads)
        print(Fore.MAGENTA + "loads and gears:", temp_list)

    def input_GTE_parameters(self) -> None:
        """Ввод параметров ГТД"""
        correct_input = False
        while not correct_input:
            self.oxidizer_var = [ox.upper() for ox in input("oxidizer: ").split()]
            for oxidizer in self.oxidizer_var:
                if oxidizer not in OXIDIZERS:
                    print(Fore.RED + "No such oxidizer: " + oxidizer)
                    print("Possible options:", OXIDIZERS)
                    correct_input = False
                    break
                correct_input = True

        correct_input = False
        while not correct_input:
            self.fuel_var = [fuel.upper() for fuel in input("fuel: ").split()]
            for fuel in self.fuel_var:
                if fuel not in FUELS:
                    print(Fore.RED + "No such fuel: " + fuel)
                    print("Possible options:", FUELS)
                    correct_input = False
                    break
                correct_input = True

        correct_input = False
        while not correct_input:
            self.H_var = [H for H in input("H [м] = ").split()]
            for H in self.H_var:
                if not isnum(H):
                    print(Fore.RED + "H must be a number!")
                    correct_input = False
                    break
                correct_input = True
        self.H_var = [float(H) for H in self.H_var]

        while not self.M_var and not self.v_var:
            correct_input = False
            while not correct_input:
                self.M_var = [M.lower() for M in input("M [] = ").split()]
                if not self.M_var:
                    break
                for M in self.M_var:
                    if not isnum(M) or float(M) < 0:
                        print(Fore.RED + "M must me a number >= 0 or empty string!")
                        correct_input = False
                        break
                    correct_input = True
            self.M_var = [float(M) for M in self.M_var]

            if not self.M_var:
                correct_input = False
                while not correct_input:
                    self.v_var = [v.lower() for v in input("v [м/с] = ").split()]
                    if not self.v_var:
                        break
                    for v in self.v_var:
                        if not isnum(v) or float(v) < 0:
                            print(Fore.RED + "v must be a number >= 0 or empty string!")
                            correct_input = False
                            break
                        correct_input = True
                self.v_var = [float(v) for v in self.v_var]

        correct_input = False
        while not correct_input:
            self.R_var = [R for R in input("R [Н] = ").split()]
            for R in self.R_var:
                if not isnum(R) or float(R) < 0:
                    print(Fore.RED + "R must be a number >= 0!")
                    correct_input = False
                    break
                correct_input = True
        self.R_var = [float(R) for R in self.R_var]

        correct_input = False
        while not correct_input:
            self.resource_var = [res for res in input("resource [ч] = ").split()]
            for res in self.resource_var:
                if not isnum(res) or float(res) < 0:
                    print(Fore.RED + "resource is num > 0!")
                    correct_input = False
                    break
                correct_input = True
        self.resource_var = [float(res) * 3600 for res in self.resource_var]

    def input_node_parameters(self) -> None:
        """Ввод параметров узлов ГТД"""
        for contour in self.scheme:
            print(
                "Parameters of nodes of "
                + Fore.MAGENTA
                + "contour "
                + to_roman(contour),
                ":",
                sep="",
            )
            for node in self.scheme[contour]:
                print("Parameters of " + Fore.YELLOW + type(node).__name__, ":", sep="")
                node.input_parameters()

    def get_mass_flow(self) -> None:
        """Расчет абсолютного массового расхода"""
        N_η = 0  # N/η деление мощности потребителя на КПД передачи как комплекс
        Lt, Lc = 0, 0  # удельные работы турбин и компрессоров [Дж/кг]
        for node in self.__path:
            if type(node) is Turbine:
                Lt += node.L * self.m[find_node_in_scheme(self.scheme, node)[0]]
                for shaft in node._shafts:
                    N, η = 0, 1
                    for tp, els in shaft.items():
                        for el in els:
                            if type(el) is Compressor:
                                Lc += (
                                    el.L
                                    * self.m[find_node_in_scheme(self.scheme, el)[0]]
                                )
                            else:
                                if type(el) is Gear:
                                    η *= el.η
                                else:
                                    N += el.N
                        N_η += N / η

        try:
            self.G[1] = N_η / (Lt - Lc)
        except ZeroDivisionError:
            pass

        if self.R != 0 and not isnan(self.R_) and self.R_ != 0:
            G1 = (self.R / self.R_) / sum(self.m.values())
            if self.G[1] < G1:
                self.G[1] = (self.R / self.R_) / sum(self.m.values())
            else:
                self.warnings[3] = {"Не достигнута заданная тяга!"}

        for contour in self.scheme:
            self.G[contour] = self.G[1] * self.m[contour]

        self.G_fuel = self.g_fuel * self.G[1]

    def check_power(self) -> bool:
        """Расчет мощностей узлов ГТД"""
        N_plus, N_minus = 0, 0

        for contour in self.scheme:
            for node in self.scheme[contour]:
                if type(node) is Turbine:
                    node.N = node.L * self.G[contour]
                    N_plus += node.N
                elif type(node) is Compressor:
                    node.N = node.L * self.G[contour]
                    N_minus += node.N
        for node in self.loads:
            if type(node) is Load or type(node) is Propeller:
                N_minus += node.N
        self.N = N_plus - N_minus
        if eps("rel", N_plus, N_minus) <= error or isnan(eps("rel", N_plus, N_minus)):
            return True
        else:
            return False

    def __calculate(self, how="all", error: float = 0.01, Niter: int = 100, **kwargs):
        for iteration in range(Niter):
            for i, node in enumerate(self.__path):
                self.__path[i].solve(
                    how=how,
                    error=error,
                    Niter=Niter,  # параметры расчета
                    scheme=self.scheme,  # схема ГТД
                    m=self.m,  # степень контурностей
                    # рабочее тело, окислитель, горючее, теплоноситель
                    substance=self.substance,
                    fuel=self.fuel,
                    H=self.H,
                    M=self.M,
                    v=self.v,  # высотно-скоростные характеристики ГТД
                    resource=self.resource,  # ресурс ГТД
                    G=self.G[1],
                )
                if self.__path[i].warnings[3]:
                    return

            self.g_fuel = 0  # относительный массовый расход горючего []
            self.R_ = 0  # суммарная удельная тяга ГТД [м/с]
            for contour in self.scheme:
                for node in self.scheme[contour]:
                    if type(node) is CombustionChamber:
                        self.g_fuel += node.g_fuel * self.m[contour]
                    if hasattr(node, "R_"):
                        self.R_ += node.R_ * self.m[contour]
            self.R_ = self.R_ / sum(self.m.values())

            self.get_mass_flow()
            if self.check_power():
                break
        else:
            self.warnings[2].add("Нарушен баланс мощностей!")
            print(
                Fore.RED
                + f"Iteration limit in module {GTE.__name__} in function {self.solve.__name__}!"
            )

        self.G_fuel_N = self.G_fuel / self.N if self.N != 0 else inf
        self.G_fuel_R = self.G_fuel / self.R if self.R != 0 else inf

        # аэрокосмический научный журнал УДК 629.7.036.34 2016г
        # сопротивление газогенератора
        m = sum(self.m.values()) - 1
        Cx_1 = (
            2.0005 * 10 ** (-2) * (m - 1) ** 0
            - 0.1879 * 10 ** (-3) * (m - 1) ** 1
            + 2.8571 * 10 ** (-5) * (m - 1) ** 2
            - 3.4667 * 10 ** (-7) * (m - 1) ** 3
        )
        # сопротивление обечайки ГТД
        Cx_ = (
            3.5606 * 10 ** (-2) * (m - 1) ** 0
            - 1.2744 * 10 ** (-3) * (m - 1) ** 1
            + 3.3429 * 10 ** (-5) * (m - 1) ** 2
            - 3.6667 * 10 ** (-7) * (m - 1) ** 3
        )
        self.Cx = Cx_1 + Cx_  # коэффициент сопротивления мотогондолы

        """for contour in self.scheme:
            for node in self.scheme[contour]:
                for priority in node.warnings:
                    self.warnings[priority] = self.warnings | node.warnings
        for node in self.loads: self.warnings = self.warnings | node.warnings"""

        """for node in self.__path: print('{}: '.format(self.find_node_in_GTE(node)) + Fore.YELLOW + type(node).__name__,
                                       ': ', node.__dict__, sep='')
        print(Fore.YELLOW + GTE.__name__, ': ', self.__dict__, sep='')"""

    def solve(self, how="all", error: float = 0.01, Niter: int = 100, **kwargs):
        """Варьирование параметров ГТД"""
        data = list()
        for gte_var in self.gte_generator():
            gte_var.__calculate(how=how, error=error, Niter=Niter)
            data.append(
                {
                    # схема ГТД
                    **dict(
                        zip(
                            list(
                                f"contour{to_roman(contour)}"
                                for contour in gte_var.scheme
                            ),
                            list(
                                "+".join(
                                    list(
                                        type(node).__name__ + f"_{node.name}"
                                        for node in v
                                    )
                                )
                                for v in gte_var.scheme.values()
                            ),
                        )
                    ),
                    # нагрузка ГТД
                    "loads": "+".join(
                        list(
                            type(node).__name__ + f"_{node.name}"
                            for node in gte_var.loads
                        )
                    ),
                    # степень контурности контуров
                    **dict(
                        zip(
                            list(
                                f"m{to_roman(contour)} []" for contour in gte_var.scheme
                            ),
                            list(gte_var.m.values()),
                        )
                    ),
                    # предупреждения
                    # 'warnings': '; '.join(gte_var.warnings),
                    "H [м]": gte_var.H,
                    "M []": gte_var.M,
                    "R [кН]": gte_var.R / 1000,
                    "resource [ч]": gte_var.resource / 3600,
                    "oxidizer": gte_var.substance,
                    "fuel": gte_var.fuel,
                    "R_ [м/с]": getattr(gte_var, "R_", None),
                    "g fuel [%]": getattr(gte_var, "g_fuel", nan) * 100,
                    **dict(
                        zip(
                            list(
                                f"G{to_roman(contour)} [кг/с]"
                                for contour in gte_var.scheme
                            ),
                            list(gte_var.G.values()),
                        )
                    ),
                    "G fuel [г/с]": getattr(gte_var, "G_fuel", nan) * 1000,
                    "G fuel/N [г/Вт/ч]": getattr(gte_var, "G_fuel_N", nan)
                    * 1000
                    * 3600,
                    "G fuel/R [г/Н/ч]": getattr(gte_var, "G_fuel_R", nan) * 1000 * 3600,
                    "Cx []": getattr(self, "Cx", None),
                    # узлы ГТД
                    **dict(
                        zip(
                            list(
                                type(node).__name__ + f"_{node.name} " + k
                                for contour in gte_var.scheme
                                for node in gte_var.scheme[contour]
                                for k, v in node.__dict__.items()
                            ),
                            list(
                                v
                                for contour in gte_var.scheme
                                for node in gte_var.scheme[contour]
                                for k, v in node.__dict__.items()
                            ),
                        )
                    ),
                    # нагрузка ГТД
                    **dict(
                        zip(
                            list(
                                type(node).__name__ + f"_{node.name} " + k
                                for node in gte_var.loads
                                for k, v in node.__dict__.items()
                            ),
                            list(
                                v
                                for node in gte_var.loads
                                for k, v in node.__dict__.items()
                            ),
                        )
                    ),
                }
            )
        pd.set_option("display.max_row", None)
        pd.set_option("display.max_columns", None)
        pd.set_option("display.width", None)
        # tqdm.pandas(desc="Converting to DataFrame")
        df = pd.DataFrame(data)  # .progress_apply(lambda x: x)
        # print(df, end='\n')
        print(df.info(memory_usage="deep"))
        export2(
            df,
            file_path="exports/" + self.name,
            file_name=self.name,
            file_extension=kwargs.get("file_type", "pkl"),
            sheet_name="Cycle",
            show_time=True,
        )
        # print(f'{Fore.YELLOW}Elapsed memory: {round(getsizeof(df) / 1024 / 1024, 3)} Mb')

    def export_gte_main(self) -> None:
        data = {  # схема ГТД
            **dict(
                zip(
                    list(f"contour{to_roman(contour)}" for contour in self.scheme),
                    list(
                        "+".join(
                            list(type(node).__name__ + f"_{node.name}" for node in v)
                        )
                        for v in self.scheme.values()
                    ),
                )
            ),
            # нагрузка ГТД
            "loads": "+".join(
                list(type(node).__name__ + f"_{node.name}" for node in self.loads)
            ),
            # степень контурности контуров
            **{
                f"m{to_roman(key)} []": [[dct[key] for dct in self.m]]
                for key in self.m[0]
            },
            "H [м]": [self.H],
            "M []": [self.M],
            "R [кН]": [[R / 1000 for R in self.R]],
            "resource [ч]": [[resource / 3600 for resource in self.resource]],
            "oxidizer": [self.substance],
            "fuel": [self.fuel],
            **dict(
                zip(
                    list(
                        type(node).__name__ + f"_{node.name} " + k
                        for contour in self.scheme
                        for node in self.scheme[contour]
                        for k, v in node.__dict__.items()
                        if type(v) is list
                    ),
                    list(
                        str(v)
                        for contour in self.scheme
                        for node in self.scheme[contour]
                        for k, v in node.__dict__.items()
                        if type(v) is list
                    ),
                )
            ),
            **dict(
                zip(
                    list(
                        type(node).__name__ + f"_{node.name} " + k
                        for node in self.loads
                        for k, v in node.__dict__.items()
                        if type(v) is list
                    ),
                    list(
                        str(v)
                        for node in self.loads
                        for k, v in node.__dict__.items()
                        if type(v) is list
                    ),
                )
            ),
        }
        pd.set_option("display.max_row", None)
        pd.set_option("display.max_columns", None)
        pd.set_option("display.width", None)
        df = pd.DataFrame(data)
        print(df.T, end="\n")
        export2(
            df,
            file_path="exports/" + self.name,
            file_name=self.name + " input",
            file_extension="xlsx",
            sheet_name="input",
            show_time=True,
            index=False,
        )


class GTE_inputs(GTE):
    """Ручной ввод параметров ГТД"""

    pass


# TODO найти библиотеку перевода текста на разные языки
def input_language():
    """Ввод языка"""
    while True:
        lang = input("Input language (en, ru): ").strip().lower()
        if lang in ("en", "eng", "english"):
            lang = "en"
            return lang
        elif lang in ("ru", "rus", "russian"):
            lang = "ru"
            return lang
        else:
            print(Fore.RED + "No such language!")


if __name__ == "__main__":
    how = "cycle"
    error = 0.01
    Niter = 8

    gte = GTE("АД")
    if 0:
        print(f"{Fore.CYAN}{gte.name}")
        gte.scheme = {
            1: [
                Inlet(),
                Compressor("a"),
                Compressor("a"),
                CombustionChamber(),
                Turbine("a"),
                Turbine("a"),
            ],
            2: [Inlet(), Compressor("a"), MixingChamber(), Nozzle()],
        }
        gte.loads = []

        gte.m = [{1: 1, 2: 1.2}]
        """[{1: 1, 2: 0.6}, {1: 1, 2: 0.7}, {1: 1, 2: 0.8}, {1: 1, 2: 0.9}, {1: 1, 2: 1.0},
                 {1: 1, 2: 1.1}, {1: 1, 2: 1.2}, {1: 1, 2: 1.3}, {1: 1, 2: 1.4}, {1: 1, 2: 1.5}]"""

        gte.H = [0]  # list(linspace(-2_000, 11_000, 13 + 1))
        gte.M = [0]  # list(linspace(0, 0.8, 8 + 1))
        gte.R = [40_000]

        gte.substance = ["AIR"]
        gte.fuel = ["KEROSENE"]
        gte.resource = [3_000 * 3600]

        gte.scheme[1][0].σ = [0.95]
        gte.scheme[1][0].g_leak = [0]

        gte.scheme[1][1].ππ = [3]
        gte.scheme[1][1].ηη = [0.89]
        gte.scheme[1][1].g_leak = [0]

        gte.scheme[1][2].ππ = list(linspace(4, 12, 8 + 1))
        gte.scheme[1][2].ηη = [0.87]
        gte.scheme[1][2].g_leak = [0.05]

        gte.scheme[1][3].T_fuel = [40 + 273.15]
        gte.scheme[1][3].η_burn = [0.99]
        gte.scheme[1][3].TT3 = [1600]  # list(linspace(1300, 1600, 6 + 1))
        gte.scheme[1][3].σ = [0.94]
        gte.scheme[1][3].T_lim = [1200]
        gte.scheme[1][3].g_leak = [0]

        gte.scheme[1][4]._shafts = [{"-": [gte.scheme[1][2]]}]
        gte.scheme[1][4].ηη = [0.92]
        gte.scheme[1][4].η_mechanical = [0.99]
        gte.scheme[1][4].T_lim = [1200]
        gte.scheme[1][4].g_leak = [0.05]

        gte.scheme[1][5]._shafts = [{"-": [gte.scheme[1][1], gte.scheme[2][1]]}]
        gte.scheme[1][5].ηη = [0.92]
        gte.scheme[1][5].η_mechanical = [0.99]
        gte.scheme[1][5].T_lim = [1200]
        gte.scheme[1][5].g_leak = [0.05]

        gte.scheme[2][0].σ = [0.95]
        gte.scheme[2][0].g_leak = [0.005]

        gte.scheme[2][1].ππ = [3]
        gte.scheme[2][1].ηη = [0.89]
        gte.scheme[2][1].g_leak = [0.05]

        gte.scheme[2][2]._nodes_add = [gte.scheme[1][5]]
        gte.scheme[2][2].g_leak = [0.05]

        gte.scheme[2][3].PP3 = [101325]
        gte.scheme[2][3].ηη = [0.92]
        gte.scheme[2][3].v_ = [0.99]
        gte.scheme[2][3].g_leak = [0.001]

        gte.validate_scheme()
        gte.export_gte_main()
        gte.solve(how="cycle", error=error, Niter=Niter, file_type="xlsx")

    gte = GTE("Jumo 004b")
    if 1:
        print(f"{Fore.CYAN}{gte}")
        gte.scheme = {
            1: [Inlet(), Compressor("a"), CombustionChamber(), Turbine("a"), Nozzle()]
        }
        gte.loads = []

        gte.m = [{1: 1}]

        gte.H = [0]
        gte.M = [0]
        gte.R = [10_000]
        gte.resource = [150 * 3600]

        gte.substance = ["AIR"]
        gte.fuel = ["KEROSENE"]

        gte.scheme[1][0].σ = [0.98]
        gte.scheme[1][0].g_leak = [0.005]

        gte.scheme[1][1].ππ = list(linspace(3, 43, 40 + 1))
        gte.scheme[1][1].ηη = [0.86]
        gte.scheme[1][1].g_leak = [0.05]

        gte.scheme[1][2].T_fuel = [40 + 273.15]
        gte.scheme[1][2].η_burn = [0.99]
        gte.scheme[1][2].TT3 = list(linspace(800, 1200, 8 + 1))
        gte.scheme[1][2].T_lim = [1000]
        gte.scheme[1][2].σ = [0.94]
        gte.scheme[1][2].g_leak = [0]

        gte.scheme[1][3]._shafts = [{"-": [gte.scheme[1][1]]}]
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
        gte.solve(how=how, error=error, Niter=Niter, file_type="xlsx")

    gte = GTE("ТВ3-117")
    if 0:
        print(f"{Fore.CYAN}{gte.name}")
        gte.scheme = {
            1: [
                Inlet(),
                Compressor("a"),
                CombustionChamber(),
                Turbine("a"),
                Turbine("a"),
                Nozzle(),
            ]
        }
        gte.loads = [Gear(), Propeller()]
        gte.loads[0].η = [0.92]
        gte.loads[1].N = [2.2 * 10**6]

        gte.m = [{1: 1}]

        gte.H = list(linspace(-2000, 6000, 8 + 1))
        gte.M = list(linspace(0, 0.6, 6 + 1))
        gte.R = [0]
        gte.resource = [3000 * 3600]

        gte.substance = ["AIR"]
        gte.fuel = ["KEROSENE"]

        gte.scheme[1][0].σ = [0.98]
        gte.scheme[1][0].g_leak = [0.005]

        gte.scheme[1][1].ππ = list(linspace(5, 30, 25 + 1))
        gte.scheme[1][1].ηη = [0.86]
        gte.scheme[1][1].g_leak = [0.05]

        gte.scheme[1][2].T_fuel = [40 + 273.15]
        gte.scheme[1][2].η_burn = [0.99]
        gte.scheme[1][2].TT3 = list(linspace(1200, 1600, 8 + 1))
        gte.scheme[1][2].σ = [0.94]
        gte.scheme[1][2].T_lim = [1000]
        gte.scheme[1][2].g_leak = [0]

        gte.scheme[1][3]._shafts = [{"-": [gte.scheme[1][1]]}]
        gte.scheme[1][3].ηη = [0.92]
        gte.scheme[1][3].η_mechanical = [0.99]
        gte.scheme[1][3].T_lim = [1200]
        gte.scheme[1][3].g_leak = [0.05]

        gte.scheme[1][4].shafts = [[gte.loads[0], gte.loads[1]]]
        gte.scheme[1][4].PP3 = [10**5]
        gte.scheme[1][4].ηη = [0.92]
        gte.scheme[1][4].η_mechanical = [0.99]
        gte.scheme[1][4].T_lim = [1200]
        gte.scheme[1][4].g_leak = [0.05]

        gte.scheme[1][5].σ = [0.98]
        gte.scheme[1][5].g_leak = [0.005]

        gte.validate_scheme()
        gte.export_gte_main()
        gte.solve(how=how, error=error, Niter=Niter, file_type="xlsx")

    gte = GTE("ТВ7-117")
    if 0:
        print(f"{Fore.CYAN}{gte.name}")
        gte.scheme = {
            1: [
                Inlet(),
                Compressor("a"),
                Compressor("r"),
                CombustionChamber(),
                Turbine("a"),
                Turbine("a"),
                Turbine("a"),
                Outlet(),
            ]
        }
        gte.loads = [Gear(), Propeller()]
        gte.loads[0].η = [0.92]
        gte.loads[1].N = [2.2 * 10**6]

        gte.m = [{1: 1}]

        gte.H = list(linspace(-2000, 6000, 8 + 1))
        gte.M = list(linspace(0, 0.6, 6 + 1))
        gte.R = [0]
        gte.resource = [3000 * 3600]

        gte.substance = ["AIR"]
        gte.fuel = ["KEROSENE"]

        gte.scheme[1][0].σ = [0.98]
        gte.scheme[1][0].g_leak = [0.005]

        gte.scheme[1][1].ππ = list(linspace(2, 6, 4 + 1))
        gte.scheme[1][1].ηη = [0.86]
        gte.scheme[1][1].g_leak = [0.05]

        gte.scheme[1][2].ππ = list(linspace(2, 6, 4 + 1))
        gte.scheme[1][2].ηη = [0.86]
        gte.scheme[1][2].g_leak = [0.05]

        gte.scheme[1][3].T_fuel = [40 + 273.15]
        gte.scheme[1][3].η_burn = [0.99]
        gte.scheme[1][3].TT3 = list(linspace(1200, 1600, 8 + 1))
        gte.scheme[1][3].σ = [0.94]
        gte.scheme[1][3].T_lim = [1000]
        gte.scheme[1][3].g_leak = [0]

        gte.scheme[1][4].shafts = [[gte.scheme[1][2]]]
        gte.scheme[1][4].ηη = [0.92]
        gte.scheme[1][4].η_mechanical = [0.99]
        gte.scheme[1][4].T_lim = [1200]
        gte.scheme[1][4].g_leak = [0.05]

        gte.scheme[1][5].shafts = [[gte.scheme[1][1]]]
        gte.scheme[1][5].ηη = [0.92]
        gte.scheme[1][5].η_mechanical = [0.99]
        gte.scheme[1][5].T_lim = [1000]
        gte.scheme[1][5].g_leak = [0.05]

        gte.scheme[1][6].shafts = [[gte.loads[0], gte.loads[1]]]
        gte.scheme[1][6].PP3 = [10**5]
        gte.scheme[1][6].ηη = [0.92]
        gte.scheme[1][6].η_mechanical = [0.99]
        gte.scheme[1][6].T_lim = [1200]
        gte.scheme[1][6].g_leak = [0.05]

        gte.scheme[1][7].σ = [0.98]
        gte.scheme[1][7].g_leak = [0.005]

        gte.validate_scheme()
        gte.export_gte_main()
        gte.solve(how=how, error=error, Niter=Niter, file_type="xlsx")

    gte = GTE("ПС-90А")
    if 0:
        print(f"{Fore.CYAN}{gte.name}")
        gte.scheme = {
            1: [
                Inlet(),
                Compressor("a"),
                Compressor("a"),
                CombustionChamber(),
                Turbine("a"),
                Turbine("a"),
                Nozzle(),
            ],
            2: [Inlet(), Compressor("a"), Nozzle()],
        }
        gte.loads = []

        gte.m = [
            {1: 1, 2: 3},
            {1: 1, 2: 3.5},
            {1: 1, 2: 4},
            {1: 1, 2: 4.5},
            {1: 1, 2: 5},
        ]

        gte.H = [0, 11000]
        gte.M = [0, 0.8]
        gte.R = [10_000, 120_000]

        gte.substance = ["AIR"]
        gte.fuel = ["KEROSENE"]
        gte.resource = [5000 * 3600]

        gte.scheme[1][0].σ = [0.98]
        gte.scheme[1][0].g_leak = [0.005]

        gte.scheme[1][1].ππ = list(linspace(1.2, 3.2, 2 + 1))
        gte.scheme[1][1].ηη = [0.88]
        gte.scheme[1][1].g_leak = [0.05]

        gte.scheme[1][2].ππ = list(linspace(5, 25, 20 + 1))
        gte.scheme[1][2].ηη = [0.87]
        gte.scheme[1][2].g_leak = [0.05]

        gte.scheme[1][3].T_fuel = [40 + 273.15]
        gte.scheme[1][3].η_burn = [0.99]
        gte.scheme[1][3].TT3 = list(linspace(1700, 1800, 12 + 1))
        gte.scheme[1][3].σ = [0.94]
        gte.scheme[1][3].T_lim = [1800]
        gte.scheme[1][3].g_leak = [0]

        gte.scheme[1][4].shafts = [[gte.scheme[1][2]]]
        gte.scheme[1][4].ηη = [0.92]
        gte.scheme[1][4].η_mechanical = [0.99]
        gte.scheme[1][4].T_lim = [1200]
        gte.scheme[1][4].g_leak = [0.05]

        gte.scheme[1][5].shafts = [{"-": [gte.scheme[1][1], gte.scheme[2][1]]}]
        gte.scheme[1][5].ηη = [0.92]
        gte.scheme[1][5].η_mechanical = [0.99]
        gte.scheme[1][5].T_lim = [1000]
        gte.scheme[1][5].g_leak = [0.05]

        gte.scheme[1][6].ππ = [1.13]
        gte.scheme[1][6].ηη = [0.96]
        gte.scheme[1][6].v_ = [0.98]
        gte.scheme[1][6].g_leak = [0.001]

        gte.scheme[2][0].σ = [0.98]
        gte.scheme[2][0].g_leak = [0.005]

        gte.scheme[2][1].ππ = [4]
        gte.scheme[2][1].ηη = [0.86]
        gte.scheme[2][1].g_leak = [0.05]

        gte.scheme[2][2].ππ = [1.13]
        gte.scheme[2][2].ηη = [0.96]
        gte.scheme[2][2].v_ = [0.98]
        gte.scheme[2][2].g_leak = [0.001]

        gte.validate_scheme()
        gte.export_gte_main()
        gte.solve(how="cycle", error=error, Niter=Niter, file_type="xlsx")

    gte = GTE("CFM-56")
    if 0:
        print(f"{Fore.CYAN}{gte.name}")
        gte.scheme = {
            1: [
                Inlet(),
                Compressor("a"),
                Compressor("a"),
                CombustionChamber(),
                Turbine("a"),
                Turbine("a"),
                Nozzle(),
            ],
            2: [Inlet(), Compressor("a"), Nozzle()],
        }
        gte.loads = []

        gte.m = [{1: 1, 2: 4}, {1: 1, 2: 5}, {1: 1, 2: 6}, {1: 1, 2: 7}, {1: 1, 2: 8}]

        gte.H = [11_000]
        gte.M = [0.8]
        gte.R = [80_000]

        gte.substance = ["AIR"]
        gte.fuel = ["KEROSENE"]
        gte.resource = [5000 * 3600]

        gte.scheme[1][0].σ = [0.98]
        gte.scheme[1][0].g_leak = [0.005]

        gte.scheme[1][1].ππ = [1.2]
        gte.scheme[1][1].ηη = [0.88]
        gte.scheme[1][1].g_leak = [0]

        gte.scheme[1][2].ππ = [12]
        gte.scheme[1][2].ηη = [0.87]
        gte.scheme[1][2].g_leak = [0.05]

        gte.scheme[1][3].T_fuel = [40 + 273.15]
        gte.scheme[1][3].η_burn = [0.99]
        gte.scheme[1][3].TT3 = list(linspace(1600, 1800, 4 + 1))
        gte.scheme[1][3].σ = [0.94]
        gte.scheme[1][3].T_lim = [1800]
        gte.scheme[1][3].g_leak = [0]

        gte.scheme[1][4].shafts = [[gte.scheme[1][2]]]
        gte.scheme[1][4].ηη = [0.92]
        gte.scheme[1][4].η_mechanical = [0.99]
        gte.scheme[1][4].T_lim = [1200]
        gte.scheme[1][4].g_leak = [0.05]

        gte.scheme[1][5].shafts = [[gte.scheme[1][1]], [gte.scheme[2][1]]]
        gte.scheme[1][5].ηη = [0.92]
        gte.scheme[1][5].η_mechanical = [0.99]
        gte.scheme[1][5].T_lim = [1000]
        gte.scheme[1][5].g_leak = [0.05]

        gte.scheme[1][6].PP3 = [66000]
        gte.scheme[1][6].ηη = [0.96]
        gte.scheme[1][6].v_ = [0.98]
        gte.scheme[1][6].g_leak = [0.001]

        gte.scheme[2][0].σ = [0.98]
        gte.scheme[2][0].g_leak = [0.005]

        gte.scheme[2][1].ππ = [2]  # list(linspace(2, 4, 4 + 1))
        gte.scheme[2][1].ηη = [0.86]
        gte.scheme[2][1].g_leak = [0.005]

        gte.scheme[2][2].ππ = [1.13]
        gte.scheme[2][2].ηη = [0.96]
        gte.scheme[2][2].v_ = [0.98]
        gte.scheme[2][2].g_leak = [0.001]

        gte.validate_scheme()
        gte.export_gte_main()
        gte.solve(how="cycle", error=error, Niter=Niter, file_type="xlsx")

    gte = GTE("RR Trent")
    if 0:
        print(f"{Fore.CYAN}{gte.name}")
        gte.scheme = {
            1: [
                Inlet(),
                Compressor("a"),
                Compressor("a"),
                Compressor("a"),
                CombustionChamber(),
                Turbine("a"),
                Turbine("a"),
                Turbine("a"),
                Nozzle(),
            ],
            2: [Inlet(), Compressor("a"), Nozzle()],
        }
        gte.loads = [Gear()]

        gte.m = [
            {1: 1, 2: 6},
            {1: 1, 2: 7},
            {1: 1, 2: 8},
            {1: 1, 2: 9},
            {1: 1, 2: 10},
            {1: 1, 2: 11},
            {1: 1, 2: 12},
        ]

        gte.loads[0].η = [0.92]

        gte.H = [11_000]
        gte.M = [0.8]
        gte.R = [160_000]

        gte.substance = ["AIR"]
        gte.fuel = ["KEROSENE"]
        gte.resource = [5000 * 3600]

        gte.scheme[1][0].σ = [0.98]
        gte.scheme[1][0].g_leak = [0.005]

        gte.scheme[1][1].ππ = list(linspace(1.2, 3.2, 2 + 1))
        gte.scheme[1][1].ηη = [0.88]
        gte.scheme[1][1].g_leak = [0.05]

        gte.scheme[1][2].ππ = list(linspace(4, 8, 4 + 1))
        gte.scheme[1][2].ηη = [0.87]
        gte.scheme[1][2].g_leak = [0.05]

        gte.scheme[1][3].ππ = list(linspace(4, 8, 4 + 1))
        gte.scheme[1][3].ηη = [0.86]
        gte.scheme[1][3].g_leak = [0.05]

        gte.scheme[1][4].T_fuel = [40 + 273.15]
        gte.scheme[1][4].η_burn = [0.99]
        gte.scheme[1][4].TT3 = list(linspace(1200, 1800, 12 + 1))
        gte.scheme[1][4].σ = [0.94]
        gte.scheme[1][4].T_lim = [1800]
        gte.scheme[1][4].g_leak = [0]

        gte.scheme[1][5].shafts = [[gte.scheme[1][3]]]
        gte.scheme[1][5].ηη = [0.92]
        gte.scheme[1][5].η_mechanical = [0.99]
        gte.scheme[1][5].T_lim = [1200]
        gte.scheme[1][5].g_leak = [0.05]

        gte.scheme[1][6].shafts = [[gte.scheme[1][2]]]
        gte.scheme[1][6].ηη = [0.92]
        gte.scheme[1][6].η_mechanical = [0.99]
        gte.scheme[1][6].T_lim = [1200]
        gte.scheme[1][6].g_leak = [0.05]

        gte.scheme[1][7].shafts = [
            [gte.loads[0], gte.scheme[1][1]],
            [gte.loads[0], gte.scheme[2][1]],
        ]
        gte.scheme[1][7].ηη = [0.92]
        gte.scheme[1][7].η_mechanical = [0.99]
        gte.scheme[1][7].T_lim = [1000]
        gte.scheme[1][7].g_leak = [0.05]

        gte.scheme[1][8].ππ = [1.13]
        gte.scheme[1][8].ηη = [0.96]
        gte.scheme[1][8].v_ = [0.98]
        gte.scheme[1][8].g_leak = [0.001]

        gte.scheme[2][0].σ = [0.98]
        gte.scheme[2][0].g_leak = [0.005]

        gte.scheme[2][1].ππ = [1.13]
        gte.scheme[2][1].ηη = [0.86]
        gte.scheme[2][1].g_leak = [0.05]

        gte.scheme[2][2].ππ = [1.13]
        gte.scheme[2][2].ηη = [0.92]
        gte.scheme[2][2].v_ = [0.99]
        gte.scheme[2][2].g_leak = [0.001]

        gte.validate_scheme()
        gte.export_gte_main()
        gte.solve(how="cycle", error=error, Niter=Niter, file_type="xlsx")

    gte = GTE("АИ-222-25")
    if 0:
        print(f"{Fore.CYAN}{gte.name}")
        gte.scheme = {
            1: [
                Inlet(),
                Compressor("a"),
                Compressor("a"),
                CombustionChamber(),
                Turbine("a"),
                Turbine("a"),
            ],
            2: [Inlet(), Compressor("a"), MixingChamber(), Nozzle()],
        }
        gte.loads = []

        gte.m = [
            {1: 1, 2: 0.8},
            {1: 1, 2: 0.9},
            {1: 1, 2: 1.0},
            {1: 1, 2: 1.1},
            {1: 1, 2: 1.2},
            {1: 1, 2: 1.3},
            {1: 1, 2: 1.4},
        ]

        gte.H = list(linspace(0, 6_000, 6 + 1))
        gte.M = list(linspace(0, 0.6, 6 + 1))
        gte.R = [25_000]

        gte.substance = ["AIR"]
        gte.fuel = ["KEROSENE"]
        gte.resource = [3000 * 3600]

        gte.scheme[1][0].σ = [0.98]
        gte.scheme[1][0].g_leak = [0.005]

        gte.scheme[1][1].ππ = [3]
        gte.scheme[1][1].ηη = [0.88]
        gte.scheme[1][1].g_leak = [0.05]

        gte.scheme[1][2].ππ = list(linspace(10, 20, 10 + 1))
        gte.scheme[1][2].ηη = [0.87]
        gte.scheme[1][2].g_leak = [0.05]

        gte.scheme[1][3].T_fuel = [40 + 273.15]
        gte.scheme[1][3].η_burn = [0.99]
        gte.scheme[1][3].TT3 = list(linspace(1300, 1600, 6 + 1))
        gte.scheme[1][3].σ = [0.94]
        gte.scheme[1][3].T_lim = [1700]
        gte.scheme[1][3].g_leak = [0]

        gte.scheme[1][4].shafts = [[gte.scheme[1][2]]]
        gte.scheme[1][4].ηη = [0.92]
        gte.scheme[1][4].η_mechanical = [0.99]
        gte.scheme[1][4].T_lim = [1200]
        gte.scheme[1][4].g_leak = [0.05]

        gte.scheme[1][5].shafts = [[gte.scheme[1][1]], [gte.scheme[2][1]]]
        gte.scheme[1][5].ηη = [0.92]
        gte.scheme[1][5].η_mechanical = [0.99]
        gte.scheme[1][5].T_lim = [1200]
        gte.scheme[1][5].g_leak = [0.05]

        gte.scheme[2][0].σ = [0.98]
        gte.scheme[2][0].g_leak = [0.005]

        gte.scheme[2][1].ππ = [2]
        gte.scheme[2][1].ηη = [0.86]
        gte.scheme[2][1].g_leak = [0.05]

        gte.scheme[2][2].node_add = [gte.scheme[1][5]]
        gte.scheme[2][2].g_leak = [0.05]

        gte.scheme[2][3].ππ = [1.13]
        gte.scheme[2][3].ηη = [0.92]
        gte.scheme[2][3].v_ = [0.99]
        gte.scheme[2][3].g_leak = [0.001]

        gte.validate_scheme()
        gte.export_gte_main()
        gte.solve(how="cycle", error=error, Niter=Niter, file_type="xlsx")

    gte = GTE("АЛ-31Ф")
    if 0:
        print(f"{Fore.CYAN}{gte.name}")
        gte.scheme = {
            1: [
                Inlet(),
                Compressor("a"),
                Compressor("a"),
                CombustionChamber(),
                Turbine("a"),
                Turbine("a"),
            ],
            2: [
                Inlet(),
                Compressor("a"),
                HeatExchanger(),
                MixingChamber(),
                CombustionChamber(),
                Nozzle(),
            ],
        }
        gte.loads = []

        gte.m = [
            {1: 1, 2: 0.2},
            {1: 1, 2: 0.4},
            {1: 1, 2: 0.6},
            {1: 1, 2: 0.8},
            {1: 1, 2: 1.0},
            {1: 1, 2: 1.2},
            {1: 1, 2: 1.4},
            {1: 1, 2: 1.6},
            {1: 1, 2: 1.8},
            {1: 1, 2: 2.0},
        ]

        gte.H = list(linspace(-2_000, 14_000, 16 + 1))
        gte.M = list(linspace(0, 3, 30 + 1))
        gte.R = [12_000]

        gte.substance = ["AIR"]
        gte.fuel = ["KEROSENE"]
        gte.resource = [3000 * 3600]

        gte.scheme[1][0].σ = [0.98]
        gte.scheme[1][0].g_leak = [0.005]

        gte.scheme[1][1].ππ = list(linspace(2, 4, 2 + 1))
        gte.scheme[1][1].ηη = [0.88]
        gte.scheme[1][1].g_leak = [0.05]

        gte.scheme[1][2].ππ = list(linspace(10, 20, 10 + 1))
        gte.scheme[1][2].ηη = [0.87]
        gte.scheme[1][2].g_leak = [0.05]

        gte.scheme[1][3].T_fuel = [40 + 273.15]
        gte.scheme[1][3].η_burn = [0.99]
        gte.scheme[1][3].TT3 = list(linspace(1300, 1700, 5 + 1))
        gte.scheme[1][3].σ = [0.94]
        gte.scheme[1][3].T_lim = [1200]
        gte.scheme[1][3].g_leak = [0]

        gte.scheme[1][4].shafts = [[gte.scheme[1][2]]]
        gte.scheme[1][4].ηη = [0.92]
        gte.scheme[1][4].η_mechanical = [0.99]
        gte.scheme[1][4].T_lim = [1200]
        gte.scheme[1][4].g_leak = [0.05]

        gte.scheme[1][5].shafts = [[gte.scheme[1][1]], [gte.scheme[2][1]]]
        gte.scheme[1][5].ηη = [0.92]
        gte.scheme[1][5].η_mechanical = [0.99]
        gte.scheme[1][5].T_lim = [1200]
        gte.scheme[1][5].g_leak = [0.05]

        gte.scheme[2][0].σ = [0.98]
        gte.scheme[2][0].g_leak = [0.005]

        gte.scheme[2][1].ππ = [2]
        gte.scheme[2][1].ηη = [0.86]
        gte.scheme[2][1].g_leak = [0.05]

        gte.scheme[2][2].node_add = [gte.scheme[1][6]]
        gte.scheme[2][2].g_leak = [0.05]

        gte.scheme[2][3].ππ = [1.13]
        gte.scheme[2][3].ηη = [0.92]
        gte.scheme[2][3].v_ = [0.99]
        gte.scheme[2][3].g_leak = [0.001]

        gte.validate_scheme()
        gte.export_gte_main()
        gte.solve(error=error, Niter=Niter, file_type="xlsx")

    gte = GTE("ГСУ")
    if 0:
        print(f"{Fore.CYAN}{gte.name}")
        gte.scheme = {
            1: [Inlet(), Compressor(), CombustionChamber(), Turbine(), Outlet()]
        }
        gte.loads = [Gear(), Load()]

        gte.loads[0].η = list(linspace(0.8, 1, 10 + 1))
        gte.loads[1].N = list(linspace(100 * 10**3, 500 * 10**3, 4 + 1))

        gte.m = [{1: 1}]

        gte.H = list(linspace(0, 6_000, 12 + 1))
        gte.M = list(linspace(0, 0.6, 6 + 1))
        gte.R = [0]

        gte.substance = ["AIR"]
        gte.fuel = ["KEROSENE"]
        gte.resource = [3000 * 3600]

        gte.scheme[1][0].σ = [0.98]
        gte.scheme[1][0].g_leak = [0.005]

        gte.scheme[1][1].ππ = list(linspace(1.2, 12.2, 11 + 1))
        gte.scheme[1][1].ηη = list(linspace(0.78, 0.92, 14 + 1))
        gte.scheme[1][1].g_leak = [0.05]

        gte.scheme[1][2].T_fuel = [5 + 273.15]
        gte.scheme[1][2].η_burn = [0.98]
        gte.scheme[1][2].TT3 = list(linspace(800, 1000, 4 + 1))
        gte.scheme[1][2].σ = [0.94]
        gte.scheme[1][2].T_lim = [1000]
        gte.scheme[1][2].g_leak = [0]

        gte.scheme[1][3].shafts = [[gte.scheme[1][1]], [gte.loads[0], gte.loads[1]]]
        gte.scheme[1][3].ηη = [0.9]  # list(linspace(0.78, 0.92, 14 + 1))
        gte.scheme[1][3].η_mechanical = [0.99]
        gte.scheme[1][3].T_lim = [1000]
        gte.scheme[1][3].g_leak = [0.05]

        gte.scheme[1][4].σ = [0.98]
        gte.scheme[1][4].g_leak = [0.005]

        gte.validate_scheme()
        gte.export_gte_main()
        gte.solve(how="cycle", error=error, Niter=Niter, file_type="pkl")

    gte = GTE("ГТД для БПЛА")
    if 0:
        print(f"{Fore.CYAN}{gte.name}")
        gte.scheme = {
            1: [
                Inlet(),
                Compressor(),
                CombustionChamber(),
                Turbine(),
                Turbine(),
                Outlet(),
            ]
        }
        gte.loads = [Gear(), Propeller()]

        gte.loads[0].η = list(linspace(0.8, 1, 10 + 1))
        gte.loads[1].N = list(linspace(100 * 10**3, 500 * 10**3, 4 + 1))

        gte.m = [{1: 1}]

        gte.H = list(linspace(0, 6_000, 12 + 1))
        gte.M = list(linspace(0, 0.6, 6 + 1))
        gte.R = [0]

        gte.substance = ["AIR"]
        gte.fuel = ["KEROSENE"]
        gte.resource = [3000 * 3600]

        gte.scheme[1][0].σ = [0.98]
        gte.scheme[1][0].g_leak = [0.005]

        gte.scheme[1][1].ππ = list(linspace(1.2, 12.2, 11 + 1))
        gte.scheme[1][1].ηη = list(linspace(0.78, 0.92, 14 + 1))
        gte.scheme[1][1].g_leak = [0.05]

        gte.scheme[1][2].T_fuel = [5 + 273.15]
        gte.scheme[1][2].η_burn = [0.98]
        gte.scheme[1][2].TT3 = list(linspace(800, 1000, 4 + 1))
        gte.scheme[1][2].σ = [0.94]
        gte.scheme[1][2].T_lim = [1000]
        gte.scheme[1][2].g_leak = [0]

        gte.scheme[1][3].shafts = [[gte.scheme[1][1]], [gte.loads[0], gte.loads[1]]]
        gte.scheme[1][3].ηη = [0.9]  # list(linspace(0.78, 0.92, 14 + 1))
        gte.scheme[1][3].η_mechanical = [0.99]
        gte.scheme[1][3].T_lim = [1000]
        gte.scheme[1][3].g_leak = [0.05]

        gte.scheme[1][4].shafts = [[gte.scheme[1][1]], [gte.loads[0], gte.loads[1]]]
        gte.scheme[1][4].ηη = [0.9]  # list(linspace(0.78, 0.92, 14 + 1))
        gte.scheme[1][4].η_mechanical = [0.99]
        gte.scheme[1][4].T_lim = [1000]
        gte.scheme[1][4].g_leak = [0.05]

        gte.scheme[1][5].σ = [0.98]
        gte.scheme[1][5].g_leak = [0.005]

        gte.validate_scheme()
        gte.export_gte_main()
        gte.solve(how="cycle", error=error, Niter=Niter, file_type="pkl")

    gte = GTE("GTE MONSTER")  # экземпляр ГТД
    if 0:
        print(f"{Fore.CYAN}{gte.name}")
        gte.scheme = {
            1: [
                Inlet(),
                Compressor(),
                Compressor(),
                Compressor(),
                CombustionChamber(),
                Turbine(),
                Turbine(),
                Turbine(),
                Nozzle(),
            ],
            2: [Inlet(), Compressor(), HeatExchanger(), Nozzle()],
            3: [Inlet(), Compressor(), HeatExchanger(), Nozzle()],
        }
        gte.loads = [
            Gear(),
            Gear(),
            Load(),
            Gear(),
            Gear(),
            Load(),
            Gear(),
            Gear(),
            Load(),
            Gear(),
            Load(),
        ]

        gte.scheme[1][0].σ = [0.98]
        gte.scheme[1][0].g_leak = [0.005]

        gte.scheme[1][1].ππ = list(linspace(1.2, 12.2, 11 + 1))
        gte.scheme[1][1].ηη = list(linspace(0.78, 0.92, 14 + 1))
        gte.scheme[1][1].g_leak = [0.05]

        gte.scheme[1][2].ππ = list(linspace(1.2, 12.2, 11 + 1))
        gte.scheme[1][2].ηη = list(linspace(0.78, 0.92, 14 + 1))
        gte.scheme[1][2].g_leak = [0.05]

        gte.scheme[1][3].ππ = list(linspace(1.2, 12.2, 11 + 1))
        gte.scheme[1][3].ηη = list(linspace(0.78, 0.92, 14 + 1))
        gte.scheme[1][3].g_leak = [0.05]

        gte.scheme[1][4].T_fuel = [5 + 273.15]
        gte.scheme[1][4].η_burn = [0.98]
        gte.scheme[1][4].TT3 = list(linspace(800, 1000, 4 + 1))
        gte.scheme[1][4].σ = [0.94]
        gte.scheme[1][4].T_lim = [1000]
        gte.scheme[1][4].g_leak = [0]

        gte.scheme[1][5].shafts = [
            {"-": [gte.scheme[1][3]], "0": [gte.loads[0]]},
            {"-": [gte.loads[2]], "0": [gte.loads[1]]},
        ]
        gte.scheme[1][5].ηη = [0.9]
        gte.scheme[1][5].η_mechanical = [0.99]
        gte.scheme[1][5].T_lim = [1000]
        gte.scheme[1][5].g_leak = [0.05]

        gte.scheme[1][6].shafts = [
            {"-": [gte.scheme[1][2]], "0": [gte.loads[3]]},
            {"-": [gte.loads[5]], "0": [gte.loads[4]]},
        ]
        gte.scheme[1][6].ηη = [0.9]
        gte.scheme[1][6].η_mechanical = [0.99]
        gte.scheme[1][6].T_lim = [1000]
        gte.scheme[1][6].g_leak = [0.05]

        gte.scheme[1][7].shafts = [
            {"-": [gte.scheme[1][1]], "0": [gte.loads[6]]},
            {"-": [gte.loads[8]], "0": [gte.loads[7]]},
        ]
        gte.scheme[1][7].ηη = [0.9]
        gte.scheme[1][7].η_mechanical = [0.99]
        gte.scheme[1][7].T_lim = [1000]
        gte.scheme[1][7].g_leak = [0.05]

        gte.scheme[1][8].shafts = [{"-": [gte.loads[10]], "0": [gte.loads[9]]}]
        gte.scheme[1][8].ηη = [0.9]
        gte.scheme[1][8].η_mechanical = [0.99]
        gte.scheme[1][8].T_lim = [1000]
        gte.scheme[1][8].g_leak = [0.05]

        gte.scheme[1][9].ππ = [1.13]
        gte.scheme[1][9].ηη = [0.92]
        gte.scheme[1][9].v_ = [0.99]
        gte.scheme[1][9].g_leak = [0.001]

        gte.scheme[2][0].σ = [0.98]
        gte.scheme[2][0].g_leak = [0.005]

        gte.scheme[2][1].ππ = list(linspace(1.2, 12.2, 11 + 1))
        gte.scheme[2][1].ηη = list(linspace(0.78, 0.92, 14 + 1))
        gte.scheme[2][1].g_leak = [0.05]

        gte.scheme[2][2].ππ = list(linspace(1.2, 12.2, 11 + 1))
        gte.scheme[2][2].ηη = list(linspace(0.78, 0.92, 14 + 1))
        gte.scheme[2][2].g_leak = [0.05]

        gte.scheme[2][3].ππ = list(linspace(1.2, 12.2, 11 + 1))
        gte.scheme[2][3].ηη = list(linspace(0.78, 0.92, 14 + 1))
        gte.scheme[2][3].g_leak = [0.05]

        gte.scheme[2][4].ππ = [1.13]
        gte.scheme[2][4].ηη = [0.92]
        gte.scheme[2][4].v_ = [0.99]
        gte.scheme[2][4].g_leak = [0.001]

        gte.validate_scheme()
        gte.export_gte_main()
        gte.solve(error=error, Niter=Niter, file_type="xlsx")

    if False:
        GTE().version()  # версия программы
        gte = GTE(gte_name)  # экземпляр ГТД

        gte.input_scheme()  # ввод схемы ГТД
        gte.input_contouring()  # ввод контурностей контуров
        gte.input_GTE_parameters()  # ввод характеристик ГТД
        gte.input_node_parameters()  # ввод характеристик узлов
        gte.export_gte_main()
        gte.variation(
            how=how,
        )
        # gte.draw_longitudinal_section_parameters()
        # gte.draw_cycle('T(S)')
        print("The End.")
