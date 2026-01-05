from typing import Any, Dict, List, Tuple

import matplotlib.pyplot as plt
from numpy import array, cos, isnan, linspace, nan, prod, radians, sin
from scipy.optimize import root
from substance import Substance
from tqdm import tqdm

try:
    from .combustion_chamber import CombustionChamber
    from .compressor import Compressor
    from .config import parameters as gtep
    from .nozzle import Nozzle
    from .turbine import Turbine
except ImportError:
    from combustion_chamber import CombustionChamber
    from compressor import Compressor
    from config import parameters as gtep
    from nozzle import Nozzle
    from turbine import Turbine


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


'''
class GTE_OLD(Variability):
    """ГТД"""

    @classmethod
    @property
    def version(cls) -> str:
        version = 8.0
        next_version = (
            "камера смешения",
            "переходный канал",
            "охлаждение турбины",
            "продолжение расчета",
            "multiprocessing",
            "ускорение расчета до 6000 [ГТД/с]",
        )
        print(f"{cls.__name__} version: {version}")
        for i, v in enumerate(next_version):
            print(cls.__name__ + " version:", int(version) + i + 1, v)
        return str(version)

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

            vars_list = root(
                self.equations,
                tuple(vars0.values()),
                xtol=xtol,
                maxfev=100 * (len(vars0) + 1),
            )
            vars = {key: vars_list[i] for i, key in enumerate(vars0.keys())}

            if log:
                print(f"variables: {vars}")
                print(f"zeros: {self.equations(list(vars.values()))}")

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
            print("Решение не найдено!")
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
'''


class GTE:
    """ГТД"""

    __slots__ = ("name", "scheme", "shafts")

    def __init__(self, name: str, scheme: Dict, shafts: Tuple) -> None:
        assert isinstance(name, str), TypeError(f"{type(name)=} must be str")
        self.name: str = name

        all_nodes = []
        assert isinstance(scheme, dict), TypeError(f"{type(scheme)=} must be dict")
        for contour, nodes in scheme.items():
            assert isinstance(contour, int), TypeError(f"{type(contour)=} must be int")
            assert isinstance(nodes, (tuple, list)), TypeError(f"{type(nodes)=} must be tuple")
            all_nodes.extend(nodes)
        self.scheme: Dict = dict(sorted(scheme.items(), key=lambda item: item[0]))  # сортировка по контурам по возрастанию

        assert isinstance(shafts, (tuple, list)), TypeError(f"{type(shafts)=} must be tuple")
        for shaft in shafts:
            for node in shaft:
                assert isinstance(node, (Compressor, Turbine)), TypeError(f"{type(node)=} must be in {Compressor, Turbine}")
                assert node in all_nodes

    def show(self, **kwargs) -> None:
        """Визуализация схемы ГТД"""
        fg = plt.figure(figsize=kwargs.get("figsize", (max(map(len, self.scheme.values())) * 2, (len(self.scheme) + 1 + 2) * 2)))
        fg.suptitle("GTE scheme", fontsize=14, fontweight="bold")
        gs = fg.add_gridspec(len(self.scheme) + 1, 1)  # строки, столбцы

        for contour in self.scheme:
            fg.add_subplot(gs[len(self.scheme) - contour, 0])
            plt.grid(True)
            plt.axis("square")
            plt.title(f"contour {contour}", fontsize=14)
            plt.xlim(0, len(self.scheme[contour]))
            plt.ylim(0, 1)
            plt.xticks(linspace(0, len(self.scheme[contour]), len(self.scheme[contour]) + 1))
            plt.yticks(linspace(0, 1, 1 + 1))

            x0 = y0 = 0.5  # center

            for i, node in enumerate(self.scheme[contour]):
                x, y = array(node.figure, dtype="float64")
                plt.plot(x + x0, y + y0, color="black", linewidth=3, label=f"{contour}.{i + 1}: {node.__class__.__name__}")
                plt.text(x0, y0, f"{contour}.{i + 1}", fontsize=12, fontweight="bold", ha="center", va="center")
                x0 += 1

        fg.add_subplot(gs[len(self.scheme), 0])
        plt.axis("off")
        plt.grid(False)
        plt.xlim(0, max(map(len, self.scheme.values())))
        plt.ylim(0, 1)
        plt.plot([0, max(map(len, self.scheme.values()))], [0.5, 0.5], color="black", linewidth=1.5, linestyle="dashdot")

        fg.legend(
            title="Specification",
            title_fontsize=14,
            alignment="center",
            loc="lower center",
            fontsize=12,
            ncols=len(self.scheme),
            frameon=True,
            framealpha=1.0,
            facecolor="white",
            edgecolor="black",
            draggable=True,
        )

        plt.show()

    def _equations(self, x: Tuple[float], args: Dict[str, Any]) -> List[float]:
        """sum(Compressor.power) = sum(Turbine.power)"""
        result = []
        for shaft in self.shafts:
            balance_power: float = 0
            for node in shaft:
                if isinstance(node, Turbine):
                    balance_power += node.power
                elif isinstance(node, Compressor):
                    balance_power -= node.power
            result.append(balance_power)
        return result

    def calculate(self, inlet: Substance, fuel: Substance, log: bool = False) -> Substance:
        """Поузловой термодинамический расчет двигателя 'в строчку'"""
        if log:
            for s in (inlet, fuel):
                print(s.name)
                for key, value in s.parameters.items():
                    print(f"{key:<25}: {value}")
                print()

        for contour in self.scheme:
            if log:
                print(f"{contour = }")

            outlet = inlet  # выход из предыдущего узла
            for i, node in enumerate(self.scheme[contour]):
                if log:
                    print(f"\t{i}: {node = }")

                if isinstance(node, Compressor):
                    outlet = node.solve(outlet, 5)["outlet"]
                elif isinstance(node, CombustionChamber):
                    outlet = node.solve(outlet, fuel)["outlet"]
                elif isinstance(node, Turbine):
                    outlet = node.solve(outlet, 5)["outlet"]
                    """elif isinstance(node, outlet):
                    outlet = node.solve(args)["outlet"]"""
                else:
                    raise Exception

                if log:
                    for key, value in outlet.parameters.items():
                        print(f"\t\t{key:<25}: {value:.4f}")

        return outlet

    def solve(self, inlet: Substance, fuel: Substance, log: bool = False) -> Dict[str, Any]:
        """Термодинамический расчет ГТД"""
        while True:
            self.calculate(inlet, fuel, log=log)


if __name__ == "__main__":
    from fixtures import air, kerosene

    c = Compressor(
        "HPC",
        characteristic={
            gtep.effeff: lambda rotation_frequency, mass_flow: 0.85,
            gtep.pipi: lambda rotation_frequency, mass_flow: 6,
        },
    )

    cc = CombustionChamber(
        "CC",
        characteristic={
            gtep.eff_burn: lambda mass_flow: 0.99,
            gtep.p_eff: lambda mass_flow: 0.95,
        },
    )

    t = Turbine(
        "HPT",
        characteristic={
            gtep.effeff: lambda rotation_frequency, mass_flow: 0.9,
            gtep.pipi: lambda rotation_frequency, mass_flow: 3,
        },
    )

    scheme = {1: (c, cc, t)}
    shafts = [(scheme[1][0], scheme[1][2])]

    gte = GTE("Jumo 004b", scheme, shafts)

    # gte.show()

    gte.calculate(air, kerosene, log=True)

    exit()

    scheme = {
        1: [Compressor(), CombustionChamber(), Turbine(), Nozzle()],
        2: [Compressor(), Nozzle()],
    }
    shafts = {1: [gte.scheme[1][1], gte.scheme[2][1], gte.scheme[1][3]]}

    gte = GTE("CFM-56", scheme, shafts)

    gte.show()

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
