from typing import Any, Dict, List, Tuple, Union

import matplotlib.pyplot as plt
from numpy import array, linspace
from scipy.optimize import root
from substance import Substance

try:
    from .config import parameters as gtep
    from .errors import TYPE_ERROR
    from .nodes.combustion_chamber.combustion_chamber import CombustionChamber
    from .nodes.compressor.compressor import Compressor
    from .nodes.node import GTENode
    from .nodes.nozzle.nozzle import Nozzle
    from .nodes.turbine.turbine import Turbine
except ImportError:
    import os
    import sys

    sys.path.insert(0, os.getcwd())

    from gte.config import parameters as gtep
    from gte.errors import TYPE_ERROR
    from gte.nodes.combustion_chamber.combustion_chamber import CombustionChamber
    from gte.nodes.compressor.compressor import Compressor
    from gte.nodes.node import GTENode
    from gte.nodes.nozzle.nozzle import Nozzle
    from gte.nodes.turbine.turbine import Turbine


class GTE:
    """ГТД"""

    NODES = (Compressor, CombustionChamber, Turbine, Nozzle)

    __slots__ = ("name", "_GTE__scheme", "shafts", "_GTE__finder", "requirements")

    def __init__(self, scheme: Tuple, name: str = "GTE") -> None:
        """Инициализация объекта ГТД"""
        self.name: str = name

        if not isinstance(scheme, (tuple, list)):
            raise TypeError(f"{type(scheme)=} must be tuple")
        for contour, nodes in enumerate(scheme):
            if not isinstance(nodes, (tuple, list)):
                raise TypeError(f"{type(nodes)=} must be tuple")
            for node in nodes:
                if not isinstance(node, self.NODES):
                    raise TypeError(f"{type(node)=} must be in {self.NODES}")
            scheme[contour] = tuple(nodes)
        self.__scheme: Tuple = tuple(scheme)

        self.shafts = []
        self.requirements = []

    def __repr__(self) -> str:
        """Описание ГТД"""
        result = f"""{self.__class__.__name__}={self.name}\nis_solvable={self.check_solvable}\nscheme:\n"""
        for contour, nodes in enumerate(self.__scheme):
            result += f"\t{contour=}: {[node.__class__.__name__ for node in nodes]}\n"
        for contour, nodes in enumerate(self.__scheme):
            for place, node in enumerate(nodes):
                result += f"\t\tnode [{contour}][{place}]: {node}\n"
                for key, value in node.parameters.items():
                    result += f"\t\t\t{key:<25}: {value:.4f}\n"

        return result

    def __setattr__(self, name, value):
        if name == "name":
            if not isinstance(value, str):
                raise TypeError(TYPE_ERROR.format(f"type(name)={type(value)}", str))
        super().__setattr__(name, value)

    def __getitem__(self, key) -> Tuple:
        """Возврат контура"""
        return self.__scheme[key]

    @property
    def scheme(self) -> Tuple[Tuple[GTENode, ...]]:
        """Схема ГТД"""
        return self.__scheme

    def add_shaft(self, *node_places) -> None:
        """Добавление вала как связи по балансу мощностей"""
        shaft: List[Union[Compressor, Turbine]] = []
        if len(node_places) == 0:
            raise ValueError("empty shaft")
        for node_place in node_places:
            if not isinstance(node_place, (list, tuple)):
                raise TypeError(f"{type(node_place)=} must be {tuple}")
            if len(node_place) != 2:
                raise ValueError(f"{len(node_place)=} must be 2")
            contour, place = node_place
            if not (0 <= contour <= len(self.__scheme) - 1):
                raise ValueError(f"{contour=} does not exist")
            if not (0 <= place <= len(self.__scheme[contour]) - 1):
                raise ValueError(f"{place=} does not exist in {contour=}")
            node = self.__scheme[contour][place]
            if not isinstance(node, (Compressor, Turbine)):
                raise ValueError(f"{type(node_place)=} must be in {Compressor, Turbine}")
            shaft.append(node)
        self.shafts.append(shaft)

        self.__finder: Dict = {}
        for shaft in self.shafts:
            for node1 in shaft:
                for contour, nodes in enumerate(self.__scheme):
                    for place, node2 in enumerate(nodes):
                        if node1 is node2:
                            self.__finder[node1] = (contour, place)

    def add_requirement(self, parameter: str, value: float, contour: int = -1, place: int = -1) -> None:
        """Добавление требований к ГТД"""
        self.requirements.append({"parameter": parameter, "value": value, "contour": contour, "place": place})

    def plot(self, **kwargs) -> plt.Figure:
        """Визуализация схемы ГТД"""
        fg = plt.figure(figsize=kwargs.get("figsize", (max(map(len, self.__scheme)) * 2, (len(self.__scheme) + 1 + 2) * 2)))
        fg.suptitle("GTE scheme", fontsize=14, fontweight="bold")
        gs = fg.add_gridspec(len(self.__scheme) + 1, 1)  # строки, столбцы

        for contour, nodes in enumerate(self.__scheme):
            fg.add_subplot(gs[len(self.__scheme) - 1 - contour, 0])
            plt.grid(True)
            plt.axis("square")
            plt.title(f"contour {contour}", fontsize=14)
            plt.xlim(0, len(self.__scheme[contour]))
            plt.ylim(0, 1)
            plt.xticks(linspace(0, len(self.__scheme[contour]), len(self.__scheme[contour]) + 1))
            plt.yticks(linspace(0, 1, 1 + 1))

            x0 = y0 = 0.5  # center

            for i, node in enumerate(nodes):
                x, y = array(node.figure, dtype="float64")
                plt.plot(x + x0, y + y0, color="black", linewidth=3, label=f"{contour}.{i}: {node.__class__.__name__}")
                plt.text(x0, y0, f"{contour}.{i}", fontsize=12, fontweight="bold", ha="center", va="center")
                x0 += 1

        fg.add_subplot(gs[len(self.__scheme), 0])
        plt.axis("off")
        plt.grid(False)
        plt.xlim(0, max(map(len, self.__scheme)))
        plt.ylim(0, 1)
        plt.plot([0, max(map(len, self.__scheme))], [0.5, 0.5], color="black", linewidth=1.5, linestyle="dashdot")

        fg.legend(
            title="Specification",
            title_fontsize=14,
            alignment="center",
            loc="lower center",
            fontsize=12,
            ncols=len(self.__scheme),
            frameon=True,
            framealpha=1.0,
            facecolor="white",
            edgecolor="black",
            draggable=True,
        )

        return fg

    @property
    def check_solvable(self) -> Tuple[bool, str]:
        x, y = 0, len(self.shafts) + len(self.requirements)
        for nodes in self.__scheme:
            for node in nodes:
                x += node.n_vars - len(node.parameters)  # требуемое - имеющееся = необходимое
        if x < y:
            return False, f"need to add {y - x} equations"
        elif x > y:
            return False, f"need to delete {x - y} equations"
        else:
            return True, ""

    def predict(self, inlet: Substance, use_ml: bool = True) -> Tuple[Dict[str, Union[int, float]]]:
        """Начальные прибличения"""
        prediction = []
        for contour, nodes in enumerate(self.__scheme):
            for place, node in enumerate(nodes):
                if node.is_solvable[0]:
                    continue
                want = node.n_vars - len(node.parameters)  # необходимое количетсво предсказаний
                for parameter in node.variables:
                    if parameter in node.parameters:
                        continue
                    while want:  # TODO
                        if parameter == gtep.pipi:
                            prediction.append({"contour": contour, "place": place, "parameter": parameter, "value": 1 / 3})
                        elif parameter == gtep.force:
                            prediction.append({"contour": contour, "place": place, "parameter": parameter, "value": 10_000})
                        elif parameter == gtep.eff_speed:
                            prediction.append({"contour": contour, "place": place, "parameter": parameter, "value": 0.98})
                        want -= 1

        return tuple(prediction)

    def calculate(self, inlet: Substance, fuel: Substance = None, verbose: bool = False) -> Tuple[Tuple[Substance, ...]]:
        """Поузловой термодинамический расчет двигателя 'в строчку'"""
        substances: List = []
        variables: List = []
        for contour, nodes in enumerate(self.__scheme):
            if verbose:
                print(f"{contour = }")

            outlet = inlet  # выход из предыдущего узла
            ss, v = [outlet], []
            for i, node in enumerate(nodes):
                if verbose:
                    print(f"\t{i}: {node = }")

                if isinstance(node, (Compressor, Turbine)):
                    var, outlet = node.calculate(outlet, node.parameters)
                elif isinstance(node, CombustionChamber):
                    var, outlet = node.calculate(outlet, fuel, node.parameters)
                elif isinstance(node, Nozzle):
                    var, outlet = node.calculate(outlet, node.parameters)
                else:
                    raise ValueError(f"no found {node=}")

                if verbose:
                    for key, value in outlet.parameters.items():
                        print(f"\t\t{key:<25}: {value:.4f}")

                ss.append(outlet)
                v.append(var)

            substances.append(tuple(ss))
            variables.append(tuple(v))

        return tuple(substances), tuple(variables)

    def _equations(self, x0: Tuple[float], args: Dict[str, Any]) -> Tuple[float, ...]:
        """sum(Compressor.power) = sum(Turbine.power)"""

        inlet, fuel = args["inlet"], args.get("fuel")
        prediction = args["prediction"]
        verbose = args.get("verbose", False)

        for i, predict in enumerate(prediction):
            contour, place, parameter = predict["contour"], predict["place"], predict["parameter"]
            self.__scheme[contour][place].parameters[parameter] = float(x0[i])

        _, _ = self.calculate(inlet, fuel, verbose=verbose)

        power_balances: List[float] = []
        # Расчет баланса мощностей по каждому валу
        for shaft in self.shafts:
            power_balance = 0.0
            for node in shaft:
                contour, place = self.__finder[node]
                power_balance += node.parameters.get(gtep.power, 0.0)
            power_balances.append(power_balance)

        return tuple(power_balances)

    def solve(self, inlet: Substance, fuel: Substance = None, verbose: bool = False):
        """Термодинамический расчет ГТД"""
        if not self.check_solvable:
            raise ArithmeticError

        predictions = self.predict(inlet)
        x0 = tuple(prediction["value"] for prediction in predictions)

        solution = root(self._equations, x0, {"inlet": inlet, "fuel": fuel, "prediction": predictions, "verbose": verbose}, method="lm")

        result: Dict = {}
        for i, prediction in enumerate(predictions):
            contour, place, parameter, value = prediction["contour"], prediction["place"], prediction["parameter"], solution.x[i]
            result[parameter] = float(value)
            del self[contour][place].parameters[parameter]  # возврат к инициализированному состоянию

        return result

    # TODO
    def to_dict(self) -> Dict:
        return {
            "name": self.name,
            "scheme": ";".join(self.__scheme for nodes in self.__scheme),
        }


if __name__ == "__main__":
    from itertools import product

    from fixtures import air as inlet
    from fixtures import kerosene as fuel
    from numpy import arange
    from tqdm import tqdm

    if False:
        gte = GTE(
            [
                (
                    Compressor({gtep.effeff: 0.85, gtep.pipi: 6}, name="HPC"),
                    CombustionChamber({gtep.eff_burn: 0.99, gtep.pipi: 0.95}, name="CC"),
                    Turbine({gtep.effeff: 0.9, gtep.pipi: 1 / 3}, name="HPT"),
                ),
            ],
            name="Simple",
        )
        test_cases = [{}]

    if True:
        inlet.parameters[gtep.m] = 80

        gte = GTE(
            [
                (
                    Compressor({gtep.effeff: 0.85, gtep.pipi: 6}, name="HPC"),
                    CombustionChamber({gtep.eff_burn: 0.99, gtep.pipi: 0.95}, name="CC"),
                    Turbine({gtep.effeff: 0.9}, name="HPT"),
                ),
            ],
            name="AI-9",
        )
        gte.add_shaft([0, 0], [0, 2])

        test_cases = []
        for c_pipi, c_effeff in product(arange(0.85, 0.95, 0.01), arange(3, 6, 1)):
            for eff_burn, cc_pipi in product(arange(0.95, 0.99, 0.01), arange(0.93, 0.97, 0.01)):
                for t_effeff in arange(0.88, 0.94, 0.01):
                    test_cases.append(
                        {
                            (0, 0): {gtep.pipi: float(c_pipi), gtep.effeff: float(c_effeff)},
                            (0, 1): {gtep.eff_burn: float(eff_burn), gtep.pipi: float(cc_pipi)},
                            (0, 2): {gtep.effeff: float(t_effeff)},
                        }
                    )

    if False:
        inlet.parameters[gtep.m] = 50

        hpc = Compressor({gtep.effeff: 0.85, gtep.pipi: 6}, name="HPC")
        cc = CombustionChamber({gtep.eff_burn: 0.99, gtep.pipi: 0.95}, name="CC")
        hpt = Turbine({gtep.effeff: 0.9}, name="HPT")
        n = Nozzle({gtep.pipi: 1 / 1.8, gtep.eff_speed: 0.98}, name="N")

        gte = GTE([(hpc, cc, hpt, n)], name="Jumo 004b")
        gte.add_shaft([0, 0], [0, 2])
        gte.add_equation(0, 3, gtep.force, 30_000)

        test_cases = []
        for c_pipi, c_effeff in product(arange(0.85, 0.95, 0.01), arange(3, 6, 1)):
            for eff_burn, cc_pipi in product(arange(0.95, 0.99, 0.01), arange(0.93, 0.97, 0.01)):
                for t_effeff in arange(0.88, 0.94, 0.01):
                    for n_pipi in arange(1 / 1.8, 1 / 1.2, 0.01):
                        test_cases.append(
                            {
                                (0, 0): {gtep.pipi: float(c_pipi), gtep.effeff: float(c_effeff)},
                                (0, 1): {gtep.eff_burn: float(eff_burn), gtep.pipi: float(cc_pipi)},
                                (0, 2): {gtep.effeff: float(t_effeff)},
                                (0, 3): {gtep.pipi: float(n_pipi)},
                            }
                        )

    for s in (inlet, fuel):
        if not s:
            continue
        print(s.name)
        for key, value in s.parameters.items():
            print(f"{key:<25}: {value}")
        print()

    gte.plot()
    # plt.show()

    for test_case in tqdm(test_cases, desc="Solving"):
        for k, v in test_case.items():
            gte[k[0]][k[1]].parameters = v

        result = gte.solve(inlet, fuel, verbose=False)

        # print(gte)

        # print(gte.to_dict())
