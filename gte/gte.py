from collections import deque
from typing import Any, Dict, List, Set, Tuple

import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.patches import Patch
from numpy import isnan
from scipy.optimize import root
from substance import Substance

try:
    from .config import EPSREL
    from .config import parameters as gtep
    from .errors import TYPE_ERROR
    from .nodes.burner.burner import Burner
    from .nodes.channel.channel import Channel
    from .nodes.joiner.joiner import Joiner
    from .nodes.node import GTENode
    from .nodes.nozzle.nozzle import Nozzle
    from .nodes.splitter.splitter import Splitter
    from .nodes.turbocompressor.rotor.rotor import Rotor
except ImportError:
    import os
    import sys

    sys.path.insert(0, os.getcwd())

    from gte.config import EPSREL
    from gte.config import parameters as gtep
    from gte.errors import TYPE_ERROR
    from gte.nodes.burner.burner import Burner
    from gte.nodes.channel.channel import Channel
    from gte.nodes.joiner.joiner import Joiner
    from gte.nodes.node import GTENode
    from gte.nodes.nozzle.nozzle import Nozzle
    from gte.nodes.splitter.splitter import Splitter
    from gte.nodes.turbocompressor.rotor.rotor import Rotor


class GTE:
    """Ориентированный ациклический граф потока рабочего тела"""

    NODES: Tuple[GTENode] = (Rotor, Burner, Channel, Nozzle, Splitter, Joiner)

    # Цвета для разных типов узлов
    NODES_COLORS = {
        "Rotor": "lightblue",
        "Burner": "red",
        "Channel": "lightgreen",
        "Nozzle": "gray",
        "Splitter": "pink",
        "Joiner": "brown",
    }

    __slots__ = (
        "name",  # имя
        "_GTE__nodes",  # все узлы
        "_GTE__in",  # входящие связи в узел
        "_GTE__out",  # исходящие связи из узла
        "_out_index",  # (from, to) -> индекс выхода
        "_splitter_counter",  # для Splitter
        "_GTE__shafts",  # валы (механическая связь)
        "requirements",  # требования
    )

    def __init__(self, name: str):
        self.name: str = name

        self.__nodes: Set[GTENode] = set()
        self.__in: Dict[GTENode, List[GTENode]] = {}
        self.__out: Dict[GTENode, List[GTENode]] = {}

        self._out_index: Dict[Tuple[GTENode, GTENode], int] = {}
        self._splitter_counter: Dict[GTENode, int] = {}

        self.__shafts: List[Tuple[GTENode]] = []

        self.requirements: List[Dict] = []

    def __repr__(self) -> str:
        """Описание ГТД"""
        return f"{self.__class__.__name__} (name={self.name}, nodes={len(self.__nodes)}, edges={sum(len(v) for v in self.__out.values())}, shafts={len(self.__shafts)}, requirements={len(self.requirements)})"

    def __setattr__(self, name, value):
        if name == "name":
            if not isinstance(value, str):
                raise TypeError(TYPE_ERROR.format(f"type(name)={type(value)}", str))
        super().__setattr__(name, value)

    def add_node(self, node: GTENode) -> None:
        """Добавление узла"""
        if not isinstance(node, GTENode):
            raise TypeError(TYPE_ERROR.format(f"{type(node)=}", GTENode))

        if node not in self.__nodes:
            self.__nodes.add(node)
            self.__in[node], self.__out[node] = [], []

    def add_edge(self, from_node: GTENode, to_node: GTENode, outlet_index: int = None) -> None:
        """Добавление связи между узлами from_node и to_node по выходу outlet_index из узла from_node"""
        if not isinstance(from_node, GTENode):
            raise TypeError(TYPE_ERROR.format(f"{type(from_node)=}", GTENode))
        if not isinstance(to_node, GTENode):
            raise TypeError(TYPE_ERROR.format(f"{type(to_node)=}", GTENode))

        if from_node not in self.__nodes:
            self.add_node(from_node)
        if to_node not in self.__nodes:
            self.add_node(to_node)

        # Для Splitter назначаем индекс выхода
        if isinstance(from_node, Splitter):
            if outlet_index is None:  # если индекс не назначен явным образом
                if from_node not in self._splitter_counter:  # если узла from_node нет в _splitter_counter
                    self._splitter_counter[from_node] = 0  # по дефолту 0 индекс выхода
                outlet_index = self._splitter_counter[from_node]
                self._splitter_counter[from_node] += 1
            self._out_index[(from_node, to_node)] = outlet_index

        self.__in[to_node].append(from_node)
        self.__out[from_node].append(to_node)

    def add_shaft(self, *nodes: GTENode) -> None:
        """Добавление вала как механической связи по балансу мощностей"""
        if len(nodes) < 2:
            raise ValueError(f"{len(nodes)=} must be >= 2")
        for node in nodes:
            if not isinstance(node, GTENode):
                raise TypeError(TYPE_ERROR.format(f"{type(node)=}", GTENode))
            if node not in self.__nodes:
                self.add_node(node)  # автоматически добавляем узел в граф
        self.__shafts.append(tuple(nodes))

    # TODO
    def add_requirement(self, parameter: str, value: float, place: int = -1) -> None:
        """Добавление требований к ГТД"""
        self.requirements.append({"parameter": parameter, "value": value, "place": place})

    def predecessors(self, node: GTENode) -> List[GTENode]:
        """Предшественники"""
        return self.__in.get(node, [])

    def successors(self, node: GTENode) -> List[GTENode]:
        """Преемники"""
        return self.__out.get(node, [])

    def scheme(self) -> Dict[str, Any]:
        """
        Возвращает три представления графа потоков ГТД:
        - adjacency_matrix: список списков (0/1) размером N x N,
        - adjacency_list: словарь {узел: [преемники]},
        - edge_list: список кортежей (from_node, to_node).
        Дополнительно возвращает ordered_nodes (список узлов в порядке, соответствующем матрице).
        """
        # Упорядочиваем узлы для фиксированного порядка в матрице
        ordered_nodes = sorted(self.__nodes, key=lambda n: n.name)

        # Сопоставляем узел -> индекс
        node_to_index = {node: i for i, node in enumerate(ordered_nodes)}
        n = len(ordered_nodes)

        # Инициализация матрицы смежности нулями
        adj_matrix = [[0] * n for _ in range(n)]

        # Построение списка смежности и списка рёбер
        adj_list = {node: [] for node in ordered_nodes}
        edge_list = []

        for from_node in self.__nodes:
            for to_node in self.__out.get(from_node, []):
                # матрица смежности
                i, j = node_to_index[from_node], node_to_index[to_node]
                adj_matrix[i][j] = 1

                # список смежности
                adj_list[from_node].append(to_node)

                # список рёбер
                edge_list.append((from_node, to_node))

        return {
            "adjacency_matrix": adj_matrix,
            "adjacency_list": adj_list,
            "edge_list": edge_list,
            "ordered_nodes": ordered_nodes,
        }

    @property
    def nodes(self) -> Set[GTENode]:
        return self.__nodes

    @property
    def shafts(self) -> List[Tuple[GTENode]]:
        return self.__shafts

    @property
    def order(self) -> List[GTENode]:
        """Топологическая сортировка (Kahn)"""
        in_degree = {node: len(self.predecessors(node)) for node in self.__nodes}
        queue = deque(n for n in self.__nodes if in_degree[n] == 0)
        order = []
        while queue:
            node = queue.popleft()
            order.append(node)
            for successor in self.successors(node):
                in_degree[successor] -= 1
                if in_degree[successor] == 0:
                    queue.append(successor)
        if len(order) != len(self.__nodes):
            raise ValueError("Cycles in the flow graph are not supported")
        return order

    @property
    def is_solvable(self) -> Tuple[bool, str]:
        # количество неизвестных
        x = sum(node.n_vars - len(node.parameters) for node in self.__nodes)
        # количество уравнений связи и требований
        y = len(self.__shafts) + len(self.requirements)

        dif = abs(x - y)  # дисбаланс неизвестных и уравнений
        if x < y:
            return False, f"need to add {dif} equations or variables"
        elif x > y:
            return False, f"need to delete {dif} equations or variables"
        else:
            return True, ""

    def plot(self) -> plt.Figure:
        """Визуализация графа потоков ГТД"""
        g = nx.DiGraph()

        # Добавляем узлы
        for node in self.__nodes:
            node_class = node.__class__.__name__
            g.add_node(node, label=f"{node_class}\n({node.name})", type=node_class)

        # Добавляем рёбра потока
        for node in self.__nodes:
            for successor in self.successors(node):
                g.add_edge(node, successor)

        # Вычисляем глубину каждого узла (расстояние от entry_points)
        depth = {node: 0 for node in self.__nodes}
        # Топологический порядок гарантирует, что предшественники обработаны раньше
        for node in self.order:
            for predecessor in self.predecessors(node):
                # Глубина узла = максимальная глубина предшественников + 1
                depth[node] = max(depth[node], depth[predecessor] + 1)

        # Группируем узлы по глубине
        nodes_by_depth = {}
        for node, d in depth.items():
            nodes_by_depth.setdefault(d, []).append(node)
        # Позиции: x = глубина, y = индекс в своей группе
        pos = {}
        for d, nodes in nodes_by_depth.items():
            for i, node in enumerate(nodes):
                # Можно развернуть вертикально: i от 0 до len(nodes)-1
                pos[node] = (d, i - (len(nodes) - 1) / 2)  # центрируем

        fg = plt.figure()

        # Цвета для узлов
        node_colors = [self.NODES_COLORS.get(g.nodes[n]["type"], "white") for n in g.nodes]

        # Рисуем узлы
        nx.draw_networkx_nodes(g, pos, node_color=node_colors, node_size=4_000, edgecolors="black")

        # Рисуем рёбра потока (стрелки для направленного графа)
        nx.draw_networkx_edges(g, pos, arrowstyle="-|>", arrowsize=20, edge_color="black", width=1.5)

        # Рисуем механические связи (валы) пунктирной линией без стрелок
        for shaft in self.__shafts:
            for i in range(len(shaft) - 1):
                u, v = shaft[i], shaft[i + 1]
                if u in pos and v in pos:
                    nx.draw_networkx_edges(
                        g,
                        pos,
                        edgelist=[(u, v)],
                        style="dashed",
                        edge_color="blue",
                        width=1.5,
                        arrows=False,
                        arrowstyle="-",
                    )

        # Рисуем метки
        labels = {n: g.nodes[n]["label"] for n in g.nodes}
        nx.draw_networkx_labels(g, pos, labels, font_size=8, font_weight="bold")

        # Легенда

        legend_elements = [Patch(facecolor=color, label=node_type) for node_type, color in self.NODES_COLORS.items()]
        plt.legend(handles=legend_elements, loc="upper center", bbox_to_anchor=(0.5, 0), title="Nodes")

        plt.title(f"{self.__class__.__name__} '{self.name}' scheme", loc="center")
        plt.axis("off")
        plt.tight_layout()

        return fg

    # TODO
    def generator(self, variables: Dict[Tuple[int, int], Dict]):
        """Генератор решаемых ГТД"""
        if not isinstance(variables, dict):
            raise TypeError(TYPE_ERROR.format(f"{variables=}", dict))

        i = 0
        yield GTE(self.__scheme, f"{i}")

    # TODO
    def predict(self, inlet: Substance, use_ml: bool = True) -> Tuple[Dict[GTENode, Dict[str, float]], Dict]:
        """Начальные прибличения"""
        prediction: Dict[GTENode, Dict[str, float]] = {}
        for node in self.__nodes:
            if node.is_solvable[0]:
                continue
            want = node.n_vars - len(node.parameters)  # необходимое количетсво предсказаний
            prediction[node] = {}
            for parameter in node.variables:
                if parameter in node.parameters:
                    continue
                while want:  # TODO
                    if parameter == gtep.pipi:
                        prediction[node][parameter] = 1 / 3
                    elif parameter == gtep.effeff:
                        prediction[node][parameter] = 1
                    elif parameter == gtep.power:
                        prediction[node][parameter] = 0
                    elif parameter == gtep.force:
                        prediction[node][parameter] = 10_000
                    elif parameter == gtep.eff_speed:
                        prediction[node][parameter] = 1
                    want -= 1

        return prediction, {}  # для requirements

    def _equations(self, x0: Tuple[float], args: Dict[str, Any]) -> Tuple[float, ...]:
        """sum(Compressor.power) = sum(Turbine.power)"""

        inlet, fuel = args["inlet"], args.get("fuel")
        prediction = args["prediction"]
        verbose = args.get("verbose", False)

        i = 0
        for node, prediction in prediction.items():
            for parameter in prediction:
                node.parameters[parameter] = float(x0[i])
                i += 1

        vars, substances = self.calculate(inlet, fuel, verbose=verbose)

        power_balances: List[float] = [0] * len(self.__shafts)
        # Расчет баланса мощностей по каждому валу
        for i, shaft in enumerate(self.__shafts):
            for node in shaft:
                power_balances[i] += vars[node][gtep.power]

        return tuple(power_balances)

    def calculate(self, inlet: Substance, fuel: Substance = None, verbose: bool = False) -> Tuple[Dict[GTENode, Dict[str, float]], Dict[GTENode, Tuple]]:
        """Расчет двигателя 'в строчку'"""
        vars, substances = {}, {}  # Кэш: узел -> (входное вещество/вещества, выходное вещество/вещества)

        def get_inputs(node: GTENode):
            """Собирает выходы всех предшественников"""
            preds = self.predecessors(node)
            if not preds:
                return inlet  # начальный узел
            if len(preds) == 1:
                return substances[preds[0]][1]  # один вход – одно вещество
            # Несколько предшественников – Joiner ожидает кортеж веществ
            return tuple(substances[p][1] for p in preds)

        for node in self.order:
            if verbose:
                print(f"{node}")
            inputs = get_inputs(node)
            if isinstance(node, Splitter):  # Splitter имеет один вход и несколько выходов
                v, outlets = node.calculate(node.parameters, inputs)
                vars[node] = v
                substances[node] = (inputs, *outlets)
            elif isinstance(node, Joiner):  # Joiner принимает несколько входов
                v, outlet = node.calculate(node.parameters, *inputs)
                vars[node] = v
                substances[node] = (inputs, outlet)
            elif isinstance(node, Burner):
                v, outlet = node.calculate(node.parameters, inputs, fuel)
                vars[node] = v
                substances[node] = (inputs, outlet)
            else:  # Rotor, Channel, Nozzle
                v, outlet = node.calculate(node.parameters, inputs)
                vars[node] = v
                substances[node] = (inputs, outlet)

        return vars, substances

    def solve(self, inlet: Substance, fuel: Substance = None, prediction: Dict[GTENode, Dict[str, float]] = None, verbose: bool = False) -> bool:
        """Термодинамический расчет ГТД"""
        is_solvable, reason = self.is_solvable
        if not is_solvable:
            raise ArithmeticError(reason)

        predictions, _ = self.predict(inlet, use_ml=True)
        x0 = tuple(parameter for parameters in predictions.values() for parameter in parameters.values())

        result = root(self._equations, x0, {"inlet": inlet, "fuel": fuel, "prediction": predictions, "verbose": verbose}, method="lm")

        return result.success

    def validate(self, inlet: Substance, fuel: Substance = None, epsrel: float = EPSREL) -> Dict[int, float]:
        """Валиация найденного решения по уравнениям _equations"""
        args = {"inlet": inlet, "fuel": fuel, "prediction": {}}

        result: Dict[int, float] = {}
        for i, null in enumerate(self._equations([], args)):
            if isnan(null) or abs(null) > epsrel:
                result[i] = null

        return result

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "nodes": self.__nodes,
            "shafts": self.__shafts,
            "requirements": self.requirements,
        }


if __name__ == "__main__":
    from fixtures import air as inlet
    from fixtures import kerosene as fuel

    c = Rotor({gtep.effeff: 0.85, gtep.pipi: 6}, name="compressor")
    b = Burner({gtep.efficiency: 0.99, gtep.pipi: 0.95}, name="burner")
    t = Rotor({gtep.effeff: 1 / 0.9}, name="hpt")
    n = Nozzle({gtep.eff_speed: 0.98, gtep.pipi: 1 / 1.8}, name="n")

    gte = GTE("test")

    gte.add_edge(c, b)
    gte.add_edge(b, t)
    gte.add_edge(t, n)

    gte.add_shaft(c, t)

    gte.plot()
    # plt.show()

    print(f"\n{gte.is_solvable=}\n")

    print(f"\n{gte.order=}\n")

    ok = gte.solve(inlet, fuel, verbose=False)

    print(f"\n{gte.validate(inlet, fuel)=}\n")

    vars, substances = gte.calculate(inlet, fuel, verbose=False)

    for node, var in vars.items():
        print(f"{node}: {var}")
        for i, s in enumerate(substances[node]):
            print(f"{i}) {s.name}")
            for p, v in s.parameters.items():
                print(f"\t{p:<25}: {v:.4f}")
        print()
