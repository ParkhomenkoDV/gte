from typing import Any, Dict, List, Tuple, Union

from numpy import isnan
from substance import Substance

try:
    from ...checks import check_mass_flow, check_pressure, check_temperature
    from ...config import EPSREL
    from ...config import parameters as gtep
    from ...errors import TYPE_ERROR
    from ...utils import call_with_kwargs
    from ..node import Node
except ImportError:
    import os
    import sys

    sys.path.insert(0, os.getcwd())

    from gte.checks import check_mass_flow, check_pressure, check_temperature
    from gte.config import EPSREL
    from gte.config import parameters as gtep
    from gte.errors import TYPE_ERROR
    from gte.nodes.node import Node
    from gte.utils import call_with_kwargs


class Joiner(Node):
    """Камера смешения"""

    variables: Tuple[str] = tuple()
    n_vars: int = 0

    __slots__ = ()  # нет новых атрибутов

    def __init__(self, parameters: Dict[str, float], name: str = "Joiner"):
        Node.__init__(self, parameters, name)

    @classmethod
    def _equations(cls, x: Tuple[float], args: Dict[str, Any]) -> Tuple[float, float]:
        """
        m_outlet = m_inlet_0 + m_inlet_1 + ...
        T*_outlet = (m_inlet_0 * hcp_inlet_0 * T*_inlet_0 + m_inlet_1 * hcp_inlet_1 * T*_inlet_1 + ...) / (m_inlet_0 + m_inlet_1 + ...)
        P*_outlet = (m_inlet_0 * P*_inlet_0 + m_inlet_1 * P*_inlet_1 + ...) / (m_inlet_0 + m_inlet_1 + ...)
        """

        inlets_got, outlet_got = args["inlets"], args["outlet"]

        _, outlet_want = cls.calculate({}, *inlets_got)

        return (
            outlet_got.parameters[gtep.m] - outlet_want.parameters[gtep.m],
            outlet_got.parameters[gtep.TT] - outlet_want.parameters[gtep.TT],
            outlet_got.parameters[gtep.PP] - outlet_want.parameters[gtep.PP],
        )

    @classmethod
    def predict(cls, parameters: Dict[str, Union[float, int]], *inlets: Substance) -> Tuple[Dict[str, float], Substance]:
        """Начальные приближения"""
        if not isinstance(parameters, dict):
            raise TypeError(TYPE_ERROR.format(f"{type(parameters)=}", dict))
        if len(parameters) != cls.n_vars:
            raise ValueError(f"{len(parameters)=} must be {cls.n_vars}")
        for parameter, value in parameters.items():
            if parameter not in cls.variables:
                raise ValueError(f"{parameter=} not in {cls.variables}")
            if not isinstance(value, (float, int)):
                raise TypeError(TYPE_ERROR.format(f"{type(value)=}", float))

        names: List[str] = []
        oxidizer, required = 0, 0  # окислитель, теоретически необходимое кодичество окислителя
        m, m_hcp, m_t_hcp, m_p = 0, 0, 0, 0
        for inlet in inlets:
            Node.validate_substance(inlet)

            names.append(inlet.name)
            hcp = call_with_kwargs(inlet.functions[gtep.hcp], **inlet.parameters)
            m += inlet.parameters[gtep.m]
            oxidizer += inlet.parameters["oxidizer"] if gtep.eo in inlet.parameters else inlet.parameters[gtep.m]
            required += inlet.parameters.get("oxidizer", 0) / inlet.parameters.get(gtep.eo, 1)  # 0 for oxidizer
            m_hcp += inlet.parameters[gtep.m] * hcp
            m_t_hcp += inlet.parameters[gtep.m] * inlet.parameters[gtep.TT] * hcp
            m_p += inlet.parameters[gtep.m] * inlet.parameters[gtep.PP]

        def gc(total_temperature, total_pressure, excess_oxidizing) -> float:
            result = 0
            for inlet in inlets:
                gc = call_with_kwargs(inlet.functions[gtep.gc], **{gtep.TT: total_temperature, gtep.PP: total_pressure, gtep.eo: excess_oxidizing})
                result += inlet.parameters[gtep.m] * gc
            return result / m

        def hcp(total_temperature, total_pressure, excess_oxidizing) -> float:
            result = 0
            for inlet in inlets:
                hcp = call_with_kwargs(inlet.functions[gtep.hcp], **{gtep.TT: total_temperature, gtep.PP: total_pressure, gtep.eo: excess_oxidizing})
                result += inlet.parameters[gtep.m] * inlet.parameters[gtep.TT] * hcp
            return result / m_hcp

        outlet = Substance(
            "+".join(names),
            {},  # TODO
            parameters={
                gtep.m: m,
                gtep.eo: oxidizer / required,
                gtep.TT: m_t_hcp / m_hcp,
                gtep.PP: m_p / m,
            },
            functions={
                gtep.gc: gc,
                gtep.hcp: hcp,
            },
        )
        if required != 0:
            outlet.parameters["oxidizer"] = oxidizer
            outlet.parameters[gtep.eo] = oxidizer / required

        vars: Dict[str, float] = {}

        return vars, outlet

    @classmethod
    def calculate(cls, parameters: Dict[str, Union[float, int]], *inlets: Substance) -> Tuple[Dict[str, float], Substance]:
        _, outlet = cls.predict(parameters, *inlets)

        outlet = Node.calculate_substance(outlet)

        return {}, outlet

    @classmethod
    def validate(cls, *substances: Substance, epsrel: float = EPSREL) -> Dict[int, float]:
        if len(substances) < 2:
            raise ValueError(f"{len(substances)=} must be >= 2")

        inlets, outlet = substances[:-1], substances[-1]

        x0 = tuple()
        args = {"inlets": inlets, "outlet": outlet}

        result: Dict[int, float] = {}
        for i, null in enumerate(cls._equations(x0, args)):
            if isnan(null) or abs(null) > epsrel:
                result[i] = null

        return result

    @classmethod
    def check_real(cls, *substances: Substance) -> str:
        if len(substances) < 2:
            raise ValueError(f"{len(substances)=} must be >= 2")

        for i, substance in enumerate(substances):
            if not check_mass_flow(substance.parameters[gtep.m]):
                return f"substance[{i}] {gtep.m} {substance.parameters[gtep.m]}"
            if not check_temperature(substance.parameters[gtep.TT]):
                return f"substance[{i}] {gtep.TT} {substance.parameters[gtep.TT]}"
            if not check_pressure(substance.parameters[gtep.PP]):
                return f"substance[{i}] {gtep.PP} {substance.parameters[gtep.PP]}"

        return ""


if __name__ == "__main__":
    from gte.fixtures import air, exhaust

    j = Joiner({}, "test")
    print(f"{j.is_solvable=}")

    _, outlet = j.calculate({}, air, exhaust)

    for k, v in outlet.parameters.items():
        print(f"{k:<25}: {v:.4f}")

    print(f"{Joiner.validate(air, exhaust, outlet) = }")
    print(f"{Joiner.check_real(air, exhaust, outlet) = }")
    print()
