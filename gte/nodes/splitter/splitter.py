from copy import deepcopy
from typing import Any, Dict, Tuple, Union

from numpy import isnan
from substance import Substance

try:
    from ...checks import check_mass_flow, check_temperature
    from ...config import EPSREL
    from ...config import parameters as gtep
    from ...errors import TYPE_ERROR
    from ..node import GTENode
except ImportError:
    import os
    import sys

    sys.path.insert(0, os.getcwd())

    from gte.checks import check_mass_flow, check_temperature
    from gte.config import EPSREL
    from gte.config import parameters as gtep
    from gte.errors import TYPE_ERROR
    from gte.nodes.node import GTENode


class Splitter(GTENode):
    """Камера отбора"""

    variables: Tuple[str] = ("splits",)
    n_vars: int = 1

    __slots__ = ()  # нет новых атрибутов

    def __init__(self, parameters: Dict[str, float], name: str = "Splitter"):
        GTENode.__init__(self, parameters, name)

    @classmethod
    def _equations(cls, x: Tuple[float], args: Dict[str, Any]) -> Tuple[float, float, float]:
        """
        total_m = sum(m)
        """

        inlet, outlets = args["inlet"], args["outlets"]

        return (inlet.parameters[gtep.m] - sum(outlet.parameters[gtep.m] for outlet in outlets),)

    @classmethod
    def predict(cls, parameters: Dict[str, Union[float, int]], inlet: Substance) -> Tuple[Dict[str, float], Substance]:
        """Начальные приближения"""
        GTENode.validate_substance(inlet)

        if not isinstance(parameters, dict):
            raise TypeError(TYPE_ERROR.format(f"{type(parameters)=}", dict))
        if len(parameters) != cls.n_vars:
            raise ValueError(f"{len(parameters)=} must be {cls.n_vars}")
        for parameter, value in parameters.items():
            if parameter not in cls.variables:
                raise ValueError(f"{parameter=} not in {cls.variables}")
            if not isinstance(value, (tuple, list)):
                raise TypeError(TYPE_ERROR.format(f"{type(value)=}", float))

        inlet = GTENode.calculate_substance(inlet)

        splits = parameters["splits"]
        if not isinstance(splits, (tuple, list)):
            raise TypeError(TYPE_ERROR.format(f"{type(splits)}", tuple))
        for split in splits:
            if not isinstance(split, (int, float)):
                raise TypeError(TYPE_ERROR.format(f"{type(split)}", float))

        total = sum(splits)  # общая масса

        outlets = []
        for split in splits:
            fraction = split / total  # нормализация

            outlet = deepcopy(inlet)
            outlet.parameters[gtep.m] *= fraction
            outlets.append(outlet)

        vars: Dict[str, float] = {}

        return vars, tuple(outlets)

    @classmethod
    def calculate(cls, parameters: Dict[str, Union[float, int]], inlet: Substance) -> Tuple[Dict[str, float], Substance]:
        _, outlets = cls.predict(parameters, inlet)

        return {}, outlets

    @classmethod
    def validate(cls, inlet: Substance, *outlets: Substance, epsrel: float = EPSREL) -> Dict[int, float]:

        args = {
            "inlet": inlet,
            "outlets": outlets,
        }

        result: Dict[int, float] = {}
        for i, null in enumerate(cls._equations([], args)):
            if isnan(null) or abs(null) > epsrel:
                result[i] = null

        return result

    @classmethod
    def check_real(cls, inlet: Substance, *outlets: Substance) -> str:
        if not check_mass_flow(inlet.parameters[gtep.m]):
            return f"inlet {gtep.m} {inlet.parameters[gtep.m]}"
        for i, outlet in enumerate(outlets):
            if not check_mass_flow(outlet.parameters[gtep.m]):
                return f"outlet[{i}] {gtep.m} {outlet.parameters[gtep.m]}"

        if not check_temperature(inlet.parameters[gtep.TT]):
            return f"inlet {gtep.TT} {inlet.parameters[gtep.TT]}"
        for i, outlet in enumerate(outlets):
            if not check_temperature(outlet.parameters[gtep.TT]):
                return f"outlet[{i}] {gtep.TT} {outlet.parameters[gtep.TT]}"

        if not check_temperature(inlet.parameters[gtep.PP]):
            return f"inlet {gtep.PP} {inlet.parameters[gtep.PP]}"
        for i, outlet in enumerate(outlets):
            if not check_temperature(outlet.parameters[gtep.PP]):
                return f"outlet[{i}] {gtep.PP} {outlet.parameters[gtep.PP]}"

        return ""


if __name__ == "__main__":
    from gte.fixtures import air as inlet

    test_cases = ({"parameters": {"splits": [1, 3, 6]}},)
    for test_case in test_cases:
        s = Splitter(test_case["parameters"], "test")
        print(f"{s.is_solvable=}")

        _, outlets = s.calculate(s.parameters, inlet)

        for i, outlet in enumerate(outlets):
            for k, v in outlet.parameters.items():
                print(f"{k:25}: {v:.4f}")
            print()

        print(f"{Splitter.validate(inlet, *outlets) = }")
        print(f"{Splitter.check_real(inlet, *outlets) = }")
        print()
