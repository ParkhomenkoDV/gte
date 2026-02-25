from abc import ABC, abstractmethod
from typing import Any, Callable, Dict, Tuple, Union

from substance import Substance
from thermodynamics import adiabatic_index, critical_sonic_velocity
from thermodynamics import parameters as tdp

try:  # Попытка относительного импорта для удаленного пакета
    from ..config import parameters as gtep
    from ..errors import SUBSTANCE_ATTRIBUTE_ERROR, TYPE_ERROR
    from ..utils import call_with_kwargs
except ImportError:  # Резервный абсолютный импорт для локального запуска
    import os
    import sys

    sys.path.insert(0, os.getcwd())

    from gte.config import parameters as gtep
    from gte.errors import SUBSTANCE_ATTRIBUTE_ERROR, TYPE_ERROR
    from gte.utils import call_with_kwargs


class GTENode(ABC):
    """Абстрактный базовый класс узла ГТД"""

    variables: Tuple[str, ...]  # переменные узла
    models: Dict[str, Any] = {}  # ML модели
    figure: Tuple[Tuple[float, ...], Tuple[float, ...]]

    __slots__ = ["name", "characteristic"]  # list to add

    def __init__(self, characteristic: Dict[str, Callable], name: str = "node") -> None:
        self.name: str = name
        self.characteristic: Dict[str, Callable] = characteristic

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}: {self.name}"

    def __setattr__(self, name, value) -> None:
        if name == "name":
            assert isinstance(value, str), TYPE_ERROR.format(f"type({name})={type(value)}", str)
        elif name == "characteristic":
            assert isinstance(value, dict), TYPE_ERROR.format(f"type({name})={type(value)}", dict)
        super().__setattr__(name, value)

    def __delattr__(self, name: str) -> None:
        if name == "name":
            self.name = self.__class__.__name__
        elif name == "characteristic":
            raise AttributeError("deleting characteristic is prohibited!")
        else:
            super().__delattr__(name)

    @staticmethod
    def validate_substance(substance: Substance) -> None:
        """Проверка обязательных параметров рабочего тела"""
        assert isinstance(substance, Substance), TypeError(f"type substance must be {Substance}")
        # validate parameters
        assert gtep.m in substance.parameters, SUBSTANCE_ATTRIBUTE_ERROR.format(substance.name, gtep.m)
        assert gtep.TT in substance.parameters, SUBSTANCE_ATTRIBUTE_ERROR.format(substance.name, gtep.TT)
        assert gtep.PP in substance.parameters, SUBSTANCE_ATTRIBUTE_ERROR.format(substance.name, gtep.PP)
        # validate functions
        tdp_keys = tdp.values()  # разрешенный список термодинамических параметров
        for name, function in substance.functions.items():
            assert name in tdp_keys, f"function '{name}' not in {tdp_keys}"
            for func_arg in function.__code__.co_varnames:
                assert func_arg in tdp_keys, f"function '{name}' has argument '{func_arg}' not in {tdp_keys}"

    @classmethod
    @abstractmethod
    def predict(cls, inlet: Substance, use_ml: bool = True) -> Dict[str, float]:
        """Начальные приближения"""
        raise NotImplementedError

    @classmethod
    def _equations(cls, x: Tuple[float], args: Dict[str, Any]) -> Tuple[float, ...]:
        """Система уравнений"""
        return tuple()

    @classmethod
    @abstractmethod
    def calculate(cls, parameters: Dict[str, Union[float, int]] = None) -> Substance:
        """Термодинамический расчет узла по уравнениям _equations"""
        # валидация входных параметров
        # расчет входных параметров
        # расчет параметров узла
        # расчет выходных параметров
        # вывод выходных параметров
        raise NotImplementedError

    @abstractmethod
    def solve(self, inlet: Substance) -> Dict[str, Any]:
        """Термодинамический расчет узла по его характеристике и уравнениям _equations"""
        raise NotImplementedError

    @staticmethod
    def calculate_substance(substance: Substance) -> Substance:
        """Расчет термодинамических параметров вещества по массе, температуре, давлению"""
        GTENode.validate_substance(substance)

        parameters = {tdp.t: substance.parameters[gtep.TT], tdp.p: substance.parameters[gtep.PP], tdp.eo: substance.parameters.get(gtep.eo)}
        substance.parameters[gtep.gc] = call_with_kwargs(substance.functions[gtep.gc], parameters)
        substance.parameters[gtep.hcp] = call_with_kwargs(substance.functions[gtep.hcp], parameters)
        substance.parameters[gtep.DD] = substance.parameters[gtep.PP] / (substance.parameters[gtep.gc] * substance.parameters[gtep.TT])
        substance.parameters[gtep.k] = adiabatic_index(substance.parameters[gtep.gc], substance.parameters[gtep.hcp])
        substance.parameters[gtep.ss_critical] = critical_sonic_velocity(substance.parameters[gtep.k], substance.parameters[gtep.gc], substance.parameters[gtep.TT])

        return substance

    @classmethod
    @abstractmethod
    def validate(cls) -> Dict[int, float]:
        """Валиация найденного решения по уравнениям _equations"""
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def check_real(cls) -> str:
        """Проверка на физичность"""
        raise NotImplementedError
