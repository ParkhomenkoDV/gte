from abc import ABC, abstractmethod
from typing import Any, Dict, Tuple, Union

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
    """Абстрактный базовый класс узла ГТД

    Attributes:
        variables: Кортеж названий переменных узла
        n_vars: Количество необходимых параметров для решения
        models: Словарь ML моделей для предсказаний
        figure: Кортеж с данными для визуализации (x, y координаты)
        name: Имя узла
        parameters: Словарь с параметрами узла
    """

    variables: Tuple[str, ...]  # переменные узла
    n_vars: int = -1  # необходимое количество параметров для решения
    models: Dict[str, Any] = {}  # ML модели
    figure: Tuple[Tuple[float, ...], Tuple[float, ...]]

    __slots__ = ("name", "parameters")

    def __init__(self, parameters: Dict[str, float], name: str = "node") -> None:
        """Инициализация объекта узла ГТД"""
        self.parameters: Dict[str, float] = parameters
        self.name: str = name

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}: {self.name}"

    def __setattr__(self, name, value) -> None:
        if name == "name":
            if not isinstance(value, str):
                raise TypeError(TYPE_ERROR.format(f"type({name})={type(value)}", str))
        elif name == "parameters":
            if not isinstance(value, dict):
                raise TypeError(TYPE_ERROR.format(f"type({name})={type(value)}", dict))
            for parameter, v in value.items():
                if parameter not in self.variables:
                    raise ValueError(f"{parameter=} not in {self.variables}")
                if not isinstance(v, (float, int)):
                    raise TypeError(TYPE_ERROR.format(f"type(value)={type(v)}", "numeric"))
        super().__setattr__(name, value)

    def __delattr__(self, name: str) -> None:
        if name == "name":
            self.name = self.__class__.__name__
        elif name == "characteristic":
            raise AttributeError("deleting characteristic is prohibited!")
        else:
            super().__delattr__(name)

    @property
    def is_solvable(self) -> Tuple[bool, str]:
        """Проверка возможности решения"""
        if len(self.parameters) < self.n_vars:
            return False, f"need to add {self.n_vars - len(self.parameters)} parameters"
        elif len(self.parameters) > self.n_vars:
            return False, f"need to delete {len(self.parameters) - self.n_vars} parameters"
        else:
            return True, ""

    @staticmethod
    def validate_substance(substance: Substance) -> None:
        """Проверка обязательных параметров рабочего тела"""
        if not isinstance(substance, Substance):
            raise TypeError(TYPE_ERROR.format(f"{type(substance)=}", Substance))
        # validate parameters
        if gtep.m not in substance.parameters:
            raise KeyError(SUBSTANCE_ATTRIBUTE_ERROR.format(substance.name, gtep.m))
        if gtep.TT not in substance.parameters:
            raise KeyError(SUBSTANCE_ATTRIBUTE_ERROR.format(substance.name, gtep.TT))
        if gtep.PP not in substance.parameters:
            raise KeyError(SUBSTANCE_ATTRIBUTE_ERROR.format(substance.name, gtep.PP))
        # validate functions
        tdp_keys = tdp.values()  # разрешенный список термодинамических параметров
        for name, function in substance.functions.items():
            if name not in tdp_keys:
                raise KeyError(f"function '{name}' not in {tdp_keys}")
            for func_arg in function.__code__.co_varnames:
                if func_arg not in tdp_keys:
                    raise KeyError(f"function '{name}' has argument '{func_arg}' not in {tdp_keys}")

    @classmethod
    @abstractmethod
    def predict(cls, inlet: Substance) -> Tuple[Dict[str, float], Substance]:
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
