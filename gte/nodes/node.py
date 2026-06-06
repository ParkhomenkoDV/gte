from abc import ABC, abstractmethod
from itertools import product
from typing import Any, Dict, Generator, Tuple, Union

from substance import Substance
from thermodynamics import adiabatic_index, critical_sonic_velocity

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
        class:
            variables: Кортеж названий переменных узла
            n_vars: Количество необходимых параметров для решения
        object:
            name: Имя узла
            parameters: Словарь с параметрами узла
    """

    variables: Tuple[str, ...]  # переменные узла

    __slots__ = ("name", "parameters")

    def __init__(self, parameters: Dict[str, float], name: str = "node") -> None:
        """Инициализация объекта узла ГТД"""
        self.parameters: Dict[str, float] = parameters
        self.name: str = name

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}={self.name}"

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
                if not isinstance(v, (float, int, tuple, list)):
                    raise TypeError(TYPE_ERROR.format(f"type(value)={type(v)}", "numeric"))
        super().__setattr__(name, value)

    def __delattr__(self, name: str) -> None:
        if name == "name":
            self.name = self.__class__.__name__
        elif name == "characteristic":
            raise AttributeError("deleting characteristic is prohibited!")
        else:
            super().__delattr__(name)

    @classmethod
    def _equations(cls, x: Tuple[float], args: Dict[str, Any]) -> Tuple[float, ...]:
        """Система уравнений"""
        return tuple()

    @property
    def n_vars(self) -> int:
        """Необходимое количество параметров для решения"""
        return len(self._equations([], {}))

    @property
    def is_solvable(self) -> Tuple[bool, str]:
        """Проверка возможности решения"""
        dif = abs(self.n_vars - len(self.parameters))
        if len(self.parameters) < self.n_vars:
            return False, f"need to add {dif} parameters"
        elif len(self.parameters) > self.n_vars:
            return False, f"need to delete {dif} parameters"
        else:
            return True, ""

    @classmethod
    def generator(cls, **parameters) -> Generator:
        """Генератор решаемых узлов"""
        if len(parameters) != cls.n_vars:
            raise ArithmeticError(f"{len(parameters)=} must be equal {cls.n_vars=}")
        for parameter, values in parameters.items():
            if parameter not in cls.variables:
                raise KeyError(f"{parameter=} not in {cls.variables}")
            if not isinstance(values, (tuple, list)):
                raise TypeError(f"{type(values)=} must be tuple")
            if len(values) == 0:
                raise ValueError(f"{len(values)=} must be > 0")

        keys = tuple(parameters.keys())
        values = tuple(parameters[key] for key in keys)
        for name, combination in enumerate(product(*values)):
            params = dict(zip(keys, combination))
            yield cls(params, f"{name}")

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
        if gtep.eo in substance.parameters:  # в случае наличия избытка окислителя
            if "oxidizer" not in substance.parameters:  # необходимо знать массовую долю окислителя
                raise KeyError(SUBSTANCE_ATTRIBUTE_ERROR.format(substance.name, "oxidizer"))  # для перерасчета в случае смешения
        # validate functions
        tdp_keys = gtep.values()  # разрешенный список термодинамических параметров
        for name, function in substance.functions.items():
            if name not in tdp_keys:
                raise KeyError(f"function '{name}' not in {tdp_keys}")
            """for func_arg in function.__code__.co_varnames:
                if func_arg not in tdp_keys:
                    raise KeyError(f"function '{name}' has argument '{func_arg}' not in {tdp_keys}")"""

    @classmethod
    @abstractmethod
    def predict(cls, parameters: Dict[str, Union[float, int]], inlet: Substance) -> Tuple[Dict[str, float], Substance]:
        """Начальные приближения"""
        raise NotImplementedError

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

        substance.parameters[gtep.gc] = call_with_kwargs(substance.functions[gtep.gc], **substance.parameters)
        substance.parameters[gtep.hcp] = call_with_kwargs(substance.functions[gtep.hcp], **substance.parameters)
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
