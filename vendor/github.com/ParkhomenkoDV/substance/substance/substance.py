from copy import deepcopy
from math import nan
from typing import Callable, Dict, Union

import numpy as np


class Substance:
    """Вещество"""

    __slots__ = (  # запрет других атрибутов + ускорение
        "name",  # имя
        "composition",  # химический состав
        "parameters",  # параметры
        "functions",  # функции
    )

    def __init__(
        self,
        name: str,
        composition: Dict[str, float] = None,
        parameters: Dict[str, Union[int, float]] = None,
        functions: Dict[str, Callable] = None,
    ) -> None:
        """
        Инициализация вещества.

        Args:
            name: Название вещества
            composition: Химический состав смеси (элемент: массовая доля)
            parameters: Физические параметры (название: значение)
            functions: Физические функции (название: функция)
        """
        self.name: str = name
        self.composition: Dict[str, float] = composition or {}
        self.parameters: Dict[str, Union[float, int]] = parameters or {}
        self.functions: Dict[str, Callable] = functions or {}

    def __validate_attribute(self, attribute: str, value: str | dict) -> str | dict:
        """Валидирование атрибутов"""
        if not isinstance(attribute, str):
            raise TypeError(f"{attribute} must be str")
        match attribute:
            case "name":
                if not isinstance(value, str):
                    raise TypeError(f"{attribute} must be a str")
                return value
            case "composition":
                if not isinstance(value, dict):
                    raise TypeError(f"{attribute} must be a dict")
                return self.__validate_composition(value)
            case "parameters":
                if not isinstance(value, dict):
                    raise TypeError(f"{attribute} must be a dict")
                return {k: self.__validate_parameter(k, v) for k, v in value.items()}
            case "functions":
                if not isinstance(value, dict):
                    raise TypeError(f"{attribute} must be a dict")
                return {k: self.__validate_function(k, v) for k, v in value.items()}
            case _:
                raise AttributeError(f"'{attribute}' not in {self.__slots__}")

    @staticmethod
    def normalize(composition: Dict[str, float]) -> Dict[str, float]:
        """Нормализация химического состава смеси"""
        mass = sum(composition.values())
        if mass == 0:
            return composition
        for element, fraction in composition.items():
            composition[element] = fraction / mass
        return composition

    def __validate_composition(self, composition: Dict[str, float]) -> Dict[str, float]:
        """Валидация смеси химического вещества"""
        for element, fraction in composition.items():
            if not isinstance(element, str):
                raise TypeError("Composition elements must be strings")
            if not isinstance(fraction, (int, float, np.number)):
                raise TypeError("Composition fractions must be numeric")
            if not (0 < fraction <= 1):
                raise ValueError("Composition values must be in (0..1]")
        composition = self.normalize(composition)
        return composition

    def __validate_parameter(self, key: str, value: Union[int, float]) -> Union[int, float]:
        """Валидация параметров"""
        if not isinstance(key, str):
            raise TypeError(f"{key} must be a str")
        if not isinstance(value, (int, float, np.number)):
            raise TypeError(f"Parameter {key}={value} must be numeric")
        return value

    def __validate_function(self, key: str, value: callable) -> Callable:
        """Валидация функций"""
        if not isinstance(key, str):
            raise TypeError(f"{key} must be a str")
        if not callable(value):
            raise TypeError(f"Function {key} value must be callable")
        return value

    def __setattr__(self, key: str, value) -> None:
        value = self.__validate_attribute(key, value)
        object.__setattr__(self, key, value)

    def __delattr__(self, key) -> None:
        raise Exception("Deleting forbidden!")

    def __deepcopy__(self, memo):
        """
        Создает глубокую копию объекта Substance.
        Обрабатывает:
        - Копирование строки name
        - Рекурсивное копирование словарей composition и parameters
        - Корректную обработку callable-объектов в parameters
        """
        new_obj = Substance.__new__(Substance)  # новый экземпляр без вызова __init__

        # Копируем простые атрибуты в memo для предотвращения циклических ссылок
        memo[id(self)] = new_obj

        # Глубокое копирование каждого атрибута
        new_obj.name = deepcopy(self.name, memo)

        # Особое внимание словарям
        new_obj.composition = {k: deepcopy(v, memo) for k, v in self.composition.items()}
        new_obj.parameters = {k: deepcopy(v, memo) for k, v in self.parameters.items()}
        new_obj.functions = {k: v for k, v in self.functions.items()}

        return new_obj

    def eq(self, other, eps: float) -> bool:
        if len(self.Parameters) != len(other.Parameters):
            return False

        for parameter, value in self.Parameters.items():
            v = other.Parameters.get(parameter)
            if v is None:
                return False
            if abs(v - value) > eps * value:
                return False

        return True

    @property
    def humidity(self) -> float:
        """Влажность"""
        h2o = self.composition.get("H2O", 0)
        total = sum(self.composition.values())
        if total == 0:
            return nan
        return h2o / total


if __name__ == "__main__":
    air = Substance(
        "air",
        composition={"N2": 0.78, "O2": 0.21, "Ar": 0.009, "CO2": 0.0004},
        parameters={
            "m": 50.0,
            "gc": 287.14,
            "TT": 300.0,
            "PP": 101325.0,
            "hcp": 1006.0,
            "k": 1.4,
            "c": 0.0,
        },
        functions={
            "gc": lambda total_temperature: 287.3,
            "hcp": lambda total_temperature: total_temperature * 1000,
        },
    )
