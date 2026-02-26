import inspect
import os
import pickle
from collections import defaultdict
from functools import lru_cache, wraps
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np
from scipy.integrate import nquad
from scipy.interpolate import LinearNDInterpolator as LNDI
from scipy.interpolate import interp1d

try:
    from .errors import TYPE_ERROR
except ImportError:
    from errors import TYPE_ERROR


@lru_cache(maxsize=256)
def get_function_signature(function: Callable):
    """Получение сингнатуры функции с помощью inspect"""
    return inspect.signature(function)


def call_with_kwargs(function: Callable, kwargs: Dict[str, Any]):
    """Вызов функции через kwargs"""
    if not callable(function):
        raise TypeError(f"{function} must be callable")
    if not isinstance(kwargs, dict):
        raise TypeError(f"type {kwargs} must be dict")

    sig = get_function_signature(function)

    arguments: Dict[str, Any] = {}
    for key, value in sig.parameters.items():
        if key in ("cls", "self"):
            continue

        # Проверяем обязательные параметры (без значений по умолчанию)
        if value.default is value.empty and key not in kwargs:
            raise ValueError(f"kwargs has not required parameter '{key}'")

        # Добавляем значение из kwargs, если оно есть
        if key in kwargs:
            arguments[key] = kwargs[key]
        # Если параметра нет в kwargs, но у него есть значение по умолчанию -
        # функция будет использовать значение по умолчанию

    return function(**arguments)


def integrate(function: Callable, **kwargs) -> Tuple[float, float]:
    """Интегрирование"""
    if not callable(function):
        TypeError(f"{function} must be callable")

    sig = get_function_signature(function)

    fixed, ranges, other_args = {}, [], {}
    for arg in sig.parameters:
        rang = kwargs.get(arg)
        if rang is None:
            raise ValueError(f"{function=} require argument '{arg}'")
        if not isinstance(rang, (tuple, list, np.ndarray)):
            raise TypeError(f"type of range must be tuple, but has {type(rang)}")
        if len(rang) != 2:
            ValueError(f"integral has 2 ranges, but has {len(rang)}")

        if rang[0] == rang[1]:
            fixed[arg] = rang[0]
        else:
            other_args[arg] = len(ranges)
            ranges.append(rang)

    if not ranges:
        return function(**fixed), 0.0

    def partial(*args):
        return function(*[fixed[p] if p in fixed else args[other_args[p]] for p in sig.parameters])

    result, abserr = nquad(partial, ranges)

    return result, abserr


def integral_average(function: Callable, **kwargs) -> Tuple[float, float]:
    """Среднее интегральное"""

    result, abserr = integrate(function, **kwargs)  # + checks

    sig = get_function_signature(function)

    devider: float = 1.0
    for arg in sig.parameters:
        rang = kwargs[arg]
        if rang[0] != rang[1]:
            devider *= rang[1] - rang[0]

    return result / devider, abserr


def load_models(*paths) -> Dict[str, Any]:
    """Загрузка моделей по переданным путям"""
    models: Dict = {}
    for path in paths:
        if not os.path.isfile(path):
            raise FileNotFoundError(f"'{path}' not found!")
        name = os.path.dirname(path)
        with open(path, "rb") as file:
            models[name] = pickle.load(file)

    return models


class Interpolator:
    """Линейное интерполирование данных"""

    __slots__ = (
        "features",  # аргументы функции
        "target",  # целевая переменная (название функции)
        "function",  # интерполированная функция
        "fill_value",  # значение вне диаппазона интерполирования
    )

    def __init__(self, datas: List[Dict[str, float]], target: str, features: Optional[List[str]] = None, fill_value: float = np.nan) -> None:
        # Валидация данных
        if not isinstance(datas, (list, tuple)):
            raise TypeError(TYPE_ERROR.format(f"{type(datas)=}", tuple))
        if len(datas) < 3:
            raise ValueError(f"{len(datas)=} must be >= 3")
        if not isinstance(target, str):
            raise TypeError(TYPE_ERROR.format(f"{type(target)=}", str))
        if features is not None:
            if not isinstance(features, (list, tuple)):
                raise TypeError(TYPE_ERROR.format(f"{type(features)=}"), list)
            if len(features) == 0:
                raise ValueError(f"{len(features)=} must be > 0")
            if not all(isinstance(f, str) for f in features):
                raise TypeError("All features must be strings")
            features_set = set(features)
        if not isinstance(fill_value, (float, int, np.number)):
            raise TypeError(TYPE_ERROR.format(f"{type(fill_value)=}", float))

        # Сбор данных и валидация структуры
        points: Dict[str, List[float]] = defaultdict(list)
        for i, data in enumerate(datas):
            if not isinstance(data, dict):
                raise TypeError(TYPE_ERROR.format(f"type(datas[{i}])", dict))
            if target not in data:
                raise KeyError(f"{target=} not in datas[{i}]")
            if features and not features_set.issubset(set(data.keys())):
                raise KeyError(f"datas[{i}].keys()={data.keys()} not contains {features=}")
            for key, value in data.items():
                if not isinstance(key, str):
                    raise TypeError(TYPE_ERROR.format(f"{type(key)=} in datas[{i}]", str))

                if features and key not in features and key != target:  # если переданы конкретные аргументы, то рассматриваем только их
                    continue

                if not isinstance(value, (float, int)):
                    raise TypeError(TYPE_ERROR.format(f"{type(value)=} datas[{i}]", float))

                points[key].append(float(value))

        # Проверка согласованности длин
        len_target: int = len(points[target])
        for key, values in points.items():
            if len(values) != len_target:
                raise ValueError(f"Inconsistent length for key '{key}': {len(values)=} != {len_target=}")

        # Определение атрибутов
        self.target: str = target
        self.features: Tuple[str, ...] = tuple(feature for feature in points if feature != target)
        if not self.features:
            raise ValueError(f"{len(self.features)=} must be > 0")
        self.fill_value: float = float(fill_value)

        x = np.array([points[k] for k in self.features]).T  # (n_samples, n_features)
        y = np.array(points[target])  # (n_samples,)

        # Удаление дубликатов точек (если есть)
        if len(x) != len(np.unique(x, axis=0)):
            unique_indices = np.unique(x, axis=0, return_index=True)[1]
            x = x[unique_indices]
            y = y[unique_indices]

        if len(self.features) == 1:
            interpolator = interp1d(x.T[0], y, kind="linear", fill_value=self.fill_value, bounds_error=False)
        else:
            interpolator = LNDI(x, y, fill_value=self.fill_value)

        @wraps(interpolator)
        def function(**kwargs) -> float:
            try:
                args = tuple(kwargs[feature] for feature in self.features)
            except KeyError as e:  # пропуск фичи
                missing = set(self.features) - set(kwargs.keys())
                raise KeyError(f"Missing required features: {missing}") from e
            except TypeError as e:  # неверный тип данных
                raise TypeError(f"Invalid type for one of the features. Expected float/int, got: {e}") from e

            if isinstance(interpolator, interp1d):
                return float(interpolator(args)[0])
            else:
                x_query = np.array([args])  # (1, n_features)
                return float(interpolator(x_query)[0])

        self.function: Callable = function

    def __repr__(self) -> str:
        return self.target

    def __call__(self, **kwargs) -> float:
        """Вызов интерполированной функции"""
        return self.function(**kwargs)


if __name__ == "__main__":

    def cp(temperature, eo):
        return temperature + eo

    print(integrate(cp, temperature=(300, 600), eo=(2, 3)), (600**2 - 300**2) / 2 + (3**2 - 2**2) / 2)
    print(integrate(cp, temperature=(300, 600), eo=(2, 2)), (600**2 - 300**2) / 2 + (2**2 - 2**2) / 2)

    print(integral_average(cp, temperature=(300, 600), eo=(0, 1)), (600**2 - 300**2) / 2 / (600 - 300) + (1**2 - 0**2) / 2 / (1 - 0))
    print(integral_average(cp, temperature=(300, 600), eo=(1, 1)), (600**2 - 300**2) / 2 / (600 - 300) + 1)

    def z(x, y):
        return float(np.sin(x) * np.exp(y))

    datas = [{"x": x, "y": y, "z": z(x, y)} for x in range(10) for y in range(10)]

    f = Interpolator(datas, "z")
