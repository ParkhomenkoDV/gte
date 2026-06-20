import os
import pickle
from collections import defaultdict
from functools import wraps
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np
from scipy.integrate import nquad
from scipy.interpolate import LinearNDInterpolator as LNDI
from scipy.interpolate import interp1d

try:
    from ..errors import TYPE_ERROR
except ImportError:
    from errors import TYPE_ERROR


class Function:
    """Функция"""

    __slots__ = ("function", "name", "args")

    def __init__(self, function: Callable, name: str = None, args: Tuple[str] = None) -> None:
        self.function: Callable = function
        self.name: str = name if name else function.__name__

        if args:
            self.args: Tuple[str] = args
        else:
            co_varnames = function.__code__.co_varnames  # количество аргументов функции
            co_argcount = function.__code__.co_argcount  # все используемые переменные функции
            self.args: Tuple[str] = co_varnames[:co_argcount]  # только аргументы функции

    def __repr__(self) -> str:
        return self.name

    def __len__(self) -> int:
        """Количество аргументов функции"""
        return len(self.args)

    def __call__(self, kwargs: Dict[str, float]) -> float:
        """Вызов функции c избыточными параметрами"""
        return self.function(**{arg: kwargs[arg] for arg in self.args})


def integrate(function: Function, **kwargs) -> Tuple[float, float]:
    """Интегрирование"""
    if not isinstance(function, Function):
        raise TypeError(TYPE_ERROR.format(f"{type(function)}"), Function)

    fixed, ranges, other_args = {}, [], {}
    for arg in function.args:
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
        return function(fixed), 0.0

    def partial(*args) -> float:
        return function({a: fixed[a] if a in fixed else args[other_args[a]] for a in function.args})

    result, abserr = nquad(partial, ranges)

    return result, abserr


def integral_average(function: Function, **kwargs) -> Tuple[float, float]:
    """Среднее интегральное"""

    result, abserr = integrate(function, **kwargs)  # + checks

    devider: float = 1.0
    for arg in function.args:
        rang = kwargs[arg]
        if rang[0] != rang[1]:
            devider *= rang[1] - rang[0]

    return result / devider, abserr


class Solvable:
    """Решаемость"""

    __slots__ = ("reason",)

    def __init__(self, reason: str) -> None:
        self.reason = reason

    def __repr__(self) -> str:
        return str(bool(self))

    def __bool__(self) -> bool:
        if self.reason == "":
            return True
        else:
            return False


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

    def hcp(temperature, eo):
        return temperature + eo

    hcp_ = Function(hcp, "test")

    print(hcp_.function)
    print(hcp_.name)
    print(hcp_.args)
    print(hcp_({"temperature": 300, "eo": 3}))

    print(integrate(hcp, temperature=(300, 600), eo=(2, 3)), (600**2 - 300**2) / 2 + (3**2 - 2**2) / 2)
    print(integrate(hcp, temperature=(300, 600), eo=(2, 2)), (600**2 - 300**2) / 2 + (2**2 - 2**2) / 2)

    print(integral_average(hcp, temperature=(300, 600), eo=(0, 1)), (600**2 - 300**2) / 2 / (600 - 300) + (1**2 - 0**2) / 2 / (1 - 0))
    print(integral_average(hcp, temperature=(300, 600), eo=(1, 1)), (600**2 - 300**2) / 2 / (600 - 300) + 1)

    def z(x, y):
        return float(np.sin(x) * np.exp(y))

    datas = [{"x": x, "y": y, "z": z(x, y)} for x in range(10) for y in range(10)]

    f = Interpolator(datas, "z")
