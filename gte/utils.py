import inspect
from functools import lru_cache
from typing import Any, Callable, Dict, Tuple

import numpy as np
from scipy.integrate import nquad


@lru_cache(maxsize=256)
def get_function_signature(function: Callable):
    """Получить сингнатуру функции с помощью inspect"""
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
    assert callable(function), TypeError(f"{function} must be callable")

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

    devider = 1.0
    for arg in sig.parameters:
        rang = kwargs[arg]
        if rang[0] != rang[1]:
            devider *= rang[1] - rang[0]

    return result / devider, abserr


if __name__ == "__main__":

    def cp(temperature, eo):
        return temperature + eo

    print(integrate(cp, temperature=(300, 600), eo=(2, 3)), (600**2 - 300**2) / 2 + (3**2 - 2**2) / 2)
    print(integrate(cp, temperature=(300, 600), eo=(2, 2)), (600**2 - 300**2) / 2 + (2**2 - 2**2) / 2)

    print(integral_average(cp, temperature=(300, 600), eo=(0, 1)), (600**2 - 300**2) / 2 / (600 - 300) + (1**2 - 0**2) / 2 / (1 - 0))
    print(integral_average(cp, temperature=(300, 600), eo=(1, 1)), (600**2 - 300**2) / 2 / (600 - 300) + 1)
