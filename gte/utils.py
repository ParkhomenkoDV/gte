import inspect
from functools import lru_cache
from typing import Any, Callable, Dict, Tuple

import numpy as np
from scipy.integrate import nquad


@lru_cache(maxsize=128)
def get_function_signature(function: Callable):
    """Получить сингнатуру функции с помощью inspect"""
    return inspect.signature(function)


def call_with_kwargs(function: Callable, kwargs: Dict[str, Any]):
    """Вызов функции через kwargs"""
    assert callable(function), TypeError(f"{function} must be callable")
    assert isinstance(kwargs, dict), TypeError(f"type {kwargs} must be dict")

    sig = get_function_signature(function)

    arguments: Dict[str, Any] = {}
    for key, value in sig.parameters.items():
        if key in ("cls", "self"):
            continue

        # Проверяем обязательные параметры (без значений по умолчанию)
        assert value.default is not value.empty or key in kwargs, ValueError(f"kwargs has not required parameter '{key}'")

        # Добавляем значение из kwargs, если оно есть
        if key in kwargs:
            arguments[key] = kwargs[key]
        # Если параметра нет в kwargs, но у него есть значение по умолчанию -
        # функция будет использовать значение по умолчанию

    return function(**arguments)


def enthalpy(heat_capacity: Callable, **kwargs) -> float:
    """Энтальпия"""
    assert callable(heat_capacity), TypeError(f"{heat_capacity} must be callable")

    partial_args, other_args = {}, {}
    arg_names = heat_capacity.__code__.co_varnames[: heat_capacity.__code__.co_argcount]
    for arg_name in arg_names:
        rang = kwargs.get(arg_name)
        assert rang is not None, Exception(f"{heat_capacity=} require argument '{arg_name}'")
        assert isinstance(rang, (tuple, list, np.ndarray)), TypeError(f"type of range must be tuple, but has {type(rang)}")
        assert len(rang) == 2, ValueError(f"integral has 2 ranges, but has {len(rang)}")

        if rang[0] == rang[1]:
            partial_args[arg_name] = rang[0]
        else:
            other_args[arg_name] = (rang[0], rang[1])

    def partial(*args):
        all_args = partial_args.copy()

        for arg_name in arg_names:
            # Если этот аргумент не фиксирован И еще не заполнен
            if arg_name not in partial_args and arg_name not in all_args:
                all_args[arg_name] = args[len(all_args) - len(partial_args)]

        # Сортируем аргументы по порядку в функции
        sorted_args = [all_args[arg] for arg in arg_names]
        return heat_capacity(*sorted_args)

    if not other_args:
        return heat_capacity(**partial_args)

    ranges = tuple(other_args.values())
    result, _ = nquad(partial, ranges)

    return result


def integral_average(function, **kwargs) -> Tuple[float, float]:
    """Среднее интегральное"""
    assert callable(function), TypeError(f"{function} must be callable")

    partial_args, other_args, devider = {}, {}, 1.0
    arg_names = function.__code__.co_varnames[: function.__code__.co_argcount]
    for arg_name in arg_names:
        rang = kwargs.get(arg_name)
        assert rang is not None, Exception(f"{function=} require argument '{arg_name}'")
        assert isinstance(rang, (tuple, list, np.ndarray)), TypeError(f"type of range must be tuple, but has {type(rang)}")
        assert len(rang) == 2, ValueError(f"integral has 2 ranges, but has {len(rang)}")

        if rang[0] == rang[1]:
            partial_args[arg_name] = rang[0]
        else:
            other_args[arg_name] = (rang[0], rang[1])
            devider *= rang[1] - rang[0]

    def partial(*args):
        all_args = partial_args.copy()

        for arg_name in arg_names:
            # Если этот аргумент не фиксирован И еще не заполнен
            if arg_name not in partial_args and arg_name not in all_args:
                all_args[arg_name] = args[len(all_args) - len(partial_args)]

        # Сортируем аргументы по порядку в функции
        sorted_args = [all_args[arg] for arg in arg_names]
        return function(*sorted_args)

    if not other_args:
        return function(**partial_args), 0.0

    ranges = tuple(other_args.values())
    result, abserr = nquad(partial, ranges)

    return result / devider, abserr


'''
from typing import Callable, Tuple, Dict, Any, Sequence
import numpy as np
from scipy.integrate import nquad
from inspect import signature, Parameter

def _validate_and_split_args(
    func: Callable,
    kwargs: Dict[str, Any]
) -> Tuple[Dict[str, float], Dict[str, Tuple[float, float]], float]:
    """
    Валидирует аргументы и разделяет их на фиксированные и интеграционные.
    Возвращает: (fixed_args, integral_ranges, volume_multiplier)
    """
    sig = signature(func)
    fixed_args = {}
    integral_ranges = {}
    volume_multiplier = 1.0

    for param_name, param in sig.parameters.items():
        if param_name not in kwargs:
            raise ValueError(f"Function {func.__name__} requires argument '{param_name}'")

        rang = kwargs[param_name]
        
        # Проверка типа диапазона
        if not isinstance(rang, (tuple, list, np.ndarray)):
            raise TypeError(f"Range for '{param_name}' must be tuple/list/ndarray, got {type(rang)}")
        if len(rang) != 2:
            raise ValueError(f!Range for '{param_name}' must have 2 elements, got {len(rang)}")

        lower, upper = float(rang[0]), float(rang[1])
        
        if lower == upper:
            fixed_args[param_name] = lower
        else:
            integral_ranges[param_name] = (lower, upper)
            volume_multiplier *= (upper - lower)

    return fixed_args, integral_ranges, volume_multiplier



def _create_partial_function(
    func: Callable,
    fixed_args: Dict[str, float],
    param_order: Sequence[str]
) -> Callable:
    """Создаёт частичную функцию с фиксированными аргументами."""
    def partial(*args):
        all_args = fixed_args.copy()
        # Заполняем аргументы по порядку
        for i, param_name in enumerate(param_order):
            if param_name not in fixed_args:
                all_args[param_name] = args[i - len(fixed_args)]
        # Вызываем с отсортированными аргументами
        return func(**{k: all_args[k] for k in param_order})
    
    return partial

def enthalpy(heat_capacity: Callable, **kwargs) -> Tuple[float, float]:
    """
    Вычисляет энтальпию через интеграл функции теплоёмкости.
    
    Args:
        heat_capacity: функция теплоёмкости (должна принимать все аргументы из kwargs)
        **kwargs: диапазоны интегрирования для каждого параметра вида param=(lower, upper)
    
    Returns:
        (значение интеграла, оценка погрешности)
    """
    if not callable(heat_capacity):
        raise TypeError(f"heat_capacity must be callable, got {type(heat_capacity)}")

    fixed_args, integral_ranges, _ = _validate_and_split_args(heat_capacity, kwargs)
    param_names = list(signature(heat_capacity).parameters.keys())

    if not integral_ranges:
        return heat_capacity(**fixed_args), 0.0

    partial_func = _create_partial_function(heat_capacity, fixed_args, param_names)
    result, abserr = nquad(partial_func, tuple(integral_ranges.values()))
    
    return result, abserr

def integral_average(function: Callable, **kwargs) -> Tuple[float, float]:
    """
    Вычисляет среднее интегральное значение функции.
    
    Args:
        function: интегрируемая функция
        **kwargs: диапазоны интегрирования для каждого параметра вида param=(lower, upper)
    
    Returns:
        (среднее значение, оценка погрешности)
    """
    if not callable(function):
        raise TypeError(f"function must be callable, got {type(function)}")

    fixed_args, integral_ranges, volume = _validate_and_split_args(function, kwargs)
    param_names = list(signature(function).parameters.keys())

    if not integral_ranges:
        return function(**fixed_args), 0.0

    partial_func = _create_partial_function(function, fixed_args, param_names)
    integral_result, abserr = nquad(partial_func, tuple(integral_ranges.values()))
    
    return integral_result / volume, abserr

'''

if __name__ == "__main__":

    def cp(temperature, eo):
        return temperature + eo

    print(integral_average(cp, temperature=(300, 600), eo=(0, 1)), (600**2 - 300**2) / 2 / (600 - 300) + (1**2 - 0**2) / 2 / (1 - 0))
    print(integral_average(cp, temperature=(300, 600), eo=(1, 1)), (600**2 - 300**2) / 2 / (600 - 300) + 1)

    print(enthalpy(cp, temperature=(300, 600), eo=(2, 3)), (600**2 - 300**2) / 2 + (3**2 - 2**2) / 2)
    print(enthalpy(cp, temperature=(300, 600), eo=(2, 2)), (600**2 - 300**2) / 2 + (2**2 - 2**2) / 2)
