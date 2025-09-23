import numpy as np
from scipy.integrate import nquad


def call_with_kwargs(function, kwargs: dict):
    """Вызов функции через kwargs"""
    assert callable(function), TypeError(f"{function} must be callable")
    assert isinstance(kwargs, dict), TypeError(f"type {kwargs} must be dict")

    arguments = {}
    for arg in function.__code__.co_varnames:
        value = kwargs.get(arg)
        assert value is not None, ValueError(f"value {value} must be not None")
        arguments[arg] = value

    return function(**arguments)


def integral_average(function, **kwargs) -> tuple[float, float]:
    """Среднее интегральное"""
    assert callable(function), TypeError(f"{function} must be callable")

    partial_args, other_args, devider = {}, {}, 1.0
    arg_names = function.__code__.co_varnames[: function.__code__.co_argcount]
    for arg_name in arg_names:
        rang = kwargs.get(arg_name)
        assert rang is not None, Exception(
            f"function {function} require argument {arg_name}"
        )
        assert isinstance(rang, (tuple, list, np.ndarray)), TypeError(
            f"type of range must be tuple, but has {type(rang)}"
        )
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


if __name__ == "__main__":

    def test_func(temperature, y):
        return temperature * y

    print(integral_average(test_func, temperature=(300, 600), y=(1, 1)))
