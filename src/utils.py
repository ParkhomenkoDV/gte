from mathematics import integral_average


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


def integral_average_fix(function, *ranges):
    """Среднеинтегральное значение при одинаковых границах интегрирования"""
    assert callable(function), TypeError(f"{function} must be callable")

    def wrapper(*ranges):
        return function(*[])


def integral_average_kwargs(function, kwargs: dict):
    """"""
    assert callable(function), TypeError(f"{function} must be callable")
    assert isinstance(kwargs, dict), TypeError(f"type {kwargs} must be dict")

    ranges = []
    for arg in function.__code__.co_varnames:
        range_ = kwargs.get(arg)
        assert range_ is not None, ValueError(f"value {range_} must be not None")
        assert isinstance(range_, (tuple, list)), TypeError(
            f"type range must be tuple, not {type(range_)}"
        )
        ranges.append(range_)

    function_ = integral_average_fix(function, *ranges)

    return integral_average(function_, *ranges)
