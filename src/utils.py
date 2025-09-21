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


def check_efficiency(efficiency: float) -> bool:
    """Проверка значения КПД"""
    return 0 <= efficiency <= 1


def check_temperature(temperature: float) -> bool:
    """Проверка на отрицательную абсолютную температуру"""
    return 0 <= temperature
