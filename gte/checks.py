import inspect
from typing import Callable, Set


def check_characteristic(function: Callable, arguments: Set[str]) -> None:
    """Проверка характеристики"""
    if not callable(function):
        raise TypeError(f"{function} must be callable")
    if not isinstance(arguments, set):
        raise TypeError(f"{type(arguments)=} must be set")

    signature = inspect.signature(function)
    keys = set(signature.parameters.keys())
    difference = keys - arguments
    if len(difference) > 0:  # в функции есть аргументы которых нет в разрешенных arguments
        raise KeyError(f"function {function} has not arguments: {difference}")


def check_mass_flow(mass_flow: float) -> bool:
    """Проверка на отрицательный массовый расход"""
    return 0 < mass_flow


def check_temperature(temperature: float) -> bool:
    """Проверка на отрицательную абсолютную температуру"""
    return 0 <= temperature


def check_pressure(pressure: float) -> bool:
    """Проверка на отрицательное давление"""
    return 0 <= pressure


def check_excess_oxidizing(excess_oxidizing: float) -> bool:
    """Проверка коэф. избытка окислителя"""
    return 0 <= excess_oxidizing


def check_efficiency(efficiency: float) -> bool:
    """Проверка значения КПД"""
    return 0 <= efficiency <= 1


def check_pressure_ratio(pressure_ratio: float) -> bool:
    """Проверка степени повышения/понижения давления"""
    return 1 <= pressure_ratio


def check_power(power: float) -> bool:
    """Проверка на отрицательную мощность"""
    return 0 <= power
