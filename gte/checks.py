def check_efficiency(efficiency: float) -> bool:
    """Проверка значения КПД"""
    return 0 <= efficiency <= 1


def check_temperature(temperature: float) -> bool:
    """Проверка на отрицательную абсолютную температуру"""
    return 0 <= temperature


def check_pressure_ratio(pressure_ratio: float) -> bool:
    """Проверка степени повышения/понижения давления"""
    return 1 <= pressure_ratio
