from mathematics import Constants

# термодинамические параметры
parameters = Constants(
    rf="rotation_frequency",  # частота вращения
    l="thermal_conductivity",  # теплопроводность
    hc="heat_capacity",  # теплоемкость
    hcp="heat_capacity_pressure",  # теплокмкость при постоянном давлении
    hcv="heat_capacity_volume",  # теплокмкость при постоянном объеме
    k="adiabatic_index",  # показатель адиабаты
    gc="gas_const",  # газовая постоянная
    eo="excess_oxidizing",  # коэффициент избытка окислителя
    m="mass",  # масса
    v="volume",  # объем
    power="power",  # мощность
    force="force",  # сила
    # статические термодинамические параметры
    T="static_temperature",  # статическая темпрература
    P="static_pressure",  # статическое давление
    D="staticdensity",  # статическая плотность
    # полные термодинамические параметры
    TT="total_temperature",  # полная температура
    PP="total_pressure",  # полное давление
    DD="total_density",  # полная плотность
    # скорости
    ss="sound_speed",  # скорость звука
    ss_critical="critical_sound_speed",  # критическая скорость звука
    c="absolute_velocity",  # абсолютная скорость
    u="portable_velocity",  # переносная скорость
    w="relative_velocity",  # относительная скорость
    # безразмерные параметры
    titi="total_temperature_ratio",  # степень повышения полной температуры
    ti="static_temperature_ratio",  # степень повышения статической температуры
    pipi="total_pressure_ratio",  # степень повышения полного давления
    pi="static_pressure_ratio",  # степень повышения статического давления
    didi="total_density_ratio",  # степень повышения полной плотности
    di="density_ratio",  # степень повышения статической плотности
    # числа
    Mach="mach_number",  # число Маха
    Nu="nusselt_number",  # число Нуссельта
    # КПД
    efficiency="efficiency",  # КПД
    effeff="total_efficiency",  # полный КПД
)

EPSREL = 0.01  # относительная ошибка
NITER = 25  # количество итераций

"""
Порядок расчета ТД параметров:
mf -> excess_oxidizing -> gas_const -> T* -> P* -> D* -> hcp -> k -> a* -> c
"""
