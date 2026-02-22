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
    # статические термодинамические параметры
    T="static_temperature",  # статическая темпрература
    P="static_pressure",  # статическое давление
    D="staticdensity",  # статическая плотность
    # полные термодинамические параметры
    TT="total_temperature",  # полная температура
    PP="total_pressure",  # полное давление
    DD="total_density",  # полная плотность
    m="mass",  # масса
    v="volume",  # объем
    # скорости
    ss="sound_speed",  # скорость звука
    ss_critical="critical_sound_speed",  # критическая скорость звука
    c="absolute_velocity",  # абсолютная скорость
    u="portable_velocity",  # переносная скорость
    w="relative_velocity",  # относительная скорость
    # безразмерные параметры
    pipi="total_pressure_ratio",  # степень повышения полного давления
    pi="static_pressure_ratio",  # степень повышения статического давления
    titi="total_temperature_ratio",  # степень повышения полной температуры
    ti="static_temperature_ratio",  # степень повышения статической температуры
    # числа
    mach="mach_number",  # число Маха
    Nu="nusselt_number",  # число Нуссельта
    # КПД
    efficiency="efficiency",  # КПД
    effeff="total_efficiency",  # полный кпд
    p_eff="pressure_efficiency",  # коэф. сохранения давления
    eff_burn="efficiency_burn",  # КПД полноты сгорания топлива
    s_eff="speed_efficiency",  # коэф. сохранения скорости
    power="power",  # мощность
)

EPSREL = 0.01  # относительная ошибка
NITER = 25  # количество итераций

"""
Порядок расчета ТД параметров:
mf -> excess_oxidizing -> gas_const -> T* -> P* -> D* -> hcp -> k -> a* -> c
"""
