from mathematics import Constants

# термодинамические параметры
parameters = Constants(
    l="thermal_conductivity",  # теплопроводность
    C="heat_capacity",  # теплоемкость
    Cp="heat_capacity_pressure",  # теплокмкость при постоянном давлении
    Cv="heat_capacity_volume",  # теплокмкость при постоянном объеме
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
    a="sound_speed",  # скорость звука
    a_critical="critical_sound_speed",  # критическая скорость звука
    c="absolute_velocity",  # абсолютная скорость
    u="portable_velocity",  # переносная скорость
    w="relative_velocity",  # относительная скорость
    mf="mass_flow",  # массовый расход
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
    peff="pressure_efficiency",  # коэф. сохранения давления
    effburn="efficiency_burn",  # КПД полноты сгорания топлива
    power="power",  # мощность
)

EPSREL = 0.01  # относительная ошибка
NITER = 25  # количество итераций
