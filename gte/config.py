from mathematics import Constants
from matplotlib import pyplot as plt

# термодинамические параметры
parameters = Constants(
    l="thermal_conductivity",  # теплопроводность
    Cp="heat_capacity_at_constant_pressure",  # теплокмкость при постоянном давлении
    Cv="heat_capacity_at_constant_volume",  # теплокмкость при постоянном объеме
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

EPSREL = 0.001  # относительная ошибка
NITER = 25  # количество итераций


def show():
    w1 = 10
    h1 = 10 * w1

    fg, ax = plt.subplots()

    plt.title("Схема ГТД", fontsize=16)
    plt.grid(False)  # сетка
    plt.axis("square")

    x, y = w1, 0

    x += 2 * w1

    ax.add_artist(plt.Rectangle((x, y), w1, 1.5 * h1, fill=False, color=(0.5, 0.5, 0.5), lw=2))  # Вх
    x += 2 * w1
    y += 2 * w1

    ax.add_artist(plt.Rectangle((x, y), w1, 1.5 * h1 - 2 * w1, fill=False, color=(0, 1, 0), lw=2))  # КНД
    x += 2 * w1
    y += 2 * w1

    ax.add_artist(plt.Rectangle((x, y), w1, 1.5 * h1 - 4 * w1, fill=False, color=(1, 1, 0), lw=2))  # КСД
    x += 2 * w1
    y += 2 * w1

    ax.add_artist(plt.Rectangle((x, y), w1, 1.5 * h1 - 6 * w1, fill=False, color=(1, 0.5, 0), lw=2))  # КВД
    x += 2 * w1
    y += 2 * w1

    x += w1
    ax.add_artist(plt.Circle((x, y), w1, fill=False, color=(1, 0, 0), lw=2))  # КС
    ax.add_artist(plt.Rectangle((x - w1 / 2, y - 2.5 * w1), w1, w1, fill=False, color=(1, 0.25, 1), lw=2))  # РВД
    ax.add_artist(plt.Rectangle((x - w1 / 2, y - 4.5 * w1), w1, w1, fill=False, color=(1, 0.25, 1), lw=2))  # РСД
    ax.add_artist(plt.Rectangle((x - w1 / 2, y - 6.5 * w1), w1, w1, fill=False, color=(1, 0.25, 1), lw=2))  # РНД
    x += w1

    x += w1
    y -= 2 * w1
    ax.add_artist(plt.Rectangle((x, y), w1, h1 - 6 * w1, fill=False, color=(1, 0.5, 0), lw=2))  # ТВД
    x += 2 * w1
    ax.add_artist(plt.Rectangle((x, y - w1 / 2), w1, w1, fill=False, color=(1, 0.25, 1), lw=2))  # РВД
    x += 2.5 * w1
    ax.add_artist(plt.Circle((x, y), w1 / 2, fill=False, color=(0, 0, 1), lw=2))  # НВД
    x += 1.5 * w1
    y -= 2 * w1

    ax.add_artist(plt.Rectangle((x, y), w1, h1 - 4 * w1, fill=False, color=(1, 1, 0), lw=2))  # ТСД
    x += 2 * w1
    ax.add_artist(plt.Rectangle((x, y - w1 / 2), w1, w1, fill=False, color=(1, 0.25, 1), lw=2))  # РСД
    x += 2.5 * w1
    ax.add_artist(plt.Circle((x, y), w1 / 2, fill=False, color=(0, 0, 1), lw=2))  # НСД
    x += 1.5 * w1
    y -= 2 * w1

    ax.add_artist(plt.Rectangle((x, y), w1, h1 - 2 * w1, fill=False, color=(0, 1, 0), lw=2))  # ТНД
    x += 2 * w1
    ax.add_artist(plt.Rectangle((x, y - w1 / 2), w1, w1, fill=False, color=(1, 0.25, 1), lw=2))  # РНД
    x += 2.5 * w1
    ax.add_artist(plt.Circle((x, y), w1 / 2, fill=False, color=(0, 0, 1), lw=2))  # ННД
    x += 1.5 * w1
    y -= 2 * w1

    ax.add_artist(plt.Rectangle((x, y), w1, h1, fill=False, color=(0.5, 0.25, 1), lw=2))  # СТ
    x += 2 * w1
    ax.add_artist(plt.Rectangle((x, y - w1 / 2), w1, w1, fill=False, color=(1, 0.25, 1), lw=2))  # РСТ
    x += 2.5 * w1
    ax.add_artist(plt.Circle((x, y), w1 / 2, fill=False, color=(0, 0, 1), lw=2))  # НСТ
    x += 1.5 * w1

    ax.add_artist(plt.Rectangle((x, y), w1, 1.5 * h1, fill=False, color=(1, 0, 0), lw=2))  # ФК
    x += 2 * w1

    ax.add_artist(plt.Rectangle((x, y), w1, 1.5 * h1, fill=False, color=(0, 0, 0.7, 0.7), lw=2))  # С
    x += 2 * w1

    ax.add_artist(plt.Rectangle((x, y), w1, 1.5 * h1, fill=False, color=(0.5, 0.5, 0.5), lw=2))  # Вых
    x += 2 * w1

    x += 2 * w1

    # коричневое сопло
    plt.plot((3 * w1, x - 3 * w1), (h1 * 1.5, h1 * 1.5), color=(0, 0, 0), lw=2)  # II контур
    plt.plot((3 * w1, x - 3 * w1), (h1, h1), color=(0, 0, 0), lw=2)  # I контур
    plt.plot((0, x), (0, 0), ls="-.", color=(0, 0, 0), lw=1.5)  # ось
    plt.plot((10 * w1, 11.5 * w1), (6 * w1, 6 * w1), (12.5 * w1, 14 * w1), (6 * w1, 6 * w1), (15 * w1, 16 * w1), (6 * w1, 6 * w1), (17 * w1, 18 * w1), (6 * w1, 6 * w1), ls="-", color=(1, 0.5, 0))
    plt.plot((8 * w1, 11.5 * w1), (4 * w1, 4 * w1), (12.5 * w1, 20 * w1), (4 * w1, 4 * w1), (21 * w1, 22 * w1), (4 * w1, 4 * w1), (23 * w1, 24 * w1), (4 * w1, 4 * w1), ls="-", color=(1, 1, 0))
    plt.plot((6 * w1, 11.5 * w1), (2 * w1, 2 * w1), (12.5 * w1, 26 * w1), (2 * w1, 2 * w1), (27 * w1, 28 * w1), (2 * w1, 2 * w1), (29 * w1, 30 * w1), (2 * w1, 2 * w1), ls="-", color=(0, 1, 0))
    plt.xlim(0, x)
    plt.ylim(-2 * h1, 2 * h1)

    plt.show()
