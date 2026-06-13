import numpy as np
from numpy import isnan, nan
from numpy import log as ln
from scipy import interpolate

try:
    from .parameters import parameters as tdp  # Попытка относительного импорта
except ImportError:
    from parameters import parameters as tdp  # Резервный абсолютный импорт


np.seterr(invalid="ignore")  # игнорирование nan ошибок

T0 = 273.15  # Абсолютный ноль температуры
GAS_CONST = 8.314_462_618_153_24  # Универсальная газовая постоянная


def sonic_velocity(k: float, gas_const: float, temperature: float) -> float:
    """Cкорость звука"""
    return (k * gas_const * temperature) ** 0.5


def critical_sonic_velocity(k: float, gas_const: float, total_temperature: float) -> float:
    """Критическая скорость звука"""
    return (2 * k / (k + 1) * gas_const * total_temperature) ** 0.5


def gdf(
    parameter: str,
    equivalent_speed: float,
    adiabatic_index: float = None,
) -> float:
    """Газодинамические функции"""
    if not isinstance(parameter, str):
        raise TypeError(f"type {parameter} must be str")
    if not isinstance(equivalent_speed, (int, float, np.number)):
        raise TypeError(f"type {equivalent_speed} must be numeric")
    parameter = parameter.upper()
    if parameter == "T":
        if not isinstance(adiabatic_index, (int, float, np.number)):
            raise TypeError(f"type {adiabatic_index} must be numeric")
        return 1 - equivalent_speed**2 * ((adiabatic_index - 1) / (adiabatic_index + 1))
    elif parameter == "P":
        return gdf("T", equivalent_speed, adiabatic_index) ** (adiabatic_index / (adiabatic_index - 1))
    elif parameter == "D":
        return gdf("T", equivalent_speed, adiabatic_index) ** (1 / (adiabatic_index - 1))
    elif parameter in ("G", "MF"):
        return ((adiabatic_index + 1) / 2) ** (1 / (adiabatic_index - 1)) * equivalent_speed * gdf("D", equivalent_speed, adiabatic_index)
    elif parameter in ("I", "MV"):
        return equivalent_speed + 1 / equivalent_speed
    else:
        raise ValueError(f'{parameter} not in ("T", "P", "D", "G", "MF", "I", "MV")')


def temperature_atmosphere_standard(height: float) -> tuple[float, str]:
    """Статическая температура стандартной атмосферы"""
    return 288.15 - 0.00651 * height if height < 11_000 else 216.65, "K"


def pressure_atmosphere_standard(height: float) -> tuple[float, str]:
    """Статическое давление стандартной атмосферы"""
    return (101_325 * (temperature_atmosphere_standard(height)[0] / 288.15) ** 5.2533 if height < 11_000 else 22_699.9 * np.exp((11_000 - height) / 6318)), "Pa"


def atmosphere_standard(height: float) -> dict[str : tuple[float, str]]:
    """Атмосфера стандартная ГОСТ 4401-81"""
    return {
        tdp.t: temperature_atmosphere_standard(height),
        tdp.p: pressure_atmosphere_standard(height),
    }


def efficiency_polytropic(process="", pipi=nan, effeff=nan, k=nan) -> float:
    """Политропический КПД"""
    if pipi == 1:
        return 1
    if process.upper().startswith("C"):  # COMPRESSION
        return ((k - 1) / k) * ln(pipi) / ln((pipi ** ((k - 1) / k) - 1) / effeff + 1)
    if process.upper().startswith("E"):  # EXTENSION
        return -(k / (k - 1)) / ln(pipi) * ln(effeff * (1 / (pipi ** ((k - 1) / k)) - 1) + 1)
    raise Exception(f'{process} not in ("C", "E")')


def chemical_formula_to_dict(formula: str) -> dict[str:int]:
    """Разбор хмической формулы поэлементно"""
    result = {}

    i = 0
    while i < len(formula):
        if i + 1 < len(formula) and formula[i + 1].islower():
            atom = formula[i : i + 2]
            i += 2
        else:
            atom = formula[i]
            i += 1

        count = 0
        while i < len(formula) and formula[i].isdigit():
            count = count * 10 + int(formula[i])
            i += 1

        if count == 0:
            count = 1

        if atom in result:
            result[atom] += count
        else:
            result[atom] = count

    return result


def adiabatic_index(gas_const: int | float, cp: int | float) -> float:
    """Показатель адиабаты"""
    if cp == gas_const:
        return nan
    return cp / (cp - gas_const)


def gas_const(substance: str) -> float:
    """Газовая постоянная (Дж/кг/К)"""
    if not isinstance(substance, str):
        raise TypeError(f"type {substance} must be str")

    if substance.upper() == "AIR":
        return 287.14
    elif substance == "N2":
        return 297
    elif substance == "NH3":
        return 488.5
    elif substance == "Ar":
        return 208.2
    elif substance == "H2":
        return 4118.2
    elif substance == "He":
        return 2078.2
    elif substance == "O2":
        return 260
    elif substance == "Kr":
        return 99.3
    elif substance == "Xe":
        return 63.4
    elif substance == "Ne":
        return 412.2
    elif substance == "CO2":
        return 189
    else:
        raise ValueError(f"{substance} not found")


def gas_const_exhaust(
    excess_oxidizing: int | float,
    stoichiometry: int | float,
    composition: dict[str:float],
) -> float:
    """Газовая постоянная продуктов сгорания (Дж/кг/К)"""
    if not isinstance(excess_oxidizing, (int, float, np.number)):
        raise TypeError(f"type {excess_oxidizing} must be numeric")
    if not isinstance(stoichiometry, (int, float, np.number)):
        raise TypeError(f"type {stoichiometry} must be numeric")
    if not isinstance(composition, dict):
        TypeError(f"type {composition} must be dict")

    H2O = composition.get("H2O", 0)  # массовая доля волы в смеси
    result = excess_oxidizing * stoichiometry
    result += 28.966 * (composition.get("H2", 0) / 4.032 + composition.get("O2", 0) / 32) * (1 - H2O)
    result += H2O * excess_oxidizing * stoichiometry * 0.60779
    result *= 29.27
    result /= 1 - H2O + excess_oxidizing * stoichiometry
    return result * 9.80665  # СИ


def gas_const_exhaust_fuel(excess_oxidizing: int | float, fuel: str) -> float:
    """Газовая постоянная продуктов сгорания (Дж/кг/К)"""
    if not isinstance(excess_oxidizing, (int, float, np.number)):
        raise TypeError(f"type {excess_oxidizing} must be numeric")
    if not isinstance(fuel, str):
        raise TypeError(f"type {fuel} must be str")
    fuel = fuel.upper()
    if fuel in ("C2H8N2", "KEROSENE"):
        return 288.1954313 + 0.691695880 / excess_oxidizing if 0 < excess_oxidizing else nan
    elif fuel in ("T1", "РЕАКТИВНОЕ ТОПЛИВО"):
        return 288.1856907 - 0.213996376 / excess_oxidizing if 0 < excess_oxidizing else nan
    elif fuel.upper() in ("TC1", "TC-1"):
        return 288.1901130 + 0.196844775 / excess_oxidizing if 0 < excess_oxidizing else nan
    elif fuel == "DIESEL":
        return 288.1792281 - 0.813445523 / excess_oxidizing if 0 < excess_oxidizing else nan
    elif fuel in ("SOLAR", "SOLAR OIL", "SOLAR_OIL"):
        return 288.1729604 - 1.432301634 / excess_oxidizing if 0 < excess_oxidizing else nan
    elif fuel == "MAZUT":
        return 288.1635509 - 2.254443192 / excess_oxidizing if 0 < excess_oxidizing else nan
    elif fuel in ("ПРИРОДНЫЙ ГАЗ", "ПРИРОДНЫЙ_ГАЗ"):
        return 290.0288864 + 12.207960640 / excess_oxidizing if 0 < excess_oxidizing else nan
    elif fuel in ("КОКСОВЫЙ ГАЗ", "КОКСОВЫЙ_ГАЗ"):
        return 288.4860344 + 20.146254880 / excess_oxidizing if 0 < excess_oxidizing else nan
    elif fuel == "BIOGAS":
        return 289.9681764 + 2.138861745 / excess_oxidizing if 0 < excess_oxidizing else nan
    else:
        raise ValueError(f"{fuel} not found")


def stoichiometry(fuel: str) -> float:
    """Стехиометрический коэффициент []"""
    if not isinstance(fuel, str):
        raise TypeError(f"type {fuel} must be str")

    fuel = fuel.upper()
    if fuel in ("C2H8N2", "KEROSENE", "T-1", "T-2", "TC-1", "TC1"):
        return 14.61
    elif fuel in ("PETROL", "GASOLINE"):
        return 14.91
    elif fuel in ("SOLAR", "SOLAR OIL", "SOLAR_OIL", "DIESEL"):
        return 14.35
    elif fuel in ("MAZUT", "Ф5", "Ф12"):
        return 13.31
    elif fuel in ("ПРИРОДНЫЙ ГАЗ", "ПРИРОДНЫЙ_ГАЗ"):
        return np.mean(
            tuple(
                {
                    "Березовский": 15.83,
                    "Войвожский": 13.69,
                    "Дашавский": 16.84,
                    "Карадагеказский": 16.51,
                    "Ленинградский": 15.96,
                    "Саратовский": 16.34,
                    "Ставропольский": 16.85,
                    "Щебеленский": 15.93,
                }.values()
            )
        )
    elif fuel in ("КОКСОВЫЙ ГАЗ", "КОКСОВЫЙ_ГАЗ"):
        return 9.908
    else:
        raise ValueError(f"{fuel} not found")


cp_kerosene = interpolate.interp1d(
    (293.15, 303.15, 313.15, 323.15, 333.15, 343.15, 353.15, 363.15, 373.15, 383.15, 393.15, 403.15, 413.15, 423.15, 433.15, 443.15, 453.15, 463.15, 473.15, 483.15, 493.15, 503.15, 513.15, 523.15, 533.15, 543.15),
    (2000, 2040, 2090, 2140, 2180, 2230, 2280, 2330, 2380, 2430, 2480, 2530, 2580, 2630, 2680, 2730, 2790, 2840, 2890, 2940, 3000, 3050, 3110, 3160, 3210, 3260),
    kind=2,
    fill_value="extrapolate",
)


def heat_capacity(substance: str, temperature) -> float:
    """Теплоемкость (Дж/кг/К)"""
    if substance.upper() in ("C2H8N2", "KEROSENE", "TC-1"):
        """Теплоемкость жидкого керосина"""
        return cp_kerosene(temperature) if -55 + T0 <= temperature <= 400 + T0 else nan  # температуры замерзания и воспламенения
    else:
        raise ValueError(f"{substance} not found")


def heat_capacity_p(substance: str, temperature: int | float | np.number) -> float:
    """Теплоемкость при постоянном давлении (Дж/кг/К)"""
    if not isinstance(substance, str):
        raise TypeError(f"type {substance} must be str")

    if substance.upper() == "AIR":
        """Теплоемкость воздуха [PTM 1677-83]"""
        t_1000 = temperature / 1000
        coefs = (0.2521923, -0.1186612, 0.3360775, -0.3073812, 0.1382207, -0.03090246, 0.002745383)
        return 4187 * sum(coef * t_1000**i for i, coef in enumerate(coefs)) if 0 <= temperature <= 3_000 else nan
    elif substance == "CO2":
        # PTM 1677-83
        t_1000 = temperature / 1000
        coefs = (0.1047056, 0.4234367, -0.3953561, 0.2249471, -0.07729786, 0.01462170, -0.001166819)
        return 4187 * sum(coef * t_1000**i for i, coef in enumerate(coefs))
    elif substance == "H2O":
        # PTM 1677-83
        t_1000 = temperature / 1000
        coefs = (0.4489375, -0.1088401, 0.4027652, -0.2638393, 0.07993751, -0.0115716, 0.0006241951)
        return 4187 * sum(coef * t_1000**i for i, coef in enumerate(coefs))
    elif substance == "O2":
        # PTM 1677-83
        t_1000 = temperature / 1000
        coefs = (0.2083632, -0.0112279, 0.2235868, -0.2732668, 0.1461334, -0.03687021, 0.003584204)
        return 4187 * sum(coef * t_1000**i for i, coef in enumerate(coefs))
    elif substance == "H2":
        # PTM 1677-83
        t_1000 = temperature / 1000
        coefs = (3.070881, 2.230734, -4.909, 5.321652, -2.756533, 0.6851526, -0.06596988)
        return 4187 * sum(coef * t_1000**i for i, coef in enumerate(coefs))
    elif substance == "N2":
        # https://www.highexpert.ru/content/gases/nitrogen.html
        return 1041 - 0.021 * temperature + 0.0003814 * temperature**2
    elif substance == "Ar":  # TODO
        return 523
    elif substance == "Ne":  # TODO
        return 1038
    else:
        raise ValueError(f"{substance} not found")


def heat_capacity_p_exhaust_eo1(
    temperature: int | float | np.number,
    composition: dict[str:float] = None,
) -> float:
    """
    Условная теплоемкость выхлопа при коэффициенте избытка воздуха = 1 (Дж/кг/К) [PTM 1677-83]

    Истинная теплоемкость выхлопа (Жд/кг/К) считается как:
    H2O = composition.get('H2O', 0)  # массовая доля волы в смеси
    result = (1 - H2O) * (условная_теплоемкость + теплоемкость_окислителя * excess_oxidizing * stoichiometry)
    result += H2O * excess_oxidizing * stoichiometry * теплоемкость_воды
    result /= (1 - H2O + excess_oxidizing * stoichiometry)
    """
    if not isinstance(temperature, (int, float, np.number)):
        raise TypeError(f"type {temperature} must be numeric")

    if composition is None:  # as default composition = {"C": 0.85, "H2": 0.15, "O2": 0, "H2O": 0}
        t_1000 = temperature / 1000
        coefs = (0.2079764, 1.211806, -1.464097, 1.291195, -0.6385396, 0.1574277, -0.01518199)  # speedup
        return 4187 * sum(coef * t_1000**i for i, coef in enumerate(coefs))
    elif isinstance(composition, dict):
        result = composition.get("C", 0) / 12.01 * (44.01 * heat_capacity_p("CO2", temperature) - 32.0 * heat_capacity_p("O2", temperature))
        result += composition.get("H2", 0) / 1.008 * (9.008 * heat_capacity_p("H2O", temperature) - 8.0 * heat_capacity_p("O2", temperature))
        result += composition.get("O2", 0) * heat_capacity_p("O2", temperature)
        return result
    else:
        raise TypeError(f"type {composition} must be dict[str:float]")


def heat_capacity_p_exhaust(
    heat_capacity_p_exhaust_eo1: float,  # условная теплоемкость выхлопных газов
    heat_capacity_p_oxidizing: float,  # теплоемкость окислителя при постоянном давлении
    heat_capacity_H20: float,  # теплоемкость пара при постоянном давлении
    excess_oxidizing: float,  # коэффициент избытка окислителя
    stoichiometry: float,  # стехиометрический коэффициент
    H2O: float,  # массовая доля волы в смеси
):
    """
    Теплоемкость выхлопных газов при постоянном давлении с учетом влажности.

    Формула расчета:
        hcp = (hcp_dry + hcp_wet) / (1 - H2O + excess_oxidizing * stoichiometry)

    где:
        - hcp_dry = (1 - H2O) * (heat_capacity_p_exhaust_eo1 + heat_capacity_p_oxidizing * excess_oxidizing * stoichiometry)
        - hcp_wet = H2O * excess_oxidizing * stoichiometry * heat_capacity_H20

    Returns:
        Теплоемкость выхлопных газов (float)

    Raises:
        ZeroDivisionError: Если знаменатель равен нулю
        ValueError: Если входные параметры имеют некорректные значения
    """
    if not (0 <= H2O <= 1):
        raise ValueError(f"{H2O=} not in [0, 1]")
    if excess_oxidizing < 0:
        raise ValueError(f"{excess_oxidizing=} must be >= 0")
    if stoichiometry <= 0:
        raise ValueError(f"{stoichiometry=} must be > 0")

    # Сухой состав
    hcp_dry = (1 - H2O) * (heat_capacity_p_exhaust_eo1 + heat_capacity_p_oxidizing * excess_oxidizing * stoichiometry)
    # Влажный состав
    hcp_wet = H2O * excess_oxidizing * stoichiometry * heat_capacity_H20
    # Перерасчет на сухой состав
    denominator = 1 - H2O + excess_oxidizing * stoichiometry

    hcp = (hcp_dry + hcp_wet) / denominator
    return hcp


def lower_heat(fuel: str) -> float:
    """Низшая теплота сгорания горючего при коэффициенте избытка окислителя = 1"""
    if not isinstance(fuel, str):
        raise TypeError(f"type {fuel} must be str")
    fuel = fuel.upper()
    if fuel in ("C2H8N2", "KEROSENE", "TC1", "TC-1", "PETROL"):
        return 0.5 * (43_600_000 + 42_700_000)
    elif fuel in ("T-6", "T-8"):
        return 42_900_000
    elif fuel == "DIESEL":
        return nan
    elif fuel in ("ПРИРОДНЫЙ ГАЗ",):
        return nan
    else:
        raise ValueError(f"'{fuel}' not found")


def dynamic_viscosity(
    substance: str,
    temperature: int | float | np.number,
    excess_oxidizing: int | float | np.number = nan,
) -> float:
    """Динамическая вязкость (Па*с)"""
    if not isinstance(substance, str):
        raise TypeError(f"type {substance} must be str")
    if not isinstance(temperature, (int, float, np.number)):
        raise TypeError(f"type {temperature} must be numeric")

    if substance.upper() == "EXHAUST":
        if isnan(excess_oxidizing):
            raise ValueError(f"{excess_oxidizing} must be numeric")
        coefs = (+0.505, +4.849, -1.333, +0.229)
        t_1000 = temperature / 1_000
        return 10 ** (-5) * (sum(coef * t_1000**i for i, coef in enumerate(coefs)) - 0.275 / excess_oxidizing)
    elif substance == "N2":
        return 16.67 / 10**6 * (temperature / T0) ** 0.68
    elif substance == "NH3":
        return 9.7 / 10**6 * (temperature / T0) ** 1.06
    elif substance == "Ar":
        return 21.08 / 10**6 * (temperature / T0) ** 0.72
    elif substance == "H2":
        return 8.36 / 10**6 * (temperature / T0) ** 0.68
    elif substance.upper() == "AIR":
        return 17.16 / 10**6 * (temperature / T0) ** 0.68
    elif substance == "He":
        return 18.64 / 10**6 * (temperature / T0) ** 0.68
    elif substance == "O2":
        return 19.42 / 10**6 * (temperature / T0) ** 0.69
    elif substance == "Kr":
        return 23.44 / 10**6 * (temperature / T0) ** 0.83
    elif substance == "Xe":
        return 21.08 / 10**6 * (temperature / T0) ** 0.89
    elif substance == "Ne":
        return 29.71 / 10**6 * (temperature / T0) ** 0.65
    elif substance == "CO2":
        return 13.65 / 10**6 * (temperature / T0) ** 0.82
    else:
        raise ValueError(f"'{substance}' not found")


def thermal_conductivity(substance: str, temperature: int | float | np.number) -> float:
    """Теплопроводность (Вт/м/К)"""
    if not isinstance(substance, str):
        raise TypeError(f"type {substance} must be str")
    if not isinstance(temperature, (int, float, np.number)):
        raise TypeError(f"type {temperature} must be numeric")

    if substance == "N2":
        return 241.9 / 10**4 * (temperature / T0) ** 0.8
    elif substance == "NH3":
        return 212 / 10**4 * (temperature / T0) ** 1.53
    elif substance == "Ar":
        return 165.1 / 10**4 * (temperature / T0) ** 0.8
    elif substance == "H2":
        return 1721.2 / 10**4 * (temperature / T0) ** 0.78
    elif substance.upper() == "AIR":
        return 244.2 / 10**4 * (temperature / T0) ** 0.82
    elif substance == "He":
        return 1425.8 / 10**4 * (temperature / T0) ** 0.73
    elif substance == "O2":
        return 245.4 / 10**4 * (temperature / T0) ** 0.87
    elif substance == "Kr":
        return 88.9 / 10**4 * (temperature / T0) ** 0.86
    elif substance == "Xe":
        return 52.3 / 10**4 * (temperature / T0) ** 0.91
    elif substance == "Ne":
        return 464 / 10**4 * (temperature / T0) ** 0.71
    elif substance == "CO2":
        return 147 / 10**4 * (temperature / T0) ** 1.23
    else:
        raise ValueError(f"'{substance}' not found")


def g_cool_ciam(temperature_input, temperature_output, temperature_lim) -> float:
    """Эмпирический относительный массовый расход на охлаждение
    по max температурам до и после КС и допустимой температуре,
    отнесенный к расходу на входе в горячую часть"""
    θ = (temperature_output - temperature_lim) / (temperature_output - temperature_input)
    g_cool = 0.059 * θ / (1 - 1.42 * θ)
    return g_cool if g_cool > 0 else 0


def g_cool_bmstu(temperature_max, temperature_lim=1000) -> float:
    """Эмпирический относительный массовый расход на охлаждение по max и допустимой температуре,
    отнесенный к расходу на входе в горячую часть"""
    coefs = (0.01, 0.09, 0.2, 0.16)
    dt = (temperature_max - temperature_lim) / 1_000
    g_cool = sum(coefs[i] * dt**i for i in range(len(coefs)))
    return g_cool if g_cool > 0 else 0
