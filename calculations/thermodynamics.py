import sys
import os
from tqdm import tqdm
from functools import lru_cache

import pandas as pd
import numpy as np
from numpy import nan, isnan, log
from scipy import interpolate, integrate
import matplotlib.pyplot as plt

import mendeleev  # ПСХЭ Менделеева

sys.path.append('D:/Programming/Python')

from decorators import timeit

from tools import isnum, av, eps
from colorama import Fore

np.seterr(invalid='ignore')  # игнорирование ошибок с nan

T0 = 273.15  # Абсолютный ноль температуры
gas_const = 8.314_462_618_153_24  # Универсальная газовая постоянная
FUELS = ('КЕРОСИН', 'KEROSENE', 'ТС-1',
         'БЕНЗИН', 'PETROL', 'GASOLINE',
         'ДИЗЕЛЬ', 'DIESEL',
         'ПРИРОДНЫЙ ГАЗ', 'ПРИРОДНЫЙ_ГАЗ', 'ПГ',
         'КОКСОВЫЙ ГАЗ', 'КОКСОВЫЙ_ГАЗ', 'КС')
OXIDIZERS = ('ВОЗДУХ', 'AIR',
             'КИСЛОРОД', 'OXYGEN', 'O2')


def GDF(what: str, λ: float = nan, k: float = nan) -> float:
    """Газодинамические функции"""
    if what in ('T', 'τ'):
        return 1 - λ ** 2 * ((k - 1) / (k + 1))
    elif what in ('P', 'π'):
        return GDF('T', λ=λ, k=k) ** (k / (k - 1))
    elif what in ('ρ', 'ε'):
        return GDF('T', λ=λ, k=k) ** (1 / (k - 1))
    elif what in ('G', 'q'):
        return ((k + 1) / 2) ** (1 / (k - 1)) * λ * GDF('ρ', λ=λ, k=k)
    elif what in ('p', 'mv'):
        return λ + 1 / λ
    else:
        raise 'No name' + str(what)


file_path = 'libraries/Атмосфера стандартная.xlsx'
try:
    EXCEL_atmosphere_standard = pd.read_excel(file_path, header=None)
    T_atmosphere_standard = interpolate.interp1d(EXCEL_atmosphere_standard[0].iloc[1:],
                                                 EXCEL_atmosphere_standard[1].iloc[1:],
                                                 kind='linear', bounds_error=False)
    P_atmosphere_standard = interpolate.interp1d(EXCEL_atmosphere_standard[0].iloc[1:],
                                                 EXCEL_atmosphere_standard[2].iloc[1:],
                                                 kind='linear', bounds_error=False)
    del EXCEL_atmosphere_standard
except IOError as exception:
    print(Fore.RED + f'file "{file_path}" did not found!' + Fore.RESET)
    T_atmosphere_standard = lambda H: 288.15 - 0.00651 * H if H < 11_000 else 216.65
    P_atmosphere_standard = lambda H: 101_325 * (T_atmosphere_standard(H) / 288.15) ** 5.2533 if H < 11_000 \
        else 22_699.9 * np.exp((11_000 - H) / 6318)


@lru_cache(maxsize=None)
def atmosphere_standard(H: int | float) -> dict[str:float]:
    """Атмосфера стандартная ГОСТ 4401-81"""
    return {'T': float(T_atmosphere_standard(H)), 'P': float(P_atmosphere_standard(H))}


'''
def Cp_air(T1, T2, P1, P2):
    """Теплоемкость воздуха"""
    if EXCEL_Cp_air.T[0].iloc[1] <= T1 <= EXCEL_Cp_air.T[0].iloc[-1] and \
            EXCEL_Cp_air.T[0].iloc[1] <= T2 <= EXCEL_Cp_air.T[0].iloc[-1] and \
            EXCEL_Cp_air[0].iloc[1] <= P1 <= EXCEL_Cp_air[0].iloc[-1] and \
            EXCEL_Cp_air[0].iloc[1] <= P2 <= EXCEL_Cp_air[0].iloc[-1]:
        return av(Cp_air_, T1, T2, P1, P2)
    else:
        print('P')
        return integrate.quad(lambda T: Cp_air__(T), T1, T2)[0] / (T2 - T1)
'''


def η_polytropic(what='', ππ=nan, ηη=nan, k=nan) -> float:
    """Политропический КПД"""
    if ππ == 1:
        return 1
    else:
        if what.upper() in ('С', 'СЖ', 'СЖАТИЕ', 'C', 'COMPRESSION'):
            return ((k - 1) / k) * log(ππ) / log((ππ ** ((k - 1) / k) - 1) / ηη + 1)
        if what.upper() in ('Р', 'РАС', 'РАСШИРЕНИЕ', 'E', 'EXTENSION'):
            return -(k / (k - 1)) / log(ππ) * log(ηη * (1 / (ππ ** ((k - 1) / k)) - 1) + 1)
        print(Fore.RED + 'No find what is', Fore.RED + what, Fore.RED + 'in η_polytropic')
        return nan


def R_gas(substance, a_ox=nan, fuel='', **kwargs) -> float:
    """Газовая постоянная [Дж/кг/К]"""
    if substance.upper() in ('AIR', 'ВОЗДУХ'):
        """Газовая постоянная воздуха"""
        return 287.14
    elif substance.upper() in ('EXHAUST', 'ВЫХЛОП') and a_ox is not nan and fuel != '':
        """Газовая постоянная продуктов сгорания"""
        if fuel.upper() in ('KEROSENE',
                            'КЕРОСИН'): return 288.1954313 + 0.691695880 / a_ox
        if fuel.upper() in ('T1',
                            'Т1', 'РЕАКТИВНОЕ ТОПЛИВО'): return 288.1856907 - 0.213996376 / a_ox
        if fuel.upper() in ('TC1', 'ТС-1',
                            'ТС1', 'TC-1'): return 288.1901130 + 0.196844775 / a_ox
        if fuel.upper() in ('DIESEL', 'DT',
                            'ДИЗЕЛЬ', 'ДТ', 'ДИЗЕЛЬНОЕ ТОПЛИВО', 'ДИЗЕЛЬНОЕ_ТОПЛИВО'):
            return 288.1792281 - 0.813445523 / a_ox
        if fuel.upper() in ('SOLAR', 'SOLAR OIL', 'SOLAR_OIL',
                            'СОЛЯРКА', 'СОЛЯРОВОЕ МАСЛО', 'СОЛЯРОВОЕ_МАСЛО'): return 288.1729604 - 1.432301634 / a_ox
        if fuel.upper() in ('MAZUT',
                            'МАЗУТ', 'Ф12'): return 288.1635509 - 2.254443192 / a_ox
        if fuel.upper() in ('ПРИРОДНЫЙ ГАЗ', 'ПРИРОДНЫЙ_ГАЗ'): return 290.0288864 + 12.207960640 / a_ox
        if fuel.upper() in ('КОКСОВЫЙ ГАЗ', 'КОКСОВЫЙ_ГАЗ'): return 288.4860344 + 20.146254880 / a_ox
        if fuel.upper() in ('BIOGAS',
                            'БИОГАЗ'): return 289.9681764 + 2.138861745 / a_ox
    else:
        print(Fore.RED + 'fuel not found!' + ' in function ' + R_gas.__name__)
        return nan


def l_stoichiometry(fuel: str) -> float:
    """Стехиометрический коэффициент []"""
    if fuel.upper() in ('KEROSENE', 'T-1', 'T-2', 'TC-1', 'ТС1',
                        'КЕРОСИН', 'Т-1', 'Т-2', 'ТС-1', 'TC1'):
        return 14.61
    elif fuel.upper() in ('PETROL', 'GASOLINE', 'БЕНЗИН'):
        return 14.91
    elif fuel.upper() in ('SOLAR', 'SOLAR OIL', 'SOLAR_OIL',
                          'СОЛЯРКА', 'СОЛЯРОВОЕ МАСЛО', 'СОЛЯРОВОЕ_МАСЛО',
                          'DIESEL',
                          'ДИЗЕЛЬ'):
        return 14.35
    elif fuel.upper() in ('MAZUT',
                          'МАЗУТ', 'Ф5', 'Ф12'):
        return 13.31
    elif fuel.upper() in ('ПРИРОДНЫЙ ГАЗ', 'ПРИРОДНЫЙ_ГАЗ'):
        return np.mean({'Березовский': 15.83,
                        'Войвожский': 13.69,
                        'Дашавский': 16.84,
                        'Карадагеказский': 16.51,
                        'Ленинградский': 15.96,
                        'Саратовский': 16.34,
                        'Ставропольский': 16.85,
                        'Щебеленский': 15.93}.values())
    elif fuel.upper() in ('КОКСОВЫЙ ГАЗ', 'КОКСОВЫЙ_ГАЗ'):
        return 9.908
    else:
        print(Fore.RED + 'Fuel not found!' + ' in function ' + l_stoichiometry.__name__)
        return nan


EXCEL_Cp_air = pd.read_excel("libraries/Теплоёмкость воздуха.xlsx", header=None)
Cp_air = interpolate.RectBivariateSpline(EXCEL_Cp_air.T[0].iloc[1:],  # T
                                         EXCEL_Cp_air[0].iloc[1:],  # P
                                         EXCEL_Cp_air.T.iloc[1:, 1:],  # Cp_air
                                         kx=1, ky=1)  # x: linear, y: linear
del EXCEL_Cp_air

EXCEL_Cp_clean_kerosene = pd.read_excel('libraries/Чистая теплоемкость керосина.xlsx', header=None)
Cp_clean_kerosene = interpolate.interp1d(EXCEL_Cp_clean_kerosene[0].iloc[1:],
                                         EXCEL_Cp_clean_kerosene[1].iloc[1:],
                                         kind='linear', fill_value='extrapolate')
del EXCEL_Cp_clean_kerosene

EXCEL_Cp_clean_diesel = pd.read_excel('libraries/Чистая теплоемкость дизеля.xlsx', header=None)
Cp_clean_diesel = interpolate.interp1d(EXCEL_Cp_clean_diesel[0].iloc[1:],
                                       EXCEL_Cp_clean_diesel[1].iloc[1:],
                                       kind='linear', bounds_error=False)
del EXCEL_Cp_clean_diesel

# TODO предупреждение об экстраполяции
EXCEL_Cp_kerosene = pd.read_excel('libraries/Теплоёмкость жидкого керосина.xlsx', header=None)
Cp_kerosene = interpolate.interp1d(EXCEL_Cp_kerosene[1].iloc[1:],
                                   EXCEL_Cp_kerosene[2].iloc[1:],
                                   kind='linear', fill_value='extrapolate')
del EXCEL_Cp_kerosene


def Cp(substance, T=nan, P=nan, a_ox=nan, fuel='', **kwargs) -> float:
    """Теплоемкость при постоянном давлении"""
    if substance.upper() in ('AIR', 'ВОЗДУХ'):
        """Теплоемкость воздуха"""
        if P is not nan and 1 == 0:  # TODO интерполяция по поверхности
            return Cp_air(T, P)[0][0]
        else:
            # PTM 1677-83
            _T = T / 1000
            coefs = (0.2521923, -0.1186612, 0.3360775, -0.3073812, 0.1382207, -0.03090246, 0.002745383)
            return 4187 * sum([coef * _T ** i for i, coef in enumerate(coefs)])

    if substance.upper() in ('ЧИСТЫЙ ВЫХЛОП', 'ЧИСТЫЙ_ВЫХЛОП', 'CLEAN_EXHAUST', 'CLEAN EXHAUST') or \
            substance.upper() in ('EXHAUST', 'ВЫХЛОП') and a_ox == 1:
        """Чистая теплоемкость"""
        if fuel.upper() in ('KEROSENE', 'КЕРОСИН', 'ТС-1', 'PETROL', 'БЕНЗИН'):
            return Cp_clean_kerosene(T)
        elif fuel.upper() in ('ДИЗЕЛЬ', 'DIESEL'):
            return Cp_clean_diesel(T)
        else:
            print(Fore.RED + 'Not found substance' + ' in function ' + Cp.__name__)
            return nan
    if substance.upper() in ('EXHAUST', 'ВЫХЛОП'):
        return ((1 + l_stoichiometry(fuel)) * Cp('EXHAUST', T=T, a_ox=1, fuel=fuel) +
                (a_ox - 1) * l_stoichiometry(fuel) * Cp('AIR', T=T)) / (1 + a_ox * l_stoichiometry(fuel))
    if substance.upper() in ('KEROSENE', 'TC-1',
                             'КЕРОСИН', 'ТС-1'):
        return Cp_kerosene(T)

    print(Fore.RED + 'Not found substance' + ' in function ' + Cp.__name__)
    return nan


# print(Cp('AIR', T=75, P=10000000))


def Qa1(fuel) -> float:
    """Низшая теплота сгорания горючего при коэффициенте избытка окислителя = 1"""
    if fuel.upper() in ('KEROSENE', 'TC1', 'TC-1', 'PETROL',
                        'КЕРОСИН', 'ТС1', 'ТС-1', 'БЕНЗИН'):
        return 0.5 * (43_600 + 42_700) * 1000
    elif fuel.upper() in ('T-6', 'T-8',
                          'Т-6', 'Т-8'):
        return 42_900 * 1000
    elif fuel.upper in ('ДИЗЕЛЬ', 'DIESEL'):
        return nan
    elif fuel.upper in ('ПРИРОДНЫЙ ГАЗ', 'ПРИРОДНЫЙ_ГАЗ'):
        return nan
    elif fuel.upper in ('КОКСОВЫЙ ГАЗ', 'КОКСОВЫЙ_ГАЗ'):
        return nan
    else:
        print(Fore.RED + 'not found fuel!' + ' in function ' + Qa1.__name__ + Fore.RESET)
        return nan


def viscosity(substance: str, T=nan, a_ox=nan) -> float:  # dynamic viscosity -> kinematic viscosity добавить
    """Динамическая вязкость"""
    if substance.upper() in ('EXHAUST', 'ВЫХЛОП'):
        return 10 ** (-5) * \
            (0.229 * (T / 1000) ** 3 - 1.333 * (T / 1000) ** 2 + 4.849 * (T / 1000) + 0.505 - 0.275 / a_ox)
    else:
        print(Fore.RED + 'Not found substance' + ' in function ' + Cp.__name__)
        return nan


def g_cool_CIAM(TT1, TT2, T_lim) -> float:
    """Эмпирический относительный массовый расход на охлаждение
    по max температурам до и после КС и допустимой температуре,
    отнесенный к расходу на входе в горячую часть"""
    θ = (TT2 - T_lim) / (TT2 - TT1)
    g_cool = 0.059 * θ / (1 - 1.42 * θ)
    return g_cool if g_cool > 0 else 0


def g_cool_BMSTU(T, T_lim=1000) -> float:
    """Эмпирический относительный массовый расход на охлаждение по max и допустимой температуре,
    отнесенный к расходу на входе в горячую часть"""
    g_cool = 0.01 + 0.09 * ((T - T_lim) / 1000) + 0.2 * ((T - T_lim) / 1000) ** 2 + 0.16 * ((T - T_lim) / 1000) ** 3
    return g_cool if g_cool > 0 else 0


def mixing_param(params: list, mass_flows: list, contourings: list, error=0.01, Niter=20) -> float:
    """Расчет параметров смешения"""
    mix_param = params[0]
    for iteration in range(Niter):
        m_p = sum([params[i] * mass_flows[i] * contourings[i] for i in range(len(params))])
        m_p /= sum([mass_flows[i] * contourings[i] for i in range(len(params))])
        if eps('rel', mix_param, m_p) <= error: return m_p
        mix_param = m_p
    else:
        print(Fore.RED + f'Limit of iteration in {mixing_param.__name__}!' + Fore.RESET)
        return nan


class Substance:
    """Вещество"""

    def __init__(self, composition: dict) -> None:
        self.__composition = dict()
        if self.validate_composition(composition):
            self.__composition = composition

    def validate_composition(self, composition: dict) -> bool:
        assert isinstance(composition, dict), 'type(composition) is dict!'
        assert all(isinstance(key, str) for key in composition.keys()), \
            'type(composition.keys()) is str!'
        assert all(isinstance(value, (int, float)) for value in composition.values()), \
            'type(composition.values()) is int or float!'
        assert all(value >= 0 for value in composition.values()), 'composition.values() >= 0!'
        return True

    def __add__(self, other):
        assert isinstance(other, Substance), 'isinstance(other, Substance)'
        m = sum(self.composition.values()) + sum(other.composition.values())
        composition = dict()
        # TODO а если элемент повторяется
        for el in self.composition.keys():
            composition[el] = self.composition[el]
        for el in other.composition.keys():
            composition[el] = other.composition[el]

        return Substance(composition)

    @staticmethod
    def formula_to_dict(formula) -> dict[str: int]:
        result = {}
        i = 0
        while i < len(formula):
            if i + 1 < len(formula) and formula[i + 1].islower():
                atom = formula[i:i + 2]
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

    @property
    def composition(self) -> dict[str:int]:
        return self.__composition

    @property
    def mol_mass(self) -> tuple[float, str]:
        """Молярная масса"""
        if hasattr(self, "_Substance__mol_mass"): return self.__mol_mass
        m = sum(self.composition.values())
        self.__mol_mass = 0

        @lru_cache(maxsize=None)  # кэширование медленнее, если не применяется!
        def get_element_mass(element: str):
            return mendeleev.element(element).mass

        for formula, fraction in self.composition.items():
            formula_dict = Substance.formula_to_dict(formula)
            for el, atoms in formula_dict.items():
                self.__mol_mass += get_element_mass(el) * atoms * fraction

        '''for formula, fraction in self.composition.items():
            for el, atoms in Substance.formula_to_dict(formula).items():
                self.__mol_mass += mendeleev.element(el).mass * atoms * fraction'''

        self.__mol_mass = self.__mol_mass / m / 1000, 'kg/mol'
        return self.__mol_mass

    @property
    def gas_const(self) -> tuple[float, str]:
        """Газовая постоянная"""
        if hasattr(self, "_Substance__gas_const"): return self.__gas_const
        self.__gas_const = gas_const / self.mol_mass[0], 'J/kg/K'
        return self.__gas_const

    @property
    def Cp(self, T, P) -> tuple[float, str]:
        """Теплоемкость при постоянном давлении"""
        Cp = 2
        return Cp, 'J/kg/K'

    @property
    def excess_oxidizing(self) -> float:
        """Коэффициент избытка окислителя"""
        return 0

    @timeit()
    def summary(self) -> dict:
        print(f'composition: {self.composition}')
        print(f'mol_mass: {self.mol_mass}')
        print(f'gas_const: {self.gas_const}')

        return {'composition': self.composition,
                'mol_mass': self.mol_mass,
                'gas_const': self.gas_const}


if __name__ == '__main__':
    if 0:
        print(Fore.YELLOW + f'testing {atmosphere_standard.__name__}' + Fore.RESET)
        for H in (-2000, 0, 4_000, 11_000, 16_000):
            print(f'H = {H}: {atmosphere_standard(H)}')

    if 0:
        print(Fore.YELLOW + f'testing {Substance.__name__}' + Fore.RESET)
        s1 = Substance({'N2': 75.5})
        s1.summary()
        s2 = Substance({'O2': 23.15})
        s2.summary()
        s3 = Substance({'Ar': 1.292})
        s3.summary()
        s4 = Substance({'Ne': 0.0014})
        s4.summary()
        s5 = Substance({'H': 0.0008})

        s = s1 + s2 + s3 + s4 + s5
        s.summary()

        s = Substance({'N2': 0.755, 'O2': 0.2315, 'Ar': 0.01292, 'Ne': 0.000014, 'H': 0.000008})
        s.summary()

        s = Substance({'H2O': 11.25})
        s.summary()

        s = Substance({'CO2': 12})
        s.summary()

        print(dir(s))

    exit()
    print(av(Cp, [400, 500], [10 ** 5, 2 * 10 ** 5]))
    print(Cp('air', T=458))
    print(Cp('air', T=458, P=2 * 10 ** 5))
