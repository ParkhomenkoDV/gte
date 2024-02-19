"""Функции, использующиеся в расчете турбины"""

import pandas as pd
from scipy import integrate
import numpy as np


def add_to_data(data_frame, data):
    """Добавляет data_frame к словарю данных data"""
    for k, v, u in zip(data_frame['Обозначение'], data_frame['Значение'], data_frame['Единицы']):
        data[k] = v


def C_pi(T, gas='Воздух сухой'):
    """Находит значение истинной массовой теплоемкости воздуха и продуктов сгорания при alpha = 1 и T"""
    a = C_pi.A[gas]
    _T = T / 1000
    res = 0
    for i in range(7):
        res = res + a[i] * _T ** i
    return res * 10 ** 3


C_pi.A = pd.read_excel('Input_Data.xlsx', sheet_name='C_pi')


def C_pm(T1, T2, gas='Воздух сухой'):
    """Находит значение средней массовой теплоемкости в диапазоне температур (T1...T2)"""
    c_pi = lambda T: C_pi(T, gas)
    if T1 != T2:
        res, err = integrate.quad(c_pi, T1, T2)
        return res / (T2 - T1)
    else:
        return c_pi(T1)


def C_pg(alpha, T1, T2, fuel='Природный газ'):
    """Находит теплоемкость при произвольном значении alpha"""
    c_pv = C_pm(T1, T2)
    c_pg1 = C_pm(T1, T2, fuel)
    l_0 = C_pg.l_0[fuel][1]
    return c_pv + (1 + l_0) / (1 + alpha * l_0) * (c_pg1 - c_pv)


C_pg.l_0 = pd.read_excel('Input_Data.xlsx', sheet_name='l_0')


def R_G(alpha, fuel):
    """Возвращает значение газовой постоянной при alpha > 1"""
    a = R_G.A[fuel]
    return a[0] + a[1] / alpha


R_G.A = pd.read_excel('Input_Data.xlsx', sheet_name='R')


def find_keys(fname, D_name):
    """Создает упорядоченный по времени создания список ключей List словаря с именем D_name"""
    List = []
    for line in open(fname + '.py', encoding='utf-8'):
        if line.lstrip().startswith(D_name + '['):
            var = line.lstrip().split("'")[1]
            if var not in List:
                List.append(var)
    return List


def U_C0_eta(rho):
    return 0.25 * rho + 39 / 80


def convert_to_vector(x1, y1, x2, y2):
    """Находит координаты вектора по двум точкам с координатами {x1, y1} и {x2, y2}"""
    return x2 - x1, y2 - y1


def find_angle(x1, y1, x2, y2):
    """Находит угол в градусах между двумя векторами a = {x1, y1} и b = {x2, y2}"""
    a_vectimes_b = x1 * x2 + y1 * y2
    mod_a = np.sqrt(x1 ** 2 + y1 ** 2)
    mod_b = np.sqrt(x2 ** 2 + y2 ** 2)
    cos_alpha = a_vectimes_b / (mod_a * mod_b)
    angle = np.rad2deg(np.arccos(cos_alpha))
    if angle > 90:
        return angle - 90
    else:
        return angle


def blade_length(l_start, l_mid, l_end, z):
    """Возвращает массив длин лопаток, распределенных по ступеням турбины z параболическим законом по трем точкам
    l_start, l_mid, l_end"""
    c = np.polyfit([0, 1, z - 1], [l_start, l_mid, l_end], 2)
    p = lambda x: c[0] * x ** 2 + c[1] * x + c[2]
    l = np.zeros(z)
    for i in range(z):
        l[i] = p(i)
    return l


def polar_to_car(rho, phi):
    """Возвращает координаты в декартовой системе координат по полярным координатам rho и phi (град)"""
    return rho * np.cos(np.deg2rad(phi)), rho * np.sin(np.deg2rad(phi))


def vectorize(x1, y1, x2, y2, name):
    """Возвращает словарь с параметрами вектора"""
    return {'name': name,
            'x1': x1, 'y1': y1,
            'x2': x2, 'y2': y2,
            'angle': f'\\alpha_{name[-1]}' if name.startswith('c') else f'\\beta_{name[-1]}',
            'len': np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2),
            'x_mid': (x2 + x1) / 2, 'y_mid': (y2 + y1) / 2}
