"""
Список литературы:

[1] = лекции Арбекова

"""

import numpy as np
from scipy import interpolate


class Bearing:
    """Подшипник"""

    def __init__(self, name: str, bearing_type: str):
        assert isinstance(name, str)
        assert isinstance(bearing_type, str)
        bearing_type = bearing_type.strip().lower()
        assert bearing_type in ('rolling', 'plain')

        self.name = name
        self.type = bearing_type  # тип подшипника

    def calc(self, tp, d, D):
        """Ориентировочные расчеты"""
        if tp == 'шариковый':
            b = np.array([0.40, 0.48]) * (D - d)  # ширина
            Dw = np.array([0.27, 0.30]) * (D - d)  # диаметр шарика
            z = round(2.9 * (D + d) / (D - d))
        elif tp == 'роликовый':
            b = np.array([0.30, 0.49]) * (D - d)  # ширина
            Dw = np.array([0.20, 0.27]) * (D - d)  # диаметр ролика
            h = np.array([0.16, 0.25]) * Dw  # высота направляющих буртов
        else:
            raise Exception('Тип подшипника неизвестен')
        return {'b': b, 'Dw': Dw}

    def solve(self, mu, rotation_frequency, d, D, l, e):
        """Радиальная нгрузка масляного слоя в подшипнике скольжения"""
        ksi = (D - d) / d
        eta = e / (D / 2 - d / 2)
        f_m = interpolate.interp1d([1.0, 1.2, 1.5], [0.85, 1.0, 1.1], kind='linear', fill_value=(np.nan, np.nan))
        m = f_m(l / d)
        C_R = m / (1 - eta) - m
        return mu * rotation_frequency * d * l * C_R / ksi ** 2, 'N'


if __name__ == '__main__':
    bearing = Bearing('1', 'rolling')
    print(bearing.calc(tp='шариковый', d=40, D=60))
    print(bearing.solve(5 * 10 ** -6, 2800, 40, 60, 40, 0.001))
