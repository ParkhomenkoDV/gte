import sys
import numpy as np
import pandas as pd
from scipy import interpolate, integrate
import matplotlib.pyplot as plt

from material import Material

sys.path.append('D:\Programming\Python')

import decorators


class Dick:
    def __init__(self, radius: tuple | list | np.ndarray, thickness: tuple | list | np.ndarray,
                 material: dict,
                 nholes=tuple(), rholes=tuple(), dholes=tuple()):

        assert type(radius) in (tuple, list, np.ndarray)
        assert type(thickness) in (tuple, list, np.ndarray)

        assert type(material) is Material

        assert type(nholes) in (tuple, list, np.ndarray)
        assert type(rholes) in (tuple, list, np.ndarray)
        assert type(dholes) in (tuple, list, np.ndarray)

        assert len(radius) == len(thickness)
        assert len(nholes) == len(rholes) == len(dholes)

        assert all(map(lambda i: type(i) in (int, float, np.float64), radius))
        assert all(map(lambda i: type(i) in (int, float, np.float64), thickness))

        assert all(map(lambda i: type(i) in (int, float, np.float64), nholes))
        assert all(map(lambda i: type(i) in (int, float, np.float64), rholes))
        assert all(map(lambda i: type(i) in (int, float, np.float64), dholes))

        assert all(map(lambda i: i >= 0, radius))
        assert all(map(lambda i: i >= 0, thickness))

        assert all(map(lambda i: i >= 0, nholes))
        assert all(map(lambda i: i >= 0, rholes))
        assert all(map(lambda i: i >= 0, dholes))

        assert all(radius[i] < radius[i + 1] for i in range(len(radius) - 1))

        self.radius, self.thickness = np.array(radius), np.array(thickness)
        self.material = material
        self.nholes = nholes

    # TODO: numpy vectorizaoin
    def __calculation1(self, rotation_frequency, pressure,
                       av_radius, av_thickness,
                       av_E, av_mu, av_tetta) -> dict[str: tuple | list | np.ndarray]:
        """1й расчет метода двух расчетов"""
        st, sr = [400 * 10 ** 6], [pressure[0]]  # -p (посадка)
        for i, (r, t, E, mu, tetta) in enumerate(zip(av_radius[1:], thickness[1:], av_E[1:], av_mu[1:], av_tetta[1:])):
            st.append(mu * (av_thickness[i - 1] / t) * sr[i - 1] +
                      (E / av_E[i - 1]) * (st[i - 1] - av_mu[i - 1] * sr[i - 1]))
            sr.append((thickness[i - 1] / t) * sr[i - 1])

            c1 = 0.5 * (st[i] + sr[i] +
                        (4 + 4 * mu) / 8 * self.material['density'] * (rotation_frequency * r) ** 2 +
                        E * tetta)
            c2 = 0.5 * r ** 2 * (st[i] - sr[i] -
                                 ((2 - 2 * mu / 8 * self.material['density'] * (
                                         rotation_frequency * r) ** 2 - E * tetta)))

        return {'σ_t': st, 'σ_r': sr}

    def __calculation2(self, rotation_frequency, pressure,
                       av_radius, av_thickness,
                       av_E, av_mu, av_tetta) -> dict[str: tuple | list | np.ndarray]:
        """2й расчет метода двух расчетов"""
        st, sr = [1], [1]
        for i, (r, t) in enumerate(zip(av_radius, av_thickness)):
            pass
        return {'σ_t': st, 'σ_r': sr}

    @decorators.timeit(6)
    def tension(self, rotation_frequency: float, temperature0: int | float,
                pressure: tuple | list, temperature: tuple | list) -> dict[str: tuple | list | np.ndarray]:
        """Метод двух расчетов"""

        average_radius = [(self.radius[i] + self.radius[i + 1]) / 2 for i in range(len(self.radius) - 1)]

        f_thickness = interpolate.interp1d(self.radius, self.thickness, kind='linear')
        avarege_thickness = [f_thickness(av_r) for av_r in average_radius]

        f_temperature = lambda r: temperature[0] + (temperature[-1] - temperature[0]) * (
                (r - self.radius[0]) / (self.radius[-1] - self.radius[0])) ** 2
        avarege_temperature = [f_temperature(av_r) for av_r in average_radius]

        avarege_E = [self.material.E(av_t) for av_t in avarege_temperature]

        avarege_mu = [self.material['mu'](av_t) for av_t in avarege_temperature]

        tetta = self.material['alpha'] * np.array([temperature[0] - temperature0, temperature[-1] - temperature0])
        f_tetta = lambda r: tetta[0] + (tetta[-1] - tetta[0]) * (
                (r - self.radius[0]) / (self.radius[-1] - self.radius[0])) ** 2
        avarege_tetta = [f_tetta(av_t) for av_t in avarege_temperature]

        # TODO: multiprocessing
        calc1 = self.__calculation1(rotation_frequency, pressure,
                                    average_radius, avarege_thickness,
                                    avarege_E, avarege_mu, avarege_tetta)
        calc2 = self.__calculation2(rotation_frequency, pressure,
                                    average_radius, avarege_thickness,
                                    avarege_E, avarege_mu, avarege_tetta)
        k = (pressure[-1] - calc1['σ_r'][-1]) / calc2['σ_r'][-1]  # коэффициент Мора

        st, sr = [], []
        return {'σ_t': st, 'σ_r': sr}

    def show(self, **kwargs) -> None:
        plt.figure(figsize=kwargs.pop('figsize', (8, 8)))
        plt.title("Disk", fontsize=14, fontweight='bold')
        plt.plot([-self.thickness[0] / 1.5, self.thickness[0] / 1.5], [0, 0],
                 color='orange', linestyle='dashdot', linewidth=1.5)
        plt.plot(self.thickness / 2, self.radius, color='black', linestyle='solid', linewidth=3)
        plt.plot(-self.thickness / 2, self.radius, color='black', linestyle='solid', linewidth=3)
        plt.plot([-self.thickness[-1] / 2, self.thickness[-1] / 2], [self.radius[-1], self.radius[-1]],
                 color='black', linestyle='solid', linewidth=3)
        if self.radius[0] > 0:
            plt.plot([-self.thickness[0] / 2, self.thickness[0] / 2], [self.radius[0], self.radius[0]],
                     color='black', linestyle='solid', linewidth=3)
        plt.grid(True)
        plt.axis('equal')
        plt.xlabel("Thickness", fontsize=12)
        plt.ylabel("Radius", fontsize=12)

        plt.show()

    def show_tension(self, **kwargs) -> None:
        fg = plt.figure(figsize=kwargs.pop('figsize', (16, 8)))
        gs = fg.add_gridspec(1, 2)  # строки, столбцы

        fg.add_subplot(gs[0, 0])
        plt.title("Disk", fontsize=14, fontweight='bold')
        plt.grid(True)
        plt.axis('equal')
        plt.xlabel("Thickness", fontsize=12)
        plt.ylabel("Radius", fontsize=12)

        fg.add_subplot(gs[0, 1])
        plt.title('Tension', fontsize=14, fontweight='bold')
        plt.grid(True)
        plt.axis('equal')
        plt.xlabel("Tension", fontsize=12)
        plt.ylabel("Radius", fontsize=12)

        plt.show()


if __name__ == "__main__":
    radius = np.array([20, 26, 30.62, 37.26, 56.94, 60.67, 72.95, 75.95, 102.41, 106.52, 109.82]) / 1000
    thickness = np.array([36, 36, 15.43, 11.27, 10, 12, 12, 8, 6, 11, 11]) / 1000

    print(pd.DataFrame({'radius': radius, 'thickness': thickness}))

    material = Material('10Х11Н20ТЗР',
                        {
                            'density': 8400,
                            'alpha': lambda t: 18 * 10 ** -6 if 400 <= t <= 800 else nan,
                            "E": [(400, 1.74 * 10 ** 11),
                                  (500, 1.66 * 10 ** 11),
                                  (600, 1.57 * 10 ** 11),
                                  (700, 1.47 * 10 ** 11),
                                  (800, 1.32 * 10 ** 11)],
                            "mu": [(400, 0.384),
                                   (500, 0.379),
                                   (600, 0.371),
                                   (700, 0.361),
                                   (800, 0.347)]
                        })

    disk = Dick(radius=radius, thickness=thickness,
                material=material)

    disk.show()
    disk.tension(rotation_frequency=2806.2, temperature0=293.15,
                 pressure=(0, 120.6 * 10 ** 6), temperature=(350, 650))
    disk.show_tension()
