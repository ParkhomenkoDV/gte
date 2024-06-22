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
    def __calculation1(self, rotation_frequency: int | float, pressure0: int | float,
                       radius: np.ndarray, thickness: np.ndarray, tetta: np.ndarray,
                       av_density, av_E, av_mu) -> dict[str: np.ndarray]:
        """1‑й расчет метода двух расчетов"""
        sigma_t, sigma_r = np.zeros((len(radius) - 1, 2)), np.zeros((len(radius) - 1, 2))
        for i in range(len(radius) - 1):
            if i == 0:
                sigma_t[i][0] = 400 * 10 ** 6
                sigma_r[i][0] = pressure0  # -p (посадка)
            else:
                sigma_t[i][0] = (av_mu[i] * (thickness[i - 1] / thickness[i]) * sigma_r[i - 1][1] +
                                 (av_E[i] / av_E[i - 1]) * (sigma_t[i - 1][1] - av_mu[i - 1] * sigma_r[i - 1][1]))
                sigma_r[i][0] = ((thickness[i - 1] / thickness[i]) * sigma_r[i - 1][1])
            c1 = (sigma_t[i][0] + sigma_r[i][0] +
                  (1 + av_mu[i]) / 2 * av_density[i] * (rotation_frequency * radius[i]) ** 2 +
                  av_E[i] * tetta[i]) / 2
            c2 = (sigma_t[i][0] - sigma_r[i][0] -
                  (1 - av_mu[i]) / 4 * av_density[i] * (rotation_frequency * radius[i]) ** 2 +
                  av_E[i] * tetta[i]) / 2 * radius[i] ** 2
            tettas = ((tetta[i] - (tetta[i + 1] - tetta[i]) / (radius[i + 1] - radius[i]) * radius[i]) *
                      (radius[i + 1] ** 2 - radius[i] ** 2) / 2 +
                      (tetta[i + 1] - tetta[i]) / (radius[i + 1] - radius[i]) *
                      (radius[i + 1] ** 3 - radius[i] ** 3) / 3)
            sigma_t[i][1] = ((c1 + c2 / radius[i + 1] ** 2 -
                              (1 + 3 * av_mu[i]) / 8 * av_density[i] * (rotation_frequency * radius[i + 1]) ** 2 +
                              (av_E[i] / radius[i + 1] ** 2) * tettas) - av_E[i] * tetta[i + 1])
            sigma_r[i][1] = ((c1 - c2 / radius[i + 1] ** 2) -
                             (3 + 1 * av_mu[i]) / 8 * av_density[i] * (rotation_frequency * radius[i + 1]) ** 2 -
                             (av_E[i] / radius[i + 1] ** 2) * tettas)

        return {'tension_t': sigma_t, 'tension_r': sigma_r}

    def __calculation2(self, radius: np.ndarray, thickness: np.ndarray, av_E, av_mu) -> dict[str: np.ndarray]:
        """2-й расчет метода двух расчетов"""
        sigma_t, sigma_r = np.zeros((len(radius) - 1, 2)), np.zeros((len(radius) - 1, 2))
        for i in range(len(radius) - 1):
            if i == 0:
                sigma_t[i][0] = 400 * 10 ** 6
                sigma_r[i][0] = sigma_t[i][0] if radius[0] == 0 else 0
            else:
                sigma_t[i][0] = (av_mu[i] * (thickness[i - 1] / thickness[i]) * sigma_r[i - 1][1] +
                                 (av_E[i] / av_E[i - 1]) * (sigma_t[i - 1][1] - av_mu[i - 1] * sigma_r[i - 1][1]))
                sigma_r[i][0] = (thickness[i - 1] / thickness[i]) * sigma_r[i - 1][1]
            c1 = (sigma_t[i][0] + sigma_r[i][0]) / 2
            c2 = (sigma_t[i][0] - sigma_r[i][0]) / 2 * radius[i] ** 2
            sigma_t[i][1] = c1 + c2 / radius[i + 1] ** 2
            sigma_r[i][1] = c1 - c2 / radius[i + 1] ** 2

        return {'tension_t': sigma_t, 'tension_r': sigma_r}

    @decorators.timeit(4)
    def tension(self, rotation_frequency: float, temperature0: int | float,
                pressure: tuple | list, temperature: tuple | list) -> dict[str: tuple | list | np.ndarray]:
        """Метод двух расчетов"""

        average_radius = np.array([(self.radius[i] + self.radius[i + 1]) / 2 for i in range(len(self.radius) - 1)])

        tetta = [self.material.alpha((temperature[0] + temperature0) / 2) * (temperature[0] - temperature0),
                 self.material.alpha((temperature[-1] + temperature0) / 2) * (temperature[-1] - temperature0)]
        # TODO: параболический закон распределения температур
        f_tetta = lambda r: (tetta[0] +
                             (tetta[-1] - tetta[0]) * ((r - self.radius[0]) / (self.radius[-1] - self.radius[0])) ** 2)
        self.tetta = np.array([f_tetta(r) for r in self.radius])
        avarege_tetta = [(self.tetta[i] + self.tetta[i + 1]) / 2 for i in range(len(self.radius) - 1)]

        # TODO: параболический закон распределения температур
        f_temperature = lambda r: (temperature[0] + (temperature[-1] - temperature[0]) *
                                   ((r - self.radius[0]) / (self.radius[-1] - self.radius[0])) ** 2)
        avarege_temperature = [f_temperature(av_r) for av_r in average_radius]

        # TODO: средняя или абс?
        avarege_density = np.array([self.material.density(av_t) for av_t in avarege_temperature])
        avarege_E = np.array([self.material.E(av_t) for av_t in avarege_temperature])
        avarege_mu = np.array([self.material.mu(av_t) for av_t in avarege_temperature])

        # TODO: multiprocessing
        calc1 = self.__calculation1(rotation_frequency, pressure[0],
                                    self.radius, self.thickness, self.tetta,
                                    avarege_density, avarege_E, avarege_mu)
        calc2 = self.__calculation2(self.radius, self.thickness, avarege_E, avarege_mu)
        k = (pressure[-1] - calc1['tension_r'][-1][-1]) / calc2['tension_r'][-1][-1]  # коэффициент Мора
        sigma_t, sigma_r = np.zeros((len(radius) - 1, 2)), np.zeros((len(radius) - 1, 2))
        for i in range(len(self.radius) - 1):
            for j in range(2):  # сечения
                sigma_t[i][j] = calc1['tension_t'][i][j] + k * calc2['tension_t'][i][j]
                sigma_r[i][j] = calc1['tension_r'][i][j] + k * calc2['tension_r'][i][j]

        av_sigma_t = np.array([sigma_t[0][0]] +
                              [(sigma_t[i][1] + sigma_t[i + 1][0]) / 2 for i in range(0, len(self.radius) - 2)] +
                              [sigma_t[-1][-1]])
        av_sigma_r = np.array([sigma_r[0][0]] +
                              [(sigma_r[i][1] + sigma_r[i + 1][0]) / 2 for i in range(0, len(self.radius) - 2)] +
                              [sigma_r[-1][-1]])
        sigma = np.array([np.sqrt(av_sigma_t[i] ** 2 - av_sigma_t[i] * av_sigma_r[i] + av_sigma_r[i] ** 2)
                          for i in range(len(radius))])

        tensions = {'tension': sigma, 'tension_t': av_sigma_t, 'tension_r': av_sigma_r}

        self.show_tension(tensions)

        return tensions

    def show_tension(self, tensions, **kwargs) -> None:
        fg = plt.figure(figsize=kwargs.pop('figsize', (16, 8)))
        gs = fg.add_gridspec(1, 2)  # строки, столбцы

        fg.add_subplot(gs[0, 0])
        plt.title("Disk", fontsize=14, fontweight='bold')

        plt.plot([-self.thickness[0] / 1.5, self.thickness[0] / 1.5], [0, 0],
                 color='orange', linestyle='dashdot', linewidth=1.5)  # ось вращения
        plt.plot(self.thickness / 2, self.radius, -self.thickness / 2, self.radius,
                 color='black', linestyle='solid', linewidth=3)
        plt.plot([-self.thickness[-1] / 2, self.thickness[-1] / 2], [self.radius[-1], self.radius[-1]],
                 color='black', linestyle='solid', linewidth=3)  # периферия
        if self.radius[0] > 0:
            plt.plot([-self.thickness[0] / 2, self.thickness[0] / 2], [self.radius[0], self.radius[0]],
                     color='black', linestyle='solid', linewidth=3)  # втулка
        # TODO: добавить прямоугольники сечений
        plt.grid(True)
        plt.axis('equal')
        plt.xlabel("Thickness", fontsize=12)
        plt.ylabel("Radius", fontsize=12)

        fg.add_subplot(gs[0, 1])
        plt.title('Tension', fontsize=14, fontweight='bold')

        plt.plot(tensions['tension'] / 10 ** 6, self.radius * 1000, label='экв',
                 color='black', linestyle='solid', linewidth=3)
        plt.plot(tensions['tension_t'] / 10 ** 6, self.radius * 1000, label='t',
                 color='red', linestyle='solid', linewidth=2)
        plt.plot(tensions['tension_r'] / 10 ** 6, self.radius * 1000, label='r',
                 color='blue', linestyle='solid', linewidth=2)

        plt.plot(tensions['tension'] / 10 ** 6, self.radius * 1000, )

        plt.legend()
        plt.grid(True)
        plt.xlabel("Tension [MPa]", fontsize=12)
        plt.ylabel("Radius [mm]", fontsize=12)

        plt.show()

    def show(self, **kwargs) -> None:
        plt.figure(figsize=kwargs.pop('figsize', (8, 8)))
        plt.title("Disk", fontsize=14, fontweight='bold')
        plt.plot([-self.thickness[0] / 1.5, self.thickness[0] / 1.5], [0, 0],
                 color='orange', linestyle='dashdot', linewidth=1.5)  # ось вращения
        plt.plot(self.thickness / 2, self.radius, -self.thickness / 2, self.radius,
                 color='black', linestyle='solid', linewidth=3)
        plt.plot([-self.thickness[-1] / 2, self.thickness[-1] / 2], [self.radius[-1], self.radius[-1]],
                 color='black', linestyle='solid', linewidth=3)  # периферия
        if self.radius[0] > 0:
            plt.plot([-self.thickness[0] / 2, self.thickness[0] / 2], [self.radius[0], self.radius[0]],
                     color='black', linestyle='solid', linewidth=3)  # втулка
        plt.grid(True)
        plt.axis('equal')
        plt.xlabel("Thickness", fontsize=12)
        plt.ylabel("Radius", fontsize=12)

        plt.show()


if __name__ == "__main__":
    radius = np.array([20, 26, 30.62, 37.26, 56.94, 60.67, 72.95, 75.95, 102.41, 106.52, 109.82]) / 1000
    thickness = np.array([36, 36, 15.43, 11.27, 10, 12, 12, 8, 6, 11, 11]) / 1000

    print(pd.DataFrame({'radius': radius, 'thickness': thickness}))

    material = Material('10Х11Н20ТЗР',
                        {
                            "density": 8400,
                            "alpha": 18 * 10 ** -6,
                            "E": interpolate.interp1d(list(range(400, 800 + 1, 100)),
                                                      np.array([1.74, 1.66, 1.57, 1.47, 1.32]) * 10 ** 11,
                                                      kind='cubic', bounds_error=False, fill_value='extrapolate'),
                            "mu": interpolate.interp1d(list(range(400, 800 + 1, 100)),
                                                       [0.384, 0.379, 0.371, 0.361, 0.347],
                                                       kind='cubic', bounds_error=False, fill_value='extrapolate')
                        })

    disk = Dick(radius=radius, thickness=thickness,
                material=material)

    # disk.show()
    print(disk.tension(rotation_frequency=2806.2, temperature0=293.15,
                       pressure=(0, 120.6 * 10 ** 6), temperature=(350, 650)))
