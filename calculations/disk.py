import numpy as np
import pandas as pd
from scipy import interpolate, integrate
import matplotlib as mpl, matplotlib.pyplot as plt

from material import Material


class Dick:

    @classmethod
    def version(cls):
        version = '2.3'
        print('Увеличить дискритизацию')
        print('Подбор толщин и радиусов для напряжения')
        return version

    def __init__(self, material: Material,
                 radius: tuple | list | np.ndarray, thickness: tuple | list | np.ndarray,
                 nholes=tuple(), rholes=tuple(), dholes=tuple()):
        # проверка на тип данных radius и thickness
        assert type(radius) in (tuple, list, np.ndarray)
        assert type(thickness) in (tuple, list, np.ndarray)
        # проверка на тип данных material
        assert type(material) is Material
        # проверка на тип данных nholes, rholes, dholes
        assert type(nholes) in (tuple, list, np.ndarray)
        assert type(rholes) in (tuple, list, np.ndarray)
        assert type(dholes) in (tuple, list, np.ndarray)
        # проверка на равенство массивов radius, thickness и nholes, rholes, dholes
        assert len(radius) == len(thickness)
        assert len(nholes) == len(rholes) == len(dholes)
        # проверка на тип данных внутри массивов radius, thickness
        assert all(map(lambda i: type(i) in (int, float, np.float64), radius))
        assert all(map(lambda i: type(i) in (int, float, np.float64), thickness))
        # проверка на тип данных внутри массивов nholes, rholes, dholes
        assert all(map(lambda i: type(i) in (int, float, np.float64), nholes))
        assert all(map(lambda i: type(i) in (int, float, np.float64), rholes))
        assert all(map(lambda i: type(i) in (int, float, np.float64), dholes))
        # проверка на адекватность данных radius, thickness
        assert all(map(lambda i: i >= 0, radius))
        assert all(map(lambda i: i >= 0, thickness))
        # проверка на адекватность данных nholes, rholes, dholes
        assert all(map(lambda i: i >= 0, nholes))
        assert all(map(lambda i: i >= 0, rholes))
        assert all(map(lambda i: i >= 0, dholes))
        # сортировка пузырьком радиусов по возрастанию вместе с соответствующими толщинами
        swapped = False
        for i in range(len(radius) - 1, 0, -1):
            for j in range(i):
                if radius[j] > radius[j + 1]:
                    radius[j], radius[j + 1] = radius[j + 1], radius[j]
                    thickness[j], thickness[j + 1] = thickness[j + 1], thickness[j]
                    swapped = True
            if swapped:
                swapped = False
            else:
                break
        # проверка на отсортированность radius по возрастанию
        assert all(radius[i] < radius[i + 1] for i in range(len(radius) - 1))
        # проверка на адекватность rholes по отношению к radius
        assert all(map(lambda i: radius[0] <= i <= radius[-1], rholes))
        # расчет и проверка расстояния между краями отверстий по окружности
        self.b = np.array([2 * np.pi * rholes[i] / nholes[i] - dholes[i] for i in range(len(nholes))])
        assert all(map(lambda i: i > 0, self.b))

        self.radius, self.thickness = np.array(radius), np.array(thickness)
        self.material = material
        self.nholes, self.rholes, self.dholes = np.array(nholes), np.array(rholes), np.array(dholes)

    @staticmethod
    def slicing(point0: tuple, point1: tuple, ndis: int) -> tuple:
        """Дробление сечений ndis раз"""
        x, y = list(), list()
        k = (point1[1] - point0[1]) / (point1[0] - point0[0]) if (point1[0] - point0[0]) != 0 else np.inf
        b = point0[1] - k * point0[0]
        delta = (point1[0] - point0[0]) / ndis
        for i in range(ndis):
            x.append(point0[0] + delta * i)
            y.append(k * x[-1] + b)
        return np.array(x), np.array(y)

    @staticmethod
    def equivalent_energy_tension(sigma_t: float | int, sigma_r: float | int) -> float:
        """Эквивалентное напряжение по энергетической теории"""
        sigma1, sigma3 = max(sigma_t, sigma_r), min(sigma_t, sigma_r)
        return np.sqrt(sigma1 ** 2 - sigma1 * sigma3 + sigma3 ** 2)

    @staticmethod
    def calculation2(rotation_frequency: int | float, pressure: tuple | list | np.ndarray,
                     radius: np.ndarray, thickness: np.ndarray, tetta: np.ndarray,
                     av_density: np.ndarray, av_E: np.ndarray, av_mu: np.ndarray) -> dict[str: np.ndarray]:
        """Метод двух расчетов"""
        f_tetta = lambda r: (tetta[0] +
                             (tetta[-1] - tetta[0]) * ((r - radius[0]) / (radius[-1] - radius[0])) ** 2)
        sigma_t = np.full((len(radius) - 1, 2), 400 * 10 ** 6)
        sigma_r = sigma_t.copy() if radius[0] == 0 else np.full((len(radius) - 1, 2), pressure[0])
        for i in range(len(radius) - 1):
            if i != 0:
                sigma_t[i][0] = (av_mu[i] * (thickness[i - 1] / thickness[i]) * sigma_r[i - 1][1] +
                                 (av_E[i] / av_E[i - 1]) * (sigma_t[i - 1][1] - av_mu[i - 1] * sigma_r[i - 1][1]))
                sigma_r[i][0] = (thickness[i - 1] / thickness[i]) * sigma_r[i - 1][1]
            c1 = (sigma_t[i][0] + sigma_r[i][0] +
                  (1 + av_mu[i]) / 2 * av_density[i] * (rotation_frequency * radius[i]) ** 2 +
                  av_E[i] * tetta[i]) / 2
            c2 = (sigma_t[i][0] - sigma_r[i][0] -
                  (1 - av_mu[i]) / 4 * av_density[i] * (rotation_frequency * radius[i]) ** 2 +
                  av_E[i] * tetta[i]) / 2 * radius[i] ** 2
            I = integrate.quad(lambda r: f_tetta(r) * r, radius[i], radius[i + 1])[0]
            sigma_t[i][1] = ((c1 + c2 / radius[i + 1] ** 2 -
                              (1 + 3 * av_mu[i]) / 8 * av_density[i] * (rotation_frequency * radius[i + 1]) ** 2 +
                              (av_E[i] / radius[i + 1] ** 2) * I) - av_E[i] * tetta[i + 1])
            sigma_r[i][1] = ((c1 - c2 / radius[i + 1] ** 2) -
                             (3 + 1 * av_mu[i]) / 8 * av_density[i] * (rotation_frequency * radius[i + 1]) ** 2 -
                             (av_E[i] / radius[i + 1] ** 2) * I)

        return {'tension_t': sigma_t, 'tension_r': sigma_r}

    def tension(self, rotation_frequency: float, temperature0: int | float,
                pressure: tuple | list | np.ndarray, temperature: tuple | list | np.ndarray,
                ndis: int = 10) -> dict[str:  np.ndarray]:
        """Расчет напряжений в диске"""

        assert type(rotation_frequency) in (float, int)
        assert type(temperature0) in (float, int) and temperature0 > 0
        assert type(pressure) in (tuple, list, np.ndarray)
        assert type(temperature) in (tuple, list, np.ndarray)
        assert all(map(lambda i: type(i) in (int, float, np.float64), pressure))
        assert all(map(lambda i: type(i) in (int, float, np.float64), temperature))
        assert len(pressure) == 2
        assert len(temperature) == 2 or len(temperature) == len(self.radius)

        tetta = [self.material.alpha((temperature[0] + temperature0) / 2) * (temperature[0] - temperature0),
                 self.material.alpha((temperature[-1] + temperature0) / 2) * (temperature[-1] - temperature0)]
        # параболический закон распределения деформаций
        f_tetta = lambda r: \
            (tetta[0] + (tetta[-1] - tetta[0]) * ((r - self.radius[0]) / (self.radius[-1] - self.radius[0])) ** 2)
        self.tetta = np.array([f_tetta(r) for r in self.radius])

        if len(temperature) == len(self.radius):
            f_temperature = interpolate.interp1d(self.radius, temperature, kind='cubic', fill_value='extrapolate')
        else:
            f_temperature = lambda r: (temperature[0] + (temperature[-1] - temperature[0]) *
                                       ((r - self.radius[0]) / (self.radius[-1] - self.radius[0])) ** 2)

        # TODO: сюда дискритезацию
        radius, thickness = list(), list()
        for i in range(len(self.radius) - 1):
            r, th = self.slicing((self.radius[i], self.thickness[i]),
                                 (self.radius[i + 1], self.thickness[i + 1]),
                                 ndis)
            radius.extend(r)
            thickness.extend(th)
        radius, thickness = np.array(radius), np.array(thickness)

        average_radius = np.array([(radius[i] + radius[i + 1]) / 2 for i in range(len(radius) - 1)])
        avarege_temperature = [f_temperature(av_r) for av_r in average_radius]
        avarege_density = np.array([self.material.density(av_t) for av_t in avarege_temperature])
        avarege_E = np.array([self.material.E(av_t) for av_t in avarege_temperature])
        avarege_mu = np.array([self.material.mu(av_t) for av_t in avarege_temperature])
        tetta = np.array([f_tetta(r) for r in radius])

        calc1 = self.calculation2(rotation_frequency, pressure,
                                  radius, thickness, tetta,
                                  avarege_density, avarege_E, avarege_mu)
        calc2 = self.calculation2(0, pressure,
                                  radius, thickness, np.zeros(len(radius)),
                                  avarege_density, avarege_E, avarege_mu)

        k = (pressure[-1] - calc1['tension_r'][-1][-1]) / calc2['tension_r'][-1][-1]  # коэффициент Мора

        sigma_t, sigma_r = np.zeros((len(radius) - 1, 2)), np.zeros((len(radius) - 1, 2))
        for i in range(len(radius) - 1):  # участки
            for j in range(2):  # сечения
                sigma_t[i][j] = calc1['tension_t'][i][j] + k * calc2['tension_t'][i][j]
                sigma_r[i][j] = calc1['tension_r'][i][j] + k * calc2['tension_r'][i][j]

        sigma_t = np.array([sigma_t[0][0]] +
                           [(sigma_t[i][1] + sigma_t[i + 1][0]) / 2 for i in range(0, len(radius) - 2)] +
                           [sigma_t[-1][-1]])
        sigma_r = np.array([sigma_r[0][0]] +
                           [(sigma_r[i][1] + sigma_r[i + 1][0]) / 2 for i in range(0, len(radius) - 2)] +
                           [sigma_r[-1][-1]])

        sigma = np.array([self.equivalent_energy_tension(sigma_t[i], sigma_r[i]) for i in range(len(radius))])

        result = {'radius': radius, 'thickness': thickness,
                  'tension': sigma, 'tension_t': sigma_t, 'tension_r': sigma_r}

        df = pd.DataFrame({'radius [mm]': result['radius'] * 1000,
                           'thickness [mm]': result['thickness'] * 1000,
                           'tension [MPa]': result['tension'] / 10 ** 6,
                           'tension_t': result['tension_t'] / 10 ** 6,
                           'tension_r': result['tension_r'] / 10 ** 6}).sort_values(by='radius [mm]', ascending=False)
        print(df)

        f_sigma_t = interpolate.interp1d(result['radius'], result['tension_t'], kind='linear')
        f_sigma_r = interpolate.interp1d(result['radius'], result['tension_r'], kind='linear')

        for i in range(len(self.nholes)):
            k = 3 - self.dholes[i] / self.b[i] - f_sigma_r(self.rholes[i]) / f_sigma_t(self.rholes[i])
            sigma_t_hole = k * f_sigma_t(self.rholes[i]) / 10 ** 6
            print(f'''
            holes: {i}
            nholes []: {self.nholes[i]}, rholes [mm]: {self.rholes[i] * 1000}, dholes [mm]: {self.dholes[i] * 1000}
            tension_t [MPa] in {(sigma_t_hole * 1.1, sigma_t_hole * 1.15)}
            ''')

        self.show_tension(result)

        return result

    def show_tension(self, tensions, **kwargs) -> None:

        radius, thickness = tensions['radius'] * 1_000, tensions['thickness'] * 1_000  # приведение к [мм]
        for key in tensions:
            if key.startswith('tension'):
                tensions[key] = tensions[key] / 10 ** 6  # приведение к [МПа]

        l, k = max(radius), 1.2
        ylim = 0 - l * (k - 1) / 2, l + l * (k - 1) / 2

        fg = plt.figure(figsize=kwargs.pop('figsize', (16, 8)))
        gs = fg.add_gridspec(1, 2)  # строки, столбцы

        fg.add_subplot(gs[0, 0])
        plt.title("Disk", fontsize=14, fontweight='bold')

        plt.plot(thickness / 2, radius, -thickness / 2, radius,
                 color='black', linestyle='solid', linewidth=3)
        plt.plot([-thickness[-1] / 2, thickness[-1] / 2], [radius[-1], radius[-1]],
                 color='black', linestyle='solid', linewidth=3)  # периферия
        if radius[0] > 0:
            plt.plot([-thickness[0] / 2, thickness[0] / 2], [radius[0]] * 2,
                     color='black', linestyle='solid', linewidth=3)  # втулка
        for i in range(len(radius) - 1):  # сечения
            av_th = (thickness[i] + thickness[i + 1]) / 2
            plt.gca().add_patch(mpl.patches.Rectangle((-av_th / 2, radius[i]),
                                                      av_th, radius[i + 1] - radius[i],
                                                      angle=0, rotation_point='xy', alpha=0.5))
        plt.plot([-thickness[0] / 1.5, thickness[0] / 1.5], [0] * 2,
                 color='orange', linestyle='dashdot', linewidth=1.5)  # ось вращения
        for i in range(len(self.nholes)):
            plt.plot([-thickness[0] / 1.5, thickness[0] / 1.5], [self.rholes[i] * 1000] * 2,
                     color='orange', linestyle='dashdot', linewidth=1.5)
            plt.plot([-thickness[0] / 1.5, thickness[0] / 1.5], [(self.rholes[i] + self.dholes[i] / 2) * 1_000] * 2,
                     [-thickness[0] / 1.5, thickness[0] / 1.5], [(self.rholes[i] - self.dholes[i] / 2) * 1_000] * 2,
                     color='black', linestyle='dashed', linewidth=1.5)

        # TODO: добавить оси направлений r и t и прорисовку НУ на 1ой картинке
        # plt.arrow (x= 4 , y= 18 , dx= 0 , dy= 5 , width= .08 )
        # plt.annotate('General direction', xy = (3.4, 17)) #add annotation
        plt.grid(True)
        plt.axis('equal')
        plt.ylim(ylim)
        plt.xlabel("Thickness [mm]", fontsize=12)
        plt.ylabel("Radius [mm]", fontsize=12)

        fg.add_subplot(gs[0, 1])
        plt.title('Tension', fontsize=14, fontweight='bold')

        plt.plot(tensions['tension'], radius, label='equivalent', color='black', linestyle='solid', linewidth=3)
        plt.plot(tensions['tension_t'], radius, label='tangential', color='red', linestyle='solid', linewidth=2)
        plt.plot(tensions['tension_r'], radius, label='radial', color='blue', linestyle='solid', linewidth=2)

        plt.legend(fontsize=12)
        plt.grid(True)
        plt.ylim(ylim)
        plt.xlabel("Tension [MPa]", fontsize=12)
        plt.ylabel("Radius [mm]", fontsize=12)

        plt.show()

    def show(self, **kwargs) -> None:
        radius, thickness = self.radius * 1_000, self.thickness * 1_000  # приведение к [мм]

        plt.figure(figsize=kwargs.pop('figsize', (8, 8)))
        plt.title("Disk", fontsize=14, fontweight='bold')

        plt.plot(thickness / 2, radius, -thickness / 2, radius,
                 color='black', linestyle='solid', linewidth=3)
        plt.plot([-thickness[-1] / 2, thickness[-1] / 2], [radius[-1], radius[-1]],
                 color='black', linestyle='solid', linewidth=3)  # периферия
        if radius[0] > 0:
            plt.plot([-thickness[0] / 2, thickness[0] / 2], [radius[0], radius[0]],
                     color='black', linestyle='solid', linewidth=3)  # втулка
        plt.plot([-thickness[0] / 1.5, thickness[0] / 1.5], [0, 0],
                 color='orange', linestyle='dashdot', linewidth=1.5)  # ось вращения
        for i in range(len(self.nholes)):
            plt.plot([-thickness[0] / 1.5, thickness[0] / 1.5], [self.rholes[i] * 1_000] * 2,
                     color='orange', linestyle='dashdot', linewidth=1.5)
            plt.plot([-thickness[0] / 1.5, thickness[0] / 1.5], [(self.rholes[i] + self.dholes[i] / 2) * 1_000] * 2,
                     [-thickness[0] / 1.5, thickness[0] / 1.5], [(self.rholes[i] - self.dholes[i] / 2) * 1_000] * 2,
                     color='black', linestyle='dashed', linewidth=1.5)

        plt.grid(True)
        plt.axis('equal')
        plt.xlabel("Thickness [mm]", fontsize=12)
        plt.ylabel("Radius [mm]", fontsize=12)

        plt.show()

    def frequency_safety_factor(self, rotation_frequency: int | float) -> tuple[float]:
        """Запас по разрушающей частоте"""
        f_thickness = interpolate.interp1d(self.radius, self.thickness, kind='cubic')  # TODO: linear!
        n = self.material.sigma_temp(0) * integrate.quad(f_thickness, self.radius[0], self.radius[-1])[0]
        n /= (pressure[-1] * self.radius[-1] * self.thickness[-1] +
              (rotation_frequency ** 2 * self.material.density(0) *
               integrate.quad(lambda r: f_thickness(r) * r ** 2, self.radius[0], self.radius[-1])[0]))
        # TODO: плотность и предел временной прочности зависит от температуры
        n = np.sqrt(n)
        return n * 0.9, n * 0.95  # действительный интервал запаса разрушающей частоты


if __name__ == "__main__":
    print(Dick.version())
    if 1:
        material = Material('10Х11Н20ТЗР',
                            {
                                "density": 8400,
                                "alpha": 18 * 10 ** -6,
                                "E": interpolate.interp1d(list(range(400, 800 + 1, 100)),
                                                          np.array([1.74, 1.66, 1.57, 1.47, 1.32]) * 10 ** 11,
                                                          kind='cubic', bounds_error=False, fill_value='extrapolate'),
                                "mu": interpolate.interp1d(list(range(400, 800 + 1, 100)),
                                                           [0.384, 0.379, 0.371, 0.361, 0.347],
                                                           kind='cubic', bounds_error=False, fill_value='extrapolate'),
                                "sigma_temp": 900 * 10 ** 6
                            })
        radius = np.array([20, 26, 30.62, 37.26, 56.94, 60.67, 72.95, 75.95, 102.41, 106.52, 109.82]) / 1000
        thickness = np.array([36, 36, 15.43, 11.27, 10, 12, 12, 8, 6, 11, 11]) / 1000
        nholes, rholes, dholes = [5], [66.8 / 1000], [6.2 / 1000]

        print(pd.DataFrame({'radius': radius, 'thickness': thickness}))

        rotation_frequency = 2806.2
        temperature0 = 293.15
        pressure = (0, 120.6 * 10 ** 6)
        temperature = (350, 650)

    if 0:
        material = Material('10Х11Н20ТЗР',
                            {
                                "density": 8000,
                                "alpha": 20 * 10 ** -6,
                                "E": 1.3 * 10 ** 11,
                                "mu": 0.33,
                                "sigma_temp": 600 * 10 ** 6
                            })
        radius = np.array(
            [0, 272, 311, 352, 388, 434, 506.5, 509, 544, 564, 584, 619, 621.5, 726, 737, 748, 763, 775, 860, 864, 868,
             880, 884, 887, 905]) / 1000
        thickness = np.array(
            [132, 132, 120, 108, 97, 83, 88.5, 154, 154, 154, 154, 88.5, 83, 90, 107, 122, 96, 83, 89, 96, 127, 96, 89,
             83, 83]) / 1000
        nholes, rholes, dholes = [10], [544 / 1000], [40 / 1000]

        print(pd.DataFrame({'radius': radius, 'thickness': thickness}))

        rotation_frequency = 3000 * 2 * np.pi / 60
        temperature0 = 620
        pressure = (0, 150 * 10 ** 6)
        temperature = (620, 800)

    if 0:
        material = Material('10Х11Н20ТЗР',
                            {
                                "density": 8200,
                                "alpha": 18 * 10 ** -6,
                                "E": 1.6 * 10 ** 11,
                                "mu": 0.33,
                                "sigma_temp": 900 * 10 ** 6
                            })
        radius = np.array(
            [84, 95.1, 107.7, 120.2, 132.8, 145.4, 157.9, 170.6, 176.6, 189.2, 203.1, 218.9, 233.3, 245.9, 258.4, 271.9,
             277.8, 283.5, 289.3, 303]) / 1000
        thickness = np.array(
            [83, 83, 91, 69, 64, 58, 58, 58, 48, 43, 39, 35.8, 31, 27.8, 29.6, 35.6, 40, 50, 62, 62]) / 1000
        nholes, rholes, dholes = [], [], []

        print(pd.DataFrame({'radius': radius, 'thickness': thickness}))

        rotation_frequency = 9150 * 2 * np.pi / 60
        temperature0 = 300
        pressure = (0, 110 * 10 ** 6)
        temperature = (800, 1050)

    disk = Dick(material=material,
                radius=radius, thickness=thickness,
                nholes=nholes, rholes=rholes, dholes=dholes)

    disk.show()
    disk.tension(rotation_frequency=rotation_frequency, temperature0=temperature0,
                 pressure=pressure, temperature=temperature, ndis=10)
    print(f'frequency_safety_factor: {disk.frequency_safety_factor(rotation_frequency)}')
