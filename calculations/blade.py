"""
Список литературы:

[1] = Елисеев Ю.С., Крымов В.В., Манушин Э.А. и др.
Конструирование и расчет на прочность турбомашин ГТ и КУ:
Учебник для студентов вузов / Под общей ред. М.И. Осипова. – М.: МГТУ, 2009

[2] = Иноземцев А.А. Динамика и прочность авиационных двигателей и энергетических установок:
учебник / А.А. Иноземцев, М.А. Нихамкин, В.Л. Сандрацкий. – М.: Машиностроение, 2008. – Т.2. – 368 с.; Т.4. – 192 с.

[3] = Костюк А.Г. Динамика и прочность турбомашин:
Учебник для вузов. – 2-е изд., перераб. и доп. – М.: Издательство МЭИ, 2000. – 480 с.

[4] = Щегляев А.В. Паровые турбины. – М.: Энергоатомиздат, 1993. – В двух книгах.

[5] = Малинин Н.Н. Прочность турбомашин. – М.: Машиностроние, 1962. – 291 с.

[6] = Колебания / В.С. Чигрин. - Пособие по лабораторному практикуму - Рыбинск: РГАТА, 2005. -20 с.
"""

import numpy as np
from numpy import nan, isnan, pi, sqrt, array, linspace, arange
from scipy import interpolate, integrate
import matplotlib.pyplot as plt

from material import Material


class Blade:
    """Лопатка/винт/лопасть"""

    @classmethod
    def version(cls):
        version = '1.3'
        print('3D построение')
        return version

    def __init__(self, material: Material, points: tuple | list | np.ndarray, amount: int = 1) -> None:
        # проверка на тип данных material
        assert isinstance(material, Material)

        assert isinstance(points, (tuple, list, np.ndarray))
        assert all(isinstance(point, (tuple, list, np.ndarray)) for point in points)
        assert all(len(point) == 3 for point in points)
        assert all('float' in type(cord).__name__ or 'int' in type(cord).__name__ for point in points for cord in point)

        assert isinstance(amount, int) and 1 <= amount

        self.material = material
        self.points = np.array(sorted(points, key=lambda p: p[2], reverse=False))
        self.amount = amount

        self.height = 0  # TODO

    def show(self, **kwargs):
        """"""
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection

        plt.figure(figsize=kwargs.pop('figsize', (9, 9)))
        ax = plt.axes(projection='3d')
        ax.axis('equal')
        *_, zs = self.points.T
        for i in np.unique(zs):
            xx, yy, zz = [], [], []
            for x, y, z in self.points:
                if i == z:
                    xx.append(x)
                    yy.append(y)
                    zz.append(z)
            vertices = [list(zip(xx, yy, zz))]
            poly = Poly3DCollection(vertices, alpha=0.8)
            ax.add_collection3d(poly)
        ax.set_title(kwargs.pop('title', 'Blade'), fontsize=14, fontweight='bold')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        plt.show()

    def solve(self, omega, *args, **kwargs):
        pass

    def natural_frequencies(self, max_k: int):
        """Частота собственных колебаний [6]"""
        self.I = self.b * self.c ** 3 / 12  # момент инерции сечения
        self.F = self.b * self.c  # площадь сечения

        if "крепление в заделке":
            k = array([1.875, 4.694] + [pi * (i - 0.5) for i in range(3, max_k + 1, 1)])
        elif "крепление шарнирное":
            k = array([0] + [pi(i - 0.75) for i in range(2, max_k + 1, 1)])
        else:
            raise

        f = self.material.E(0) * self.I / (self.material.density(0) * np.mean(self.F))
        f = sqrt(f)
        f *= k ** 2 / (2 * pi * self.h ** 2)
        return f, '1/s'

    def campbell_diagram(self, max_rotation_frequency: int, k=arange(1, 11, 1), **kwargs):
        """Диаграмма Кэмпбелла [6]"""

        rotation_frequency = arange(0, max_rotation_frequency + 1, 1) * (2 * pi)  # перевод из рад/c в 1/c=об/c=Гц
        # динамическая частота колебаний вращающегося диска
        # self.R[0] = радиус корня
        f = sqrt(self.natural_frequencies(max(k))[0] ** 2 + (1 + self.R[0] / self.h) * rotation_frequency ** 2)
        resonance = set()  # резонансные частоты [1/с]

        plt.figure(figsize=kwargs.pop('figsize', (8, 8)))
        plt.title('Campbell diagram', fontsize=14, fontweight='bold')
        for k_ in k:
            plt.plot([0, rotation_frequency[-1]], [0, rotation_frequency[-1] * k_],
                     color='orange', linestyle='solid', linewidth=1)
            plt.text(rotation_frequency[-1], rotation_frequency[-1] * k_, f'k{k_}',
                     fontsize=12, ha='left', va='center')
            if k_ ** 2 - (1 + self.R[0] / self.h) >= 0:
                x0 = f[0] / sqrt(k_ ** 2 - (1 + self.R[0] / self.h))  # f
                if not isnan(x0) and x0 <= rotation_frequency[-1]:
                    resonance.add(round(x0, 6))
                    plt.scatter(x0, k_ * x0, color='red')
        plt.plot(rotation_frequency, f, color='black', linestyle='solid', linewidth=2, label='f')
        plt.xlabel('Frequency [1/s]', fontsize=12)
        plt.ylabel('Frequency [1/s]', fontsize=12)
        plt.grid(True)
        plt.legend(fontsize=12)
        plt.show()

        return sorted(list(resonance), reverse=False), '1/s'


if __name__ == '__main__':
    print(Blade.version())

    blades = list()

    if 1:
        material = Material('ЖС-40', {'density': 8600})

        points = [
            (0.00, +0.000, 0.00),
            (0.20, +0.005, 0.00),
            (0.50, +0.010, 0.00),
            (1.00, +0.000, 0.00),
            (0.50, -0.010, 0.00),
            (0.00, +0.000, 0.00),

            (0.00, +0.000, 0.10),
            (0.20, +0.005, 0.10),
            (0.50, +0.010, 0.10),
            (1.00, +0.000, 0.10),
            (0.50, -0.010, 0.10),
            (0.00, +0.000, 0.10),

            (0.00, +0.000, 0.20),
            (0.20, +0.005, 0.20),
            (0.50, +0.010, 0.20),
            (1.00, +0.000, 0.20),
            (0.50, -0.010, 0.20),
            (0.00, +0.000, 0.20),
        ]

        blade = Blade(material=material, points=points)
        blades.append(blade)

    for blade in blades:
        blade.show()
