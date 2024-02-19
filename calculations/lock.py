from math import nan, pi, cos, sin, tan, radians, sqrt
from numpy import linspace
import matplotlib.pyplot as plt
import pandas as pd
from tools import cot

df_dove_tail = pd.read_excel('ОСТ 1 11031-81.xlsx', header=None)  # ласточкин хвост
df_christmas_tree = pd.read_excel('ОСТ 1 10975-81.xlsx', header=None)  # елка


class Lock:
    """Замок"""

    def __init__(self):
        self.Z = 1  # количество замков

    color_geom = 'orange'
    color_disk = 'blue'
    color_blade = 'red'
    color_info = 'black'


class LockDoveTail(Lock):
    """Замок типа ласточкин хвост"""

    def get_params(self, typesize: int) -> None:
        """Получение параметров замка из БД"""
        print(df_dove_tail)
        params = df_dove_tail.iloc[typesize]
        print(params)
        params = params.to_dict()
        print(params)

    def draw(self, typesize=None):
        plt.figure(figsize=(8, 8))  # размер окна в дюймах
        plt.grid(True)  # сетка
        plt.axis('equal')
        plt.xlim(-(R - H) / 2, (R - H) / 2)
        plt.ylim(R - H, R)
        plt.axis('on')  # СК

        plt.plot([0, 0], [(R - H) * 0.95, R * 1.05], color=self.color_geom, ls='-.')  # ось симметрии

        betta = pi / 2 - a / 2

        # точка пересечения периферии и поверхности лопатки
        x0 = cot(a / 2) / (cot(a / 2) ** 2 + 1) * \
             (H - R + sqrt(2 * H * R * tan(a / 2) ** 2 - H ** 2 * tan(a / 2) ** 2 + R ** 2))

        plt.plot(list(R * cos(alpha) + 0 for alpha in linspace(0, 2 * pi, 360)),
                 list(R * sin(alpha) + 0 for alpha in linspace(0, 2 * pi, 360)),
                 color='black', ls='-')
        plt.plot(list((R - H) * cos(alpha) + 0 for alpha in linspace(0, 2 * pi, 360)),
                 list((R - H) * sin(alpha) + 0 for alpha in linspace(0, 2 * pi, 360)),
                 color='black', ls='-')
        plt.plot([0, l / 2], [R - H + 0.2 / 1000, R - H + 0.2 / 1000], color='black', ls='-', linewidth=2)
        plt.plot([l / 2, l / 2], [R - H + 0.2 / 1000, R - H], color='black', ls='-', linewidth=2)

        plt.plot([0, B / 2], [R - H, R - H], color='black', ls='-', linewidth=1)
        plt.plot(s([0, B / 2]), [R - H, R - H], color='black', ls='-', linewidth=1)

        # основная поверхность
        plt.plot([B / 2, x0], [R - H, cot(a / 2) * x0 + (R - H)],
                 color='black', ls='-', linewidth=1)
        plt.plot(s([B / 2, x0]), [R - H, cot(a / 2) * x0 + (R - H)],
                 color='black', ls='-', linewidth=1)

        # линия
        plt.plot([0, b / 2], [R - H + h, R - H + h], color='black', ls='-', linewidth=1)
        plt.plot(s([0, b / 2]), [R - H + h, R - H + h], color='black', ls='-', linewidth=1)

        plt.show()


class LockChristmasTree(Lock):
    """Замок типа елка"""

    def draw(self):
        pass


def s(x):
    return list(-i for i in x)


if __name__ == '__main__':
    print(df_dove_tail)
    print(df_christmas_tree)
    LockDoveTail().get_params(1)
    #LockDoveTail.draw(typesize=1)

    R = 150 / 1000
    H = 20 / 1000
    l = 2 / 1000
    h = 5 / 1000
    b = 10 / 1000
    B = 30 / 1000
    a = radians(60)
