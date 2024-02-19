import scipy as sp
import numpy as np
import time


def f1(x, y):
    coef = 1
    return coef * (y ** 2 + x ** 2)


# дискретность разбития данных
discr = 30

# дискретные значения x и y
X = np.linspace(0, 100, discr)
Y = np.linspace(0, 100, discr)

# сетка данных
xg, yg = np.meshgrid(X, Y, indexing='ij')
data = f1(xg, yg)

# вывод данных
for x in X:
    for y in Y:
        print(f'{f1(x, y):8.1f}', end=' ')
    print()

# интерполируемая функция
interp_func = sp.interpolate.RegularGridInterpolator((X, Y), data)

# интегрирование
start = time.monotonic()
print(f'Результат интегрирования обычной функции: {sp.integrate.dblquad(f1, 10, 20, 50, 60)[0]:.2f}. '
      f'Потрачено {time.monotonic() - start:.2f} секунд')

start = time.monotonic()
print(
    f'Результат интегрирования интерполированных данных: {sp.integrate.dblquad(lambda x, y: interp_func((x, y)), 10, 20, 50, 60)[0]:.2f}. '
    f'Потрачено {time.monotonic() - start:.2f} секунд')
