'''import numpy
import openpyxl
import scipy

ORIGIN = 1
N = 5

print(10**-3)


def f(t, x):
    return numpy.exp(-x * t) / t ** N


# Теплоемкость воздуха
def Cp_air(PP, TT):
    wb = openpyxl.reader.excel.load_workbook(filename="Теплоёмкость воздуха.xlsx",
                                             read_only=True,
                                             data_only=True)
    wb.active = 1 - ORIGIN  # первый лист
    sheet = wb.active

    P = []
    T = []
    Cp = []
    Cpp = []

    for row in range(2, 30 + ORIGIN):
        P.append(sheet[row][1 - ORIGIN].value)
    for col in range(2 - ORIGIN, 58 + ORIGIN):
        T.append(sheet[1][col].value)
    for col in range(2 - ORIGIN, 58 + ORIGIN):
        for row in range(2, 30 + ORIGIN):
            Cpp.append(sheet[row][col].value)
        Cp.append(Cpp)
        Cpp = []
    del Cpp
    return scipy.interpolate.interp2d(P, T, Cp, kind='cubic')(PP, TT)[0]


# (0.20000000000002294, 1.2239614263187945e-08)

#print(Cp_air(10 ** 6, 500))
#print(Cp_air(2 * 10 ** 6, 600))
#print(f(45,7))
print(scipy.integrate.nquad(f, [[500, 600], [10 ** 6, 2 * 10 ** 6]]))'''

import random
import time
import numpy as np


from multiprocessing import Pool, cpu_count

print(cpu_count())
total_iterations = 101

# Разделение общего количества итераций между процессами
iterations_per_process = total_iterations // cpu_count()
print(iterations_per_process)

# Создание пула процессов
pool = Pool(cpu_count())
# Запуск процессов с использованием метода map. Каждому процессу передается диапазон итераций
pool.map(process_function, [(i + 1) * iterations_per_process for i in range(cpu_count())])
# Закрытие пула процессов
pool.close()
pool.join()
