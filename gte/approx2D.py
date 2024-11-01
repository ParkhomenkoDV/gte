from numpy import array, polyfit, polyval
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Исходные данные
x = array([1, 2, 3, 4, 5])
y = array([1, 4, 9, 16, 25])


# ----------------------------------------------------------------

# Определение функции для аппроксимации
def func(x, a, b, c): return a * x ** 2 + b * x + c


popt, pcov = curve_fit(func, x, y)  # Аппроксимация данных
for i in range(len(popt)): print('p' + str(len(popt) - 1 - i), round(popt[i], 8))
# print('pcov:\n', pcov)
'''Массив popt содержит оптимальные значения параметров функции, 
которые были найдены в результате аппроксимации данных. 

Массив pcov содержит ковариационную матрицу, 
которая показывает, насколько точно определены значения параметров. 
Если диагональные элементы этой матрицы близки к нулю, 
то это означает, что значения параметров определены с высокой точностью. 
Если же диагональные элементы большие, 
то это может указывать на неопределенность в значениях параметров.'''

# График исходных данных и аппроксимационной кривой
plt.plot(x, y, 'o')
plt.plot(x, func(x, *popt), '-')
plt.show()

# ----------------------------------------------------------------


p = polyfit(x, y, 2)  # коэффициенты полинома
for i in range(len(p)): print('p' + str(len(p) - 1 - i), round(p[i], 8))

# График исходных данных и аппроксимационной кривой
plt.plot(x, y, 'o')
plt.plot(x, polyval(p, x), '-')
plt.show()
