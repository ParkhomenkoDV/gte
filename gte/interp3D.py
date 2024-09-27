from numpy import linspace, meshgrid, random, nan, sin, cos, pi, array
import pandas as pd
from scipy.interpolate import RectBivariateSpline, interp2d

EXCEL_Cp_air = pd.read_excel('Теплоёмкость воздуха.xlsx', header=None)
# print(EXCEL_Cp_air)

'''
Cp_air = interp2d(EXCEL_Cp_air.T[0].iloc[1:],  # T
                  EXCEL_Cp_air[0].iloc[1:],  # P
                  EXCEL_Cp_air.iloc[1:, 1:],  # Cp_Air
                  kind='linear', fill_value=nan)  # None
'''

# Создание случайных точек на поверхности
x = array(EXCEL_Cp_air.T[0].iloc[1:])
print(x, '\n', len(x), '\n')
y = array(EXCEL_Cp_air[0].iloc[1:])
print(y, '\n', len(y), '\n')
z = array(EXCEL_Cp_air.iloc[1:, 1:]).T
print(z, '\n', len(z), len(z[0]), '\n')

# Создание сетки для интерполяции
xi, yi = meshgrid(linspace(0, 1, 50), linspace(0, 1, 40))
print(xi, '\n', yi, '\n')

# Визуализация результатов
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z, c='r', marker='o')
ax.plot_surface(xi, yi, F(xi, yi))
plt.show()

# create some sample data
x = asarray(EXCEL_Cp_air.T[0].iloc[1:])
y = asarray(EXCEL_Cp_air[0].iloc[1:])
z = asmatrix(EXCEL_Cp_air.iloc[1:, 1:]).T
print(len(z))
print(len(random.rand(29, 10)))

# create the spline function
spline = RectBivariateSpline(x, y, z)

print(spline(485, 2 * 100_000)[0][0])
