"""
Список литературы:

[1] = Марочник сталей и сплавов.
2-е изд., исправл. и доп. / Зубченко А.С., Колосков М.М., Каширский Ю.В. и др. Под ред. А.С. Зубченко.
М.: Машиностроение, 2003. 784 с.
"""

from colorama import Fore
import pandas as pd
import numpy as np
from numpy import array, nan, isnan, sqrt, arange, linspace
# from pint import UnitRegistry  # СИ
from scipy import interpolate
import matplotlib.pyplot as plt

T0 = 273.15  # абсолютный температурный ноль
M = 10 ** 6  # приставка Мега

hardness = pd.read_excel('libraries/Hardness.xlsx').drop(['d10mm'], axis=1)  # [1, c.784]


class Material:
    __PARAMETERS = {'density': 'плотность [кг/м^3]',
                    'alpha': 'коэффициент линейного расширения [1/К]',
                    'E': 'модуль Юнга I рода [Па]',
                    'G': 'модуль (сдвига) Юнга II рода [Па]',
                    'mu': 'коэффициент Пуассона []',
                    'sigma_t': 'предел текучести [Па]',
                    'sigma_s': 'предел прочности [Па]',
                    'conductivity': 'теплопроводность [Вт/м/К]',
                    'heat_capacity': 'теплоемкость [Дж/кг/К]',
                    'KCU': 'ударная вязкость [Дж/м^2]',
                    'HB': 'твердость по Бринеллю [Па]',
                    'HRC': 'твердость по Роквеллу [Па]',
                    'HV': 'твердость по Виккерсу [Па]',
                    '*': 'частный параметр [...]'}

    @classmethod
    def __version__(cls):
        version = '3.0'
        print(version)

    @classmethod
    def help(cls):
        print(Fore.CYAN + 'Material parameters:' + Fore.RESET)
        for k, v in Material.__PARAMETERS.items():
            print('\t' + f'{k}: {v}')
        print(Fore.RED + 'type value must be int, float, array with shape (-1,2) or callable(int | float)' + Fore.RESET)

    def __init__(self, name: str, parameters: dict, composition=None, reference='', kind: int = 1, fill_value=nan):
        assert isinstance(name, str)
        self.__name = name

        assert isinstance(kind, int) and 0 <= kind <= 3
        assert isnan(fill_value) or fill_value == 'extrapolate'

        assert isinstance(parameters, dict)
        assert all(isinstance(el, str) for el in parameters.keys())  # все ключи - стоки

        for parameter, value in parameters.items():  # есть возможность создавать свои свойства
            if isinstance(value, (int, float)):  # const
                setattr(self, parameter, interpolate.interp1d((273.15,), (value,), kind=0,  # для константы надо 0
                                                              bounds_error=False, fill_value='extrapolate'))
            elif isinstance(value, (tuple, list, np.ndarray)):  # таблица значений
                value = array(value).T
                assert len(value.shape) == 2 and value.shape[0] == 2 and value.shape[1] > 3
                setattr(self, parameter, interpolate.interp1d(value[0], value[1], kind=kind,
                                                              bounds_error=False, fill_value=fill_value))
            elif callable(value):  # функция
                try:
                    value(273.15)  # проверка на вызов от численного значения
                    setattr(self, parameter, value)
                except Exception:
                    print(f'parameter "{parameter}" has not callable value!')
            else:
                raise Exception('type of values parameters is in (int, float) or callable(int, float)')

        if composition is None:
            self.composition = composition
        elif isinstance(composition, dict):
            assert all(isinstance(el, str) for el in composition.keys())
            assert all(isinstance(el, (float, int)) for el in composition.values())
            assert all(0 <= el <= 1 for el in composition.values())
            self.composition = composition
        else:
            raise ValueError('type(composition) must be dict!')

        assert isinstance(reference, str)
        self.reference = reference

    @property
    def name(self) -> str:
        return self.__name

    @name.setter
    def name(self, name: str) -> None:
        assert isinstance(name, str)
        self.__name = name

    @name.deleter
    def name(self) -> None:
        raise

    @staticmethod
    def jung(**kwargs) -> float:
        """Модуль Юнга I и II рода"""
        mu = kwargs.pop('mu', None)
        E = kwargs.pop('E', None)
        G = kwargs.pop('G', None)
        assert isinstance(mu, (float, int, np.number)) and 0 < mu
        if isinstance(E, (float, int, np.number)):
            assert 0 < E
            return E / (2 * (mu + 1))
        elif isinstance(G, (float, int, np.number)):
            assert 0 < G
            return 2 * G * (mu + 1)
        else:
            raise Exception('isinstance(E, (float, int, np.number)) or isinstance(G, (float, int, np.number))')

    @staticmethod
    def hardness(h: str, value: float | int | np.number) -> dict[str:float]:
        """Перевод твердости"""
        assert h in hardness.columns, hardness.columns
        assert isinstance(value, (float, int, np.number)) and 0 <= value
        group = hardness.groupby(h).mean()  # замена дубликатов их средним значением
        result = dict()
        for column in group.columns:
            temp = group.dropna(subset=column)  # удаление nan
            func = interpolate.interp1d(temp.index, temp[column], kind=3, fill_value=nan, bounds_error=False)
            result[column] = float(func(value))
        return result

    def show(self, temperature: tuple | list | np.ndarray, **kwargs) -> None:
        assert isinstance(temperature, (tuple, list, np.ndarray))
        parameters = [k for k, v in self.__dict__.items() if callable(v)]

        fg = plt.figure(figsize=kwargs.pop("figsize", (4 * len(parameters), 8)))
        fg.suptitle(self.__name, fontsize=16, fontweight='bold')
        gs = fg.add_gridspec(1, len(parameters))  # строки, столбцы

        for i, param in enumerate(parameters):
            if not hasattr(self, param): continue
            xy = array([(t, getattr(self, param)(t))
                        for t in linspace(temperature[0], temperature[-1], 1_000, endpoint=True)
                        if not isnan(getattr(self, param)(t))])
            fg.add_subplot(gs[0, i])
            plt.grid(True)
            plt.xlim(temperature[0], temperature[-1]),
            plt.xticks(temperature)
            plt.xlabel('temperature', fontsize=12)
            plt.ylabel(param, fontsize=12)
            plt.plot(*xy.T)

        plt.show()


references = (
    '''Справочник по конструкционным материалам:
    Справочник / Б.Н. Арзамасов, Т.В. Соловьева, С.А. Герасимов и др.;
    Под ред. Б.Н. Арзамасова, Т.В. Соловьевой.
    - М.: Изд-во МГТУ им Н.Э. Баумана, 2006. с.: ил.''',
    '''
    '''
)

materials = list()

materials.append(Material('ХН70МВТЮБ',
                          {
                              'sigma_s': array((array((20, 600, 700, 800, 850, 900)) + T0,
                                                array((1060, 980, 930, 720, 600, 380)) * 10 ** 6)).T,
                              'sigma_t': array((array((20, 600, 700, 800, 850, 900)) + T0,
                                                array((560, 550, 530, 450, 400, 220)) * 10 ** 6)).T,
                              'KCU': array((array((700, 750, 800, 850)) + 273.15,
                                            array((0.8, 0.7, 0.6, 0.7)) * 10 ** 6)).T,
                              'sigma_100': array((array((650, 700, 800, 850)) + T0,
                                                  array((620, 480, 250, 180)) * 10 ** 6)).T,
                              'sigma_200': array((array((650, 700, 800, 850)) + T0,
                                                  array((600, 420, 230, 230)) * 10 ** 6)).T, },
                          reference=references[0] + ', c. 412-413'))
materials.append(Material('ХН80ТБЮ',
                          {
                              'sigma_s': array((array((29, 500, 600, 630, 650, 700)) + T0,
                                                array((960, 1000, 830, 790, 700, 680)) * M)).T,
                              'sigma_t': array((array((29, 500, 600, 630, 650, 700)) + T0,
                                                array((650, 610, 600, 600, 550, 500)))).T,
                              'KCU': interpolate.interp1d(array((29, 650, 700)) + T0, array((0.7, 1.0, 1.2)) * M,
                                                          kind=2, bounds_error=False, fill_value='extrapolate'),
                              'sigma_1000': interpolate.interp1d(array((650, 700)) + T0, array((450, 280)) * M,
                                                                 kind=1, bounds_error=False, fill_value=nan),
                              'sigma_5000': interpolate.interp1d(array((650, 700)) + T0, array((320, 220)) * M,
                                                                 kind=1, bounds_error=False, fill_value=nan),
                              'sigma_10000': interpolate.interp1d(array((650, 700)) + T0, array((280, 170)) * M,
                                                                  kind=1, bounds_error=False, fill_value=nan), },
                          reference=references[0] + ', c. 413'))
materials.append(Material('ХН70ВМТЮ',
                          {
                              'sigma_s': array((array((20, 500, 600, 650, 700, 750, 800, 900, 950, 1000)) + T0,
                                                array((1030, 1020, 970, 990, 890, 710, 570, 300, 140, 80)) * M)).T,
                              'sigma_t': array((array((20, 500, 600, 650, 700, 750, 800, 900, 950, 1000)) + T0,
                                                array((670, 640, 600, 600, 580, 580, 500, 280, 120, 70)) * M)).T,
                              'KCU': array((array((20, 500, 600, 650, 700, 750, 800)) + T0,
                                            array((0.8, 0.9, 0.9, 0.8, 0.9, 0.85, 1.05)) * M)).T,
                              'sigma_100': interpolate.interp1d(array((750, 800)) + T0, array((360, 240)) * M,
                                                                kind=1, fill_value=nan, bounds_error=False),
                              'sigma_1000': interpolate.interp1d(array((600, 650, 700, 750, 800)) + T0,
                                                                 array((650, 550, 310, 250, 175)) * M,
                                                                 kind=2, fill_value=nan, bounds_error=False),
                              'sigma_2000': interpolate.interp1d(array((600, 650, 700, 750)) + T0,
                                                                 array((600, 400, 270, 200)) * M,
                                                                 kind=3, fill_value=nan, bounds_error=False),
                              'sigma_3000': interpolate.interp1d(array((650, 750, 800)) + T0,
                                                                 array((470, 215, 145)) * M,
                                                                 kind=2, fill_value=nan, bounds_error=False),
                              'sigma_5000': interpolate.interp1d(array((600, 650, 750, 800)) + T0,
                                                                 array((560, 440, 185, 130)) * M,
                                                                 kind=3, fill_value=nan, bounds_error=False),
                              'sigma_10000': interpolate.interp1d(array((600, 650, 750, 800)) + T0,
                                                                  array((530, 385, 170, 125)) * M,
                                                                  kind=3, fill_value=nan, bounds_error=False),
                              'sigma_20000': interpolate.interp1d(array((600, 650, 750)) + T0,
                                                                  array((500, 340, 190)) * M,
                                                                  kind=2, fill_value=nan, bounds_error=False)},
                          reference=references[0] + ', c. 414-415'))


def test():
    """Тестирование"""
    Material.__version__()
    Material.help()

    material = Material('test',  # тестируемый материал
                        {
                            "density": 8400,  # int
                            "alpha": interpolate.interp1d((400, 600, 800),
                                                          array((18, 18, 18)) * 10 ** -6,
                                                          kind=1, bounds_error=False, fill_value='extrapolate'),
                            "E": interpolate.interp1d(arange(400, 800 + 1, 100),
                                                      array([1.74, 1.66, 1.57, 1.47, 1.32]) * 10 ** 11,
                                                      kind=3, bounds_error=False, fill_value=nan),
                            "mu": interpolate.interp1d(arange(400, 800 + 1, 100),
                                                       (0.384, 0.379, 0.371, 0.361, 0.347),
                                                       kind=3, bounds_error=False, fill_value='extrapolate'),
                            "heat_capacity": lambda t: 4200,  # lambda
                            "conductivity": ((0, 16), (100, 18), (200, 19), (400, 19.5)),  # tuple
                            "HV": array(((0, 16), (100, 18), (200, 19), (400, 19.5))),  # array
                            "smth": 3.1415  # float
                        })
    materials.insert(0, material)

    temperature, t = arange(200, 1_200 + 1, 50), 700
    for material in materials:
        print(Fore.MAGENTA + material.name + Fore.RESET)
        for k, v in material.__dict__.items():
            if callable(v):
                print('\t' + f'{k}({t}): {v(t)}')
        material.show(temperature)

    for column in hardness.columns:
        for h in range(0, 1_000 + 1, 10):
            print(f'"{column}": {h}, {Material.hardness(column, h)}')


if __name__ == "__main__":
    import cProfile

    cProfile.run('test()', sort='cumtime')
