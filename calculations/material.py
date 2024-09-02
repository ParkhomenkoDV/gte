from colorama import Fore
import numpy as np
from numpy import array, nan, isnan, sqrt, arange, linspace
# from pint import UnitRegistry  # СИ
from scipy import interpolate
import matplotlib.pyplot as plt


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
                    'HRC': 'твердость по Роквеллу [Па]'}

    @classmethod
    def help(cls):
        print(Fore.MAGENTA + 'Material parameters:' + Fore.RESET)
        for k, v in Material.__PARAMETERS.items():
            print('\t' + f'{k}: {v}')
        print(Fore.RED + 'type value must be int, float, array with shape (-1,2) or callable(int | float)' + Fore.RESET)

    def __init__(self, name: str, parameters: dict, composition=None, reference=''):
        assert isinstance(name, str)
        self.__name = name

        assert isinstance(parameters, dict)
        assert all(isinstance(el, str) for el in parameters.keys())  # все ключи - стоки

        for parameter, value in parameters.items():  # есть возможность создавать свои свойства
            if isinstance(value, (int, float)):  # const
                setattr(self, parameter, interpolate.interp1d((273.15,), (value,), kind=0,
                                                              bounds_error=False, fill_value='extrapolate'))
            elif isinstance(value, (tuple, list, np.ndarray)):  # таблица значений
                value = array(value).T
                assert len(value.shape) == 2 and value.shape[0] == 2 and value.shape[1] > 3
                setattr(self, parameter, interpolate.interp1d(value[0], value[1], kind=1,
                                                              bounds_error=False, fill_value='extrapolate'))
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
    def E2G(E: int | float, mu: int | float) -> float:
        """Модуль сдвига Юнга II рода"""
        assert isinstance(E, (int, float)) and isinstance(mu, (int, float))
        return E / (2 * (mu + 1))

    @staticmethod
    def G2E(G: int | float, mu: int | float) -> float:
        """Модуль Юнга I рода"""
        assert isinstance(G, (int, float)) and isinstance(mu, (int, float))
        return 2 * G * (mu + 1)

    @staticmethod
    def HB2HRC(HB: int | float) -> float:
        """Перевод твердости по Бринеллю в твердость по Роквеллу"""
        return nan  # TODO

    @staticmethod
    def HRC2HB(HRC: int | float) -> float:
        """Перевод твердости по Роквеллу в твердость по Бринеллю"""
        return nan  # TODO

    def show(self, temperature, **kwargs) -> None:
        assert isinstance(temperature, (tuple, list, np.ndarray))

        fg = plt.figure(figsize=kwargs.pop("figsize", (8, 8)))
        fg.suptitle(self.__name, fontsize=16, fontweight='bold')
        gs = fg.add_gridspec(1, len(Material.__PARAMETERS))  # строки, столбцы

        for i, param in enumerate(Material.__PARAMETERS):
            if not hasattr(self, param): continue
            x, y = [], []
            for t in temperature:
                if not isnan(getattr(self, param)(t)):
                    x.append(t)
                    y.append(getattr(self, param)(t))
            fg.add_subplot(gs[0, i])
            plt.grid(True)
            plt.xlim(temperature[0], temperature[-1]),
            plt.xticks(temperature)
            plt.xlabel('Temperature', fontsize=12)
            plt.ylabel(param, fontsize=12)
            plt.plot(x, y)

        plt.show()


references = (
    '''Справочник по конструкционным материалам:
    Справочник / Б.Н. Арзамасов, Т.В. Соловьева, С.А. Герасимов и др.;
    Под ред. Б.Н. Арзамасова, Т.В. Соловьевой.
    - М.: Изд-во МГТУ им Н.Э. Баумана, 2006. с.: ил.''',)

materials = list()

materials.append(Material('ХН70МВТЮБ',
                          {
                              'sigma_': array((array((20, 600, 700, 800, 850, 900)) + 273.15,
                                               array((1060, 980, 930, 720, 600, 380)) * 10 ** 6)).T,

                          },
                          reference=references[0] + 'c.412'))


def test():
    """Тестирование"""
    Material.help()

    material = Material('10Х11Н20ТЗР',
                        {
                            "density": 8400,
                            "alpha": interpolate.interp1d((400, 600, 800),
                                                          array((18, 18, 18)) * 10 ** -6,
                                                          kind=1, bounds_error=False, fill_value='extrapolate'),
                            "E": interpolate.interp1d(arange(400, 800 + 1, 100),
                                                      array([1.74, 1.66, 1.57, 1.47, 1.32]) * 10 ** 11,
                                                      kind=3, bounds_error=False, fill_value=nan),
                            "mu": interpolate.interp1d(arange(400, 800 + 1, 100),
                                                       (0.384, 0.379, 0.371, 0.361, 0.347),
                                                       kind=3, bounds_error=False, fill_value='extrapolate'),
                            "heat_capacity": lambda t: 4200,
                            "conductivity": ((0, 16), (100, 18), (200, 19), (400, 19.5)),
                            "smth": 3.1415
                        })
    materials.insert(0, material)

    for material in materials:
        print(material.name)
        print(material.density(500))
        print(material.alpha(500))
        print(material.E(500))
        print(material.heat_capacity(20))
        print(material.conductivity(20))
        print(material.smth(20))
        print(material.__dict__)
        material.show(arange(200, 1_000 + 1, 100))


if __name__ == "__main__":
    import cProfile

    cProfile.run('test()', sort='cumtime')
