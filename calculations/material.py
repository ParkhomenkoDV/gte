import numpy as np
from numpy import array, full, nan, isnan, pi, sqrt, arange, linspace
# from pint import UnitRegistry  # СИ
from scipy import interpolate
import matplotlib.pyplot as plt


class Material:
    def __init__(self, name: str, parameters: dict, composition=None, reference=''):
        assert isinstance(name, str)
        self.name = name

        assert isinstance(parameters, dict)
        assert all(map(lambda k: isinstance(k, str), parameters.keys()))
        assert all(map(lambda v: isinstance(v, (int, float)) or callable(v), parameters.values()))

        density = parameters.pop("density", nan)
        if type(density) in (int, float):
            self.density = lambda T: density
        else:
            self.density = density

        alpha = parameters.pop("alpha", nan)
        if type(alpha) in (int, float):
            self.alpha = lambda T: alpha
        else:
            self.alpha = alpha

        E = parameters.pop("E", nan)
        if type(E) in (int, float):
            self.E = interpolate.interp1d(list(range(0, 500, 100)), [E] * 5,
                                          kind='cubic', fill_value='extrapolate')
        else:
            self.E = E

        mu = parameters.pop("mu", nan)
        if type(mu) in (int, float):
            self.mu = interpolate.interp1d(list(range(0, 500, 100)), [mu] * 5,
                                           kind='cubic', fill_value='extrapolate')
        else:
            self.mu = mu

        sigma_temp = parameters.pop("sigma_temp", nan)
        if type(sigma_temp) in (int, float):
            self.sigma_temp = interpolate.interp1d(list(range(0, 500, 100)), [sigma_temp] * 5,
                                                   kind='cubic', fill_value='extrapolate')
        else:
            self.sigma_temp = sigma_temp

    @staticmethod
    def G(E: int | float, mu: int | float) -> float:
        """Модуль сдвига Юнга II рода"""
        assert isinstance(E, (int, float)) and isinstance(mu, (int, float))
        return E / (2 * (mu + 1))

    @staticmethod
    def E(G: int | float, mu: int | float) -> float:
        """Модуль Юнга I рода"""
        assert isinstance(G, (int, float)) and isinstance(mu, (int, float))
        return 2 * G * (mu + 1)

    def show(self, **kwargs) -> None:
        fg = plt.figure(figsize=kwargs.pop("figsize", (8, 8)))
        fg.suptitle(self.name, fontsize=16, fontweight='bold')
        gs = fg.add_gridspec(1, 4)  # строки, столбцы

        T = list(range(200, 2_000 + 1, 50))

        for i, param in enumerate(('density', 'alpha', 'E', 'mu')):
            x, y = [], []
            for t in T:
                if not isnan(getattr(self, param)(t)):
                    x.append(t)
                    y.append(getattr(self, param)(t))
            fg.add_subplot(gs[0, i])
            plt.grid(True)
            plt.xlim(T[0], T[-1]),
            plt.xticks(T)
            plt.xlabel('Temperature [K]', fontsize=12)
            plt.ylabel(param, fontsize=12)
            plt.plot(x, y)

        plt.show()


def test():
    """Тестирование"""
    material = Material('10Х11Н20ТЗР',
                        {
                            "density": 8400,
                            "alpha": interpolate.interp1d([400, 600, 800],
                                                          [18 * 10 ** -6, 18 * 10 ** -6, 18 * 10 ** -6],
                                                          kind='linear', bounds_error=False, fill_value='extrapolate'),
                            "E": interpolate.interp1d(list(range(400, 800 + 1, 100)),
                                                      np.array([1.74, 1.66, 1.57, 1.47, 1.32]) * 10 ** 11,
                                                      kind='cubic', bounds_error=False, fill_value=np.nan),
                            "mu": interpolate.interp1d(list(range(400, 800 + 1, 100)),
                                                       [0.384, 0.379, 0.371, 0.361, 0.347],
                                                       kind='cubic', bounds_error=False, fill_value='extrapolate')
                        })
    print(material.density(500))
    print(material.alpha(500))
    print(material.E(500))
    material.show()


if __name__ == "__main__":
    import cProfile

    cProfile.run('test()', sort='cumtime')
