import numpy as np
from numpy import nan, isnan
from scipy import interpolate
import matplotlib.pyplot as plt


class Material:
    def __init__(self, name: str, parameters: dict):
        assert type(name) == str
        self.name = name

        assert type(parameters) is dict
        assert all(map(lambda k: type(k) is str, parameters.keys()))
        assert all(map(lambda v: type(v) in (int, float, tuple, list, np.ndarray) or callable(v), parameters.values()))
        for key, value in parameters.items():
            if type(value) in (int, float) or callable(value):
                setattr(self, key, value)
            elif type(value) in (tuple, list, np.ndarray):
                assert all(map(lambda i: type(i) in (tuple, list, np.ndarray), value))
                assert len(value) >= 3
                assert all(map(lambda i: len(i) == 2, value))
                # TODO: подумать над экстраполяцией
                setattr(self, key, interpolate.interp1d([v[0] for v in value], [v[1] for v in value],
                                                        kind='cubic', fill_value="extrapolate"))
            else:
                raise

    @staticmethod
    def G(E: int | float, mu: int | float) -> float:
        """Модуль Юнга II рода"""
        return E / (2 * (mu + 1))

    def show(self, **kwargs) -> None:
        plt.figure(figsize=kwargs.pop("figsize", (8, 8)))
        plt.title(self.name, fontsize=14, fontweight="bold")

        plt.show()


if __name__ == "__main__":
    material = Material('10Х11Н20ТЗР',
                        {
                            "density": 8400,
                            "alpha": lambda t: 18 * 10 ** -6 if 400 <= t <= 800 else nan,
                            "E": [(400, 1.74 * 10 ** 11),
                                  (500, 1.66 * 10 ** 11),
                                  (600, 1.57 * 10 ** 11),
                                  (700, 1.47 * 10 ** 11),
                                  (800, 1.32 * 10 ** 11)],
                            "mu": [(400, 0.384),
                                   (500, 0.379),
                                   (600, 0.371),
                                   (700, 0.361),
                                   (800, 0.347)]})
    material.show()
