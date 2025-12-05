from typing import Dict

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RBFInterpolator

try:
    from .config import parameters as gtep
except ImportError:
    from config import parameters as gtep


class Characteristic:
    """Характеристика лопаточной машины"""

    KERNELS = ("linear",)

    __slots__ = (gtep.rf, gtep.mf, gtep.effeff, gtep.pipi, "__efficiency", "__pi")

    def __init__(self, data, kernel: str = "linear", smoothing: float = 0.0):
        setattr(self, gtep.rf, [])
        setattr(self, gtep.mf, [])
        setattr(self, gtep.effeff, [])
        setattr(self, gtep.pipi, [])

        # Валидация данных
        assert isinstance(data, (tuple, list)), TypeError(f"{type(data)=} must be tuple")
        assert len(data) >= 4, ValueError(f"{len(data)=} must be >= 4")
        for d in data:
            assert isinstance(d, dict), TypeError(f"{type(d)=} must be dict")
            for k, value in d.items():
                assert k in (gtep.rf, gtep.mf, gtep.effeff, gtep.pipi)
                assert isinstance(value, (float, int)), TypeError(f"{type(value)=} must be float")
                assert value > 0, ValueError(f"{value=} must be > 0")

            getattr(self, gtep.rf).append(d[gtep.rf])
            getattr(self, gtep.mf).append(d[gtep.mf])
            getattr(self, gtep.effeff).append(d[gtep.effeff])
            getattr(self, gtep.pipi).append(d[gtep.pipi])

        assert kernel in self.KERNELS, ValueError(f"{kernel=} not in {self.KERNELS}")
        assert isinstance(smoothing, float), TypeError(f"{type(smoothing)=} must be float")

        points = np.column_stack((getattr(self, gtep.rf), getattr(self, gtep.mf)))

        self.__efficiency = RBFInterpolator(points, getattr(self, gtep.effeff), kernel=kernel, smoothing=smoothing)
        self.__pi = RBFInterpolator(points, getattr(self, gtep.pipi), kernel=kernel, smoothing=smoothing)

    def efficiency(self, rf: float, mf: float) -> float:
        """КПД"""
        assert isinstance(rf, (float, int)), TypeError(f"{type(rf)=} must be float")
        assert isinstance(mf, (float, int)), TypeError(f"{type(mf)=} must be float")
        return float(self.__efficiency([(rf, mf)])[0])

    def pi(self, rf: float, mf: float) -> float:
        """Степень повышения давления"""
        assert isinstance(rf, (float, int)), TypeError(f"{type(rf)=} must be float")
        assert isinstance(mf, (float, int)), TypeError(f"{type(mf)=} must be float")
        return float(self.__pi([(rf, mf)])[0])

    def show(self, figsize=(4, 8)):
        """Визуализация характеристики лопаточной машины"""
        fg = plt.figure(figsize=figsize)
        gs = fg.add_gridspec(2, 1)  # строки, столбцы
        plt.suptitle("Characteristic", fontsize=16)

        # efficiency
        fg.add_subplot(gs[0, 0])
        plt.title("efficiency")
        plt.xlabel(gtep.rf, fontsize=12)
        plt.ylabel(gtep.mf, fontsize=12)
        z = [[self.efficiency(rf, mf) for mf in getattr(self, gtep.mf)] for rf in getattr(self, gtep.rf)]
        heatmap = plt.imshow(z, cmap="bwr")
        cbar = plt.colorbar(heatmap)

        # pi
        fg.add_subplot(gs[1, 0])
        plt.title("pipi")
        plt.xlabel(gtep.rf, fontsize=12)
        plt.ylabel(gtep.mf, fontsize=12)
        z = [[self.pi(rf, mf) for mf in getattr(self, gtep.mf)] for rf in getattr(self, gtep.rf)]
        heatmap = plt.imshow(z)
        cbar = plt.colorbar(heatmap)

        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    data = (
        # ГУ
        {gtep.rf: 1.05, gtep.mf: 0.98, gtep.effeff: 1.02, gtep.pipi: 1.085},
        {gtep.rf: 1.05, gtep.mf: 1.13, gtep.effeff: 0.90, gtep.pipi: 0.970},
        {gtep.rf: 0.80, gtep.mf: 0.84, gtep.effeff: 0.90, gtep.pipi: 0.71},
        {gtep.rf: 0.80, gtep.mf: 0.59, gtep.effeff: 0.90, gtep.pipi: 0.82},
    )

    character = Characteristic(data)

    for k in character.__slots__:
        if k.startswith("_"):
            continue
        print(k, getattr(character, k))

    for d in data:
        print("eff:", character.efficiency(d[gtep.rf], d[gtep.mf]))
        print("pi: ", character.pi(d[gtep.rf], d[gtep.mf]))

    character.show()
