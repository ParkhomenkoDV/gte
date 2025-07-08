from math import nan, prod
import numpy as np


class Node:
    """Базовый класс узла"""

    def __init__(self):
        self.init()

    def init(self) -> None:
        """Инициализация параметров"""
        self.ox_inlet, self.ox_outlet = nan, nan
        self.R_gas_inlet, self.R_gas_outlet = nan, nan
        self.k_inlet, self.k_outlet = nan, nan
        self.Cp_inlet, self.Cp_outlet = nan, nan
        self.G_inlet, self.G_outlet = nan, nan
        self.g_inlet, self.g_outlet = nan, nan
        self.g_leak = nan

        # полные параметры
        self.TT_inlet, self.TT_outlet = nan, nan
        self.PP_inlet, self.PP_outlet = nan, nan
        self.RoRo_inlet, self.RoRo_outlet = nan, nan

        # статические параметры
        self.T_inlet, self.T_outlet = nan, nan
        self.P_inlet, self.P_outlet = nan, nan
        self.Ro_inlet, self.Ro_outlet = nan, nan
        self.Mc_inlet, self.Mc_outlet = nan, nan
        self.c_inlet, self.c_outlet = nan, nan

        self.warnings = set()  # предупреждения

    def get_variability(self) -> int:
        """Максимальное количество комбинаций варьируемых параметров"""
        return prod(
            [
                len(value)
                for key, value in self.__dict__.items()
                if isinstance(value, (list, tuple, np.array))
                and len(value)
                and not key.startswith("_")
            ]
        )

    @property
    def inlet_parameters(self) -> dict[str:float]:
        return {k: v for k, v in self.__dict__.items() if k.endswith("inlet")}
    
    @property
    def outlet_parameters(self) -> dict[str:float]:
        return {k: v for k, v in self.__dict__.items() if k.endswith("outlet")}

    '''def get_inlet_parameters(self) -> None:
        """Расчет параметров перед"""
        if scheme:
            if hasattr(scheme[c][n - 1], 'a_ox3'): self.a_ox1 = scheme[c][n - 1].a_ox3
            self.TT1 = scheme[c][n - 1].TT3
            self.PP1 = scheme[c][n - 1].PP3
            self.g1 = scheme[c][n - 1].g3
        assert hasattr(self, 'TT1'), f'{type(self).__name__} object has no attribute TT1!'
        assert hasattr(self, 'PP1'), f'{type(self).__name__} object has no attribute PP1!'
        if not hasattr(self, 'g1'): self.g1 = 1
        self.R_gas1 = R_gas(self.substance, a_ox=getattr(self, 'a_ox1', None), fuel=fuel)
        self.ρρ1 = self.PP1 / (self.R_gas1 * self.TT1)
        self.Cp1 = Cp(self.substance, T=self.TT1, P=self.PP1, a_ox=getattr(self, 'a_ox1', None), fuel=fuel)
        self.k1 = self.Cp1 / (self.Cp1 - self.R_gas1)'''


if __name__ == "__main__":
    bn = Node()
    print(bn.__dict__)
    print(bn.inlet_parameters)
    print(bn.outlet_parameters)

