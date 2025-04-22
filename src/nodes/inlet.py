import sys
from numpy import nan, sqrt, prod
from thermodynamics import atmosphere_standard, Cp, R_gas
from colorama import Fore

sys.path.append('D:/Programming/Python')

from tools import isnum, eps


class Inlet:
    """Входное устройство"""

    def __init__(self, name='Inlet'):
        self.name = name
        self.warnings = {0: set(), 1: set()}  # предупреждения

    def __str__(self):
        return self.name

    def get_variability(self) -> int:
        """Максимальное количество комбинаций варьируемых параметров"""
        return prod([len(value) for key, value in self.__dict__.items()
                     if type(value) is list and len(value) and not key.startswith('_')])

    def set_combination(self, combination, inlet_main) -> None:
        """Установка комбинации"""
        varible_params = [key for key, value in inlet_main.__dict__.items()
                          if type(value) is list and len(value) and not key.startswith('_')]
        positions = [0] * len(varible_params)

        for i in range(combination):
            for j, param in enumerate(varible_params):
                if positions[j] == len(getattr(inlet_main, varible_params[j])) - 1:
                    positions[j] = 0
                else:
                    positions[j] += 1
                    continue

        for j, param in enumerate(varible_params):
            setattr(self, varible_params[j], getattr(inlet_main, varible_params[j])[positions[j]])

    def input_parameters(self) -> list:

        correct_input = False
        while not correct_input:
            self.σ_var = [sigma for sigma in input('σ [] = ').strip().split()]
            if not self.σ_var:
                break
            for sigma in self.σ_var:
                if not isnum(sigma) or float(sigma) < 0 or float(sigma) > 1:
                    print(Fore.RED + 'σ must be a number in [0..1] or empty string!')
                    correct_input = False
                    break
                correct_input = True
        self.σ_var = [float(sigma) for sigma in self.σ_var]

        correct_input = False
        while not correct_input:
            self.g_leak_var = [gleak for gleak in input('g утечки [] = ').split()]
            for gleak in self.g_leak_var:
                if not isnum(gleak) or float(gleak) < 0 or float(gleak) > 1:
                    print(Fore.RED + 'g_leak must be a number in [0..1]!')
                    correct_input = False
                    break
                correct_input = True
        self.g_leak_var = [float(gleak) for gleak in self.g_leak_var]

    @staticmethod
    def σ_inlet(M: float) -> float:
        """Коэффициент сохранения полного давления на входе"""
        if M >= 2:  # работает только для ПВРД!!!
            k = (0.7345454545, 0.4873659674, 0.3040559441, 0.05421911422, 0.003263403263)
            return k[0] + k[1] * M - k[2] * M ** 2 + k[3] * M ** 3 - k[4] * M ** 4
        else:
            return nan

    def get_inlet_parameters(self, **kwargs) -> None:
        """Расчет параметров перед"""
        # Невозмущенные параметры
        self.R_gas1 = R_gas(self.substance, a_ox=getattr(self, 'a_ox1', None), fuel=fuel)
        if not hasattr(self, 'T1') and self.substance.strip().upper() == 'AIR':
            self.T1 = atmosphere_standard(kwargs.get('H', None))['T']
        if not hasattr(self, 'P1') and self.substance.strip().upper() == 'AIR':
            self.P1 = atmosphere_standard(kwargs.get('H', None))['P']
        assert self.T1 and self.P1, f'{type(self).__name__} object has no attributes T1 and P1!'
        self.ρ1 = self.P1 / (self.R_gas1 * self.T1)
        self.g1 = 1
        self.Cp1 = Cp(self.substance, T=self.T1, P=self.P1, a_ox=getattr(self, 'a_ox1', None), fuel=fuel)
        self.k1 = self.Cp1 / (self.Cp1 - self.R_gas1)
        self.Mc1 = 0
        self.c1 = 0
        self.TT1 = self.T1
        self.PP1 = self.P1
        self.ρρ1 = self.ρ1

    def get_outlet_parameters(self, how='all', error=1 / 100, Niter=100, **kwargs) -> None:
        """Расчет параметров после"""
        self.R_gas3 = R_gas(self.substance, a_ox=getattr(self, 'a_ox1', None), fuel=fuel)

        self.Mc3 = kwargs.get('M', None)
        self.c3 = kwargs.get('v', None)
        assert self.Mc3 is not None or self.c3 is not None, \
            f'{Inlet.__name__}.{Inlet.solve.__name__}() method has no arguments M or v!'
        if self.Mc3: self.Mc3 = self.c3 / sqrt(self.k1 * self.R_gas1 * self.T1)
        if self.c3: self.c3 = self.Mc3 * sqrt(self.k1 * self.R_gas1 * self.T1)

        if not hasattr(self, 'σ'): self.σ = self.σ_inlet(self.Mc3)

        self.k3 = self.k1  # нулевое приближение
        for iteration in range(Niter):
            k2 = 0.5 * (self.k1 + self.k3)
            self.TT3 = self.T1 * (1 + (k2 - 1) / 2 * self.Mc3 ** 2)
            self.PP3 = self.P1 * (1 + (k2 - 1) / 2 * self.Mc3 ** 2) ** (k2 / (k2 - 1)) * self.σ
            self.Cp3 = Cp(self.substance, T=self.TT3, P=self.PP3, a_ox=getattr(self, 'a_ox1', None), fuel=fuel)
            if abs(eps('rel', self.Cp3 / (self.Cp3 - self.R_gas3), self.k3)) <= error: break
            self.k3 = self.Cp3 / (self.Cp3 - self.R_gas3)
        else:
            print(f'{Fore.RED}Iteration limit'
                  f'in class {Inlet.__name__} in method {Inlet.get_outlet_parameters.__name__}!')
        self.ρρ3 = self.PP3 / (self.R_gas3 * self.TT3)
        assert hasattr(self, 'g_leak'), f'{type(self).__name__} object has no attribute g_leak!'
        self.g3 = self.g1 - self.g_leak

    def __calculate(self, how='all', error: float = 1 / 100, Niter: int = 100, **kwargs) -> None:
        global fuel

        self.substance = kwargs.get('substance', '')
        assert self.substance, f'{type(self).__name__} object has no attribute substance!'
        fuel = kwargs.get('fuel', '')  # горючее

        self.get_inlet_parameters(**kwargs)
        self.get_outlet_parameters(how=how, error=error, Niter=Niter, **kwargs)

    def solve(self, how='all', error: float = 1 / 100, Niter: int = 100, **kwargs) -> None:
        self.__calculate(how=how, error=error, Niter=Niter, **kwargs)


if __name__ == '__main__':
    inlet = Inlet()

    inlet.σ = 0.98
    inlet.g_leak = 0.005

    inlet.solve(how='cycle', substance='AIR', H=0, M=0)
    for k, v in inlet.__dict__.items(): print(k, '=', v)
