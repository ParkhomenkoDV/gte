from numpy import prod
from thermodynamics import η_polytropic, Cp, R_gas
from tools import isnum, eps
from colorama import Fore


def find_node_in_scheme(scheme, node2find) -> tuple:
    """Поиск положения узла в схеме"""
    for contour in scheme:
        for i, node in enumerate(scheme[contour]):
            if scheme[contour][i] is node2find: return contour, i


class Compressor:
    """Компрессор"""

    def __init__(self, compressor_type, name='Compressor'):
        self.name = name
        self.type = compressor_type  # тип компрессора
        self.warnings = {0: set(), 1: set(), 2: set(), 3: set()}  # предупреждения

    def __str__(self) -> str:
        return self.name

    def get_variability(self) -> int:
        """Максимальное количество комбинаций варьируемых параметров"""
        return prod([len(value) for key, value in self.__dict__.items()
                     if type(value) is list and len(value) and not key.startswith('_')])

    def set_combination(self, combination, compressor_main) -> None:
        """Установка комбинации"""
        varible_params = [key for key, value in compressor_main.__dict__.items()
                          if type(value) is list and len(value) and not key.startswith('_')]
        positions = [0] * len(varible_params)

        for i in range(combination):
            for j, param in enumerate(varible_params):
                if positions[j] == len(getattr(compressor_main, varible_params[j])) - 1:
                    positions[j] = 0
                else:
                    positions[j] += 1
                    continue
        for j, param in enumerate(varible_params):
            setattr(self, varible_params[j], getattr(compressor_main, varible_params[j])[positions[j]])

    def input_parameters(self) -> list:
        correct_input = False
        while not correct_input:
            self.ππ_var = [pp for pp in input('π* [] = ').split()]
            for pp in self.ππ_var:
                if not isnum(pp) or float(pp) < 1:
                    print(Fore.RED + 'π* must be a number >= 1!')
                    correct_input = False
                    break
                correct_input = True
        self.ππ_var = [float(pp) for pp in self.ππ_var]

        correct_input = False
        while not correct_input:
            self.ηη_var = [eta for eta in input('η* [] = ').split()]
            for eta in self.ηη_var:
                if not isnum(eta) or float(eta) < 0 or float(eta) > 1:
                    print(f'{Fore.RED}η* must be a number in [0..1]!')
                    correct_input = False
                    break
                correct_input = True
        self.ηη_var = [float(eta) for eta in self.ηη_var]

        correct_input = False
        while not correct_input:
            self.g_leak_var = [gleak for gleak in input('g утечки [] = ').split()]
            for gleak in self.g_leak_var:
                if not isnum(gleak) or float(gleak) < 0 or float(gleak) > 1:
                    print(f'{Fore.RED}g_leak is num in [0..1]!')
                    correct_input = False
                    break
                correct_input = True
        self.g_leak_var = [float(gleak) for gleak in self.g_leak_var]

    def get_inlet_parameters(self) -> None:
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
        self.k1 = self.Cp1 / (self.Cp1 - self.R_gas1)

    def get_outlet_parameters(self, how='all', error=1 / 100, Niter=100) -> None:
        """Расчет параметров после"""
        if hasattr(self, 'a_ox1'): self.a_ox3 = self.a_ox1
        self.R_gas3 = R_gas(self.substance, a_ox=getattr(self, 'a_ox1', None), fuel=fuel)
        assert hasattr(self, 'ππ') or hasattr(self, 'PP3'), f'{type(self).__name__} object has no attributes ππ or PP3!'
        if not hasattr(self, 'ππ'): self.ππ = self.PP3 / self.PP1
        self.PP3 = self.PP1 * self.ππ
        self.k3 = self.k1  # нулевое приближение
        for iteration in range(Niter):
            k2 = 0.5 * (self.k1 + self.k3)
            assert hasattr(self, 'ηη'), f'{type(self).__name__} object has no attributes ηη!'
            self.TT3 = self.TT1 * (1 + (self.ππ ** ((k2 - 1) / k2) - 1) / self.ηη)
            self.Cp3 = Cp(self.substance, T=self.TT3, P=self.PP3, a_ox=getattr(self, 'a_ox1', None), fuel=fuel)
            if abs(eps('rel', self.Cp3 / (self.Cp3 - self.R_gas3), self.k3)) <= error: break
            self.k3 = self.Cp3 / (self.Cp3 - self.R_gas3)
        else:
            print(f'{Fore.RED}Iteration limit'
                  f'in class {Compressor.__name__} in method {Compressor.get_outlet_parameters.__name__}!')
        self.ρρ3 = self.PP3 / (self.R_gas3 * self.TT3)
        assert hasattr(self, 'g_leak'), f'{type(self).__name__} object has no attribute g_leak!'
        self.g3 = self.g1 - self.g_leak
        self.ηn = η_polytropic(what='C', ππ=self.ππ, ηη=self.ηη, k=k2)
        self.L = self.Cp3 * self.TT3 - self.Cp1 * self.TT1

    def __calculate(self, how='all', error: float = 1 / 100, Niter: int = 100, **kwargs) -> None:

        global c, n, scheme, fuel
        c = n = None
        scheme = kwargs.get('scheme', {})

        if scheme:
            c, n = find_node_in_scheme(scheme, self)
            self.substance = scheme[c][n - 1].substance
        else:
            self.substance = kwargs.get('substance', '')  # рабочее тело
        assert self.substance, f'{type(self).__name__} object has no attribute substance!'
        fuel = kwargs.get('fuel', '')  # горючее

        self.get_inlet_parameters()
        self.get_outlet_parameters(how=how, error=error, Niter=Niter)

    def solve(self, how='all', error=0.01, Niter=100, **kwargs) -> None:
        self.__calculate(how=how, error=error, Niter=Niter, **kwargs)


if __name__ == '__main__':
    compressor = Compressor('a')

    compressor.TT1 = 300
    compressor.PP1 = 115666

    compressor.ηη = 0.86
    compressor.ππ = 6
    compressor.g_leak = 0.05

    compressor.solve(substance='AIR')
    for k, v in compressor.__dict__.items(): print(k, '=', v)
