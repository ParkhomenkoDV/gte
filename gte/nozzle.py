from numpy import nan, inf, isnan, sqrt, prod
from thermodynamics import T0, Cp, R_gas, η_polytropic
from tools import isnum, eps
from colorama import Fore


def find_node_in_scheme(scheme, node2find) -> tuple:
    """Поиск положения узла в схеме"""
    for contour in scheme:
        for i, node in enumerate(scheme[contour]):
            if scheme[contour][i] is node2find: return contour, i


class Nozzle:
    """Сопло"""

    def __init__(self, name='Nozzle'):
        self.name = name
        self.warnings = {0: set(), 1: set(), 2: set(), 3: set()}  # предупреждения

    def get_variability(self) -> int:
        """Максимальное количество комбинаций варьируемых параметров"""
        return prod([len(value) for key, value in self.__dict__.items()
                     if type(value) is list and len(value) and not key.startswith('_')])

    def set_combination(self, combination, nozzle_main) -> None:
        """Установка комбинации"""
        varible_params = [key for key, value in nozzle_main.__dict__.items()
                          if type(value) is list and len(value) and not key.startswith('_')]
        positions = [0] * len(varible_params)

        for i in range(combination):
            for j, param in enumerate(varible_params):
                if positions[j] == len(getattr(nozzle_main, varible_params[j])) - 1:
                    positions[j] = 0
                else:
                    positions[j] += 1
                    continue
        for j, param in enumerate(varible_params):
            setattr(self, varible_params[j], getattr(nozzle_main, varible_params[j])[positions[j]])

    def input_parameters(self) -> list:

        while True:  # главные характеристики
            if not self.ππ_var or not self.PP3_var:
                correct_input = False
                while not correct_input:
                    self.ππ_var = [pp for pp in input('π* [] = ').split()]
                    if not self.ππ_var:
                        break
                    for pp in self.ππ_var:
                        if not isnum(pp) or float(pp) < 1:
                            print(Fore.RED + 'π* must be a number >= 1 or empty string!')
                            correct_input = False
                            break
                        correct_input = True
                self.ππ_var = [float(pp) for pp in self.ππ_var]

            if not self.ππ_var:
                correct_input = False
                while not correct_input:
                    self.PP3_var = [PP3 for PP3 in input('P*3 [Па] = ').split()]
                    if not self.PP3_var:
                        break
                    for PP3 in self.PP3_var:
                        if isnum(self.PP3) and 0 <= float(self.PP3):
                            print(Fore.RED + 'P*3 must be a number >= 0 or empty string!')
                            correct_input = False
                            break
                        correct_input = True
                self.PP3_var = [float(PP3) for PP3 in self.PP3_var]

            if self.ππ_var or self.PP3_var:
                break

        correct_input = False
        while not correct_input:
            self.ηη_var = [eta for eta in input('η* [] = ').split()]
            for eta in self.ηη_var:
                if not isnum(eta) or float(eta) < 0 or float(eta) > 1:
                    print(Fore.RED + 'η* must be a number in [0..1]!')
                    correct_input = False
                    break
                correct_input = True
        self.ηη_var = [float(eta) for eta in self.ηη_var]

        correct_input = False
        while not correct_input:
            self.v__var = [v for v in input('v_ [] = ').split()]
            for v in self.v__var:
                if not isnum(v) or float(v) < 0 or float(v) > 1:
                    print(Fore.RED + 'v_ must be a number in [0..1]!')
                    correct_input = False
                    break
                correct_input = True
        self.v__var = [float(v) for v in self.v__var]

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

    def get_outlet_parameters(self, error=0.01, Niter=100, **kwargs):
        """Расчет параметров после"""
        if hasattr(self, 'a_ox1'): self.a_ox3 = self.a_ox1
        self.R_gas3 = R_gas(self.substance, a_ox=getattr(self, 'a_ox1', None), fuel=fuel)
        self.TT3 = self.TT1

        assert hasattr(self, 'ππ') or hasattr(self, 'PP3'), f'{type(self).__name__} object has no attributes ππ and PP3'
        if hasattr(self, 'PP3'): self.ππ = self.PP1 / self.PP3

        if self.ππ < 1: self.warnings[3].add('π* < 1!'); return

        self.PP3 = self.PP1 / self.ππ
        self.ρρ3 = self.PP3 / (self.R_gas3 * self.TT3)
        assert hasattr(self, 'g_leak'), f'{type(self).__name__} object has no attribute g_leak!'
        self.g3 = self.g1 - self.g_leak
        self.Cp3 = Cp(self.substance, T=self.TT3, a_ox=getattr(self, 'a_ox3', None), fuel=fuel)
        self.k3 = self.Cp3 / (self.Cp3 - self.R_gas3)

        R_gas2 = self.R_gas1
        Cp2 = 0.5 * (self.Cp1 + self.Cp3)  # TODO решить через интеграл
        k2 = Cp2 / (Cp2 - R_gas2)
        assert hasattr(self, 'ηη'), f'{type(self).__name__} object has no attribute ηη!'
        self.ηn = η_polytropic(what='E', ππ=self.ππ, ηη=self.ηη, k=k2)

        self.ππ_max = ((k2 + 1) / 2) ** (k2 / (k2 - 1))
        if self.ππ < self.ππ_max:
            self.c3 = self.v_ * sqrt(2 * Cp2 * self.TT3 * (1 - self.ππ ** ((1 - k2) / k2)))
            self.T3 = self.TT3 - (self.k3 - 1) / (self.k3 * self.R_gas3) * self.c3 ** 2 / 2
            self.R_ = self.c3 * self.g3 - scheme[c][0].c3 * scheme[c][0].g3
        else:
            self.c3 = self.v_ * sqrt(2 * Cp2 * self.TT3)
            self.T3 = self.TT3 - (self.k3 - 1) / (self.k3 * self.R_gas3) * self.c3 ** 2 / 2
            self.R_ = self.c3 * self.g3 - scheme[c][0].c3 * scheme[c][0].g3
            self.R_ += self.R_gas3 * self.T3 / self.c3 * (1 - self.PP3 * self.ππ_max / self.PP1) * self.g3

    def __calculate(self, how='all', error: float = 1 / 100, Niter: int = 100, **kwargs) -> None:

        global c, n, scheme, fuel, m
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

    def solve(self, how='all', error: float = 0.01, Niter: int = 100, **kwargs) -> None:
        self.__calculate(how=how, error=error, Niter=Niter, **kwargs)


if __name__ == '__main__':
    n = Nozzle()

    n.TT1 = 1500
    n.PP1 = 600_000
    n.a_ox1 = 3

    n.solve(substance='EXHAUST', fuel='КЕРОСИН')
    for k, v in n.__dict__.items(): print(k, '=', v)
