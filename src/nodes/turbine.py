from numpy import nan, inf, isnan, prod
from thermodynamics import Cp, R_gas, η_polytropic, g_cool_BMSTU, mixing_param
from tools import isnum, eps
from colorama import Fore


def find_node_in_scheme(scheme, node2find) -> tuple:
    """Поиск положения узла в схеме"""
    for contour in scheme:
        for i, node in enumerate(scheme[contour]):
            if scheme[contour][i] is node2find: return contour, i


class Turbine:
    """Турбина"""

    def __init__(self, turbine_type, name='Turbine'):
        self.name = name
        self.type = turbine_type  # тип турбины
        self.coolers = []  # охладители

        self.g_cool = nan  # относительный расход охлаждения
        self.PP_cool = nan  # полное давление охлаждения
        self.TT_cool = nan  # полная температура охлаждения

        self.η_mechanical = nan  # механический КПД

        self.warnings = {0: set(), 1: set(), 2: set(), 3: set()}  # предупреждения

    def __str__(self):
        return self.name

    def get_variability(self) -> int:
        """Максимальное количество комбинаций варьируемых параметров"""
        return prod([len(value) for key, value in self.__dict__.items()
                     if type(value) is list and len(value) and not key.startswith('_')])

    def set_combination(self, combination, turbine_main) -> None:
        """Установка комбинации"""
        varible_params = [key for key, value in turbine_main.__dict__.items()
                          if type(value) is list and len(value) and not key.startswith('_')]
        positions = [0] * len(varible_params)

        for i in range(combination):
            for j, param in enumerate(varible_params):
                if positions[j] == len(getattr(turbine_main, varible_params[j])) - 1:
                    positions[j] = 0
                else:
                    positions[j] += 1
                    continue
        for j, param in enumerate(varible_params):
            setattr(self, varible_params[j], getattr(turbine_main, varible_params[j])[positions[j]])

    def input_parameters(self) -> list:
        correct_input = False
        while not correct_input:
            self.PP3_var = [PP3.lower() for PP3 in input('P*3 [Па] = ').split()]
            if not self.PP3_var:
                break
            for PP3 in self.PP3_var:
                if not isnum(PP3) or float(PP3) < 0:
                    print(Fore.RED + 'P*3 must be a number >= 0 or empty string!')
                    correct_input = False
                    break
                correct_input = True
        self.PP3_var = [float(PP3) for PP3 in self.PP3_var]

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
            self.T_lim_var = [Tlim for Tlim in input('T доп [К] = ').split()]
            for tlim in self.T_lim_var:
                if not isnum(tlim) or float(tlim) <= 0:
                    print(Fore.RED + 'T_lim must be a number > 0!')
                    correct_input = False
                    break
                correct_input = True
        self.T_lim_var = [float(tlim) for tlim in self.T_lim_var]

        correct_input = False
        while not correct_input:
            self.η_mechanical_var = [eta for eta in input('η_мех = ').split()]
            for eta in self.η_mechanical_var:
                if not isnum(eta) or float(eta) < 0 or float(eta) > 1:
                    print(Fore.RED + 'η_мех must be a number in [0..1]!')
                    correct_input = False
                    break
                correct_input = True
        self.η_mechanical_var = [float(eta) for eta in self.η_mechanical_var]

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
        if hasattr(self, 'a_ox1'): self.a_ox3 = self.a_ox1 * (1 + self.g_cool)
        self.R_gas3 = R_gas(self.substance, a_ox=getattr(self, 'a_ox1', None), fuel=fuel)

        assert hasattr(self, 'ηη'), f'{type(self).__name__} object has no attribute ηη!'
        self.k3 = self.k1  # нулевое приближение
        if self.L:
            for iteration in range(Niter):
                k2 = 0.5 * (self.k1 + self.k3)
                Cp2 = k2 / (k2 - 1) * R_gas(self.substance, a_ox=self.a_ox3, fuel=fuel)
                assert hasattr(self, 'ηη'), f'{type(self).__name__} object has no attribute ηη!'
                self.ππ = (1 - self.L / (Cp2 * self.TT1 * self.ηη)) ** (k2 / (1 - k2))
                self.TT3 = self.TT1 - self.L / Cp2
                # if self.g_cool: self.TT3 = mixing_param([self.TT3, *TT_cools], [], [])
                self.Cp3 = Cp(self.substance, T=self.TT3, a_ox=self.a_ox3, fuel=fuel)
                if isnan(k2): self.warnings[3].add(f'T*3 < 0 [К]! Не хватает теплоперепада в турбине_{self}!'); return
                if abs(eps('rel', self.Cp3 / (self.Cp3 - self.R_gas3), self.k3)) <= error: break
                self.k3 = self.Cp3 / (self.Cp3 - self.R_gas3)
            else:
                print(Fore.RED + f'Iteration limit in module {Turbine.__name__} in function {self.solve.__name__} L!')
            self.PP3 = self.PP1 / self.ππ
        elif hasattr(self, 'PP3'):
            self.ππ = self.PP1 / self.PP3
            for iteration in range(Niter):
                k2 = 0.5 * (self.k1 + self.k3)
                self.TT3 = self.TT1 * (1 - (1 - self.ππ ** ((1 - k2) / k2)) * self.ηη)
                self.Cp3 = Cp(self.substance, T=self.TT3, a_ox=self.a_ox3, fuel=fuel)
                if abs(eps('rel', self.Cp3 / (self.Cp3 - self.R_gas3), self.k3)) <= error: break
                self.k3 = self.Cp3 / (self.Cp3 - self.R_gas3)
            else:
                print(Fore.RED + f'Iteration limit in module {Turbine.__name__} in function {self.solve.__name__} PP3!')
            self.L = self.Cp1 * self.TT1 - self.Cp3 * self.TT3
        k2 = 0.5 * (self.k1 + self.k3)
        self.ηn = η_polytropic(what='E', ππ=self.ππ, ηη=self.ηη, k=k2)
        self.ρρ3 = self.PP3 / (self.R_gas3 * self.TT3)
        assert hasattr(self, 'g_leak'), f'{type(self).__name__} object has no attribute g_leak!'
        self.g3 = self.g1 - self.g_leak

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

        if self._shafts:
            m = kwargs.get('m', {})
            assert m, f'{type(self).__name__} object has no attribute m!'
            L_C = N = 0
            η = 1
            for shaft in self._shafts:
                for tp, els in shaft.items():
                    if tp == '-':
                        for el in els:
                            if hasattr(el, 'L'):
                                L_C += el.L * m[find_node_in_scheme(scheme, el)[0]]
                            elif hasattr(el, 'N'):
                                N += shaft[-1].N
                            else:
                                raise '!'
                    elif tp == '0':
                        for el in els:
                            if hasattr(el, 'η'):
                                η *= el.η
                    else:
                        raise '!'
            self.L = L_C / η / self.η_mechanical  # удельная работа от компрессоров TODO сделать КПД мех атрибутом вала
            G = kwargs.get('G', inf)
            if isnan(G) or G != 0: G = inf
            L_N = (N / (G / m[
                find_node_in_scheme(scheme, self)[0]])) / η / self.η_mechanical  # удельная работа от нагрузок
            self.L += L_N  # удельная работа от нагрузок
            assert hasattr(self, 'L') or L_N, f'{type(self).__name__} object has no attributes L and/or N!'
            self.L /= m[find_node_in_scheme(scheme, self)[0]]
        assert hasattr(self, 'L') or (hasattr(self, 'N') and hasattr(self, 'G')), \
            f'{type(self).__name__} object has no attributes L and/or (N and G)'

        g_leaks = self.g_leak
        g_fuels = 0
        if scheme:
            for node in range(n):
                g_leaks += scheme[c][node].g_leak  # суммарные утечки до турбины
                g_fuels += scheme[c][node].g_fuel if hasattr(scheme[c][node], 'g_fuel') else 0
        assert hasattr(self, 'T_lim'), f'{type(self).__name__} object has no attribute T_lim!'
        self.g_cool = g_cool_BMSTU(self.TT1, T_lim=self.T_lim)
        self.g_cool = self.g_cool * (1 - g_leaks) / (1 + self.g_cool - g_fuels)

        '''if self.g_cool:
            for node in self.coolers:
                TT_cools = [TT for TT in node.TT3]
                PP_cools = [PP for PP in node.PP3]
                ρρ_cools = [ρρ for ρρ in node.ρρ3]'''

        self.get_outlet_parameters(how=how, error=error, Niter=Niter)

    def solve(self, how='all', error: float = 0.01, Niter: int = 100, **kwargs) -> None:
        self.__calculate(how=how, error=error, Niter=Niter, **kwargs)


if __name__ == '__main__':
    t = Turbine('a')

    t.TT1 = 1600
    t.PP1 = 1_200_000

    t.L = 250000
    t.ηη = 0.91
    t.η_mechanical = 0.99
    t.T_lim = 1200
    t.g_leak = 0.05
    t.a_ox1 = 3

    t.solve(substance='EXHAUST', fuel='KEROSENE')
    for k, v in t.__dict__.items(): print(k, '=', v)
