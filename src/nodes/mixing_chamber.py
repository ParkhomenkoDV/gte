from numpy import nan, inf, isnan, prod
from thermodynamics import Cp, R_gas
from tools import isnum
from colorama import Fore


def find_node_in_scheme(scheme, node2find) -> tuple:
    """Поиск положения узла в схеме"""
    for contour in scheme:
        for i, node in enumerate(scheme[contour]):
            if scheme[contour][i] is node2find: return contour, i


class MixingChamber:
    """Камера смешения"""

    def __init__(self, name='MixingChamber'):
        self.name = name
        self.warnings = {0: set(), 1: set(), 2: set(), 3: set()}  # предупреждения

        self.nodes_add = []  # узел подмешивания
        self.TT_add = []  # полная температура подмешивания
        self.PP_add = []  # полное давление подмешивания
        self.ρρ_add = []  # полная плотность подмешивания
        self.g_add = []  # относительный массовый расход подмешивания
        self.Cp_add = []

    def __str__(self) -> str:
        return self.name

    def get_variability(self):
        """Максимальное количество комбинаций варьируемых параметров"""
        return prod([len(value) for key, value in self.__dict__.items()
                     if type(value) is list and len(value) and not key.startswith('_')])

    def set_combination(self, combination, mixingchamber_main):
        """Установка комбинации"""
        varible_params = [key for key, value in mixingchamber_main.__dict__.items()
                          if type(value) is list and len(value) and not key.startswith('_')]
        positions = [0] * len(varible_params)

        for i in range(combination):
            for j, param in enumerate(varible_params):
                if positions[j] == len(getattr(mixingchamber_main, varible_params[j])) - 1:
                    positions[j] = 0
                else:
                    positions[j] += 1
                    continue
        for j, param in enumerate(varible_params):
            setattr(self, varible_params[j], getattr(mixingchamber_main, varible_params[j])[positions[j]])

    def input_parameters(self):
        while True:
            n_add = input('nodes_add []: ').strip()
            if isnum(n_add, type_num='int') and 1 <= int(n_add):
                self.node_add = [None] * int(n_add)
                break
            else:
                print(Fore.RED + 'nodes_add is int num >= 1!')

        for i, _ in enumerate(self.nodes_add):
            while True:
                TT_add = input('T*_add []').strip()
                if isnum(TT_add) and 0 < float(TT_add):
                    self.TT_add.append(float(TT_add))
                    break
                else:
                    print(Fore.RED + 'T*_add is float num > 0!')

            while True:
                PP_add = input('P*_add []').strip()
                if isnum(PP_add) and 0 < float(PP_add):
                    self.PP_add.append(float(PP_add))
                    break
                else:
                    print(Fore.RED + 'P*_add is float num > 0!')

            while True:
                ρρ_add = input('ρ*_add []').strip()
                if isnum(ρρ_add) and 0 < float(ρρ_add):
                    self.ρρ_add.append(float(ρρ_add))
                    break
                else:
                    print(Fore.RED + 'ρ*_add is float num > 0!')

            while True:
                g_add = input('g_add []').strip()
                if isnum(g_add) and 0 < float(g_add):
                    self.g_add.append(float(g_add))
                    break
                else:
                    print(Fore.RED + 'g_add is float num > 0!')

    def get_inlet_parameters(self, *args, **kwargs):
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

    def get_add_parameters(self, *args, **kwargs):
        if self._nodes_add:
            for node in self._nodes_add:
                self.TT_add.append(node.TT3)
                self.PP_add.append(node.PP3)
                self.ρρ_add.append(node.ρρ3)
                self.g_add.append(node.g3)
                self.Cp_add.append(node.Cp3)
                self.substance = node.substance if node.substance == 'EXHAUST' else self.substance
        else:
            self.nodes_add = [None] * len(self.TT_add)

    def get_outlet_parameters(self, *args, **kwargs):
        """Расчет параметров после"""
        m = kwargs.get('m', {})
        m_self = 1 if not m else m[c]
        self.TT3 = self.Cp1 * self.TT1 * self.g1 * m_self
        self.PP3 = self.PP1 * self.g1 * m_self
        self.ρρ3 = self.ρρ1 * self.g1 * m_self
        self.g3 = (self.g1 - self.g_leak) * m_self
        gm = 0  # sum of complex g*m
        cpgm = 0
        for i, node in enumerate(self.nodes_add):
            m_node = 1 if not m else m[find_node_in_scheme(scheme, node)[0]]
            self.TT3 += self.Cp_add[i] * self.TT_add[i] * self.g_add[i] * m_node
            self.PP3 += self.PP_add[i] * self.g_add[i] * m_node
            self.ρρ3 += self.ρρ_add[i] * self.g_add[i] * m_node
            self.g3 += self.g_add[i] * m_node
            gm += self.g_add[i] * m_node
            cpgm += self.Cp_add[i] * self.g_add[i] * m_node
        self.TT3 /= self.g1 * m_self + cpgm
        self.PP3 /= self.g1 * m_self + gm
        self.ρρ3 /= self.g1 * m_self + gm
        self.g3 /= self.g1 * m_self  # относительный расход отнесенный ко входу в соответствующий контур

        self.R_gas3 = self.PP3 / (self.ρρ3 * self.TT3)

        if scheme:
            gmm = 0
            for node_add in self.nodes_add:
                c_node_add, n_node_add = find_node_in_scheme(scheme, node_add)

                for i, node in enumerate(scheme[c_node_add][:n_node_add]):
                    gmm = 0
                    if i == n_node_add - 1:
                        gmm += self.g_add * (m[c_node_add] / m_self)  # complex of sum g_add*(m_from/m_to)

            self.a_ox3 = getattr(self, 'a_ox1', 1) * (1 + gmm)
            self.Cp3 = Cp(self.substance, T=self.TT1, P=self.PP1, a_ox=self.a_ox3, fuel=fuel)

    def __calculate(self, how='all', **kwargs) -> None:

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
        assert fuel, f'{MixingChamber.__name__}.{MixingChamber.solve.__name__}() method has no arguments fuel!'

        self.get_inlet_parameters()
        self.get_add_parameters()
        self.get_outlet_parameters(how=how)

    def solve(self, how='all', error=0.01, Niter=100, **kwargs):
        self.__calculate(how=how, error=error, Niter=Niter, **kwargs)


if __name__ == '__main__':
    mc = MixingChamber()

    mc.TT_add = [700]
    mc.PP_add = [2 * 10 ** 5]
    mc.ρρ_add = [4]
    mc.g_add = [1]

    mc.g_leak = 0.05

    mc.solve(substance='AIR', TT3=400, PP3=1.8 * 10 ** 5)
    for k, v in mc.__dict__.items(): print(k, '=', v)
