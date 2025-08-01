from copy import deepcopy

from node import GTENode
from numpy import nan
from substance import Substance


class Inlet(GTENode):
    """Входное устройство"""

    def __init__(self, name: str = "Inlet"):
        GTENode.__init__(self, name=name)

        self.loss_pressure = nan

    def set_combination(self, combination, inlet_main) -> None:
        """Установка комбинации"""
        varible_params = [
            key
            for key, value in inlet_main.__dict__.items()
            if type(value) is list and len(value) and not key.startswith("_")
        ]
        positions = [0] * len(varible_params)

        for i in range(combination):
            for j, param in enumerate(varible_params):
                if positions[j] == len(getattr(inlet_main, varible_params[j])) - 1:
                    positions[j] = 0
                else:
                    positions[j] += 1
                    continue

        for j, param in enumerate(varible_params):
            setattr(
                self,
                varible_params[j],
                getattr(inlet_main, varible_params[j])[positions[j]],
            )

    @staticmethod
    def get_loss_pressure(mach: float) -> float:
        """Коэффициент сохранения полного давления на входе"""
        if mach >= 2:  # работает только для ПВРД!
            coefs = (
                0.7345454545,
                0.4873659674,
                -0.3040559441,
                0.05421911422,
                -0.003263403263,
            )
            return sum(coefs[i] * mach**i for i in range(len(coefs)))
        else:
            return nan

    '''
    def get_outlet_parameters(self, error=1 / 100, Niter=100, **kwargs) -> None:
        """Расчет параметров после"""

        self.Mc3 = kwargs.get("M", None)
        self.c3 = kwargs.get("v", None)
        assert self.Mc3 is not None or self.c3 is not None, (
            f"{Inlet.__name__}.{Inlet.solve.__name__}() method has no arguments M or v!"
        )
        if self.Mc3:
            self.Mc3 = self.c3 / sqrt(self.k1 * self.R_gas1 * self.T1)
        if self.c3:
            self.c3 = self.Mc3 * sqrt(self.k1 * self.R_gas1 * self.T1)

        if not hasattr(self, "σ"):
            self.σ = self.loss_pressure(self.Mc3)

        self.k3 = self.k1  # нулевое приближение
        for iteration in range(Niter):
            k2 = 0.5 * (self.k1 + self.k3)
            self.TT3 = self.T1 * (1 + (k2 - 1) / 2 * self.Mc3**2)
            self.PP3 = (
                self.P1 * (1 + (k2 - 1) / 2 * self.Mc3**2) ** (k2 / (k2 - 1)) * self.σ
            )
            self.Cp3 = Cp(
                self.substance,
                T=self.TT3,
                P=self.PP3,
                a_ox=getattr(self, "a_ox1", None),
                fuel=fuel,
            )
            if abs(eps("rel", self.Cp3 / (self.Cp3 - self.R_gas3), self.k3)) <= error:
                break
            self.k3 = self.Cp3 / (self.Cp3 - self.R_gas3)
        else:
            print(
                f"{Fore.RED}Iteration limit"
                f"in class {Inlet.__name__} in method {Inlet.get_outlet_parameters.__name__}!"
            )
        self.ρρ3 = self.PP3 / (self.R_gas3 * self.TT3)
        assert hasattr(self, "g_leak"), (
            f"{type(self).__name__} object has no attribute g_leak!"
        )
        self.g3 = self.g1 - self.g_leak
        '''

    def calculate(self, substance_inlet: Substance) -> Substance:
        self.substance_inlet = deepcopy(substance_inlet)
        self.substance_outlet = deepcopy(self.substance_inlet)

        return self.substance_outlet


if __name__ == "__main__":
    s = Substance("air", parameters={"T": 300, "P": 101_325, "mach": 0})

    inlet = Inlet()
    inlet.calculate(s)
    print(inlet.summary)
    print(inlet.__dict__)

    exit()

    inlet.Mc = 0
    inlet.sigma = 0.98
    inlet.g_leak = 0.005
    print(inlet.__dict__)

    inlet.solve(how="cycle", substance="AIR", H=0, M=0)
    for k, v in inlet.__dict__.items():
        print(k, "=", v)
