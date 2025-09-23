from copy import deepcopy

from numpy import isclose, isnan, nan
from scipy.optimize import fsolve
from substance import Substance
from thermodynamics import adiabatic_index

from src.checks import check_efficiency, check_temperature
from src.config import EPSREL
from src.config import parameters as gtep
from src.nodes.node import GTENode
from src.utils import call_with_kwargs, integral_average


class Compressor(GTENode):
    """Компрессор"""

    def __init__(self, name="Compressor"):
        GTENode.__init__(self, name=name)

        setattr(self, gtep.pipi, nan)
        setattr(self, gtep.effeff, nan)
        setattr(self, gtep.power, nan)

    @property
    def variables(self):
        return {
            gtep.pipi: getattr(self, gtep.pipi),
            gtep.effeff: getattr(self, gtep.effeff),
            gtep.power: getattr(self, gtep.power),
        }

    @property
    def __x0(self) -> dict[str:float]:
        """Начальные приближения"""
        x0 = {
            "outlet_" + gtep.TT: self.inlet.parameters[gtep.TT],
            "outlet_" + gtep.PP: self.inlet.parameters[gtep.PP],
        }

        pipi = getattr(self, gtep.pipi)
        effeff = getattr(self, gtep.effeff)
        power = getattr(self, gtep.power)

        if isnan(power) and not isnan(pipi) and not isnan(effeff):
            x0[gtep.power] = 20 * 10**6  # TODO: model
        elif isnan(effeff) and not isnan(pipi) and not isnan(power):
            x0[gtep.effeff] = 0.8  # TODO: model
        elif isnan(pipi) and not isnan(effeff) and not isnan(power):
            x0[gtep.pipi] = 6  # TODO: model
        elif not isnan(pipi) and not isnan(effeff) and not isnan(power):
            return x0
        else:
            raise "недоопределено"

        return x0

    def equations(self, x: tuple, args: dict) -> tuple:
        """Уравнения"""  # вызыввается в цикле!
        self.outlet.parameters[gtep.TT] = x[0]
        self.outlet.parameters[gtep.PP] = x[1]
        if gtep.power not in args:
            setattr(self, gtep.power, x[2])
        elif gtep.effeff not in args:
            setattr(self, gtep.effeff, x[2])
        elif gtep.pipi not in args:
            setattr(self, gtep.pipi, x[2])
        elif gtep.pipi in args and gtep.effeff in args and gtep.power in args:
            pass

        TT_i = self.inlet.parameters[gtep.TT]
        PP_i = self.inlet.parameters[gtep.PP]
        f_gc = self.inlet.functions[gtep.gc]
        f_Cp = self.inlet.functions[gtep.Cp]

        mf = (self.inlet.parameters[gtep.mf] + self.outlet.parameters[gtep.mf]) / 2
        gc = integral_average(
            f_gc,
            **{
                gtep.TT: (TT_i, self.outlet.parameters[gtep.TT]),
                gtep.PP: (PP_i, self.outlet.parameters[gtep.PP]),
            },
        )[0]
        Cp = integral_average(
            f_Cp,
            **{
                gtep.TT: (TT_i, self.outlet.parameters[gtep.TT]),
                gtep.PP: (PP_i, self.outlet.parameters[gtep.PP]),
            },
        )[0]
        k = adiabatic_index(gc, Cp)

        return (
            getattr(self, gtep.power) - mf * Cp * (self.outlet.parameters[gtep.TT] - TT_i),
            self.outlet.parameters[gtep.TT] - TT_i * (1 + (getattr(self, gtep.pipi) ** ((k - 1) / k) - 1) / getattr(self, gtep.effeff)),
            getattr(self, gtep.pipi) - self.outlet.parameters[gtep.PP] / PP_i,
        )

    def validate(self, epsrel: float = EPSREL) -> bool:
        """Проверка найденного решения"""
        x0 = (
            self.outlet.parameters[gtep.TT],
            self.outlet.parameters[gtep.PP],
            1,  #
        )
        args = {
            gtep.pipi: getattr(self, gtep.pipi),
            gtep.effeff: getattr(self, gtep.effeff),
            gtep.power: getattr(self, gtep.power),
        }

        for null in self.equations(x0, args):
            print(f"{null:.6f}")
        return all(isclose(null, 0, rtol=epsrel) for null in self.equations(x0, args))

    def calculate(self, substance_inlet: Substance, x0=None) -> Substance:
        count_variables = sum((1 if isnan(v) else 0 for v in self.variables.values()))
        if count_variables < 1:
            raise ArithmeticError("система переопределена")
        elif count_variables > 1:
            raise ArithmeticError("система недоопределена")

        GTENode.validate_substance(self, substance_inlet)
        self.inlet = deepcopy(substance_inlet)
        self.outlet = deepcopy(self.inlet)

        self.outlet.parameters[gtep.mf] = self.inlet.parameters[gtep.mf] - self.mass_flow_leak

        pipi = getattr(self, gtep.pipi)
        effeff = getattr(self, gtep.effeff)
        power = getattr(self, gtep.power)

        args = {}
        if not isnan(effeff) and not isnan(pipi) and isnan(power):
            args.update({gtep.effeff: effeff, gtep.pipi: pipi})
        elif not isnan(pipi) and not isnan(power) and isnan(effeff):
            args.update({gtep.pipi: pipi, gtep.power: power})
        elif not isnan(power) and not isnan(effeff) and isnan(pipi):
            args.update({gtep.power: power, gtep.effeff: effeff})
        else:
            raise f"{x0=}"

        fsolve(self.equations, tuple(self.__x0.values()), args)

        self.outlet.parameters[gtep.gc] = call_with_kwargs(self.outlet.functions[gtep.gc], self.outlet.parameters)
        self.outlet.parameters[gtep.Cp] = call_with_kwargs(self.outlet.functions[gtep.Cp], self.outlet.parameters)
        self.outlet.parameters[gtep.DD] = self.outlet.parameters[gtep.PP] / (self.outlet.parameters[gtep.gc] * self.outlet.parameters[gtep.TT])
        self.outlet.parameters[gtep.k] = adiabatic_index(self.outlet.parameters[gtep.gc], self.outlet.parameters[gtep.Cp])

        return self.outlet

    @property
    def is_real(self):
        checks = (
            check_efficiency(getattr(self, gtep.effeff)),
            check_temperature(self.outlet.parameters[gtep.TT]),
        )
        return all(checks)


if __name__ == "__main__":
    for k, v in gtep.items():
        print(f"{k:<10}: {v}")

    substance_inlet = Substance(
        "air",
        parameters={
            gtep.gc: 287,
            gtep.TT: 300,
            gtep.PP: 101_325,
            gtep.mf: 100,
            gtep.Cp: 1006,
            gtep.k: 1.4,
            gtep.c: 0,
        },
        functions={
            gtep.gc: lambda total_temperature: 287,
            gtep.Cp: lambda total_temperature: 1006,
        },
    )

    compressor = Compressor()
    compressor.summary

    setattr(compressor, gtep.pipi, 6)
    setattr(compressor, gtep.effeff, 0.87)
    compressor.mass_flow_leak = 0.03

    compressor.calculate(substance_inlet)

    compressor.summary

    print(f"{compressor.validate() = }")
    print(f"{compressor.is_real = }")
