from copy import deepcopy

import numpy as np
from mathematics import integral_average
from numpy import isclose, isnan, nan
from scipy import integrate
from scipy.optimize import fsolve
from substance import Substance
from thermodynamics import adiabatic_index

from src.config import EPSREL
from src.config import parameters as gtep
from src.errors import EFFICIENCY_ERROR
from src.nodes.node import GTENode
from src.utils import check_efficiency

# TODO подправить integral_average так, чтобы считалось


def average_integral(f, *borders) -> float:
    """Среднеинтегральное значение"""
    for border in borders:
        assert isinstance(border, (tuple, list))
        for b in border:
            assert isinstance(b, (int, float, np.number)), f"{type(b)}"

    if borders[0][0] == borders[0][1]:
        return f(borders[0][0])
    else:
        return integrate.quad(f, borders[0][0], borders[0][1])[0] / (
            borders[0][1] - borders[0][0]
        )


class Compressor(GTENode):
    """Компрессор"""

    def __init__(self, name="Compressor"):
        GTENode.__init__(self, name=name)

        self.pipi = nan
        self.effeff = nan
        self.power = nan

        """
        setattr(self, gtep.pipi, nan)
        setattr(self, gtep.effeff, nan)
        setattr(self, gtep.power, nan)
        """

    """
    def __setattr__(self, name, value):
        if name == "pipi":
            assert isinstance(value, (float, int))
            assert 1 <= value
        elif name == gtep.effeff:
            assert isinstance(value, (float, int))
            assert check_efficiency(value)
        elif name == gtep.power:
            assert isinstance(value, (float, int))
            assert 0 <= value
        return super().__setattr__(name, value)
    """

    @property
    def variables(self):
        return {gtep.pipi: self.pipi, gtep.effeff: self.effeff, gtep.power: self.power}

    @property
    def __x0(self) -> dict[str:float]:
        """Начальные приближения"""
        x0 = {
            "outlet_" + gtep.TT: self.inlet.parameters[gtep.TT],
            "outlet_" + gtep.PP: self.inlet.parameters[gtep.PP],
        }

        if isnan(self.power) and not isnan(self.pipi) and not isnan(self.effeff):
            x0[gtep.power] = 20 * 10**6  # TODO: model
        elif isnan(self.effeff) and not isnan(self.pipi) and not isnan(self.power):
            x0[gtep.effeff] = 0.8  # TODO: model
        elif isnan(self.pipi) and not isnan(self.effeff) and not isnan(self.power):
            x0[gtep.pipi] = 6  # TODO: model
        elif not isnan(self.pipi) and not isnan(self.effeff) and not isnan(self.power):
            return x0
        else:
            raise "недоопределено"

        return x0

    def equations(self, x: tuple, args: dict) -> tuple:
        """Уравнения"""  # вызыввается в цикле!
        self.outlet.parameters[gtep.TT] = x[0]
        self.outlet.parameters[gtep.PP] = x[1]
        if gtep.power not in args:
            self.power = x[2]
        elif gtep.effeff not in args:
            self.effeff = x[2]
        elif gtep.pipi not in args:
            self.pipi = x[2]
        elif gtep.pipi in args and gtep.effeff in args and gtep.power in args:
            pass

        TT_i = self.inlet.parameters[gtep.TT]
        PP_i = self.inlet.parameters[gtep.PP]
        f_gc = self.inlet.functions[gtep.gc]
        f_Cp = self.inlet.functions[gtep.Cp]

        mf = (self.inlet.parameters[gtep.mf] + self.outlet.parameters[gtep.mf]) / 2
        gc = average_integral(f_gc, (TT_i, self.outlet.parameters[gtep.TT]))
        Cp = average_integral(f_Cp, (TT_i, self.outlet.parameters[gtep.TT]))
        k = adiabatic_index(gc, Cp)

        return (
            self.power - mf * Cp * (self.outlet.parameters[gtep.TT] - TT_i),
            self.outlet.parameters[gtep.TT]
            - TT_i * (1 + (self.pipi ** ((k - 1) / k) - 1) / self.effeff),
            self.pipi - self.outlet.parameters[gtep.PP] / PP_i,
        )

    def validate(self, epsrel: float = EPSREL) -> bool:
        """Проверка найденного решения"""
        x0 = (
            self.outlet.parameters[gtep.TT],
            self.outlet.parameters[gtep.PP],
            1,
        )
        args = {gtep.pipi: self.pipi, gtep.effeff: self.effeff, gtep.power: self.power}

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

        self.outlet.parameters[gtep.mf] = (
            self.inlet.parameters[gtep.mf] - self.mass_flow_leak
        )

        args = {}
        if not isnan(self.effeff) and not isnan(self.pipi) and isnan(self.power):
            args.update({gtep.effeff: self.effeff, gtep.pipi: self.pipi})
        elif not isnan(self.pipi) and not isnan(self.power) and isnan(self.effeff):
            args.update({gtep.pipi: self.pipi, gtep.power: self.power})
        elif not isnan(self.power) and not isnan(self.effeff) and isnan(self.pipi):
            args.update({gtep.power: self.power, gtep.effeff: self.effeff})
        else:
            raise f"{x0=}"

        fsolve(self.equations, tuple(self.__x0.values()), args)

        assert check_efficiency(self.effeff), AssertionError(
            EFFICIENCY_ERROR.format(self.effeff)
        )

        self.outlet.parameters[gtep.gc] = self.outlet.functions[gtep.gc](
            total_temperature=self.outlet.parameters[gtep.TT]
        )
        self.outlet.parameters[gtep.Cp] = self.outlet.functions[gtep.Cp](
            total_temperature=self.outlet.parameters[gtep.TT]
        )
        self.outlet.parameters[gtep.DD] = self.outlet.parameters[gtep.PP] / (
            self.outlet.parameters[gtep.gc] * self.outlet.parameters[gtep.TT]
        )
        self.outlet.parameters[gtep.k] = adiabatic_index(
            self.outlet.parameters[gtep.gc], self.outlet.parameters[gtep.Cp]
        )

        return self.outlet


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
            gtep.c: 100,
        },
        functions={
            gtep.gc: lambda total_temperature: 287,
            gtep.Cp: lambda total_temperature: 1006,
        },
    )

    compressor = Compressor()
    compressor.summary

    compressor.effeff = 0.87
    compressor.pipi = 6
    compressor.mass_flow_leak = 0.03

    compressor.calculate(substance_inlet)

    compressor.summary

    print(f"{compressor.validate() = }")
