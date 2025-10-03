from copy import deepcopy

from numpy import isnan, nan
from scipy.optimize import root
from substance import Substance
from thermodynamics import adiabatic_index, gas_const, heat_capacity_at_constant_pressure

from src.checks import check_efficiency, check_temperature
from src.config import EPSREL, substance_mixing
from src.config import parameters as gtep
from src.nodes.node import GTENode
from src.utils import call_with_kwargs, integral_average


class Turbine(GTENode):
    """Турбина"""

    def __init__(self, name="Turbine"):
        GTENode.__init__(self, name=name)

        setattr(self, gtep.pipi, nan)
        setattr(self, gtep.effeff, nan)
        setattr(self, gtep.power, nan)

    @property
    def variables(self) -> dict[str:float]:
        return {
            gtep.pipi: getattr(self, gtep.pipi),
            gtep.effeff: getattr(self, gtep.effeff),
            gtep.power: getattr(self, gtep.power),
        }

    @property
    def _x0(self) -> dict[str:float]:
        """Начальные приближения"""
        x0 = {
            f"outlet_{gtep.TT}": self.inlet.parameters[gtep.TT] * 0.7,
            f"outlet_{gtep.PP}": self.inlet.parameters[gtep.PP] / 2.5,
        }
        for k, v in self.variables.items():
            if not isnan(v):
                continue
            if k == gtep.pipi:
                x0[k] = 3  # TODO: model or formula
            elif k == gtep.effeff:
                x0[k] = 1.0
            elif k == gtep.power:
                x0[k] = 24 * 10**6  # TODO: model or formula
        return x0

    def equations(self, x, args: dict) -> tuple:
        """Уравнения"""
        self.outlet.parameters[gtep.TT] = x[0]
        self.outlet.parameters[gtep.PP] = x[1]
        idx = 2
        for k in self.variables:
            if k not in args:
                setattr(self, k, x[idx])
                idx += 1

        eo_i = self.inlet.parameters[gtep.eo]
        TT_i = self.inlet.parameters[gtep.TT]
        PP_i = self.inlet.parameters[gtep.PP]

        eo_o = self.outlet.parameters[gtep.eo]

        mf = (self.inlet.parameters[gtep.mf] + self.outlet.parameters[gtep.mf]) / 2
        gc, _ = integral_average(
            self.inlet.functions[gtep.gc],
            **{
                gtep.eo: (eo_i, eo_o),
                gtep.TT: (TT_i, self.outlet.parameters[gtep.TT]),
                gtep.PP: (PP_i, self.outlet.parameters[gtep.PP]),
            },
        )
        Cp, _ = integral_average(
            self.inlet.functions[gtep.Cp],
            **{
                gtep.eo: (eo_i, eo_o),
                gtep.TT: (TT_i, self.outlet.parameters[gtep.TT]),
                gtep.PP: (PP_i, self.outlet.parameters[gtep.PP]),
            },
        )
        k = adiabatic_index(gc, Cp)

        return (
            getattr(self, gtep.power) - mf * Cp * (TT_i - self.outlet.parameters[gtep.TT]),
            self.outlet.parameters[gtep.TT] - TT_i * (1 - getattr(self, gtep.effeff) * (1 - getattr(self, gtep.pipi) ** ((1 - k) / k))),
            getattr(self, gtep.pipi) - PP_i / self.outlet.parameters[gtep.PP],
        )

    def calculate(self, substance_inlet: Substance, substance_mixing: Substance = substance_mixing, x0: dict = None) -> Substance:
        GTENode.validate_substance(self, substance_inlet)
        self.inlet = deepcopy(substance_inlet)
        GTENode.validate_substance(self, substance_mixing)
        self.mixing = deepcopy(substance_mixing)
        self.outlet.name = self.inlet.name
        self.outlet.functions = self.inlet.functions

        self.outlet.parameters[gtep.mf] = self.inlet.parameters[gtep.mf] + self.mixing.parameters[gtep.mf] - self.mass_flow_leak
        self.outlet.parameters[gtep.eo] = self.inlet.parameters[gtep.eo]  # TODO посчитать через массу!

        if x0 is None:
            x0 = tuple(self._x0.values())
        else:
            assert isinstance(x0, dict), TypeError(f"type x0 must be dict, but has {type(x0)}")
            for k, v in self._x0.items():
                if k not in x0:
                    x0[k] = v
        args = {k: v for k, v in self.variables.items() if not isnan(v)}
        count_variables = sum(1 if k not in args else 0 for k in self.variables)
        count_equations = len(self.equations(x0, args)) - 2  # outlet_TT, outlet_PP

        if count_variables < count_equations:
            raise ArithmeticError(f"{count_variables=} < {count_equations=}")
        elif count_variables > count_equations:
            raise ArithmeticError(f"{count_variables=} > {count_equations=}")

        root(self.equations, x0, args, method="lm")

        self.outlet.parameters[gtep.gc] = call_with_kwargs(self.outlet.functions[gtep.gc], self.outlet.parameters)
        self.outlet.parameters[gtep.Cp] = call_with_kwargs(self.outlet.functions[gtep.Cp], self.outlet.parameters)
        self.outlet.parameters[gtep.DD] = self.outlet.parameters[gtep.PP] / (self.outlet.parameters[gtep.gc] * self.outlet.parameters[gtep.TT])
        self.outlet.parameters[gtep.k] = adiabatic_index(self.outlet.parameters[gtep.gc], self.outlet.parameters[gtep.Cp])

        return self.outlet

    def validate(self, epsrel: float = EPSREL) -> bool:
        """Проверка найденного решения"""
        x0 = (self.outlet.parameters[gtep.TT], self.outlet.parameters[gtep.PP])
        args = {k: v for k, v in self.variables.items() if not isnan(v)}

        result = True
        for i, null in enumerate(self.equations(x0, args)):
            if abs(null) > epsrel:
                result = False
                print(f"{i}: {null:.6f}")

        return result

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
        "exhaust",
        parameters={
            gtep.mf: 51,
            gtep.eo: 3,
            gtep.gc: 287,
            gtep.TT: 1600,
            gtep.PP: 101_325 * 23,
            gtep.Cp: 1206,
            gtep.k: 1.33,
            gtep.c: 100,
        },
        functions={
            gtep.gc: lambda excess_oxidizing: gas_const("EXHAUST", excess_oxidizing, fuel="kerosene"),
            gtep.Cp: lambda total_temperature, excess_oxidizing: heat_capacity_at_constant_pressure("EXHAUST", total_temperature, excess_oxidizing, fuel="kerosene"),
        },
    )

    t = Turbine()
    t.summary

    setattr(t, gtep.power, 32 * 10**6)
    setattr(t, gtep.effeff, 0.9)

    t.calculate(substance_inlet)

    t.summary

    print(f"{t.validate(0.0001) = }")
    print(f"{t.is_real = }")
