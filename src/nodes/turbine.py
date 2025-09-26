from copy import deepcopy

from mathematics import eps
from numpy import isinf, isnan, nan
from scipy.optimize import fsolve
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
            f"outlet_{gtep.TT}": self.inlet.parameters[gtep.TT],
            f"outlet_{gtep.PP}": self.inlet.parameters[gtep.PP],
        }
        for k, v in self.variables.items():
            if not isnan(v):
                continue
            if k == gtep.pipi:
                x0[k] = 4  # TODO: model or formula
            elif k == gtep.effeff:
                x0[k] = 1.0
            elif k == gtep.power:
                x0[k] = 20 * 10**6  # TODO: model or formula
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

        TT_i = self.inlet.parameters[gtep.TT]
        PP_i = self.inlet.parameters[gtep.PP]
        eo_i = self.inlet.parameters[gtep.eo]

        mf = (self.inlet.parameters[gtep.mf] + self.outlet.parameters[gtep.mf]) / 2
        gc, _ = integral_average(
            self.inlet.functions[gtep.gc],
            **{
                gtep.TT: (TT_i, self.outlet.parameters[gtep.TT]),
                gtep.PP: (PP_i, self.outlet.parameters[gtep.PP]),
                gtep.eo: (eo_i, self.outlet.parameters[gtep.eo]),
            },
        )
        Cp, _ = integral_average(
            self.inlet.functions[gtep.Cp],
            **{
                gtep.TT: (TT_i, self.outlet.parameters[gtep.TT]),
                gtep.PP: (PP_i, self.outlet.parameters[gtep.PP]),
                gtep.eo: (eo_i, self.outlet.parameters[gtep.eo]),
            },
        )
        k = adiabatic_index(gc, Cp)

        return (
            getattr(self, gtep.power) - mf * Cp * (self.outlet.parameters[gtep.TT] - TT_i),
            self.outlet.parameters[gtep.TT] - TT_i * (1 - (1 - getattr(self, gtep.pipi) ** ((1 - k) / k) - 1) * getattr(self, gtep.effeff)),
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
        self.outlet.parameters[gtep.eo] = 1

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

        fsolve(self.equations, x0, args)

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
            epsilon = eps("rel", null, 0)
            if epsilon > epsrel and not isinf(epsilon):
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
            gtep.eo: 1.4,
            gtep.gc: 287,
            gtep.TT: 1700,
            gtep.PP: 101_325 * 5,
            gtep.Cp: 1206,
            gtep.k: 1.33,
            gtep.c: 100,
        },
        functions={
            gtep.gc: lambda excess_oxidizing: gas_const("EXHAUST", excess_oxidizing, fuel="kerosene"),
            gtep.Cp: lambda total_temperature, excess_oxidizing: heat_capacity_at_constant_pressure("AIR", total_temperature, excess_oxidizing, fuel="kerosene"),
        },
    )

    t = Turbine()
    t.summary

    setattr(t, gtep.power, 25 * 10**6)
    setattr(t, gtep.effeff, 0.9)

    t.summary

    t.calculate(substance_inlet)

    for k, v in t.__dict__.items():
        print(k, "=", v)
