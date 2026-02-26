import os
from typing import Any, Callable, Dict, List, Tuple, Union

import matplotlib.pyplot as plt
from numpy import arange, cos, isnan, linspace, nan, radians, sin
from numpy.typing import ArrayLike
from scipy.optimize import root
from substance import Substance
from thermodynamics import T0, heat_capacity_p, heat_capacity_p_exhaust
from thermodynamics import parameters as tdp

try:
    from ...checks import check_characteristic, check_efficiency, check_mass_flow, check_temperature
    from ...config import EPSREL
    from ...config import parameters as gtep
    from ...errors import TYPE_ERROR
    from ...utils import integrate
    from ..node import GTENode
except ImportError:
    import os
    import sys

    sys.path.insert(0, os.getcwd())

    from gte.checks import check_characteristic, check_efficiency, check_mass_flow, check_temperature
    from gte.config import EPSREL
    from gte.config import parameters as gtep
    from gte.errors import TYPE_ERROR
    from gte.nodes.node import GTENode
    from gte.utils import integrate


class CombustionChamber(GTENode):
    """Камера сгорания"""

    variables: Tuple[str, str] = (gtep.eff_burn, gtep.p_eff)
    models: Dict[str, Any] = {}
    figure: Tuple[Tuple[float, ...], Tuple[float, ...]] = (
        tuple(0.4 * cos(alpha) for alpha in linspace(0, radians(360), 360, endpoint=True)),
        tuple(0.4 * sin(alpha) for alpha in linspace(0, radians(360), 360, endpoint=True)),
    )

    @classmethod
    def __validate_fuel(cls, fuel: Substance) -> None:
        """Проверка параметров горючего"""
        if len(fuel.composition) == 0:
            raise ValueError(f"{fuel.composition=} must be not empty")

        stoichiometry = fuel.parameters.get("stoichiometry")
        if stoichiometry is None:
            raise KeyError("fuel has not parameter 'stoichiometry'")
        if not isinstance(stoichiometry, (int, float)):
            raise TypeError(TYPE_ERROR.format(f"{type(stoichiometry)=}", float))

        lower_heat = fuel.parameters.get("lower_heat")
        if lower_heat is None:
            raise KeyError("fuel has not parameter 'lower_heat'")
        if not isinstance(lower_heat, (int, float)):
            raise TypeError(f"{type(lower_heat)=}", float)

        if fuel.functions.get(gtep.gc) is None:
            raise KeyError(f"fuel has not function '{gtep.gc}'")

    def __init__(self, efficiency_burn: Callable, pressure_efficiency: Callable, name="CombustionChamber"):
        """Инициализация объекта камеры сгорания"""

        for function in (efficiency_burn, pressure_efficiency):
            check_characteristic(function, {gtep.m})

        GTENode.__init__(self, {gtep.eff_burn: efficiency_burn, gtep.p_eff: pressure_efficiency}, name)

    def plot_characteristic(
        self,
        mass_flow: Union[Tuple[float], List[float], ArrayLike],
        figsize: Tuple[int, int] = (8, 10),
    ) -> plt.Figure:
        fg = plt.figure(figsize=figsize)
        gs = fg.add_gridspec(2, 1)  # строки, столбцы

        for i, (name, func) in enumerate(self.characteristic.items()):
            ax = fg.add_subplot(gs[i, 0])
            ax.axis("equal")
            ax.set_xlabel(gtep.m, fontsize=12)
            ax.set_ylabel(name, fontsize=12)
            ax.grid()

            y = [func(**{gtep.m: m}) for m in mass_flow]
            ax.plot(mass_flow, y)

        fg.tight_layout()
        return fg

    def solve(self, inlet: Substance, fuel: Substance) -> Dict[str, Any]:
        eff_burn = self.characteristic[gtep.eff_burn](inlet.parameters[gtep.m])
        p_eff = self.characteristic[gtep.p_eff](inlet.parameters[gtep.m])

        outlet = self.calculate(inlet, fuel, parameters={gtep.eff_burn: eff_burn, gtep.p_eff: p_eff})

        return {"outlet": outlet}

    @classmethod
    def predict(cls, inlet: Substance, fuel: Substance, parameters: Dict[str, float | int], use_ml: bool = True) -> Dict[str, float]:
        """Начальные приближения"""
        GTENode.validate_substance(inlet)
        GTENode.validate_substance(fuel)
        cls.__validate_fuel(fuel)

        if not isinstance(parameters, dict):
            raise TypeError(TYPE_ERROR.format(f"{type(parameters)=}", dict))
        assert len(parameters) == 2, f"{len(parameters)=} must be 2"
        for key, value in parameters.items():
            assert key in (gtep.eff_burn, gtep.p_eff)
            assert isinstance(value, (float, int)), TypeError(f"{type(value)=} must be numeric")

        if not isinstance(use_ml, bool):
            raise TypeError(TYPE_ERROR.format(f"{type(use_ml)=}", bool))

        prediction = {
            f"outlet_{gtep.TT}": inlet.parameters[gtep.TT],  # TODO: model or formula
            f"outlet_{gtep.PP}": inlet.parameters[gtep.PP] * parameters[gtep.p_eff],
        }
        return prediction

    @classmethod
    def _equations(cls, x: Tuple[float], args: Dict[str, Any]) -> Tuple[float, float]:
        """
        (m_i * enthalpy_i) + m_f * (Q*eff_burn + enthalpy_f) = (m_i + m_f) * enthalpy_o
        p_eff = P*_outlet / P*_inlet
        """
        eff_burn, p_eff = args.get(gtep.eff_burn), args.get(gtep.p_eff)
        outlet_TT, outlet_PP = x[0], x[1]

        inlet, fuel, outlet = args["inlet"], args["fuel"], args["outlet"]

        return (
            (
                outlet.parameters[gtep.m]
                * integrate(
                    outlet.functions[gtep.hcp],
                    **{
                        tdp.t: (T0 + 15, outlet_TT),
                        tdp.p: (101_325, outlet_PP),
                        tdp.eo: (outlet.parameters[gtep.eo], outlet.parameters[gtep.eo]),
                    },
                )[0]
            )
            - (inlet.parameters[gtep.m] * inlet.parameters.get("enthalpy", nan))
            - fuel.parameters[gtep.m] * (fuel.parameters.get("enthalpy", nan) + eff_burn * fuel.parameters.get("lower_heat", nan)),
            outlet_PP - inlet.parameters[gtep.PP] * p_eff,
        )

    @classmethod
    def calculate(cls, inlet: Substance, fuel: Substance, parameters: Dict[str, float | int]) -> Substance:
        prediction: Dict[str, float] = cls.predict(inlet, fuel, parameters, use_ml=False)

        outlet = Substance("exhaust")
        outlet.functions[gtep.gc] = fuel.functions[gtep.gc]

        stoichiometry = fuel.parameters["stoichiometry"]
        H2O = fuel.composition.get("H2O", 0)  # массовая доля волы в смеси
        outlet.functions[gtep.hcp] = lambda temperature, excess_oxidizing: (
            ((1 - H2O) * (heat_capacity_p_exhaust(temperature, fuel.composition) + inlet.functions[gtep.hcp](temperature) * excess_oxidizing * stoichiometry) + H2O * excess_oxidizing * stoichiometry * heat_capacity_p("H2O", temperature))
            / (1 - H2O + excess_oxidizing * stoichiometry)
        )

        outlet.parameters[gtep.m] = inlet.parameters[gtep.m] + fuel.parameters[gtep.m]
        outlet.parameters[gtep.eo] = inlet.parameters[gtep.m] / fuel.parameters[gtep.m] / fuel.parameters["stoichiometry"]

        inlet.parameters["enthalpy"], _ = integrate(inlet.functions[gtep.hcp], **{tdp.t: (T0 + 15, inlet.parameters[gtep.TT]), tdp.p: (101325, inlet.parameters[gtep.PP])})
        fuel.parameters["enthalpy"], _ = integrate(fuel.functions[gtep.hc], **{tdp.t: (T0 + 15, fuel.parameters[gtep.TT])})

        args: Dict[str, Any] = {"inlet": inlet, "fuel": fuel, "outlet": outlet, **parameters}  # НУ
        x0 = tuple(prediction.values())

        result = root(cls._equations, x0, args, method="lm")

        outlet.parameters[gtep.TT], outlet.parameters[gtep.PP] = float(result.x[0]), float(result.x[1])

        outlet = GTENode.calculate_substance(outlet)

        return outlet

    @classmethod
    def validate(cls, inlet: Substance, fuel: Substance, outlet: Substance, epsrel: float = EPSREL) -> Dict[int, float]:
        eff_burn = cls.efficiency_burn(inlet, fuel, outlet)
        p_eff = cls.pressure_efficiency(inlet, fuel, outlet)

        x0 = (outlet.parameters[gtep.TT], outlet.parameters[gtep.PP])
        args = {"inlet": inlet, "fuel": fuel, "outlet": outlet, gtep.eff_burn: eff_burn, gtep.p_eff: p_eff}

        result: Dict[int, float] = {}
        for i, null in enumerate(cls._equations(x0, args)):
            if isnan(null) or abs(null) > epsrel:
                result[i] = null

        return result

    @classmethod
    def check_real(cls, inlet: Substance, fuel: Substance, outlet: Substance) -> str:
        if not check_mass_flow(inlet.parameters[gtep.m]):
            return f"inlet {gtep.m} {inlet.parameters[gtep.m]}"
        if not check_mass_flow(fuel.parameters[gtep.m]):
            return f"fuel {gtep.m} {fuel.parameters[gtep.m]}"
        if not check_mass_flow(outlet.parameters[gtep.m]):
            return f"outlet {gtep.m} {outlet.parameters[gtep.m]}"

        if not check_temperature(inlet.parameters[gtep.TT]):
            return f"inlet {gtep.TT} {inlet.parameters[gtep.TT]}"
        if not check_temperature(fuel.parameters[gtep.TT]):
            return f"fuel {gtep.TT} {fuel.parameters[gtep.TT]}"
        if not check_temperature(outlet.parameters[gtep.TT]):
            return f"outlet {gtep.TT} {outlet.parameters[gtep.TT]}"

        if not (inlet.parameters[gtep.TT] <= outlet.parameters[gtep.TT]):
            return f"{inlet.parameters[gtep.TT] <= outlet.parameters[gtep.TT]}"

        if not (0 <= outlet.parameters[gtep.eo]):
            return f"{outlet.parameters[gtep.eo]}"

        eff_burn = cls.efficiency_burn(inlet, fuel, outlet)
        if not check_efficiency(eff_burn):
            return f"{eff_burn}"

        return ""

    @classmethod
    def efficiency_burn(cls, inlet: Substance, fuel: Substance, outlet: Substance) -> float:
        """КПД полноты сгорания топлива"""
        GTENode.validate_substance(inlet)
        GTENode.validate_substance(fuel)
        cls.__validate_fuel(fuel)

        inlet_enthalpy, _ = integrate(inlet.functions[gtep.hcp], **{tdp.t: (T0 + 15, inlet.parameters[gtep.TT]), tdp.p: (101325, inlet.parameters[gtep.PP])})
        fuel_enthalpy, _ = integrate(fuel.functions[gtep.hc], **{tdp.t: (T0 + 15, fuel.parameters[gtep.TT])})
        outlet_enthalpy, _ = integrate(
            outlet.functions[gtep.hcp],
            **{
                tdp.t: (T0 + 15, outlet.parameters[gtep.TT]),
                tdp.p: (101_325, outlet.parameters[gtep.PP]),
                tdp.eo: (outlet.parameters[gtep.eo], outlet.parameters[gtep.eo]),
            },
        )

        return (outlet.parameters[gtep.m] * outlet_enthalpy - inlet.parameters[gtep.m] * inlet_enthalpy - fuel.parameters[gtep.m] * fuel_enthalpy) / (fuel.parameters[gtep.m] * fuel.parameters["lower_heat"])

    @classmethod
    def pressure_efficiency(cls, inlet: Substance, fuel: Substance, outlet: Substance) -> float:
        """Коэф. сохранения давления"""
        GTENode.validate_substance(inlet)
        GTENode.validate_substance(fuel)
        cls.__validate_fuel(fuel)

        return outlet.parameters[gtep.PP] / inlet.parameters[gtep.PP]


if __name__ == "__main__":
    from gte.fixtures import air as inlet
    from gte.fixtures import kerosene as fuel

    inlet.parameters[gtep.TT] = 600
    inlet.parameters[gtep.PP] = 101325 * 6

    outlet = CombustionChamber.calculate(inlet, fuel, {gtep.eff_burn: 0.99, gtep.p_eff: 0.95})

    for k, v in outlet.parameters.items():
        print(f"{k:<40}: {v}")

    print(f"{CombustionChamber.validate(inlet, fuel, outlet) = }")
    print(f"{CombustionChamber.check_real(inlet, fuel, outlet) = }")
    print()

    cc = CombustionChamber(
        efficiency_burn=lambda mass: 0.99,
        pressure_efficiency=lambda mass: 0.95,
    )

    cc.plot_characteristic(mass_flow=arange(0.5, 1.1, 0.05))
    plt.show()

    result = cc.solve(inlet, fuel)
    print(result)
