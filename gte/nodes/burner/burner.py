from typing import Any, Dict, Tuple, Union

from numpy import isnan, nan
from scipy.optimize import root
from substance import Substance
from thermodynamics import T0, heat_capacity_p, heat_capacity_p_exhaust, heat_capacity_p_exhaust_eo1

try:
    from ...checks import check_efficiency, check_mass_flow, check_temperature
    from ...config import EPSREL
    from ...config import parameters as gtep
    from ...errors import TYPE_ERROR
    from ...utils.utils import Function, integrate
    from ..node import Node
except ImportError:
    import os
    import sys

    sys.path.insert(0, os.getcwd())

    from gte.checks import check_efficiency, check_mass_flow, check_temperature
    from gte.config import EPSREL
    from gte.config import parameters as gtep
    from gte.errors import TYPE_ERROR
    from gte.nodes.node import Node
    from gte.utils.utils import Function, integrate


class Burner(Node):
    """Камера сгорания"""

    variables: Tuple[str, str] = (gtep.efficiency, gtep.pipi)
    n_vars: int = 2

    __slots__ = ()  # нет новых атрибутов

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

    def __init__(self, parameters: Dict[str, float], name="Burner"):
        Node.__init__(self, parameters, name)

    @classmethod
    def _equations(cls, x: Tuple[float], args: Dict[str, Any]) -> Tuple[float, float]:
        """
        (m_i * enthalpy_i) + m_f * (Q*efficiency + enthalpy_f) = (m_i + m_f) * enthalpy_o
        pipi = P*_outlet / P*_inlet
        """
        efficiency, pipi = args.get(gtep.efficiency), args.get(gtep.pipi)
        outlet_TT, outlet_PP = x[0], x[1]

        inlet, fuel, outlet = args["inlet"], args["fuel"], args["outlet"]

        return (
            (
                outlet.parameters[gtep.m]
                * integrate(
                    outlet.functions[gtep.hcp],
                    **{
                        gtep.TT: (T0 + 15, outlet_TT),
                        gtep.PP: (101_325, outlet_PP),
                        gtep.eo: (outlet.parameters[gtep.eo], outlet.parameters[gtep.eo]),
                    },
                )[0]
            )
            - (inlet.parameters[gtep.m] * inlet.parameters.get("enthalpy", nan))
            - fuel.parameters[gtep.m] * (fuel.parameters.get("enthalpy", nan) + efficiency * fuel.parameters.get("lower_heat", nan)),
            outlet_PP - inlet.parameters[gtep.PP] * pipi,
        )

    @classmethod
    def predict(cls, parameters: Dict[str, Union[float, int]], inlet: Substance, fuel: Substance) -> Tuple[Dict[str, float], Substance]:
        """Начальные приближения"""
        Node.validate_substance(inlet)
        Node.validate_substance(fuel)
        cls.__validate_fuel(fuel)

        if not isinstance(parameters, dict):
            raise TypeError(TYPE_ERROR.format(f"{type(parameters)=}", dict))
        if len(parameters) != 2:
            raise ValueError(f"{len(parameters)=} must be {cls.n_vars}")
        for parameter, value in parameters.items():
            assert parameter in cls.variables
            assert isinstance(value, (float, int)), TypeError(f"{type(value)=} must be numeric")

        hcp_i: float = inlet.functions[gtep.hcp](inlet.parameters)
        hc_f: float = fuel.functions[gtep.hc](fuel.parameters)

        T15 = T0 + 15  # начальная температура измерения теплоемкости энтальпии

        outlet = Substance(
            "outlet",
            parameters={
                gtep.TT: T15
                + (inlet.parameters[gtep.m] * hcp_i * (inlet.parameters[gtep.TT] - T15) + fuel.parameters[gtep.m] * (hc_f * (fuel.parameters[gtep.TT] - T15) + fuel.parameters["lower_heat"] * parameters[gtep.efficiency]))
                / (inlet.parameters[gtep.m] + fuel.parameters[gtep.m])
                / hcp_i,  # TODO: hcp_exhaust
                gtep.PP: inlet.parameters[gtep.PP] * parameters[gtep.pipi],
            },
        )

        return parameters, outlet

    @classmethod
    def calculate(cls, parameters: Dict[str, float | int], inlet: Substance, fuel: Substance) -> Tuple[Dict[str, float], Substance]:
        _, outlet_ = cls.predict(parameters, inlet, fuel)

        outlet = Substance(
            "exhaust",
            parameters={
                gtep.m: inlet.parameters[gtep.m] + fuel.parameters[gtep.m],
                gtep.PP: inlet.parameters[gtep.PP] * parameters[gtep.pipi],
            },
            functions={
                gtep.gc: fuel.functions[gtep.gc],
            },
        )

        H2O = fuel.composition.get("H2O", 0)  # массовая доля волы в смеси
        outlet.functions[gtep.hcp] = Function(
            lambda total_temperature, excess_oxidizing: heat_capacity_p_exhaust(
                heat_capacity_p_exhaust_eo1(total_temperature, fuel.composition),
                inlet.functions[gtep.hcp]({gtep.eo: excess_oxidizing, gtep.TT: total_temperature, gtep.PP: outlet.parameters[gtep.PP]}),
                heat_capacity_p("H2O", total_temperature),
                excess_oxidizing,
                fuel.parameters["stoichiometry"],
                H2O,
            ),
            name=gtep.hcp,
            args=(gtep.TT, gtep.eo),
        )

        if gtep.eo in inlet.parameters:  # exhaust
            outlet.parameters["oxidizer"] = inlet.parameters["oxidizer"]
            outlet.parameters[gtep.eo] = outlet.parameters["oxidizer"] / (inlet.parameters["oxidizer"] / inlet.parameters[gtep.eo] + fuel.parameters[gtep.m] * fuel.parameters["stoichiometry"])
        else:  # oxidizer
            outlet.parameters["oxidizer"] = inlet.parameters[gtep.m]
            outlet.parameters[gtep.eo] = outlet.parameters["oxidizer"] / (fuel.parameters[gtep.m] * fuel.parameters["stoichiometry"])

        inlet.parameters["enthalpy"], _ = integrate(inlet.functions[gtep.hcp], **{gtep.TT: (T0 + 15, inlet.parameters[gtep.TT]), gtep.PP: (101325, inlet.parameters[gtep.PP]), gtep.eo: (1, inlet.parameters.get(gtep.eo, 1))})
        fuel.parameters["enthalpy"], _ = integrate(fuel.functions[gtep.hc], **{gtep.TT: (T0 + 15, fuel.parameters[gtep.TT])})

        args: Dict[str, Any] = {"inlet": inlet, "fuel": fuel, "outlet": outlet, **parameters}  # НУ
        x0 = [outlet_.parameters[gtep.TT], outlet_.parameters[gtep.PP]]

        result = root(cls._equations, x0, args, method="lm")

        outlet.parameters[gtep.TT], outlet.parameters[gtep.PP] = float(result.x[0]), float(result.x[1])

        outlet = Node.calculate_substance(outlet)

        efficiency = parameters.get(gtep.efficiency, cls.efficiency(inlet, fuel, outlet))
        pipi = parameters.get(gtep.pipi, cls.total_pressure_ratio(inlet, fuel, outlet))

        return {gtep.efficiency: efficiency, gtep.pipi: pipi}, outlet

    @classmethod
    def validate(cls, inlet: Substance, fuel: Substance, outlet: Substance, epsrel: float = EPSREL) -> Dict[int, float]:
        efficiency = cls.efficiency(inlet, fuel, outlet)
        pipi = cls.total_pressure_ratio(inlet, fuel, outlet)

        x0 = (outlet.parameters[gtep.TT], outlet.parameters[gtep.PP])
        args = {"inlet": inlet, "fuel": fuel, "outlet": outlet, gtep.efficiency: efficiency, gtep.pipi: pipi}

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

        efficiency = cls.efficiency(inlet, fuel, outlet)
        if not check_efficiency(efficiency):
            return f"{efficiency}"

        return ""

    @classmethod
    def efficiency(cls, inlet: Substance, fuel: Substance, outlet: Substance) -> float:
        """КПД полноты сгорания топлива"""
        Node.validate_substance(inlet)
        Node.validate_substance(fuel)
        cls.__validate_fuel(fuel)

        inlet_enthalpy, _ = integrate(inlet.functions[gtep.hcp], **{gtep.eo: (inlet.parameters.get(gtep.eo, 1), inlet.parameters.get(gtep.eo, 1)), gtep.TT: (T0 + 15, inlet.parameters[gtep.TT]), gtep.PP: (101325, inlet.parameters[gtep.PP])})
        fuel_enthalpy, _ = integrate(fuel.functions[gtep.hc], **{gtep.TT: (T0 + 15, fuel.parameters[gtep.TT])})
        outlet_enthalpy, _ = integrate(
            outlet.functions[gtep.hcp],
            **{
                gtep.eo: (outlet.parameters[gtep.eo], outlet.parameters[gtep.eo]),
                gtep.TT: (T0 + 15, outlet.parameters[gtep.TT]),
                gtep.PP: (101_325, outlet.parameters[gtep.PP]),
            },
        )

        return (outlet.parameters[gtep.m] * outlet_enthalpy - inlet.parameters[gtep.m] * inlet_enthalpy - fuel.parameters[gtep.m] * fuel_enthalpy) / (fuel.parameters[gtep.m] * fuel.parameters["lower_heat"])

    @classmethod
    def total_pressure_ratio(cls, inlet: Substance, fuel: Substance, outlet: Substance) -> float:
        """Коэф. сохранения давления"""
        Node.validate_substance(inlet)
        Node.validate_substance(fuel)
        cls.__validate_fuel(fuel)

        return outlet.parameters[gtep.PP] / inlet.parameters[gtep.PP]


if __name__ == "__main__":
    from gte.fixtures import air as inlet
    from gte.fixtures import kerosene as fuel

    inlet.parameters[gtep.TT] = 600
    inlet.parameters[gtep.PP] = 101325 * 6

    cc = Burner({gtep.efficiency: 0.99, gtep.pipi: 0.95}, name="test")
    print(f"{cc.is_solvable=}")

    vars, outlet = cc.calculate(cc.parameters, inlet, fuel)

    for k, v in outlet.parameters.items():
        print(f"{k:<25}: {v:.4f}")
    print(vars)

    print(f"{Burner.validate(inlet, fuel, outlet) = }")
    print(f"{Burner.check_real(inlet, fuel, outlet) = }")
    print()
