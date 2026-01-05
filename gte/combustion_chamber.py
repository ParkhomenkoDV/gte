import os
import pickle
from typing import Any, Callable, Dict, Tuple

from numpy import cos, isnan, linspace, nan, radians, sin
from scipy.optimize import root
from substance import Substance
from thermodynamics import T0, heat_capacity_p, heat_capacity_p_exhaust
from thermodynamics import parameters as tdp

try:
    from .checks import check_efficiency, check_mass_flow, check_temperature
    from .config import EPSREL
    from .config import parameters as gtep
    from .node import GTENode
    from .utils import integrate
except ImportError:
    from checks import check_efficiency, check_mass_flow, check_temperature
    from config import EPSREL
    from config import parameters as gtep
    from node import GTENode
    from utils import integrate


models = {}
for model in (gtep.TT, gtep.PP, gtep.pipi, gtep.effeff, gtep.power):
    path = f"gte/models/compressor_{model}.pkl"
    if os.path.exists(path):
        with open(path, "rb") as file:
            models[model] = pickle.load(file)
    else:
        print(f"'{path}' not found!")


class CombustionChamber(GTENode):
    """Камера сгорания"""

    variables = (gtep.eff_burn, gtep.p_eff)
    models: Dict[str, Any] = models
    figure: Tuple[Tuple[float, ...], Tuple[float, ...]] = (
        tuple(0.4 * cos(alpha) for alpha in linspace(0, radians(360), 360)),
        tuple(0.4 * sin(alpha) for alpha in linspace(0, radians(360), 360)),
    )

    @classmethod
    def __validate_fuel(cls, fuel: Substance) -> None:
        """Проверка параметров горючего"""
        assert len(fuel.composition) > 0, ValueError(f"{fuel.composition = }")

        stoichiometry = fuel.parameters.get("stoichiometry")
        assert stoichiometry is not None, KeyError("fuel has not parameter 'stoichiometry'")
        assert isinstance(stoichiometry, (int, float)), TypeError(f"{type(stoichiometry)=} must be numeric")

        lower_heat = fuel.parameters.get("lower_heat")
        assert lower_heat is not None, KeyError("fuel has not parameter 'lower_heat'")
        assert isinstance(lower_heat, (int, float)), TypeError(f"{type(lower_heat)=} must be numeric")

        assert fuel.functions.get(gtep.gc) is not None, KeyError(f"fuel has not function '{gtep.gc}'")

    def __init__(self, name="CombustionChamber", characteristic: Dict[str, Callable] = None):
        GTENode.__init__(self, name, characteristic)

        assert isinstance(characteristic, dict), TypeError(f"{type(characteristic)=} must be dict")
        assert gtep.eff_burn in characteristic, KeyError(f"{gtep.eff_burn} not in {characteristic=}")
        assert gtep.p_eff in characteristic, KeyError(f"{gtep.p_eff} not in {characteristic=}")
        eff_burn, p_eff = characteristic[gtep.eff_burn], characteristic[gtep.p_eff]
        # TODO partial
        self.characteristic: Dict[str, Callable] = {gtep.eff_burn: eff_burn, gtep.p_eff: p_eff}

    def solve(self, inlet: Substance, fuel: Substance) -> Dict[str, Any]:
        eff_burn = self.characteristic[gtep.eff_burn](inlet.parameters[gtep.mf])
        p_eff = self.characteristic[gtep.p_eff](inlet.parameters[gtep.mf])

        outlet = self.calculate(inlet, fuel, parameters={gtep.eff_burn: eff_burn, gtep.p_eff: p_eff})

        return {"outlet": outlet}

    @classmethod
    def predict(cls, inlet: Substance, fuel: Substance, parameters: Dict[str, float | int], use_ml: bool = True) -> Dict[str, float]:
        """Начальные приближения"""
        GTENode.validate_substance(inlet)
        GTENode.validate_substance(fuel)
        cls.__validate_fuel(fuel)

        assert isinstance(parameters, dict), TypeError(f"{type(parameters)=} must be dict")
        assert len(parameters) == 2, ArithmeticError(f"{len(parameters)=} must be 2")
        for key, value in parameters.items():
            assert key in (gtep.eff_burn, gtep.p_eff)
            assert isinstance(value, (float, int)), TypeError(f"{type(value)=} must be numeric")

        assert isinstance(use_ml, bool), TypeError(f"{type(use_ml)=} must be bool")

        prediction = {
            f"outlet_{gtep.TT}": inlet.parameters[gtep.TT],  # TODO: model or formula
            f"outlet_{gtep.PP}": inlet.parameters[gtep.PP] * parameters[gtep.p_eff],
        }
        return prediction

    @classmethod
    def _equations(cls, x: Tuple[float], args: Dict[str, Any]) -> Tuple[float, float]:
        """
        (mf_i * enthalpy_i) + mf_f * (Q*eff_burn + enthalpy_f) = (mf_i + mf_f) * enthalpy_o
        p_eff = P*_outlet / P*_inlet
        """
        eff_burn, p_eff = args.get(gtep.eff_burn), args.get(gtep.p_eff)
        outlet_TT, outlet_PP = x[0], x[1]

        inlet, fuel, outlet = args["inlet"], args["fuel"], args["outlet"]

        return (
            (
                outlet.parameters[gtep.mf]
                * integrate(
                    outlet.functions[gtep.hcp],
                    **{
                        tdp.t: (T0 + 15, outlet_TT),
                        tdp.p: (101_325, outlet_PP),
                        tdp.eo: (outlet.parameters[gtep.eo], outlet.parameters[gtep.eo]),
                    },
                )[0]
            )
            - (inlet.parameters[gtep.mf] * inlet.parameters.get("enthalpy", nan))
            - fuel.parameters[gtep.mf] * (fuel.parameters.get("enthalpy", nan) + eff_burn * fuel.parameters.get("lower_heat", nan)),
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
            (1 - H2O) * (heat_capacity_p_exhaust(temperature, fuel.composition) + inlet.functions[gtep.hcp](temperature) * excess_oxidizing * stoichiometry) + H2O * excess_oxidizing * stoichiometry * heat_capacity_p("H2O", temperature)
        ) / (1 - H2O + excess_oxidizing * stoichiometry)

        outlet.parameters[gtep.mf] = inlet.parameters[gtep.mf] + fuel.parameters[gtep.mf]
        outlet.parameters[gtep.eo] = inlet.parameters[gtep.mf] / fuel.parameters[gtep.mf] / fuel.parameters["stoichiometry"]

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
        if not check_mass_flow(inlet.parameters[gtep.mf]):
            return f"inlet {gtep.mf} {inlet.parameters[gtep.mf]}"
        if not check_mass_flow(fuel.parameters[gtep.mf]):
            return f"fuel {gtep.mf} {fuel.parameters[gtep.mf]}"
        if not check_mass_flow(outlet.parameters[gtep.mf]):
            return f"outlet {gtep.mf} {outlet.parameters[gtep.mf]}"

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

        return (outlet.parameters[gtep.mf] * outlet_enthalpy - inlet.parameters[gtep.mf] * inlet_enthalpy - fuel.parameters[gtep.mf] * fuel_enthalpy) / (fuel.parameters[gtep.mf] * fuel.parameters["lower_heat"])

    @classmethod
    def pressure_efficiency(cls, inlet: Substance, fuel: Substance, outlet: Substance) -> float:
        """Коэф. сохранения давления"""
        GTENode.validate_substance(inlet)
        GTENode.validate_substance(fuel)
        cls.__validate_fuel(fuel)

        return outlet.parameters[gtep.PP] / inlet.parameters[gtep.PP]


if __name__ == "__main__":
    from fixtures import air, kerosene

    air.parameters[gtep.TT] = 600
    air.parameters[gtep.PP] = 101325 * 6

    outlet = CombustionChamber.calculate(air, kerosene, {gtep.eff_burn: 0.99, gtep.p_eff: 0.95})

    for k, v in outlet.parameters.items():
        print(f"{k:<40}: {v}")

    print(f"{CombustionChamber.validate(air, kerosene, outlet) = }")
    print(f"{CombustionChamber.check_real(air, kerosene, outlet) = }")
    print()
