import os
import pickle
from typing import Any, Callable, Dict, Tuple

from numpy import inf, isnan, nan
from scipy.optimize import root
from substance import Substance
from thermodynamics import parameters as tdp

try:
    from .checks import check_efficiency, check_mass_flow, check_temperature
    from .config import EPSREL
    from .config import parameters as gtep
    from .node import GTENode

except ImportError:
    from checks import check_efficiency, check_mass_flow, check_temperature
    from config import EPSREL
    from config import parameters as gtep
    from node import GTENode


'''
class Nozzle:
    def get_outlet_parameters(self, error=0.01, Niter=100, **kwargs):
        """Расчет параметров после"""

        self.R_gas3 = R_gas(self.substance, a_ox=getattr(self, "a_ox1", None), fuel=fuel)
        self.TT3 = self.TT1

        assert hasattr(self, "ππ") or hasattr(self, "PP3"), f"{type(self).__name__} object has no attributes ππ and PP3"
        if hasattr(self, "PP3"):
            self.ππ = self.PP1 / self.PP3

        if self.ππ < 1:
            self.warnings[3].add("π* < 1!")
            return

        self.PP3 = self.PP1 / self.ππ
        self.ρρ3 = self.PP3 / (self.R_gas3 * self.TT3)
        assert hasattr(self, "g_leak"), f"{type(self).__name__} object has no attribute g_leak!"
        self.g3 = self.g1 - self.g_leak
        self.Cp3 = Cp(self.substance, T=self.TT3, a_ox=getattr(self, "a_ox3", None), fuel=fuel)
        self.k3 = self.Cp3 / (self.Cp3 - self.R_gas3)

        R_gas2 = self.R_gas1
        Cp2 = 0.5 * (self.Cp1 + self.Cp3)  # TODO решить через интеграл
        k2 = Cp2 / (Cp2 - R_gas2)

        self.ππ_max = ((k2 + 1) / 2) ** (k2 / (k2 - 1))
        if self.ππ < self.ππ_max:
            self.c3 = self.v_ * sqrt(2 * Cp2 * self.TT3 * (1 - self.ππ ** ((1 - k2) / k2)))
            self.T3 = self.TT3 - (self.k3 - 1) / (self.k3 * self.R_gas3) * self.c3**2 / 2
            self.R_ = self.c3 * self.g3 - scheme[c][0].c3 * scheme[c][0].g3
        else:
            self.c3 = self.v_ * sqrt(2 * Cp2 * self.TT3)
            self.T3 = self.TT3 - (self.k3 - 1) / (self.k3 * self.R_gas3) * self.c3**2 / 2
            self.R_ = self.c3 * self.g3 - scheme[c][0].c3 * scheme[c][0].g3
            self.R_ += self.R_gas3 * self.T3 / self.c3 * (1 - self.PP3 * self.ππ_max / self.PP1) * self.g3
'''


models = {}
for model in (gtep.TT, gtep.PP, gtep.pipi, gtep.effeff, gtep.power):
    path = f"gte/models/compressor_{model}.pkl"
    if os.path.exists(path):
        with open(path, "rb") as file:
            models[model] = pickle.load(file)
    else:
        print(f"'{path}' not found!")


class Nozzle(GTENode):
    """Выходное устройство"""

    variables = (gtep.p_eff, gtep.s_eff)
    models: Dict[str, Any] = models
    figure: Tuple[Tuple[float, ...], Tuple[float, ...]] = (
        (+0.4, -0.4, -0.4, +0.4),
        (+0.4, +0.4, -0.4, -0.4),
    )

    def __init__(self, name="Outlet", characteristic: Dict[str, Callable] = None):
        GTENode.__init__(self, name, characteristic)

        assert isinstance(characteristic, dict), TypeError(f"{type(characteristic)=} must be dict")
        assert gtep.p_eff in characteristic, KeyError(f"{gtep.p_eff} not in {characteristic=}")
        assert gtep.s_eff in characteristic, KeyError(f"{gtep.s_eff} not in {characteristic=}")
        p_eff, s_eff = characteristic[gtep.p_eff], characteristic[gtep.s_eff]
        # TODO partial
        self.characteristic: Dict[str, Callable] = {gtep.p_eff: p_eff, gtep.s_eff: s_eff}

    def solve(self, inlet: Substance) -> Dict[str, Any]:
        p_eff = self.characteristic[gtep.p_eff](inlet.parameters[gtep.mf])
        s_eff = self.characteristic[gtep.s_eff](inlet.parameters[gtep.mf])

        outlet = self.calculate(inlet, parameters={gtep.p_eff: p_eff, gtep.s_eff: s_eff})

        return {"outlet": outlet}

    @classmethod
    def predict(cls, inlet: Substance, parameters: Dict[str, float | int], use_ml: bool = True) -> Dict[str, float]:
        """Начальные приближения"""
        GTENode.validate_substance(inlet)

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
        p_eff, s_eff = args.get(gtep.p_eff), args.get(gtep.s_eff)
        outlet_TT, outlet_PP = x[0], x[1]

        inlet, outlet = args["inlet"], args["outlet"]

        return (outlet_PP - inlet.parameters[gtep.PP] * p_eff,)

    @classmethod
    def calculate(cls, inlet: Substance, parameters: Dict[str, float | int]) -> Substance:
        prediction: Dict[str, float] = cls.predict(inlet, parameters, use_ml=False)

        outlet = Substance(
            inlet.name,
            inlet.composition,
            parameters={
                gtep.mf: inlet.parameters[gtep.mf],
                gtep.eo: inlet.parameters.get(gtep.eo),
            },
            functions=inlet.functions,
        )

        args: Dict[str, Any] = {"inlet": inlet, "outlet": outlet, **parameters}  # НУ
        x0 = tuple(prediction.values())

        result = root(cls._equations, x0, args, method="lm")

        outlet.parameters[gtep.TT], outlet.parameters[gtep.PP] = float(result.x[0]), float(result.x[1])

        outlet = GTENode.calculate_substance(outlet)

        return outlet

    @classmethod
    def validate(cls, inlet: Substance, outlet: Substance, epsrel: float = EPSREL) -> Dict[int, float]:
        eff_burn = cls.efficiency_burn(inlet, outlet)
        p_eff = cls.pressure_efficiency(inlet, outlet)

        x0 = (outlet.parameters[gtep.TT], outlet.parameters[gtep.PP])
        args = {"inlet": inlet, "outlet": outlet, gtep.eff_burn: eff_burn, gtep.p_eff: p_eff}

        result: Dict[int, float] = {}
        for i, null in enumerate(cls._equations(x0, args)):
            if isnan(null) or abs(null) > epsrel:
                result[i] = null

        return result

    @classmethod
    def check_real(cls, inlet: Substance, outlet: Substance) -> str:
        if not check_mass_flow(inlet.parameters[gtep.mf]):
            return f"inlet {gtep.mf} {inlet.parameters[gtep.mf]}"
        if not check_mass_flow(outlet.parameters[gtep.mf]):
            return f"outlet {gtep.mf} {outlet.parameters[gtep.mf]}"

        if not check_temperature(inlet.parameters[gtep.TT]):
            return f"inlet {gtep.TT} {inlet.parameters[gtep.TT]}"
        if not check_temperature(outlet.parameters[gtep.TT]):
            return f"outlet {gtep.TT} {outlet.parameters[gtep.TT]}"

        if not (0 <= outlet.parameters[gtep.eo]):
            return f"{outlet.parameters[gtep.eo]}"

        eff_burn = cls.efficiency_burn(inlet, outlet)
        if not check_efficiency(eff_burn):
            return f"{eff_burn}"

        return ""


if __name__ == "__main__":
    from fixtures import exhaust

    exhaust.parameters[gtep.TT] = 1200
    exhaust.parameters[gtep.PP] = 101325 * 2

    outlet = Nozzle.calculate(exhaust, {gtep.p_eff: 0.98, gtep.s_eff: 0.99})

    for k, v in outlet.parameters.items():
        print(f"{k:<40}: {v}")

    print(f"{Nozzle.validate(exhaust,  outlet) = }")
    print(f"{Nozzle.check_real(exhaust,  outlet) = }")
    print()
