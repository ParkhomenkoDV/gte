from copy import deepcopy

from mathematics import eps
from node import GTENode
from numpy import nan
from scipy import integrate
from substance import Substance
from thermodynamics import adiabatic_index, efficiency_polytropic

from src.errors import ITERATION_LIMIT
from src.parameters import parameters as gtep


class Compressor(GTENode):
    """Компрессор"""

    def __init__(self, name="Compressor"):
        GTENode.__init__(self, name=name)

        self.pipi = nan
        self.eff = nan
        self.power = nan

    def calculate(
        self,
        substance_inlet: Substance,
        epsrel: float = 0.01,
        niter: int = 10,
        **kwargs,
    ) -> Substance:
        GTENode.validate_substance(self, substance_inlet)
        self.inlet = deepcopy(substance_inlet)
        self.outlet = deepcopy(self.inlet)

        self.outlet.parameters[gtep.PP] = self.inlet.parameters[gtep.PP] / self.pipi

        for _ in range(niter):
            k = 0.5 * (self.inlet.parameters[gtep.k] + self.outlet.parameters[gtep.k])
            self.outlet.parameters[gtep.TT] = self.inlet.parameters[gtep.TT] * (
                1 + (self.pipi ** ((k - 1) / k) - 1) / self.eff
            )

            self.outlet.parameters[gtep.Cp] = self.outlet.functions[gtep.Cp](
                T=self.outlet.parameters[gtep.TT], P=self.outlet.parameters[gtep.PP]
            )
            k_outlet = adiabatic_index(
                self.outlet.parameters[gtep.gc],
                self.outlet.parameters[gtep.Cp],
            )
            if abs(eps("rel", self.outlet.parameters[gtep.k], k_outlet)) <= epsrel:
                self.outlet.parameters[gtep.k] = k_outlet
                break
            self.outlet.parameters[gtep.k] = k_outlet
        else:
            raise AssertionError(ITERATION_LIMIT.format(self.name))

        self.outlet.parameters[gtep.DD] = self.outlet.parameters[gtep.PP] / (
            self.outlet.parameters[gtep.gc] * self.outlet.parameters[gtep.TT]
        )
        self.outlet.parameters[gtep.mf] = (
            1 - self.mass_flow_leak
        ) * self.inlet.parameters[gtep.mf]

        self.ηn = efficiency_polytropic("C", pipi=self.pipi, effeff=self.eff, k=k)
        self.power = (
            integrate.quad(
                self.outlet.functions[gtep.Cp],
                self.inlet.parameters[gtep.TT],
                self.outlet.parameters[gtep.TT],
            )[0]
            * (self.inlet.parameters[gtep.mf] + self.outlet.parameters[gtep.mf])
            / 2
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
        },
        functions={
            gtep.Cp: lambda T: 1006,
        },
    )

    compressor = Compressor()
    compressor.pipi = 6
    compressor.eff = 0.86
    compressor.mass_flow_leak = 0.03

    compressor.calculate(substance_inlet)

    for k, v in compressor.summary.items():
        print(f"{k:<40}: {v}")
