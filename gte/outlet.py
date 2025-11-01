from numpy import inf, nan
from substance import Substance
from thermodynamics import adiabatic_index, gas_const, heat_capacity_at_constant_pressure

try:
    from .checks import check_efficiency, check_temperature
    from .config import EPSREL
    from .config import parameters as gtep
    from .node import GTENode
    from .utils import call_with_kwargs, integral_average
except ImportError:
    from checks import check_efficiency, check_temperature
    from config import EPSREL
    from config import parameters as gtep
    from node import GTENode
    from utils import call_with_kwargs, integral_average


class Outlet(GTENode):
    """Выходное устройство"""

    __slots__ = (gtep.peff,)

    def __init__(self, name="Outlet"):
        GTENode.__init__(self, name=name)

        setattr(self, gtep.peff, nan)

    @property
    def variables(self) -> dict[str:float]:
        return {
            gtep.pipi: getattr(self, gtep.pipi),
            gtep.peff: getattr(self, gtep.peff),
        }

    def calculate(self, substance_inlet: Substance, x0=None):
        return self.outlet


if __name__ == "__main__":
    from colorama import Fore

    for k, v in gtep.items():
        print(f"{k:<10}: {v}")

    outlet = Outlet()
    outlet.σ = 0.98

    outlet.TT_inlet = 1500
    outlet.PP_inlet = 600_000
    outlet.a_ox = 3
    print(outlet.solve(substance="EXHAUST", fuel="КЕРОСИН"))
