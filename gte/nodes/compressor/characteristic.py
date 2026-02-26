from typing import Callable, Dict, Tuple

from numpy import nan

try:
    from ...config import parameters as gtep
    from ...utils import Interpolator
except ImportError:
    import os
    import sys

    sys.path.insert(0, os.getcwd())

    from gte.config import parameters as gtep
    from gte.utils import Interpolator


# TODO: add mode
characteristic_13_25: Tuple[Dict[str, float]] = (
    #
    {gtep.rf: 0.800, gtep.m: 0.590, gtep.effeff: 0.90, gtep.pipi: 0.820},
    {gtep.rf: 0.800, gtep.m: 0.640, gtep.effeff: 0.95, gtep.pipi: 0.800},
    {gtep.rf: 0.800, gtep.m: 0.760, gtep.effeff: 0.95, gtep.pipi: 0.760},
    {gtep.rf: 0.800, gtep.m: 0.835, gtep.effeff: 0.90, gtep.pipi: 0.710},
    #
    {gtep.rf: 0.825, gtep.m: 0.610, gtep.effeff: 0.90, gtep.pipi: 0.860},
    {gtep.rf: 0.825, gtep.m: 0.640, gtep.effeff: 0.95, gtep.pipi: 0.850},
    {gtep.rf: 0.825, gtep.m: 0.685, gtep.effeff: 1.00, gtep.pipi: 0.827},
    {gtep.rf: 0.825, gtep.m: 0.748, gtep.effeff: 1.00, gtep.pipi: 0.805},
    {gtep.rf: 0.825, gtep.m: 0.810, gtep.effeff: 0.95, gtep.pipi: 0.775},
    {gtep.rf: 0.825, gtep.m: 0.880, gtep.effeff: 0.90, gtep.pipi: 0.740},
    #
    {gtep.rf: 0.850, gtep.m: 0.630, gtep.effeff: 0.90, gtep.pipi: 0.895},
    {gtep.rf: 0.850, gtep.m: 0.660, gtep.effeff: 0.95, gtep.pipi: 0.885},
    {gtep.rf: 0.850, gtep.m: 0.690, gtep.effeff: 1.00, gtep.pipi: 0.875},
    {gtep.rf: 0.850, gtep.m: 0.800, gtep.effeff: 1.00, gtep.pipi: 0.830},
    {gtep.rf: 0.850, gtep.m: 0.855, gtep.effeff: 0.95, gtep.pipi: 0.800},
    {gtep.rf: 0.850, gtep.m: 0.920, gtep.effeff: 0.90, gtep.pipi: 0.760},
    #
    {gtep.rf: 0.900, gtep.m: 0.700, gtep.effeff: 0.95, gtep.pipi: 0.950},
    {gtep.rf: 0.900, gtep.m: 0.740, gtep.effeff: 1.00, gtep.pipi: 0.940},
    {gtep.rf: 0.900, gtep.m: 0.790, gtep.effeff: 1.05, gtep.pipi: 0.920},
    {gtep.rf: 0.900, gtep.m: 0.825, gtep.effeff: 1.05, gtep.pipi: 0.910},
    {gtep.rf: 0.900, gtep.m: 0.875, gtep.effeff: 1.00, gtep.pipi: 0.880},
    {gtep.rf: 0.900, gtep.m: 0.925, gtep.effeff: 0.95, gtep.pipi: 0.855},
    {gtep.rf: 0.900, gtep.m: 0.980, gtep.effeff: 0.90, gtep.pipi: 0.820},
    #
    {gtep.rf: 0.950, gtep.m: 0.800, gtep.effeff: 1.00, gtep.pipi: 0.990},
    {gtep.rf: 0.950, gtep.m: 0.840, gtep.effeff: 1.05, gtep.pipi: 0.975},
    {gtep.rf: 0.950, gtep.m: 0.900, gtep.effeff: 1.05, gtep.pipi: 0.960},
    {gtep.rf: 0.950, gtep.m: 0.940, gtep.effeff: 1.00, gtep.pipi: 0.940},
    {gtep.rf: 0.950, gtep.m: 0.990, gtep.effeff: 0.95, gtep.pipi: 0.910},
    #
    {gtep.rf: 1.000, gtep.m: 0.925, gtep.effeff: 1.05, gtep.pipi: 1.040},
    {gtep.rf: 1.000, gtep.m: 0.965, gtep.effeff: 1.05, gtep.pipi: 1.020},
    {gtep.rf: 1.000, gtep.m: 1.000, gtep.effeff: 1.00, gtep.pipi: 0.990},
    {gtep.rf: 1.000, gtep.m: 1.050, gtep.effeff: 0.95, gtep.pipi: 0.960},
    #
    {gtep.rf: 1.025, gtep.m: 0.960, gtep.effeff: 1.05, gtep.pipi: 1.060},
    {gtep.rf: 1.025, gtep.m: 0.975, gtep.effeff: 1.05, gtep.pipi: 1.045},
    {gtep.rf: 1.025, gtep.m: 1.025, gtep.effeff: 1.00, gtep.pipi: 1.020},
    {gtep.rf: 1.025, gtep.m: 1.070, gtep.effeff: 0.95, gtep.pipi: 0.980},
    #
    {gtep.rf: 1.050, gtep.m: 1.050, gtep.effeff: 1.00, gtep.pipi: 1.040},
    {gtep.rf: 1.050, gtep.m: 1.090, gtep.effeff: 0.95, gtep.pipi: 1.000},
)

characteristics: Dict[str, Dict[str, Callable]] = {
    "1.3..2.5": {
        gtep.effeff: Interpolator(characteristic_13_25, gtep.effeff, features=(gtep.rf, gtep.m), fill_value=nan),
        gtep.pipi: Interpolator(characteristic_13_25, gtep.pipi, features=(gtep.rf, gtep.m), fill_value=nan),
    },
    # TODO: add more
}

if __name__ == "__main__":
    print(characteristics)
