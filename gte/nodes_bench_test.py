import numpy as np
import pytest
from combustion_chambler import CombustionChamber
from compressor import Compressor
from config import parameters as gtep
from fixtures import air, exhaust, kerosene
from turbine import Turbine


@pytest.mark.benchmark
class TestCompressor:
    """Бенчмарки Compressor"""

    def test_init(self, benchmark):
        def benchfunc():
            return Compressor()

        node = benchmark(benchfunc)
        assert isinstance(node, Compressor)

    @pytest.mark.parametrize(
        "node, kwargs",
        [
            (Compressor(), {gtep.pipi: 6.0, gtep.effeff: 0.85}),
            (Compressor(), {gtep.effeff: 0.85, gtep.power: 24 * 10**6}),
            (Compressor(), {gtep.pipi: 6.0, gtep.power: 24 * 10**6}),
        ],
    )
    def test_calculate(self, benchmark, node, kwargs):
        def benchfunc(node, kwargs):
            for k, v in kwargs.items():
                setattr(node, k, v)

            node.calculate(air)

            for k in node.variables:
                setattr(node, k, np.nan)

        benchmark(benchfunc, node, kwargs)


@pytest.mark.benchmark
class TestCombustionChamber:
    """Бенчмарки CombustionChamber"""

    def test_init(self, benchmark):
        def benchfunc():
            return CombustionChamber()

        node = benchmark(benchfunc)
        assert isinstance(node, CombustionChamber)

    @pytest.mark.parametrize(
        "node, kwargs",
        [
            (CombustionChamber(), {gtep.peff: 0.95, "efficiency_burn": 0.99}),
        ],
    )
    def test_calculate(self, benchmark, node, kwargs):
        def benchfunc(node, kwargs):
            for k, v in kwargs.items():
                setattr(node, k, v)

            node.calculate(air, kerosene)

            for k in node.variables:
                setattr(node, k, np.nan)

        benchmark(benchfunc, node, kwargs)


@pytest.mark.benchmark
class TestTurbine:
    """Бенчмарки Turbine"""

    def test_init(self, benchmark):
        def benchfunc():
            return Turbine()

        node = benchmark(benchfunc)
        assert isinstance(node, Turbine)

    @pytest.mark.parametrize(
        "node, kwargs",
        [
            (Turbine(), {gtep.pipi: 4.0, gtep.effeff: 0.9}),
            (Turbine(), {gtep.effeff: 0.9, gtep.power: 24 * 10**6}),
            (Turbine(), {gtep.pipi: 4.0, gtep.power: 24 * 10**6}),
        ],
    )
    def test_calculate(self, benchmark, node, kwargs):
        def benchfunc(node, kwargs):
            for k, v in kwargs.items():
                setattr(node, k, v)

            node.calculate(exhaust)

            for k in node.variables:
                setattr(node, k, np.nan)

        benchmark(benchfunc, node, kwargs)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s", "-x", "--benchmark-only", "--benchmark-min-rounds=10"])
