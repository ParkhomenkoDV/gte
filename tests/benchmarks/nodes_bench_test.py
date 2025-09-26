import numpy as np
import pytest
from substance import Substance
from thermodynamics import T0, gas_const, heat_capacity_at_constant_pressure

from src.config import parameters as gtep
from src.nodes import CombustionChamber, Compressor


@pytest.fixture
def air():
    """Воздух для тестов"""
    return Substance(
        "air",
        parameters={
            gtep.mf: 100.0,
            gtep.gc: 287.14,
            gtep.TT: 300.0,
            gtep.PP: 101325.0,
            gtep.Cp: 1006.0,
            gtep.k: 1.4,
            gtep.c: 0.0,
        },
        functions={
            gtep.gc: lambda total_temperature: gas_const("air"),
            gtep.Cp: lambda total_temperature: heat_capacity_at_constant_pressure("air", total_temperature),
        },
    )


@pytest.fixture
def fuel():
    """Горючее для тестов"""
    return Substance(
        "kerosene",
        parameters={
            gtep.mf: 1.0,
            gtep.gc: 287.0,
            gtep.TT: T0 + 40,
            gtep.PP: 101325.0,
            gtep.Cp: 1006.0,
            gtep.k: 1.4,
            gtep.c: 0.0,
        },
        functions={
            gtep.Cp: lambda total_temperature: 200,
        },
    )


@pytest.mark.benchmark
class TestNodes:
    """Бенчмарки производительности компрессора"""

    @pytest.mark.parametrize(
        "Node, name",
        [
            # Compressor
            (Compressor, ""),
            (Compressor, "test"),
            # CombustionChamber
            (CombustionChamber, ""),
            (CombustionChamber, "test"),
        ],
    )
    def test_init(self, benchmark, Node, name):
        def benchfunc(name):
            return Node(name)

        node = benchmark(benchfunc, name)
        assert isinstance(node, Node)
        assert hasattr(node, "inlet")
        assert hasattr(node, "outlet")

    @pytest.mark.parametrize(
        "node, kwargs",
        [
            # Compressor
            (Compressor(), {gtep.pipi: 6.0, gtep.effeff: 0.85}),
            (Compressor(), {gtep.effeff: 0.85, gtep.power: 24 * 10**6}),
            (Compressor(), {gtep.pipi: 6.0, gtep.power: 24 * 10**6}),
            # CombustionChamber
            (CombustionChamber(), {gtep.peff: 0.95, "efficiency_burn": 0.99}),
        ],
    )
    def test_calculate(self, benchmark, node, kwargs, air, fuel):
        def benchfunc(node, kwargs):
            for k, v in kwargs.items():
                setattr(node, k, v)

            if isinstance(node, CombustionChamber):
                node.calculate(air, fuel)
            else:
                node.calculate(air)

            for k in node.variables:
                setattr(node, k, np.nan)

        benchmark(benchfunc, node, kwargs)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s", "--benchmark-only", "--benchmark-min-rounds=10"])
