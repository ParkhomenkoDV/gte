import numpy as np
import pytest
from substance import Substance
from thermodynamics import gas_const, heat_capacity_at_constant_pressure

from src.config import parameters as gtep
from src.nodes import CombustionChamber, Compressor


@pytest.fixture
def air():
    """Воздух для тестов"""
    return Substance(
        "air",
        parameters={
            gtep.gc: 287.14,
            gtep.TT: 300.0,
            gtep.PP: 101325.0,
            gtep.mf: 100.0,
            gtep.Cp: 1006.0,
            gtep.k: 1.4,
            gtep.c: 0.0,
        },
        functions={
            gtep.gc: lambda total_temperature: gas_const("air"),
            gtep.Cp: lambda total_temperature: heat_capacity_at_constant_pressure("air", total_temperature),
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
            (Compressor(), {gtep.pipi: 6.0, gtep.effeff: 0.85}),
            (Compressor(), {gtep.effeff: 0.85, gtep.power: 24 * 10**6}),
            (Compressor(), {gtep.pipi: 6.0, gtep.power: 24 * 10**6}),
        ],
    )
    def test_calculate(self, benchmark, node, kwargs, air):
        def benchfunc(node, kwargs):
            for k, v in kwargs.items():
                setattr(node, k, v)

            node.calculate(air)

            for k in node.variables:
                setattr(node, k, np.nan)

        benchmark(benchfunc, node, kwargs)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s", "--benchmark-only", "--benchmark-min-rounds=10"])
