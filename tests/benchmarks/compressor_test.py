import pytest
from substance import Substance
from thermodynamics import gas_const, heat_capacity_at_constant_pressure

from src.config import parameters as gtep
from src.nodes.compressor import Compressor


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
class TestCompressorBenchmarks:
    """Бенчмарки производительности компрессора"""

    def benchmark_compressor_calculation(self, benchmark):
        compressor = Compressor()
        compressor.pressure_ratio = 5.0
        compressor.efficiency = 0.87
        substance = air()

        result = benchmark(compressor.calculate, substance)
        assert result is not None

    def benchmark_multiple_calculations(self, benchmark):
        def setup():
            compressor = Compressor()
            substances = [air(temperature=300 + i * 100) for i in range(10)]
            return (compressor, substances)

        def calculate_all(compressor, substances):
            return [compressor.calculate(sub) for sub in substances]

        benchmark.pedantic(calculate_all, setup=setup, rounds=100)


if __name__ == "__main__":
    pytest.main([__file__, "--benchmark-only"])
