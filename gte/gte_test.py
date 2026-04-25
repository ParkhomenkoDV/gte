import pytest
from substance import Substance

try:
    from .config import EPSREL
    from .config import parameters as gtep
    from .fixtures import air, exhaust, kerosene
    from .gte import GTE
    from .nodes.combustion_chamber.combustion_chamber import CombustionChamber
    from .nodes.nozzle.nozzle import Nozzle
    from .nodes.turbocompressor.compressor.compressor import Compressor
    from .nodes.turbocompressor.turbine.turbine import Turbine
except ImportError:
    import os
    import sys

    sys.path.insert(0, os.getcwd())

    from gte.config import EPSREL
    from gte.config import parameters as gtep
    from gte.fixtures import air, exhaust, kerosene
    from gte.gte import GTE
    from gte.nodes.combustion_chamber.combustion_chamber import CombustionChamber
    from gte.nodes.nozzle import Nozzle
    from gte.nodes.turbocompressor.compressor.compressor import Compressor
    from gte.nodes.turbocompressor.turbine.turbine import Turbine


class TestGTE:
    """Тесты для класса GTE"""

    @pytest.mark.parametrize(
        "scheme, name",
        [
            ([[Compressor({}), CombustionChamber({}), Turbine({})]], "AI-9"),
            ([[Compressor({}), CombustionChamber({}), Turbine({}), Nozzle({})]], "Jumo 004"),
        ],
    )
    def test_init(self, scheme, name):
        """Тест инициализации ГТД"""
        gte = GTE(scheme, name)
        assert gte.name == name

    @pytest.mark.parametrize(
        "scheme, name",
        [
            ([[Compressor({}), CombustionChamber({}), Turbine({})]], "AI-9"),
            ([[Compressor({}), CombustionChamber({}), Turbine({}), Nozzle({})]], "Jumo 004"),
        ],
    )
    @pytest.mark.benchmark
    def test_gte_init(self, benchmark, scheme, name):
        """Бенчмарк инициализации ГТД"""

        def benchfunc(scheme, name):
            GTE(scheme, name)

        benchmark(benchfunc, scheme, name)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s", "-x"])
