import pytest
from substance import Substance

try:
    from .config import EPSREL
    from .config import parameters as gtep
    from .fixtures import air, exhaust, kerosene
    from .gte import GTE
    from .nodes.burner.burner import Burner
    from .nodes.channel.channel import Channel
    from .nodes.nozzle.nozzle import Nozzle
    from .nodes.turbocompressor.rotor.rotor import Rotor
    from .nodes.turbocompressor.turbine.turbine import Turbine
except ImportError:
    import os
    import sys

    sys.path.insert(0, os.getcwd())

    from gte.config import EPSREL
    from gte.config import parameters as gtep
    from gte.fixtures import air, exhaust, kerosene
    from gte.gte import GTE
    from gte.nodes.burner.burner import Burner
    from gte.nodes.channel.channel import Channel
    from gte.nodes.nozzle import Nozzle
    from gte.nodes.turbocompressor.rotor.rotor import Rotor
    from gte.nodes.turbocompressor.turbine.turbine import Turbine


class TestGTE:
    """Тесты для класса GTE"""

    @pytest.mark.parametrize(
        "scheme, name",
        [
            ([[Rotor({}), Burner({}), Turbine({})]], "AI-9"),
            ([[Rotor({}), Burner({}), Turbine({}), Turbine({})]], "TV3-117"),
            ([[Rotor({}), Burner({}), Turbine({}), Nozzle({})]], "Jumo 004"),
            (
                [
                    [Rotor({}), Rotor({}), Burner({}), Turbine({}), Turbine({}), Nozzle({})],
                    [Rotor({}), Channel({})],
                ],
                "AL-31F",
            ),
        ],
    )
    def test_init(self, scheme, name):
        """Тест инициализации ГТД"""
        gte = GTE(scheme, name)
        assert gte.name == name

    @pytest.mark.parametrize(
        "scheme, name",
        [
            ([[Rotor({}), Burner({}), Turbine({})]], "AI-9"),
            ([[Rotor({}), Burner({}), Turbine({}), Turbine({})]], "TV3-117"),
            ([[Rotor({}), Burner({}), Turbine({}), Nozzle({})]], "Jumo 004"),
            (
                [
                    [Rotor({}), Rotor({}), Burner({}), Turbine({}), Turbine({}), Nozzle({})],
                    [Rotor({}), Channel({})],
                ],
                "AL-31F",
            ),
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
