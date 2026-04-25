import pytest
from substance import Substance

try:
    from ..config import EPSREL
    from ..config import parameters as gtep
    from ..fixtures import air, exhaust, kerosene
    from .combustion_chamber.combustion_chamber import CombustionChamber
    from .nozzle.nozzle import Nozzle
    from .turbocompressor.compressor.compressor import Compressor
    from .turbocompressor.turbine.turbine import Turbine
except ImportError:
    import os
    import sys

    sys.path.insert(0, os.getcwd())

    from gte.config import EPSREL
    from gte.config import parameters as gtep
    from gte.fixtures import air, exhaust, kerosene
    from gte.nodes.combustion_chamber.combustion_chamber import CombustionChamber
    from gte.nodes.nozzle import Nozzle
    from gte.nodes.turbocompressor.compressor.compressor import Compressor
    from gte.nodes.turbocompressor.turbine.turbine import Turbine


class TestNode:
    """Тесты абстрактного класса GTENode"""

    @pytest.mark.parametrize("Node", [Compressor, CombustionChamber, Turbine, Nozzle])
    def test_name(self, Node):
        node = Node({})  # empty parameters
        assert node.name == Node.__name__  # default
        node.name = "Test"
        assert node.name == "Test"  # reset
        node = Node({}, name="Node")
        assert node.name == "Node"  # __init__


@pytest.fixture
def compressor():
    """Создает экземпляр компрессора"""
    return Compressor({}, "Test")


class TestCompressor:
    """Тесты для класса Compressor"""

    def test_init(self, compressor):
        """Тест инициализации компрессора"""
        assert compressor.name == "Test"

    @pytest.mark.benchmark
    def test_compressor_init(self, benchmark):
        """Бенчмарк инициализации компрессора"""

        def benchfunc():
            return Compressor({}, name="Bench")

        benchmark(benchfunc)

    @pytest.mark.parametrize(
        "inlet, parameters, expected_parameters, expected_outlet",
        [
            (
                air,
                {gtep.power: 12 * 10**6, gtep.effeff: 0.85},  # pipi
                {gtep.pipi: 6.098},
                Substance("outlet", parameters={gtep.TT: 539.245, gtep.PP: 617_911.9}),
            ),
            (
                air,
                {gtep.pipi: 6.0, gtep.power: 12 * 10**6},  # effeff
                {gtep.effeff: 0.8402},
                Substance("outlet", parameters={gtep.TT: 539.245, gtep.PP: 607_950.0}),
            ),
            (
                air,
                {gtep.pipi: 6.0, gtep.effeff: 0.85},  # power
                {gtep.power: 11_862_134},
                Substance("outlet", parameters={gtep.TT: 536.5, gtep.PP: 607_950.0}),
            ),
        ],
    )
    def test_predict(self, compressor, inlet, parameters, expected_parameters, expected_outlet):
        """Тест предсказания компресоора"""
        vars, outlet = compressor.predict(inlet, parameters)
        for k, v in expected_parameters.items():
            assert vars[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")
        for k, v in expected_outlet.parameters.items():
            assert outlet.parameters[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")

    @pytest.mark.parametrize(
        "inlet, parameters",
        [
            (
                air,
                {gtep.power: 12 * 10**6, gtep.effeff: 0.85},  # pipi
            ),
            (
                air,
                {gtep.pipi: 6.0, gtep.power: 12 * 10**6},  # effeff
            ),
            (
                air,
                {gtep.pipi: 6.0, gtep.effeff: 0.85},  # power
            ),
        ],
    )
    @pytest.mark.benchmark
    def test_compressor_predict(self, benchmark, compressor, inlet, parameters):
        """Бенчмарк предсказания компрессора"""

        def benchfunc(compressor, inlet, parameters):
            compressor.predict(inlet, parameters)

        benchmark(benchfunc, compressor, inlet, parameters)

    @pytest.mark.parametrize(
        "inlet, parameters, expected_parameters, expected_outlet",
        [
            (
                air,
                {gtep.power: 12 * 10**6, gtep.effeff: 0.85},  # pipi
                {gtep.pipi: 6.098},
                Substance("outlet", parameters={gtep.TT: 539.245, gtep.PP: 617_911.9}),
            ),
            (
                air,
                {gtep.pipi: 6.0, gtep.power: 12 * 10**6},  # effeff
                {gtep.effeff: 0.8402},
                Substance("outlet", parameters={gtep.TT: 539.245, gtep.PP: 607_950.0}),
            ),
            (
                air,
                {gtep.pipi: 6.0, gtep.effeff: 0.85},  # power
                {gtep.power: 11_816_242},
                Substance("outlet", parameters={gtep.TT: 532.3, gtep.PP: 607_950.0}),
            ),
        ],
    )
    def test_calculate(self, compressor, inlet, parameters, expected_parameters, expected_outlet):
        """Тест расчета компресоора"""
        vars, outlet = compressor.calculate(inlet, parameters)
        for k, v in expected_parameters.items():
            assert vars[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")
        for k, v in expected_outlet.parameters.items():
            assert outlet.parameters[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")

    @pytest.mark.parametrize(
        "inlet, parameters",
        [
            (
                air,
                {gtep.power: 12 * 10**6, gtep.effeff: 0.85},  # pipi
            ),
            (
                air,
                {gtep.pipi: 6.0, gtep.power: 12 * 10**6},  # effeff
            ),
            (
                air,
                {gtep.pipi: 6.0, gtep.effeff: 0.85},  # power
            ),
        ],
    )
    @pytest.mark.benchmark
    def test_compressor_calculate(self, benchmark, compressor, inlet, parameters):
        """Бенчмарк предсказания компрессора"""

        def benchfunc(compressor, inlet, parameters):
            compressor.calculate(inlet, parameters)

        benchmark(benchfunc, compressor, inlet, parameters)


@pytest.fixture
def combustion_chamber():
    """Создает экземпляр камеры сгорания"""
    return CombustionChamber({}, "Test")


class TestCombustionChamber:
    """Тесты для класса CombustionChamber"""

    def test_init(self, combustion_chamber):
        """Тест инициализации камеры сгорания"""
        assert combustion_chamber.name == "Test"

    @pytest.mark.benchmark
    def test_cc_init(self, benchmark):
        """Бенчмарк инициализации камеры сгорания"""

        def benchfunc():
            return CombustionChamber({}, "Bench")

        benchmark(benchfunc)

    @pytest.mark.parametrize(
        "inlet, fuel, parameters, expected_outlet",
        [
            (
                air,
                kerosene,
                {gtep.eff_burn: 0.99, gtep.pipi: 0.95},
                Substance("outlet", parameters={gtep.TT: 1137.5, gtep.PP: 96_258}),
            ),
        ],
    )
    def test_predict(self, combustion_chamber, inlet, fuel, parameters, expected_outlet):
        """Тест предсказания камеры сгорания"""
        vars, outlet = combustion_chamber.predict(inlet, fuel, parameters)
        for k, v in parameters.items():
            assert vars[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")
        for k, v in expected_outlet.parameters.items():
            assert outlet.parameters[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")

    @pytest.mark.parametrize(
        "inlet, fuel, parameters",
        [
            (
                air,
                kerosene,
                {gtep.eff_burn: 0.99, gtep.pipi: 0.95},
            ),
        ],
    )
    @pytest.mark.benchmark
    def test_cc_predict(self, benchmark, combustion_chamber, inlet, fuel, parameters):
        """Бенчмарк предсказания камеры сгорания"""

        def benchfunc(combustion_chamber, inlet, fuel, parameters):
            combustion_chamber.predict(inlet, fuel, parameters)

        benchmark(benchfunc, combustion_chamber, inlet, fuel, parameters)

    @pytest.mark.parametrize(
        "inlet, fuel, parameters, expected_outlet",
        [
            (
                air,
                kerosene,
                {gtep.eff_burn: 0.99, gtep.pipi: 0.95},
                Substance("outlet", parameters={gtep.TT: 1079, gtep.PP: 96_258}),
            ),
        ],
    )
    def test_calculate(self, combustion_chamber, inlet, fuel, parameters, expected_outlet):
        """Тест предсказания камеры сгорания"""
        vars, outlet = combustion_chamber.calculate(inlet, fuel, parameters)
        for k, v in parameters.items():
            assert vars[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")
        for k, v in expected_outlet.parameters.items():
            assert outlet.parameters[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")

    @pytest.mark.parametrize(
        "inlet, fuel, parameters",
        [
            (
                air,
                kerosene,
                {gtep.eff_burn: 0.99, gtep.pipi: 0.95},
            ),
        ],
    )
    @pytest.mark.benchmark
    def test_cc_calculate(self, benchmark, combustion_chamber, inlet, fuel, parameters):
        """Бенчмарк расчета камеры сгорания"""

        def benchfunc(combustion_chamber, inlet, fuel, parameters):
            combustion_chamber.calculate(inlet, fuel, parameters)

        benchmark(benchfunc, combustion_chamber, inlet, fuel, parameters)


@pytest.fixture
def turbine():
    """Создает экземпляр турбины"""
    return Turbine({}, "Test")


class TestTurbine:
    """Тесты для класса Turbine"""

    def test_init(self, turbine):
        """Тест инициализации турбины"""
        assert turbine.name == "Test"

    @pytest.mark.benchmark
    def test_turbine_init(self, benchmark):
        """Бенчмарк инициализации турбины"""

        def benchfunc():
            return Turbine({}, name="Bench")

        benchmark(benchfunc)

    @pytest.mark.parametrize(
        "inlet, parameters, expected_parameters, expected_outlet",
        [
            (
                exhaust,
                {gtep.power: -24 * 10**6, gtep.effeff: 0.9},
                {gtep.pipi: 1 / 4.16},
                Substance("outlet", parameters={gtep.TT: 1103.7, gtep.PP: 292837.9}),
            ),
            (
                exhaust,
                {gtep.pipi: 1 / 4.0, gtep.power: -24 * 10**6},
                {gtep.effeff: 0.9203},
                Substance("outlet", parameters={gtep.TT: 1103.7, gtep.PP: 303975.0}),
            ),
            (
                exhaust,
                {gtep.pipi: 1 / 4.0, gtep.effeff: 0.9},
                {gtep.power: -23_468_473},
                Substance("outlet", parameters={gtep.TT: 1121.8, gtep.PP: 303975.0}),
            ),
        ],
    )
    def test_predict(self, turbine, inlet, parameters, expected_parameters, expected_outlet):
        """Тест предсказания турбины"""
        vars, outlet = turbine.predict(inlet, parameters)
        for k, v in expected_parameters.items():
            assert vars[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")
        for k, v in expected_outlet.parameters.items():
            assert outlet.parameters[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")

    @pytest.mark.parametrize(
        "inlet, parameters",
        [
            (
                exhaust,
                {gtep.power: 24 * 10**6, gtep.effeff: 0.9},
            ),
            (
                exhaust,
                {gtep.pipi: 4.0, gtep.power: 24 * 10**6},
            ),
            (
                exhaust,
                {gtep.pipi: 4.0, gtep.effeff: 0.9},
            ),
        ],
    )
    @pytest.mark.benchmark
    def test_turbine_predict(self, benchmark, turbine, inlet, parameters):
        def benchfunc(turbine, inlet, parameters):
            turbine.predict(inlet, parameters)

        benchmark(benchfunc, turbine, inlet, parameters)

    @pytest.mark.parametrize(
        "inlet, parameters, expected_parameters, expected_outlet",
        [
            (
                exhaust,
                {gtep.power: -24 * 10**6, gtep.effeff: 0.9},
                {gtep.pipi: 1 / 4.16},
                Substance("outlet", parameters={gtep.TT: 1103.7, gtep.PP: 292837.9}),
            ),
            (
                exhaust,
                {gtep.pipi: 1 / 4.0, gtep.power: -24 * 10**6},
                {gtep.effeff: 0.9232},
                Substance("outlet", parameters={gtep.TT: 1103.7, gtep.PP: 303975.0}),
            ),
            (
                exhaust,
                {gtep.pipi: 1 / 4.0, gtep.effeff: 0.9},
                {gtep.power: -23_396_982},
                Substance("outlet", parameters={gtep.TT: 1121.8, gtep.PP: 303975.0}),
            ),
        ],
    )
    def test_calculate(self, turbine, inlet, parameters, expected_parameters, expected_outlet):
        """Тест предсказания турбины"""
        vars, outlet = turbine.calculate(inlet, parameters)
        for k, v in expected_parameters.items():
            assert vars[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")
        for k, v in expected_outlet.parameters.items():
            assert outlet.parameters[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")

    @pytest.mark.parametrize(
        "inlet, parameters",
        [
            (
                exhaust,
                {gtep.power: 24 * 10**6, gtep.effeff: 0.9},
            ),
            (
                exhaust,
                {gtep.pipi: 4.0, gtep.power: 24 * 10**6},
            ),
            (
                exhaust,
                {gtep.pipi: 4.0, gtep.effeff: 0.9},
            ),
        ],
    )
    @pytest.mark.benchmark
    def test_turbine_calculate(self, benchmark, turbine, inlet, parameters):
        def benchfunc(turbine, inlet, parameters):
            turbine.calculate(inlet, parameters)

        benchmark(benchfunc, turbine, inlet, parameters)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s", "-x"])
