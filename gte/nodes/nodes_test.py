import pytest
from substance import Substance

try:
    from ..config import EPSREL
    from ..config import parameters as gtep
    from ..fixtures import air, exhaust, kerosene
    from .burner.burner import Burner
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
    from gte.nodes.burner.burner import Burner
    from gte.nodes.nozzle import Nozzle
    from gte.nodes.turbocompressor.compressor.compressor import Compressor
    from gte.nodes.turbocompressor.turbine.turbine import Turbine


class TestNode:
    """Тесты абстрактного класса GTENode"""

    @pytest.mark.parametrize("Node", [Compressor, Burner, Turbine, Nozzle])
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
        "parameters, inlet, expected_parameters, expected_outlet",
        [
            (
                {gtep.power: 12 * 10**6, gtep.effeff: 0.85},  # pipi
                air,
                {gtep.pipi: 6.098},
                Substance("outlet", parameters={gtep.TT: 539.245, gtep.PP: 617_911.9}),
            ),
            (
                {gtep.pipi: 6.0, gtep.power: 12 * 10**6},  # effeff
                air,
                {gtep.effeff: 0.8402},
                Substance("outlet", parameters={gtep.TT: 539.245, gtep.PP: 607_950.0}),
            ),
            (
                {gtep.pipi: 6.0, gtep.effeff: 0.85},  # power
                air,
                {gtep.power: 11_862_134},
                Substance("outlet", parameters={gtep.TT: 536.5, gtep.PP: 607_950.0}),
            ),
        ],
    )
    def test_predict(self, compressor, parameters, inlet, expected_parameters, expected_outlet):
        """Тест предсказания компресоора"""
        vars, outlet = compressor.predict(parameters, inlet)
        for k, v in expected_parameters.items():
            assert vars[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")
        for k, v in expected_outlet.parameters.items():
            assert outlet.parameters[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")

    @pytest.mark.parametrize(
        "parameters, inlet",
        [
            (
                {gtep.power: 12 * 10**6, gtep.effeff: 0.85},  # pipi
                air,
            ),
            (
                {gtep.pipi: 6.0, gtep.power: 12 * 10**6},  # effeff
                air,
            ),
            (
                {gtep.pipi: 6.0, gtep.effeff: 0.85},  # power
                air,
            ),
        ],
    )
    @pytest.mark.benchmark
    def test_compressor_predict(self, benchmark, compressor, parameters, inlet):
        """Бенчмарк предсказания компрессора"""

        def benchfunc(compressor, parameters, inlet):
            compressor.predict(parameters, inlet)

        benchmark(benchfunc, compressor, parameters, inlet)

    @pytest.mark.parametrize(
        "parameters, inlet, expected_parameters, expected_outlet",
        [
            (
                {gtep.power: 12 * 10**6, gtep.effeff: 0.85},  # pipi
                air,
                {gtep.pipi: 6.098},
                Substance("outlet", parameters={gtep.TT: 539.245, gtep.PP: 617_911.9}),
            ),
            (
                {gtep.pipi: 6.0, gtep.power: 12 * 10**6},  # effeff
                air,
                {gtep.effeff: 0.8402},
                Substance("outlet", parameters={gtep.TT: 539.245, gtep.PP: 607_950.0}),
            ),
            (
                {gtep.pipi: 6.0, gtep.effeff: 0.85},  # power
                air,
                {gtep.power: 11_816_242},
                Substance("outlet", parameters={gtep.TT: 532.3, gtep.PP: 607_950.0}),
            ),
        ],
    )
    def test_calculate(self, compressor, parameters, inlet, expected_parameters, expected_outlet):
        """Тест расчета компресоора"""
        vars, outlet = compressor.calculate(parameters, inlet)
        for k, v in expected_parameters.items():
            assert vars[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")
        for k, v in expected_outlet.parameters.items():
            assert outlet.parameters[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")

    @pytest.mark.parametrize(
        "parameters, inlet",
        [
            (
                {gtep.power: 12 * 10**6, gtep.effeff: 0.85},  # pipi
                air,
            ),
            (
                {gtep.pipi: 6.0, gtep.power: 12 * 10**6},  # effeff
                air,
            ),
            (
                {gtep.pipi: 6.0, gtep.effeff: 0.85},  # power
                air,
            ),
        ],
    )
    @pytest.mark.benchmark
    def test_compressor_calculate(self, benchmark, compressor, parameters, inlet):
        """Бенчмарк предсказания компрессора"""

        def benchfunc(compressor, parameters, inlet):
            compressor.calculate(parameters, inlet)

        benchmark(benchfunc, compressor, parameters, inlet)


@pytest.fixture
def burner():
    """Создает экземпляр камеры сгорания"""
    return Burner({}, "Test")


class TestBurner:
    """Тесты для класса Burner"""

    def test_init(self, burner):
        """Тест инициализации камеры сгорания"""
        assert burner.name == "Test"

    @pytest.mark.benchmark
    def test_burner_init(self, benchmark):
        """Бенчмарк инициализации камеры сгорания"""

        def benchfunc():
            return Burner({}, "Bench")

        benchmark(benchfunc)

    @pytest.mark.parametrize(
        "parameters, inlet, fuel, expected_outlet",
        [
            (
                {gtep.eff_burn: 0.99, gtep.pipi: 0.95},
                air,
                kerosene,
                Substance("outlet", parameters={gtep.TT: 1137.5, gtep.PP: 96_258}),
            ),
        ],
    )
    def test_predict(self, burner, parameters, inlet, fuel, expected_outlet):
        """Тест предсказания камеры сгорания"""
        vars, outlet = burner.predict(parameters, inlet, fuel)
        for k, v in parameters.items():
            assert vars[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")
        for k, v in expected_outlet.parameters.items():
            assert outlet.parameters[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")

    @pytest.mark.parametrize(
        "parameters, inlet, fuel",
        [
            (
                {gtep.eff_burn: 0.99, gtep.pipi: 0.95},
                air,
                kerosene,
            ),
        ],
    )
    @pytest.mark.benchmark
    def test_burner_predict(self, benchmark, burner, parameters, inlet, fuel):
        """Бенчмарк предсказания камеры сгорания"""

        def benchfunc(burner, parameters, inlet, fuel):
            burner.predict(parameters, inlet, fuel)

        benchmark(benchfunc, burner, parameters, inlet, fuel)

    @pytest.mark.parametrize(
        "parameters, inlet, fuel, expected_outlet",
        [
            (
                {gtep.eff_burn: 0.99, gtep.pipi: 0.95},
                air,
                kerosene,
                Substance("outlet", parameters={gtep.TT: 1079, gtep.PP: 96_258}),
            ),
        ],
    )
    def test_calculate(self, burner, parameters, inlet, fuel, expected_outlet):
        """Тест предсказания камеры сгорания"""
        vars, outlet = burner.calculate(parameters, inlet, fuel)
        for k, v in parameters.items():
            assert vars[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")
        for k, v in expected_outlet.parameters.items():
            assert outlet.parameters[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")

    @pytest.mark.parametrize(
        "parameters, inlet, fuel",
        [
            (
                {gtep.eff_burn: 0.99, gtep.pipi: 0.95},
                air,
                kerosene,
            ),
        ],
    )
    @pytest.mark.benchmark
    def test_burner_calculate(self, benchmark, burner, parameters, inlet, fuel):
        """Бенчмарк расчета камеры сгорания"""

        def benchfunc(burner, parameters, inlet, fuel):
            burner.calculate(parameters, inlet, fuel)

        benchmark(benchfunc, burner, parameters, inlet, fuel)


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
        "parameters, inlet, expected_parameters, expected_outlet",
        [
            (
                {gtep.power: -24 * 10**6, gtep.effeff: 0.9},
                exhaust,
                {gtep.pipi: 1 / 4.16},
                Substance("outlet", parameters={gtep.TT: 1103.7, gtep.PP: 292837.9}),
            ),
            (
                {gtep.pipi: 1 / 4.0, gtep.power: -24 * 10**6},
                exhaust,
                {gtep.effeff: 0.9203},
                Substance("outlet", parameters={gtep.TT: 1103.7, gtep.PP: 303975.0}),
            ),
            (
                {gtep.pipi: 1 / 4.0, gtep.effeff: 0.9},
                exhaust,
                {gtep.power: -23_468_473},
                Substance("outlet", parameters={gtep.TT: 1121.8, gtep.PP: 303975.0}),
            ),
        ],
    )
    def test_predict(self, turbine, parameters, inlet, expected_parameters, expected_outlet):
        """Тест предсказания турбины"""
        vars, outlet = turbine.predict(parameters, inlet)
        for k, v in expected_parameters.items():
            assert vars[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")
        for k, v in expected_outlet.parameters.items():
            assert outlet.parameters[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")

    @pytest.mark.parametrize(
        "parameters, inlet",
        [
            (
                {gtep.power: 24 * 10**6, gtep.effeff: 0.9},
                exhaust,
            ),
            (
                {gtep.pipi: 4.0, gtep.power: 24 * 10**6},
                exhaust,
            ),
            (
                {gtep.pipi: 4.0, gtep.effeff: 0.9},
                exhaust,
            ),
        ],
    )
    @pytest.mark.benchmark
    def test_turbine_predict(self, benchmark, turbine, parameters, inlet):
        def benchfunc(turbine, parameters, inlet):
            turbine.predict(parameters, inlet)

        benchmark(benchfunc, turbine, parameters, inlet)

    @pytest.mark.parametrize(
        "parameters, inlet, expected_parameters, expected_outlet",
        [
            (
                {gtep.power: -24 * 10**6, gtep.effeff: 0.9},
                exhaust,
                {gtep.pipi: 1 / 4.16},
                Substance("outlet", parameters={gtep.TT: 1103.7, gtep.PP: 292837.9}),
            ),
            (
                {gtep.pipi: 1 / 4.0, gtep.power: -24 * 10**6},
                exhaust,
                {gtep.effeff: 0.9232},
                Substance("outlet", parameters={gtep.TT: 1103.7, gtep.PP: 303975.0}),
            ),
            (
                {gtep.pipi: 1 / 4.0, gtep.effeff: 0.9},
                exhaust,
                {gtep.power: -23_396_982},
                Substance("outlet", parameters={gtep.TT: 1121.8, gtep.PP: 303975.0}),
            ),
        ],
    )
    def test_calculate(self, turbine, parameters, inlet, expected_parameters, expected_outlet):
        """Тест предсказания турбины"""
        vars, outlet = turbine.calculate(parameters, inlet)
        for k, v in expected_parameters.items():
            assert vars[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")
        for k, v in expected_outlet.parameters.items():
            assert outlet.parameters[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")

    @pytest.mark.parametrize(
        "parameters, inlet",
        [
            (
                {gtep.power: 24 * 10**6, gtep.effeff: 0.9},
                exhaust,
            ),
            (
                {gtep.pipi: 4.0, gtep.power: 24 * 10**6},
                exhaust,
            ),
            (
                {gtep.pipi: 4.0, gtep.effeff: 0.9},
                exhaust,
            ),
        ],
    )
    @pytest.mark.benchmark
    def test_turbine_calculate(self, benchmark, turbine, parameters, inlet):
        def benchfunc(turbine, parameters, inlet):
            turbine.calculate(parameters, inlet)

        benchmark(benchfunc, turbine, parameters, inlet)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s", "-x"])
