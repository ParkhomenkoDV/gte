import pytest
from substance import Substance

try:
    from ..config import EPSREL
    from ..config import parameters as gtep
    from ..fixtures import air, exhaust, kerosene
    from .burner.burner import Burner
    from .channel.channel import Channel
    from .joiner.joiner import Joiner
    from .nozzle.nozzle import Nozzle
    from .splitter.splitter import Splitter
    from .turbocompressor.rotor.rotor import Rotor
except ImportError:
    import os
    import sys

    sys.path.insert(0, os.getcwd())

    from gte.config import EPSREL
    from gte.config import parameters as gtep
    from gte.fixtures import air, exhaust, kerosene
    from gte.nodes.burner.burner import Burner
    from gte.nodes.channel.channel import Channel
    from gte.nodes.joiner.joiner import Joiner
    from gte.nodes.nozzle import Nozzle
    from gte.nodes.splitter.splitter import Splitter
    from gte.nodes.turbocompressor.rotor.rotor import Rotor


class TestNode:
    """Тесты абстрактного класса GTENode"""

    @pytest.mark.parametrize("Node", [Rotor, Burner, Nozzle])
    def test_name(self, Node):
        node = Node({})  # empty parameters
        assert node.name == Node.__name__  # default
        node.name = "Test"
        assert node.name == "Test"  # reset
        node = Node({}, name="Node")
        assert node.name == "Node"  # __init__


@pytest.fixture
def rotor():
    """Создает экземпляр ротора"""
    return Rotor({}, "Test")


class TestRotor:
    """Тесты для класса Rotor"""

    def test_init(self, rotor):
        """Тест инициализации"""
        assert rotor.name == "Test"

    @pytest.mark.benchmark
    def test_rotor_init(self, benchmark):
        """Бенчмарк инициализации"""

        def benchfunc():
            return Rotor({}, name="Bench")

        benchmark(benchfunc)

    @pytest.mark.parametrize(
        "parameters, inlet, expected_parameters, expected_outlet",
        [
            # compressor
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
            # turbine
            (
                {gtep.power: -24 * 10**6, gtep.effeff: 1 / 0.9},
                exhaust,
                {gtep.pipi: 1 / 4.16},
                Substance("outlet", parameters={gtep.TT: 1103.7, gtep.PP: 292837.9}),
            ),
            (
                {gtep.pipi: 1 / 4.0, gtep.power: -24 * 10**6},
                exhaust,
                {gtep.effeff: 1 / 0.9203},
                Substance("outlet", parameters={gtep.TT: 1103.7, gtep.PP: 303975.0}),
            ),
            (
                {gtep.pipi: 1 / 4.0, gtep.effeff: 1 / 0.9},
                exhaust,
                {gtep.power: -23_468_473},
                Substance("outlet", parameters={gtep.TT: 1121.8, gtep.PP: 303975.0}),
            ),
        ],
    )
    def test_predict(self, rotor, parameters, inlet, expected_parameters, expected_outlet):
        """Тест предсказания компресоора"""
        vars, outlet = rotor.predict(parameters, inlet)
        for k, v in expected_parameters.items():
            assert vars[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")
        for k, v in expected_outlet.parameters.items():
            assert outlet.parameters[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")

    @pytest.mark.parametrize(
        "parameters, inlet",
        [
            # compressor
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
            # turbine
            (
                {gtep.power: -24 * 10**6, gtep.effeff: 1 / 0.9},
                exhaust,
            ),
            (
                {gtep.pipi: 1 / 4.0, gtep.power: -24 * 10**6},
                exhaust,
            ),
            (
                {gtep.pipi: 1 / 4.0, gtep.effeff: 1 / 0.9},
                exhaust,
            ),
        ],
    )
    @pytest.mark.benchmark
    def test_rotor_predict(self, benchmark, rotor, parameters, inlet):
        """Бенчмарк предсказания компрессора"""

        def benchfunc(rotor, parameters, inlet):
            rotor.predict(parameters, inlet)

        benchmark(benchfunc, rotor, parameters, inlet)

    @pytest.mark.parametrize(
        "parameters, inlet, expected_parameters, expected_outlet",
        [
            # compressor
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
            # turbine
            (
                {gtep.power: -24 * 10**6, gtep.effeff: 1 / 0.9},
                exhaust,
                {gtep.pipi: 1 / 4.16},
                Substance("outlet", parameters={gtep.TT: 1103.7, gtep.PP: 292837.9}),
            ),
            (
                {gtep.pipi: 1 / 4.0, gtep.power: -24 * 10**6},
                exhaust,
                {gtep.effeff: 1 / 0.9232},
                Substance("outlet", parameters={gtep.TT: 1103.7, gtep.PP: 303975.0}),
            ),
            (
                {gtep.pipi: 1 / 4.0, gtep.effeff: 1 / 0.9},
                exhaust,
                {gtep.power: -23_396_982},
                Substance("outlet", parameters={gtep.TT: 1121.8, gtep.PP: 303975.0}),
            ),
        ],
    )
    def test_calculate(self, rotor, parameters, inlet, expected_parameters, expected_outlet):
        """Тест расчета компресоора"""
        vars, outlet = rotor.calculate(parameters, inlet)
        for k, v in expected_parameters.items():
            assert vars[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")
        for k, v in expected_outlet.parameters.items():
            assert outlet.parameters[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")

    @pytest.mark.parametrize(
        "parameters, inlet",
        [
            # compressor
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
            # turbine
            (
                {gtep.power: -24 * 10**6, gtep.effeff: 1 / 0.9},
                exhaust,
            ),
            (
                {gtep.pipi: 1 / 4.0, gtep.power: -24 * 10**6},
                exhaust,
            ),
            (
                {gtep.pipi: 1 / 4.0, gtep.effeff: 1 / 0.9},
                exhaust,
            ),
        ],
    )
    @pytest.mark.benchmark
    def test_rotor_calculate(self, benchmark, rotor, parameters, inlet):
        """Бенчмарк предсказания компрессора"""

        def benchfunc(rotor, parameters, inlet):
            rotor.calculate(parameters, inlet)

        benchmark(benchfunc, rotor, parameters, inlet)


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
def nozzle():
    """Создает экземпляр сопла"""
    return Nozzle({}, "Test")


class TestNozzle:
    """Тесты для класса Nozzle"""

    def test_init(self, nozzle):
        """Тест инициализации"""
        assert nozzle.name == "Test"

    @pytest.mark.benchmark
    def test_nozzle_init(self, benchmark):
        """Бенчмарк инициализации"""

        def benchfunc():
            return Nozzle({}, name="Bench")

        benchmark(benchfunc)

    @pytest.mark.parametrize(
        "parameters, inlet, expected_parameters, expected_outlet",
        [
            (
                {gtep.eff_speed: 0.98, gtep.pipi: 1 / 1.8},
                exhaust,
                {gtep.force: 34_439},
                Substance("outlet", parameters={gtep.TT: exhaust.parameters[gtep.TT], gtep.PP: exhaust.parameters[gtep.PP] / 1.8}),
            ),
            (
                {gtep.eff_speed: 0.98, gtep.force: 34_439},
                exhaust,
                {gtep.pipi: 1 / 1.8},
                Substance("outlet", parameters={gtep.TT: exhaust.parameters[gtep.TT], gtep.PP: exhaust.parameters[gtep.PP] / 1.8}),
            ),
            (
                {gtep.pipi: 1 / 1.8, gtep.force: 34_439},
                exhaust,
                {gtep.eff_speed: 0.98},
                Substance("outlet", parameters={gtep.TT: exhaust.parameters[gtep.TT], gtep.PP: exhaust.parameters[gtep.PP] / 1.8}),
            ),
        ],
    )
    def test_predict(self, nozzle, parameters, inlet, expected_parameters, expected_outlet):
        """Тест предсказания"""
        vars, outlet = nozzle.predict(parameters, inlet)
        for k, v in expected_parameters.items():
            assert vars[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")
        for k, v in expected_outlet.parameters.items():
            assert outlet.parameters[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")

    @pytest.mark.parametrize(
        "parameters, inlet",
        [
            (
                {gtep.eff_speed: 0.98, gtep.pipi: 1 / 1.8},
                exhaust,
            ),
            (
                {gtep.eff_speed: 0.98, gtep.force: 34_439},
                exhaust,
            ),
            (
                {gtep.pipi: 1 / 1.8, gtep.force: 34_439},
                exhaust,
            ),
        ],
    )
    @pytest.mark.benchmark
    def test_nozzle_predict(self, benchmark, nozzle, parameters, inlet):
        def benchfunc(nozzle, parameters, inlet):
            nozzle.predict(parameters, inlet)

        benchmark(benchfunc, nozzle, parameters, inlet)

    @pytest.mark.parametrize(
        "parameters, inlet, expected_parameters, expected_outlet",
        [
            (
                {gtep.eff_speed: 0.98, gtep.pipi: 1 / 1.8},
                exhaust,
                {gtep.force: 34_439},
                Substance("outlet", parameters={gtep.TT: exhaust.parameters[gtep.TT], gtep.PP: exhaust.parameters[gtep.PP] / 1.8}),
            ),
            (
                {gtep.eff_speed: 0.98, gtep.force: 34_439},
                exhaust,
                {gtep.pipi: 1 / 1.8},
                Substance("outlet", parameters={gtep.TT: exhaust.parameters[gtep.TT], gtep.PP: exhaust.parameters[gtep.PP] / 1.8}),
            ),
            (
                {gtep.pipi: 1 / 1.8, gtep.force: 34_439},
                exhaust,
                {gtep.eff_speed: 0.98},
                Substance("outlet", parameters={gtep.TT: exhaust.parameters[gtep.TT], gtep.PP: exhaust.parameters[gtep.PP] / 1.8}),
            ),
        ],
    )
    def test_calculate(self, nozzle, parameters, inlet, expected_parameters, expected_outlet):
        """Тест предсказания"""
        vars, outlet = nozzle.calculate(parameters, inlet)
        for k, v in expected_parameters.items():
            assert vars[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")
        for k, v in expected_outlet.parameters.items():
            assert outlet.parameters[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")

    @pytest.mark.parametrize(
        "parameters, inlet",
        [
            (
                {gtep.eff_speed: 0.98, gtep.pipi: 1 / 1.8},
                exhaust,
            ),
            (
                {gtep.eff_speed: 0.98, gtep.force: 34_439},
                exhaust,
            ),
            (
                {gtep.pipi: 1 / 1.8, gtep.force: 34_439},
                exhaust,
            ),
        ],
    )
    @pytest.mark.benchmark
    def test_nozzle_calculate(self, benchmark, nozzle, parameters, inlet):
        def benchfunc(nozzle, parameters, inlet):
            nozzle.calculate(parameters, inlet)

        benchmark(benchfunc, nozzle, parameters, inlet)


@pytest.fixture
def channel():
    """Создает экземпляр канала"""
    return Channel({}, "Test")


class TestChannel:
    """Тесты для класса Channel"""

    def test_init(self, channel):
        """Тест инициализации"""
        assert channel.name == "Test"

    @pytest.mark.benchmark
    def test_channel_init(self, benchmark):
        """Бенчмарк инициализации"""

        def benchfunc():
            return Channel({}, name="Bench")

        benchmark(benchfunc)

    @pytest.mark.parametrize(
        "parameters, inlet, expected_outlet",
        [
            (
                {gtep.titi: 1.05, gtep.pipi: 0.95},
                air,
                Substance("outlet", parameters={gtep.TT: air.parameters[gtep.TT] * 1.05, gtep.PP: air.parameters[gtep.PP] * 0.95}),
            ),
        ],
    )
    def test_predict(self, channel, parameters, inlet, expected_outlet):
        """Тест предсказания"""
        _, outlet = channel.predict(parameters, inlet)

        for k, v in expected_outlet.parameters.items():
            assert outlet.parameters[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")

    @pytest.mark.parametrize(
        "parameters, inlet",
        [
            (
                {gtep.titi: 1.05, gtep.pipi: 0.95},
                air,
            ),
        ],
    )
    @pytest.mark.benchmark
    def test_channel_predict(self, benchmark, channel, parameters, inlet):
        def benchfunc(channel, parameters, inlet):
            channel.predict(parameters, inlet)

        benchmark(benchfunc, channel, parameters, inlet)

    @pytest.mark.parametrize(
        "parameters, inlet, expected_outlet",
        [
            (
                {gtep.titi: 1.05, gtep.pipi: 0.95},
                air,
                Substance("outlet", parameters={gtep.TT: air.parameters[gtep.TT] * 1.05, gtep.PP: air.parameters[gtep.PP] * 0.95}),
            )
        ],
    )
    def test_calculate(self, channel, parameters, inlet, expected_outlet):
        """Тест предсказания"""
        _, outlet = channel.calculate(parameters, inlet)
        for k, v in expected_outlet.parameters.items():
            assert outlet.parameters[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")

    @pytest.mark.parametrize(
        "parameters, inlet",
        [
            (
                {gtep.titi: 1.05, gtep.pipi: 0.95},
                air,
            ),
        ],
    )
    @pytest.mark.benchmark
    def test_channel_calculate(self, benchmark, channel, parameters, inlet):
        def benchfunc(channel, parameters, inlet):
            channel.calculate(parameters, inlet)

        benchmark(benchfunc, channel, parameters, inlet)


@pytest.fixture
def splitter():
    """Создает экземпляр камеры отбора"""
    return Splitter({}, "Test")


class TestSplitter:
    """Тесты для класса Splitter"""

    def test_init(self, splitter):
        """Тест инициализации"""
        assert splitter.name == "Test"

    @pytest.mark.benchmark
    def test_splitter_init(self, benchmark):
        """Бенчмарк инициализации"""

        def benchfunc():
            return Splitter({}, name="Bench")

        benchmark(benchfunc)

    @pytest.mark.parametrize(
        "parameters, inlet, expected_outlets",
        [
            (
                {"splits": (1, 3, 6)},
                air,
                (
                    Substance("air", parameters={p for p in air.parameters if p != gtep.m}.update({gtep.m: 5}), functions=air.functions),
                    Substance("air", parameters={p for p in air.parameters if p != gtep.m}.update({gtep.m: 15}), functions=air.functions),
                    Substance("air", parameters={p for p in air.parameters if p != gtep.m}.update({gtep.m: 30}), functions=air.functions),
                ),
            ),
        ],
    )
    def test_predict(self, splitter, parameters, inlet, expected_outlets):
        """Тест предсказания"""
        _, outlets = splitter.predict(parameters, inlet)

        for i, outlet in enumerate(outlets):
            for k, v in expected_outlets[i].parameters.items():
                assert outlet.parameters[k] == pytest.approx(v, rel=EPSREL), AssertionError(f"{k=} {v=}")

    @pytest.mark.parametrize(
        "parameters, inlet",
        [
            (
                {"splits": (1, 3, 6)},
                air,
            ),
        ],
    )
    @pytest.mark.benchmark
    def test_splitter_predict(self, benchmark, splitter, parameters, inlet):
        def benchfunc(splitter, parameters, inlet):
            splitter.predict(parameters, inlet)

        benchmark(benchfunc, splitter, parameters, inlet)


@pytest.fixture
def joiner():
    """Создает экземпляр камеры смешения"""
    return Joiner({}, "Test")


class TestJoiner:
    """Тесты для класса Joiner"""

    def test_init(self, joiner):
        """Тест инициализации"""
        assert joiner.name == "Test"

    @pytest.mark.benchmark
    def test_joiner_init(self, benchmark):
        """Бенчмарк инициализации"""

        def benchfunc():
            return Joiner({}, name="Bench")

        benchmark(benchfunc)

    @pytest.mark.parametrize(
        "inlets, expected_outlet",
        [
            (
                (air, exhaust),
                Substance(
                    f"{air.name}+{exhaust.name}",
                    parameters={
                        gtep.m: air.parameters[gtep.m] + exhaust.parameters[gtep.m],
                        gtep.TT: (air.parameters[gtep.m] * air.parameters[gtep.TT] * air.parameters[gtep.hcp] + exhaust.parameters[gtep.m] * exhaust.parameters[gtep.TT] * exhaust.parameters[gtep.hcp])
                        / (air.parameters[gtep.m] * air.parameters[gtep.hcp] + exhaust.parameters[gtep.m] * exhaust.parameters[gtep.hcp]),
                        gtep.PP: (air.parameters[gtep.m] * air.parameters[gtep.PP] + exhaust.parameters[gtep.m] * exhaust.parameters[gtep.PP]) / (air.parameters[gtep.m] + exhaust.parameters[gtep.m]),
                    },
                ),
            ),
        ],
    )
    def test_predict(self, joiner, inlets, expected_outlet):
        """Тест предсказания"""
        _, outlet = joiner.predict({}, *inlets)

        for k, v in expected_outlet.parameters.items():
            assert outlet.parameters[k] == pytest.approx(v, rel=EPSREL * 2), AssertionError(f"{k=} {v=}")

    @pytest.mark.parametrize(
        "inlets",
        [
            (air, exhaust),
        ],
    )
    @pytest.mark.benchmark
    def test_joiner_predict(self, benchmark, joiner, inlets):
        def benchfunc(joiner, inlets):
            joiner.predict({}, *inlets)

        benchmark(benchfunc, joiner, inlets)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s", "-x"])
