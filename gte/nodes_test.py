from collections import namedtuple

import numpy as np
import pytest
from combustion_chambler import CombustionChamber
from compressor import Compressor
from config import EPSREL
from config import parameters as gtep
from fixtures import air, exhaust, kerosene
from numpy import isnan, nan
from substance import Substance
from turbine import Turbine


class TestNode:
    """Тесты абстрактного класса GTENode"""

    @pytest.mark.parametrize("Node", [Compressor, CombustionChamber, Turbine])
    def test_init(self, Node):
        node = Node()
        assert isinstance(node, Node)
        assert isinstance(node.name, str)
        assert getattr(node, "leak") == 0

    @pytest.mark.parametrize("Node", [Compressor, CombustionChamber, Turbine])
    def test_name(self, Node):
        node = Node()
        assert node.name == Node.__name__  # default
        node.name = "C"
        assert node.name == "C"  # reset
        node = Node("123")
        assert node.name == "123"  # __init__

    @pytest.mark.parametrize("Node", [Compressor, CombustionChamber, Turbine])
    def test_leak(self, Node):
        node = Node()
        assert node.leak == 0  # default
        node.leak = 0.5
        assert node.leak == 0.5  # reset


@pytest.fixture
def compressor():
    """Создает экземпляр компрессора"""
    return Compressor("TestCompressor")


class TestCompressor:
    """Тесты для класса Compressor"""

    def test_init(self, compressor):
        """Тест инициализации компрессора"""
        assert compressor.name == "TestCompressor"
        assert isnan(getattr(compressor, gtep.pipi))
        assert isnan(getattr(compressor, gtep.effeff))
        assert isnan(getattr(compressor, gtep.power))
        assert compressor.leak == 0.0

    @pytest.mark.benchmark
    def test_compressor_init(self, benchmark):
        """Бенчмарк инициализации компрессора"""

        def benchfunc():
            return Compressor()

        benchmark(benchfunc)

    @pytest.mark.parametrize(
        "inlet_parameters, parameters, use_ml, expected",
        [
            # pipi
            ({gtep.mf: 100, gtep.TT: 300, gtep.PP: 101_325}, {gtep.power: 24 * 10**6, gtep.effeff: 0.85}, False, {"outlet_total_pressure": 101325, "outlet_total_temperature": 300, gtep.pipi: 6}),
            ({gtep.mf: 100, gtep.TT: 300, gtep.PP: 101_325}, {gtep.power: 24 * 10**6, gtep.effeff: 0.85}, True, {"outlet_total_pressure": 110000, "outlet_total_temperature": 1076, gtep.pipi: 1.1}),
            # effeff
            ({gtep.mf: 100, gtep.TT: 300, gtep.PP: 101_325}, {gtep.pipi: 6.0, gtep.power: 24 * 10**6}, False, {"outlet_total_pressure": 101325, "outlet_total_temperature": 300, gtep.effeff: 1}),
            ({gtep.mf: 100, gtep.TT: 300, gtep.PP: 101_325}, {gtep.pipi: 6.0, gtep.power: 24 * 10**6}, True, {"outlet_total_pressure": 1210000, "outlet_total_temperature": 1954, gtep.effeff: 0.85}),
            # power
            ({gtep.mf: 100, gtep.TT: 300, gtep.PP: 101_325}, {gtep.pipi: 6.0, gtep.effeff: 0.85}, False, {"outlet_total_pressure": 101325, "outlet_total_temperature": 300, gtep.power: 30_180_000}),
            ({gtep.mf: 100, gtep.TT: 300, gtep.PP: 101_325}, {gtep.pipi: 6.0, gtep.effeff: 0.85}, True, {"outlet_total_pressure": 110000, "outlet_total_temperature": 1076, gtep.power: 3_124_532}),
        ],
    )
    def test_predict(self, compressor, inlet_parameters, parameters, use_ml, expected):
        """Тест предсказания компресоора"""
        inlet = Substance("inlet")
        inlet.parameters = inlet_parameters
        for param, value in parameters.items():
            setattr(compressor, param, value)
        prediction = compressor.predict(inlet, use_ml=use_ml)
        for k, v in prediction.items():
            assert v == pytest.approx(expected[k], rel=0.2), AssertionError(f"{k=} {v=}")

    @pytest.mark.parametrize(
        "inlet_parameters, parameters, use_ml",
        [
            # pipi
            ({gtep.mf: 100, gtep.TT: 300, gtep.PP: 101_325}, {gtep.power: 24 * 10**6, gtep.effeff: 0.85}, False),
            ({gtep.mf: 100, gtep.TT: 300, gtep.PP: 101_325}, {gtep.power: 24 * 10**6, gtep.effeff: 0.85}, True),
            # effeff
            ({gtep.mf: 100, gtep.TT: 300, gtep.PP: 101_325}, {gtep.pipi: 6.0, gtep.power: 24 * 10**6}, False),
            ({gtep.mf: 100, gtep.TT: 300, gtep.PP: 101_325}, {gtep.pipi: 6.0, gtep.power: 24 * 10**6}, True),
            # power
            ({gtep.mf: 100, gtep.TT: 300, gtep.PP: 101_325}, {gtep.pipi: 6.0, gtep.effeff: 0.85}, False),
            ({gtep.mf: 100, gtep.TT: 300, gtep.PP: 101_325}, {gtep.pipi: 6.0, gtep.effeff: 0.85}, True),
        ],
    )
    @pytest.mark.benchmark
    def test_compressor_predict(self, benchmark, compressor, inlet_parameters, parameters, use_ml):
        def benchfunc(compressor, inlet_parameters, parameters, use_ml):
            inlet = Substance("inlet")
            inlet.parameters = inlet_parameters
            for param, value in parameters.items():
                setattr(compressor, param, value)
            compressor.predict(inlet, use_ml=use_ml)

        benchmark(benchfunc, compressor, inlet_parameters, parameters, use_ml)

    expected = namedtuple("expected", ["key", "value"])

    @pytest.mark.parametrize(
        "pipi, effeff, power, leak, error, expected",
        [
            (6.0, 0.85, nan, 0, False, expected(gtep.power, 11_816_242)),
            (6.0, nan, 12 * 10**6, 0, False, expected(gtep.effeff, 0.8369)),
            (nan, 0.85, 12 * 10**6, 0, False, expected(gtep.pipi, 6.1329)),
            # error
            (nan, nan, 24 * 10**6, 0, True, expected("", 0)),
            (6.0, nan, nan, 0, True, expected("", 0)),
            (nan, nan, nan, 0, True, expected("", 0)),
        ],
    )
    def test_solve(self, compressor, pipi, effeff, power, leak, error, expected):
        """Интеграционный тест расчета"""
        setattr(compressor, gtep.pipi, pipi)
        setattr(compressor, gtep.effeff, effeff)
        setattr(compressor, gtep.power, power)
        compressor.leak = leak

        if error:
            with pytest.raises(Exception):
                compressor.solve(air)
        else:
            outlet = compressor.solve(air)

            assert isinstance(outlet, Substance)
            # Проверяем основные параметры
            assert gtep.mf in outlet.parameters
            assert gtep.TT in outlet.parameters
            assert gtep.PP in outlet.parameters

            assert getattr(compressor, expected.key) == pytest.approx(expected.value, rel=EPSREL)

    @pytest.mark.parametrize(
        "kwargs",
        [
            {gtep.pipi: 6.0, gtep.effeff: 0.85},
            {gtep.effeff: 0.85, gtep.power: 24 * 10**6},
            {gtep.pipi: 6.0, gtep.power: 24 * 10**6},
        ],
    )
    @pytest.mark.benchmark
    def test_compressor_solve(self, benchmark, compressor, kwargs):
        def benchfunc(kwargs):
            for k, v in kwargs.items():
                setattr(compressor, k, v)

            compressor.solve(air)

            for k in compressor.variables:
                setattr(compressor, k, np.nan)

        benchmark(benchfunc, kwargs)

    @pytest.mark.parametrize(
        "effeff, TT, pipi, expected",
        [
            (0.75, 300, 6, True),  # нормальный КПД
            (0.0, 300, 6, True),  # нулевой КПД
            (1.5, 300, 6, False),  # КПД > 1
            (-0.1, 300, 6, False),  # отрицательный КПД
        ],
    )
    def test_is_real(self, compressor, effeff, TT, pipi, expected):
        """Тест на реальность"""
        setattr(compressor, gtep.effeff, effeff)
        setattr(compressor, gtep.pipi, pipi)
        compressor.outlet = Substance("air")
        compressor.outlet.parameters[gtep.TT] = TT

        assert compressor.is_real == expected

    def test_error(self, compressor):
        """Тест обработки ошибок"""
        # Тест с некорректным веществом
        with pytest.raises(Exception):
            compressor.solve(None)

        # Тест с веществом без необходимых параметров
        invalid_substance = Substance("invalid", parameters={}, functions={})
        with pytest.raises(Exception):
            compressor.solve(invalid_substance)


@pytest.fixture
def cc():
    """Создает экземпляр камеры сгорания"""
    return CombustionChamber("TestCombustionChamber")


class TestCombustionChamber:
    """Тесты для класса CombustionChamber"""

    def test_init(self, cc):
        """Тест инициализации камеры сгорания"""
        assert cc.name == "TestCombustionChamber"
        assert isnan(getattr(cc, "efficiency_burn"))
        assert isnan(getattr(cc, gtep.peff))
        assert cc.leak == 0.0

    @pytest.mark.benchmark
    def test_cc_init(self, benchmark):
        """Бенчмарк инициализации камеры сгорания"""

        def benchfunc():
            return CombustionChamber()

        benchmark(benchfunc)

    @pytest.mark.parametrize(
        "peff, efficiency_burn, leak, error, expected",
        [
            (0.95, 0.98, 0, False, {gtep.TT: 2332.05636}),
            # error
            (nan, nan, 0, True, {}),
            (6.0, nan, 0, True, {}),
            (nan, nan, 0, True, {}),
        ],
    )
    def test_solve(self, cc, peff, efficiency_burn, leak, error, expected):
        """Интеграционный тест расчета"""
        setattr(cc, gtep.peff, peff)
        setattr(cc, "efficiency_burn", efficiency_burn)
        cc.leak = leak

        if error:
            with pytest.raises(Exception):
                cc.solve(air, kerosene)
        else:
            outlet = cc.solve(air, kerosene)

            assert isinstance(outlet, Substance)
            # Проверяем основные параметры
            assert gtep.mf in outlet.parameters
            assert gtep.TT in outlet.parameters
            assert gtep.PP in outlet.parameters

            for k, v in expected.items():
                assert outlet.parameters[k] == pytest.approx(v, rel=EPSREL)

    @pytest.mark.parametrize(
        "node, kwargs",
        [
            (CombustionChamber(), {gtep.peff: 0.95, "efficiency_burn": 0.99}),
        ],
    )
    @pytest.mark.benchmark
    def test_cc_solve(self, benchmark, node, kwargs):
        def benchfunc(node, kwargs):
            for k, v in kwargs.items():
                setattr(node, k, v)

            node.solve(air, kerosene)

            for k in node.variables:
                setattr(node, k, np.nan)

        benchmark(benchfunc, node, kwargs)


@pytest.fixture
def turbine():
    """Создает экземпляр турбины"""
    return Turbine("TestTurbine")


class TestTurbine:
    """Тесты для класса Turbine"""

    def test_init(self, turbine):
        """Тест инициализации компрессора"""
        assert turbine.name == "TestTurbine"
        assert isnan(getattr(turbine, gtep.pipi))
        assert isnan(getattr(turbine, gtep.effeff))
        assert isnan(getattr(turbine, gtep.power))
        assert turbine.leak == 0.0

    @pytest.mark.benchmark
    def test_turbine_init(self, benchmark):
        """Бенчмарк инициализации компрессора"""

        def benchfunc():
            return Turbine()

        benchmark(benchfunc)

    expected = namedtuple("expected", ["key", "value"])

    @pytest.mark.parametrize(
        "pipi, effeff, power, leak, error, expected",
        [
            (4.0, 0.9, nan, 0, False, expected(gtep.power, 23_396_982)),
            (4.0, nan, 24 * 10**6, 0, False, expected(gtep.effeff, 0.9232)),
            (nan, 0.9, 24 * 10**6, 0, False, expected(gtep.pipi, 4.175)),
            # error
            (nan, nan, 24 * 10**6, 0, True, expected("", 0)),
            (6.0, nan, nan, 0, True, expected("", 0)),
            (nan, nan, nan, 0, True, expected("", 0)),
        ],
    )
    def test_solve(self, turbine, pipi, effeff, power, leak, error, expected):
        """Интеграционный тест расчета"""
        setattr(turbine, gtep.pipi, pipi)
        setattr(turbine, gtep.effeff, effeff)
        setattr(turbine, gtep.power, power)
        turbine.leak = leak

        if error:
            with pytest.raises(Exception):
                turbine.solve(exhaust)
        else:
            outlet = turbine.solve(exhaust)

            assert isinstance(outlet, Substance)
            # Проверяем основные параметры
            assert gtep.mf in outlet.parameters
            assert gtep.TT in outlet.parameters
            assert gtep.PP in outlet.parameters

            assert getattr(turbine, expected.key) == pytest.approx(expected.value, rel=EPSREL)

    @pytest.mark.parametrize(
        "node, kwargs",
        [
            (Turbine(), {gtep.pipi: 4.0, gtep.effeff: 0.9}),
            (Turbine(), {gtep.effeff: 0.9, gtep.power: 24 * 10**6}),
            (Turbine(), {gtep.pipi: 4.0, gtep.power: 24 * 10**6}),
        ],
    )
    @pytest.mark.benchmark
    def test_turbine_solve(self, benchmark, node, kwargs):
        def benchfunc(node, kwargs):
            for k, v in kwargs.items():
                setattr(node, k, v)

            node.solve(exhaust)

            for k in node.variables:
                setattr(node, k, np.nan)

        benchmark(benchfunc, node, kwargs)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s", "-x"])
