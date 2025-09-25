from collections import namedtuple

import numpy as np
import pytest
from numpy import isnan, nan
from substance import Substance
from thermodynamics import gas_const, heat_capacity_at_constant_pressure

from src.config import EPSREL
from src.config import parameters as gtep
from src.nodes import CombustionChamber, Compressor


@pytest.fixture
def substance_air():
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


@pytest.fixture
def substance_fuel():
    """Горючее для тестов"""
    return Substance(
        "kerosene",
        parameters={
            gtep.gc: 287.0,
            gtep.TT: 300.0,
            gtep.PP: 101325.0,
            gtep.mf: 1.0,
            gtep.Cp: 1006.0,
            gtep.k: 1.4,
            gtep.c: 0.0,
        },
        functions={
            gtep.gc: lambda total_temperature: 287.0,
            gtep.Cp: lambda total_temperature, total_pressure: 1006.0,
        },
    )


@pytest.fixture
def substance_exhaust():
    """Выхлопные газы для тестов"""
    return Substance(
        "exhaust",
        parameters={
            gtep.gc: 287.0,
            gtep.TT: 300.0,
            gtep.PP: 101325.0,
            gtep.mf: 1.0,
            gtep.Cp: 1006.0,
            gtep.k: 1.4,
            gtep.c: 0.0,
        },
        functions={
            gtep.gc: lambda excess_oxidizing: gas_const("air", excess_oxidizing, fuel="kerosene"),
            gtep.Cp: lambda total_temperature, excess_oxidizing: heat_capacity_at_constant_pressure("exhaust", total_temperature, excess_oxidizing),
        },
    )


class TestNode:
    @pytest.mark.parametrize("Node", [Compressor, CombustionChamber])
    def test_init(self, Node):
        node = Node()
        assert isinstance(node, Node)
        assert getattr(node, "mass_flow_leak") == 0
        assert hasattr(node, "inlet")
        assert hasattr(node, "outlet")
        assert isinstance(node.inlet, Substance)
        assert isinstance(node.outlet, Substance)

    @pytest.mark.parametrize("Node", [Compressor, CombustionChamber])
    def test_name(self, Node):
        node = Node()
        assert node.name == Node.__name__
        node.name = "C"
        assert node.name == "C"
        node = Node("123")
        assert node.name == "123"


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
        assert compressor.mass_flow_leak == 0.0

    @pytest.mark.parametrize(
        "pipi,effeff,power,expected_count",
        [
            (nan, 0.8, nan, 2),  # недоопределено
            (6.0, 0.8, nan, 1),  # определяется power
            (6.0, nan, 1e6, 1),  # определяется effeff
            (nan, 0.8, 1e6, 1),  # определяется pipi
            (6.0, 0.8, 1e6, 0),  # все определено
        ],
    )
    def test_variables_nan_count(self, compressor, pipi, effeff, power, expected_count):
        """Тест подсчета неопределенных переменных"""
        setattr(compressor, gtep.pipi, pipi)
        setattr(compressor, gtep.effeff, effeff)
        setattr(compressor, gtep.power, power)

        count = sum(1 if np.isnan(v) else 0 for v in compressor.variables.values())
        assert count == expected_count

    expected = namedtuple("expected", ["key", "value"])

    @pytest.mark.parametrize(
        "pipi, effeff, power, mf_leak, error, expected",
        [
            (6.0, 0.85, nan, 0, False, expected(gtep.power, 23_632_485)),
            (6.0, nan, 24 * 10**6, 0, False, expected(gtep.effeff, 0.8369)),
            (nan, 0.85, 24 * 10**6, 0, False, expected(gtep.pipi, 6.1329)),
            # error
            (nan, nan, 24 * 10**6, 0, True, expected("", 0)),
            (6.0, nan, nan, 0, True, expected("", 0)),
            (nan, nan, nan * 10**6, 0, True, expected("", 0)),
        ],
    )
    def test_calculate_integration(self, compressor, substance_air, pipi, effeff, power, mf_leak, error, expected):
        """Интеграционный тест расчета"""
        setattr(compressor, gtep.pipi, pipi)
        setattr(compressor, gtep.effeff, effeff)
        setattr(compressor, gtep.power, power)
        compressor.mass_flow_leak = mf_leak

        if error:
            with pytest.raises(Exception):
                compressor.calculate(substance_air)
        else:
            outlet = compressor.calculate(substance_air)

            assert isinstance(outlet, Substance)
            # Проверяем основные параметры
            assert gtep.mf in outlet.parameters
            assert gtep.TT in outlet.parameters
            assert gtep.PP in outlet.parameters

            assert getattr(compressor, expected.key) == pytest.approx(expected.value, rel=EPSREL)

    @pytest.mark.parametrize(
        "effeff, TT, expected",
        [
            (0.75, 300, True),  # нормальный КПД
            (0.0, 300, True),  # нулевой КПД
            (1.5, 300, False),  # КПД > 1
            (-0.1, 300, False),  # отрицательный КПД
        ],
    )
    def test_is_real(self, compressor, effeff, TT, expected):
        """Тест на реальность"""
        setattr(compressor, gtep.effeff, effeff)
        compressor.outlet.parameters[gtep.TT] = TT

        assert compressor.is_real == expected

    def test_error(self, compressor):
        """Тест обработки ошибок"""
        # Тест с некорректным веществом
        with pytest.raises(Exception):
            compressor.calculate(None)

        # Тест с веществом без необходимых параметров
        invalid_substance = Substance("invalid", parameters={}, functions={})
        with pytest.raises(Exception):
            compressor.calculate(invalid_substance)


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
        assert getattr(cc, "total_pressure_loss") == 0
        assert cc.mass_flow_leak == 0.0

        assert hasattr(cc, "fuel")
        assert isinstance(cc.fuel, Substance)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
