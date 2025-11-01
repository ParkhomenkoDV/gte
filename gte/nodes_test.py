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
        "pipi, effeff, power, expected_count",
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
            (6.0, 0.85, nan, 0, False, expected(gtep.power, 11_816_242)),
            (6.0, nan, 12 * 10**6, 0, False, expected(gtep.effeff, 0.8369)),
            (nan, 0.85, 12 * 10**6, 0, False, expected(gtep.pipi, 6.1329)),
            # error
            (nan, nan, 24 * 10**6, 0, True, expected("", 0)),
            (6.0, nan, nan, 0, True, expected("", 0)),
            (nan, nan, nan, 0, True, expected("", 0)),
        ],
    )
    def test_calculate_integration(self, compressor, pipi, effeff, power, mf_leak, error, expected):
        """Интеграционный тест расчета"""
        setattr(compressor, gtep.pipi, pipi)
        setattr(compressor, gtep.effeff, effeff)
        setattr(compressor, gtep.power, power)
        compressor.mass_flow_leak = mf_leak

        if error:
            with pytest.raises(Exception):
                compressor.calculate(air)
        else:
            outlet = compressor.calculate(air)

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
        assert isnan(getattr(cc, gtep.peff))
        assert cc.mass_flow_leak == 0.0

    @pytest.mark.parametrize(
        "peff, efficiency_burn, mf_leak, error, expected",
        [
            (0.95, 0.98, 0, False, {gtep.TT: 1167}),
            # error
            (nan, nan, 0, True, {}),
            (6.0, nan, 0, True, {}),
            (nan, nan, 0, True, {}),
        ],
    )
    def test_calculate_integration(self, cc, air, fuel, peff, efficiency_burn, mf_leak, error, expected):
        """Интеграционный тест расчета"""
        setattr(cc, gtep.peff, peff)
        setattr(cc, "efficiency_burn", efficiency_burn)
        cc.mass_flow_leak = mf_leak

        if error:
            with pytest.raises(Exception):
                cc.calculate(air, fuel)
        else:
            outlet = cc.calculate(air, fuel)

            assert isinstance(outlet, Substance)
            # Проверяем основные параметры
            assert gtep.mf in outlet.parameters
            assert gtep.TT in outlet.parameters
            assert gtep.PP in outlet.parameters

            for k, v in expected.items():
                assert outlet.parameters[k] == pytest.approx(v, rel=EPSREL)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
