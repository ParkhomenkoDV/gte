import numpy as np
import pytest
from numpy import isnan, nan
from substance import Substance

from src.config import EPSREL
from src.config import parameters as gtep
from src.nodes import Compressor


@pytest.fixture
def substance_air():
    """Базовое вещество для тестов"""
    return Substance(
        "air",
        parameters={
            gtep.gc: 287.0,
            gtep.TT: 300.0,
            gtep.PP: 101325.0,
            gtep.mf: 100.0,
            gtep.Cp: 1006.0,
            gtep.k: 1.4,
            gtep.c: 0.0,
        },
        functions={
            gtep.gc: lambda total_temperature: 287.0,
            gtep.Cp: lambda total_temperature, total_pressure: 1006.0,
        },
    )


class TestNode:
    @pytest.mark.parametrize("Node", [Compressor])
    def test_init(self, Node):
        node = Node()
        assert isinstance(node, Node)

    @pytest.mark.parametrize("Node", [Compressor])
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

    def test_initialization(self, compressor):
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

    @pytest.mark.parametrize(
        "pipi, effeff, power, mf_leak, error, expected",
        [
            (6.0, 0.85, nan, 0, False, [gtep.power, 23_690_849]),
            (6.0, nan, 24 * 10**6, 0, False, [gtep.effeff, 0.84]),
            (nan, 0.85, 24 * 10**6, 0, False, [gtep.pipi, 6.1]),
            # error
            (nan, nan, 24 * 10**6, 0, True, []),
            (6.0, nan, nan, 0, True, []),
            (nan, nan, nan * 10**6, 0, True, []),
        ],
    )
    def test_calculate_integration(
        self, compressor, substance_air, pipi, effeff, power, mf_leak, error, expected
    ):
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

            assert getattr(compressor, expected[0]) == pytest.approx(
                expected[1], rel=EPSREL
            )

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

    def test_error_handling(self, compressor):
        """Тест обработки ошибок"""
        # Тест с некорректным веществом
        with pytest.raises(Exception):
            compressor.calculate(None)

        # Тест с веществом без необходимых параметров
        invalid_substance = Substance("invalid", parameters={}, functions={})
        with pytest.raises(Exception):
            compressor.calculate(invalid_substance)
