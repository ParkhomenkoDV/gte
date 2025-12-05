import random
from typing import Tuple

import pytest

try:
    from .characteristic import Characteristic
    from .config import parameters as gtep
except ImportError:
    from characteristic import Characteristic
    from config import parameters as gtep


@pytest.fixture
def data(n: int = 100) -> Tuple:
    return tuple(
        {
            gtep.rf: 1.0 + random.random() / 10,
            gtep.mf: 1.0 + random.random() / 10,
            gtep.effeff: 0.9 + random.random() / 100,
            gtep.pipi: random.random(),
        }
        for _ in range(n)
    )


class TestCharacteristic:
    """Тесты для класса Characteristic"""

    def test_init(self, data):
        """Тест инициализации характеристики"""
        characteristic = Characteristic(data)
        assert isinstance(characteristic, Characteristic)
        assert hasattr(Characteristic, gtep.rf)
        assert hasattr(Characteristic, gtep.mf)
        assert hasattr(Characteristic, gtep.effeff)
        assert hasattr(Characteristic, gtep.pipi)

    @pytest.mark.benchmark
    def test_characteristic_init(self, benchmark, data):
        def benchfunc():
            return Characteristic(data)

        benchmark(benchfunc)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s", "-x"])
