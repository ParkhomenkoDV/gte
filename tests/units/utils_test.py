import numpy as np
import pytest

from src.utils import call_with_kwargs, integral_average


class Test_call_with_kwargs:
    @staticmethod
    def f(a, b):
        return a + b

    def test_types(self):
        assert call_with_kwargs(self.f, {"a": 2, "b": 3}) == 5
        assert call_with_kwargs(self.f, {"a": 2, "b": 3, "c": 4}) == 5

        # has not function
        with pytest.raises((AssertionError, TypeError)):
            call_with_kwargs(2, 3)

        # type(kwargs) is not dict
        with pytest.raises((AssertionError, TypeError)):
            call_with_kwargs(self.f, (2, 3))

        # нехватка аргументов
        with pytest.raises((AssertionError, TypeError)):
            call_with_kwargs(self.f, {"a": 2, "c": 4})


class TestIntegralAverage:
    """Тесты для функции integral_average"""

    ABSERR = 1e-9
    RELERR = 0.000_1

    @staticmethod
    def f0_1(x):
        return 5.0

    @pytest.mark.parametrize(
        "kwargs, expected",
        [
            ({"x": (0, 1)}, 5.0),
            ({"x": (1, 1)}, 5.0),
        ],
    )
    def test_f0_1(self, kwargs, expected):
        result, error = integral_average(self.f0_1, **kwargs)
        assert result == pytest.approx(expected, rel=self.RELERR)
        assert error < self.ABSERR

    @staticmethod
    def f0_2(x, y):
        return 5.0

    @pytest.mark.parametrize(
        "kwargs, expected",
        [
            ({"x": (0, 1), "y": (0, 2)}, 5.0),
            ({"x": (0, 1), "y": (2, 2)}, 5.0),
            ({"x": (1, 1), "y": (0, 2)}, 5.0),
            ({"x": (1, 1), "y": (2, 2)}, 5.0),
        ],
    )
    def test_f0_2(self, kwargs, expected):
        result, error = integral_average(self.f0_2, **kwargs)
        assert result == pytest.approx(expected, rel=self.RELERR)
        assert error < self.ABSERR

    @staticmethod
    def f0_3(x, y, z):
        return 5.0

    @pytest.mark.parametrize(
        "kwargs, expected",
        [
            ({"x": (0, 1), "y": (0, 2), "z": (0, 3)}, 5.0),
            ({"x": (1, 1), "y": (0, 2), "z": (0, 3)}, 5.0),
            ({"x": (0, 1), "y": (2, 2), "z": (0, 3)}, 5.0),
            ({"x": (0, 1), "y": (0, 2), "z": (3, 3)}, 5.0),
            ({"x": (1, 1), "y": (2, 2), "z": (0, 3)}, 5.0),
            ({"x": (1, 1), "y": (0, 2), "z": (3, 3)}, 5.0),
            ({"x": (0, 1), "y": (2, 2), "z": (3, 3)}, 5.0),
            ({"x": (1, 1), "y": (2, 2), "z": (3, 3)}, 5.0),
        ],
    )
    def test_f0_3(self, kwargs, expected):
        result, error = integral_average(self.f0_3, **kwargs)
        assert result == pytest.approx(expected, rel=self.RELERR)
        assert error < self.ABSERR

    @staticmethod
    def f1_1(x):
        return x

    @pytest.mark.parametrize(
        "kwargs, expected",
        [
            ({"x": (0, 2)}, (2**2) / 2 / (2 - 0)),
            ({"x": (2, 2)}, 2.0),
        ],
    )
    def test_f1_1(self, kwargs, expected):
        result, error = integral_average(self.f1_1, **kwargs)
        assert result == pytest.approx(expected, rel=self.RELERR)
        assert error < self.ABSERR

    @staticmethod
    def f1_2(x, y):
        # x**2/2 * y + y**2/2 * x
        return x + y

    @pytest.mark.parametrize(
        "kwargs, expected",
        [
            (
                {"x": (0, 1), "y": (0, 2)},
                ((1**2 - 0**2) / 2 * (2 - 0) + (2**2 - 0**2) / 2 * (1 - 0))
                / (1 - 0)
                / (2 - 0),
            ),
            (
                {"x": (1, 1), "y": (0, 2)},
                (2 + (2**2 - 0**2) / 2) / (2 - 0),
            ),
            (
                {"x": (0, 1), "y": (2, 2)},
                (2 + (1**2 / 2)) / (1 - 0),
            ),
            (
                {"x": (1, 1), "y": (2, 2)},
                1 + 2,
            ),
        ],
    )
    def test_f1_2(self, kwargs, expected):
        result, error = integral_average(self.f1_2, **kwargs)
        assert result == pytest.approx(expected, rel=self.RELERR)
        assert error < self.ABSERR

    @staticmethod
    def f1_3(x, y, z):
        return x + y + z

    @pytest.mark.parametrize(
        "kwargs, expected",
        [
            (
                {"x": (0, 1), "y": (0, 2), "z": (0, 3)},
                (
                    (1**2 / 2 * (2 - 0) * (3 - 0))
                    + (2**2 / 2 * (1 - 0) * (3 - 0))
                    + (3**2 / 2 * (1 - 0) * (2 - 0))
                )
                / (1 - 0)
                / (2 - 0)
                / (3 - 0),
            )
        ],
    )
    def test_f1_3(self, kwargs, expected):
        result, error = integral_average(self.f1_3, **kwargs)
        assert result == pytest.approx(expected, rel=self.RELERR)
        assert error < self.ABSERR

    @staticmethod
    def f2_1(x):
        return x**2

    @pytest.mark.parametrize(
        "kwargs, expected",
        [
            ({"x": (0, 1)}, 1**3 / 3),
            ({"x": (1, 1)}, 1),
        ],
    )
    def test_f2_1(self, kwargs, expected):
        result, error = integral_average(self.f2_1, **kwargs)
        assert result == pytest.approx(expected, rel=self.RELERR)
        assert error < self.ABSERR

    @staticmethod
    def exp(x):
        return np.exp(x)

    @pytest.mark.parametrize(
        "kwargs, expected",
        [
            ({"x": (0, 1)}, (np.exp(1) - 1) / (1 - 0)),
        ],
    )
    def test_exp(self, kwargs, expected):
        result, error = integral_average(self.exp, **kwargs)
        assert result == pytest.approx(expected, rel=self.RELERR)
        assert error < self.ABSERR

    @staticmethod
    def sin(x):
        return np.sin(x)

    @pytest.mark.parametrize(
        "kwargs, expected",
        [
            ({"x": (0, np.pi)}, (1 - np.cos(np.pi)) / (np.pi - 0)),
        ],
    )
    def test_trigonometric_function(self, kwargs, expected):
        result, error = integral_average(self.sin, **kwargs)
        assert result == pytest.approx(expected, rel=self.RELERR)
        assert error < self.ABSERR

    @staticmethod
    def discontinuous(x):
        return 1.0 if x < 0.5 else 2.0

    @pytest.mark.parametrize(
        "kwargs, expected",
        [
            ({"x": (0, 1)}, 1.5),
        ],
    )
    def test_edge_cases(self, kwargs, expected):
        result, error = integral_average(self.discontinuous, x=(0, 1))
        assert result == pytest.approx(expected, rel=self.RELERR)
        assert error < self.ABSERR

    def test_error_cases(self):
        """Тесты обработки ошибок"""

        def test_func(x):
            return x

        # Не callable объект
        with pytest.raises((AssertionError, TypeError)):
            integral_average("not a function", x=(0, 1))

        # Отсутствующий аргумент
        with pytest.raises((AssertionError, Exception), match="require argument"):
            integral_average(test_func, y=(0, 1))  # x отсутствует

        # Неправильный тип диапазона
        with pytest.raises((AssertionError, TypeError)):
            integral_average(test_func, x="not a range")

        # Диапазон неправильной длины
        with pytest.raises((AssertionError, ValueError)):
            integral_average(test_func, x=(0, 1, 2))
