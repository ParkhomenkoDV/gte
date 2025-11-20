import numpy as np
import pytest

try:
    from .utils import call_with_kwargs, integral_average, integrate
except ImportError:
    from utils import call_with_kwargs, integral_average, integrate


class Test_call_with_kwargs:
    """Тесты функции call_with_kwargs"""

    @staticmethod
    def f0():
        return "success"

    @staticmethod
    def f(name: str, age: int, active: bool = True):
        return f"{name} is {age} years old and active: {active}"

    class TestClass:
        def method(self, x, y):
            return x**y

    @pytest.mark.parametrize(
        "function, kwargs, expected",
        [
            (f0.__func__, {}, "success"),
            (f0.__func__, {"name": "Tom", "age": 21, "active": True}, "success"),
            (f.__func__, {"name": "Tom", "age": 21}, "Tom is 21 years old and active: True"),
            (f.__func__, {"name": "Bob", "age": 24, "active": False, "sex": "male"}, "Bob is 24 years old and active: False"),
            (lambda x, y: x**y, {"x": 2, "y": 3}, 8),
            (TestClass().method, {"x": 2, "y": 3}, 8),
        ],
    )
    def test_function(self, function, kwargs, expected):
        assert call_with_kwargs(function, kwargs) == expected

    @pytest.mark.benchmark
    def test_call_with_kwargs(self, benchmark):
        """Бенчмарк функции call_with_kwargs"""

        def benchfunc(kwargs):
            call_with_kwargs(self.f, kwargs)

        benchmark(benchfunc, {"name": "Bob", "age": 24, "active": False, "sex": "male"})

    def test_error(self):
        # has not function
        with pytest.raises((AssertionError, TypeError)):
            call_with_kwargs(2, 3)

        # type(kwargs) is not dict
        with pytest.raises((AssertionError, TypeError)):
            call_with_kwargs(self.f, (2, 3))

        # нехватка обязательных аргументов
        with pytest.raises((AssertionError, ValueError)):
            call_with_kwargs(self.f, {"name": "Alice"})


class TestIntegrate:
    """Тесты для функции integrate()"""

    ABSERR = 1e-2
    RELERR = 0.000_1

    @staticmethod
    def f0_1(T):
        return 1000.0

    @staticmethod
    def f1_1(T):
        return 800.0 + 0.5 * T

    @staticmethod
    def f2_1(T):
        return 900.0 + 0.1 * T + 0.001 * T**2

    @staticmethod
    def f1_2(T, P):
        return 1000.0 + 0.2 * T + 0.0001 * P

    @staticmethod
    def f1_3(T, P, V):
        return 950.0 + 0.15 * T + 0.0002 * P + 0.001 * V

    @pytest.mark.parametrize(
        "function, kwargs, expected",
        [
            (f0_1.__func__, {"T": (300, 500)}, 1000.0 * (500 - 300)),
            (f1_1.__func__, {"T": (300, 500)}, 800 * (500 - 300) + 0.5 * (500**2 - 300**2) / 2),
            (f2_1.__func__, {"T": (300, 500)}, 900 * (500 - 300) + 0.1 * (500**2 - 300**2) / 2 + 0.001 * (500**3 - 300**3) / 3),
            (f1_2.__func__, {"T": (300, 500), "P": (100_000, 200_000)}, 1000 * (500 - 300) * (200000 - 100000) + 0.2 * (500**2 - 300**2) / 2 * (200000 - 100000) + 0.0001 * (500 - 300) * (200000**2 - 100000**2) / 2),
            (f1_3.__func__, {"T": (300, 400), "P": (100_000, 200_000), "V": (0.3, 0.4)}, 1_032_500_350),
        ],
    )
    def test_function(self, function, kwargs, expected):
        result, error = integrate(function, **kwargs)
        assert result == pytest.approx(expected, rel=self.RELERR)
        assert error < self.ABSERR

    @pytest.mark.benchmark
    def test_integrate(self, benchmark):
        """Бенчмарк для функции integrate()"""

        def f(x):
            return ((np.sin(x) ** (2 * x)) / np.exp(np.arccos(x**0.5))) ** np.log(abs(x) + 1)

        def benchfunc(ranges):
            integrate(f, **ranges)

        benchmark(benchfunc, {"x": (-100, 100)})

    def test_reverse_integration_limits(self):
        """Тест с обратными пределами интегрирования"""
        result_forward, error_forward = integrate(self.f1_1, T=(300, 500))
        result_reverse, error_reverse = integrate(self.f1_1, T=(500, 300))
        assert abs(result_forward + result_reverse) < 1e-10
        assert error_forward < self.ABSERR
        assert error_reverse < self.ABSERR

    def test_zero_range(self):
        """Тест с нулевым диапазоном"""
        result, error = integrate(self.f0_1, T=(400, 400))
        assert result == pytest.approx(1000.0, rel=self.RELERR)
        assert error < self.ABSERR

    def test_error_cases(self):
        """Тесты обработки ошибок"""

        # Не callable объект
        with pytest.raises((AssertionError, TypeError)):
            integrate("not a function", T=(300, 500))

        # Отсутствующий аргумент
        with pytest.raises((AssertionError, Exception), match="require argument"):
            integrate(self.f1_1)  # T отсутствует

        # Неправильный тип диапазона
        with pytest.raises((AssertionError, TypeError)):
            integrate(self.f1_1, T="not a range")

        # Диапазон неправильной длины
        with pytest.raises((AssertionError, ValueError)):
            integrate(self.f1_1, T=(300, 400, 500))

        # Лишние аргументы (должны игнорироваться)
        result, error1 = integrate(self.f1_1, T=(300, 500), extra_arg=(100, 200))
        expected, error2 = integrate(self.f1_1, T=(300, 500))
        assert result == pytest.approx(expected, rel=self.RELERR)
        assert error1 == error2

    def test_lambda_functions(self):
        """Тест с lambda-функциями"""
        result, error = integrate(lambda T: 1000.0, T=(300, 500))
        assert result == pytest.approx(1000.0 * 200, rel=self.RELERR)
        assert error < self.ABSERR

        result, error = integrate(lambda T, P: 1000.0 + 0.1 * T, T=(300, 500), P=(100_000, 200_000))
        assert result == pytest.approx((1000.0 + 0.1 * 400) * 200 * 100000, rel=self.RELERR)
        assert error < self.ABSERR

    def test_edge_cases(self):
        """Тесты граничных случаев"""

        # Очень большой диапазон
        result, error = integrate(self.f0_1, T=(0, 10_000))
        assert result == pytest.approx(1000.0 * 10_000, rel=1e-6)
        assert error < self.ABSERR

        # Отрицательный диапазон
        result, error = integrate(self.f1_1, T=(-100, 100))
        expected = 800.0 * 200 + 0.5 * (10_000 - 10_000) / 2  # T² симметрично
        assert result == pytest.approx(expected, rel=1e-6)
        assert error < self.ABSERR

        # Очень маленький диапазон
        result, error = integrate(self.f0_1, T=(400, 400 + 1e-10))
        assert abs(result) < self.ABSERR
        assert error < self.ABSERR


class TestIntegralAverage:
    """Тесты для функции integral_average()"""

    ABSERR = 1e-9
    RELERR = 0.000_1

    @staticmethod
    def f0_1(x):
        return 5.0

    @staticmethod
    def f0_2(x, y):
        return 5.0

    @staticmethod
    def f0_3(x, y, z):
        return 5.0

    @staticmethod
    def f1_1(x):
        return x

    @staticmethod
    def f1_2(x, y):
        # x**2 / 2 * y + y**2 / 2 * x
        return x + y

    @staticmethod
    def f1_3(x, y, z):
        # (x**2 + y**2 + z**2) / 2
        return x + y + z

    @staticmethod
    def f2_1(x):
        # 2/3 * x**3
        return x**2

    @staticmethod
    def exp(x):
        return np.exp(x)

    @staticmethod
    def sin(x):
        return np.sin(x)

    @staticmethod
    def discontinuous(x):
        return 1.0 if x < 0.5 else 2.0

    @pytest.mark.parametrize(
        "function, kwargs, expected",
        [
            # f0_1
            (lambda x: 5.0, {"x": (0, 1)}, 5.0),
            (f0_1.__func__, {"x": (0, 1)}, 5.0),
            (f0_1.__func__, {"x": (1, 1)}, 5.0),
            # f0_2
            (lambda x, y: 5.0, {"x": (0, 1), "y": (0, 2)}, 5.0),
            (f0_2.__func__, {"x": (0, 1), "y": (0, 2)}, 5.0),
            (f0_2.__func__, {"x": (0, 1), "y": (2, 2)}, 5.0),
            (f0_2.__func__, {"x": (1, 1), "y": (0, 2)}, 5.0),
            (f0_2.__func__, {"x": (1, 1), "y": (2, 2)}, 5.0),
            # f0_3
            (lambda x, y, z: 5.0, {"x": (0, 1), "y": (0, 2), "z": (0, 3)}, 5.0),
            (f0_3.__func__, {"x": (0, 1), "y": (0, 2), "z": (0, 3)}, 5.0),
            (f0_3.__func__, {"x": (1, 1), "y": (0, 2), "z": (0, 3)}, 5.0),
            (f0_3.__func__, {"x": (0, 1), "y": (2, 2), "z": (0, 3)}, 5.0),
            (f0_3.__func__, {"x": (0, 1), "y": (0, 2), "z": (3, 3)}, 5.0),
            (f0_3.__func__, {"x": (1, 1), "y": (2, 2), "z": (0, 3)}, 5.0),
            (f0_3.__func__, {"x": (1, 1), "y": (0, 2), "z": (3, 3)}, 5.0),
            (f0_3.__func__, {"x": (0, 1), "y": (2, 2), "z": (3, 3)}, 5.0),
            (f0_3.__func__, {"x": (1, 1), "y": (2, 2), "z": (3, 3)}, 5.0),
            # f1_1
            (lambda x: x, {"x": (0, 2)}, (2**2) / 2 / (2 - 0)),
            (f1_1.__func__, {"x": (0, 2)}, (2**2) / 2 / (2 - 0)),
            (f1_1.__func__, {"x": (2, 2)}, 2.0),
            # f1_2
            (lambda x, y: x + y, {"x": (0, 1), "y": (0, 2)}, ((1**2 - 0**2) / 2 * (2 - 0) + (2**2 - 0**2) / 2 * (1 - 0)) / (1 - 0) / (2 - 0)),
            (f1_2.__func__, {"x": (0, 1), "y": (0, 2)}, ((1**2 - 0**2) / 2 * (2 - 0) + (2**2 - 0**2) / 2 * (1 - 0)) / (1 - 0) / (2 - 0)),
            (f1_2.__func__, {"x": (1, 1), "y": (0, 2)}, (1 * (2 - 0) + (2**2 - 0**2) / 2) / (2 - 0)),
            (f1_2.__func__, {"x": (0, 1), "y": (2, 2)}, (2 * (1 - 0) + (1**2 - 0**2) / 2) / (1 - 0)),
            (f1_2.__func__, {"x": (1, 1), "y": (2, 2)}, 1 + 2),
            # f1_3
            (lambda x, y, z: x + y + z, {"x": (0, 1), "y": (0, 2), "z": (0, 3)}, ((1**2 / 2 * (2 - 0) * (3 - 0)) + (2**2 / 2 * (1 - 0) * (3 - 0)) + (3**2 / 2 * (1 - 0) * (2 - 0))) / (1 - 0) / (2 - 0) / (3 - 0)),
            (f1_3.__func__, {"x": (0, 1), "y": (0, 2), "z": (0, 3)}, ((1**2 / 2 * (2 - 0) * (3 - 0)) + (2**2 / 2 * (1 - 0) * (3 - 0)) + (3**2 / 2 * (1 - 0) * (2 - 0))) / (1 - 0) / (2 - 0) / (3 - 0)),
            (f1_3.__func__, {"x": (1, 1), "y": (0, 2), "z": (0, 3)}, 1 + (0 + 2) / 2 + (0 + 3) / 2),
            (f1_3.__func__, {"x": (0, 1), "y": (2, 2), "z": (0, 3)}, (0 + 1) / 2 + 2 + (0 + 3) / 2),
            (f1_3.__func__, {"x": (0, 1), "y": (0, 2), "z": (3, 3)}, (0 + 1) / 2 + (0 + 2) / 2 + 3),
            (f1_3.__func__, {"x": (0, 1), "y": (2, 2), "z": (3, 3)}, (0 + 1) / 2 + 2 + 3),
            (f1_3.__func__, {"x": (1, 1), "y": (2, 2), "z": (0, 3)}, 1 + 2 + (0 + 3) / 2),
            (f1_3.__func__, {"x": (1, 1), "y": (0, 2), "z": (3, 3)}, 1 + (0 + 2) / 2 + 3),
            (f1_3.__func__, {"x": (1, 1), "y": (2, 2), "z": (3, 3)}, 1 + 2 + 3),
            # f2_1
            (lambda x: x**2, {"x": (0, 1)}, (1**3 - 0**3) / 3),
            (f2_1.__func__, {"x": (0, 1)}, (1**3 - 0**3) / 3),
            (f2_1.__func__, {"x": (1, 1)}, 1),
            # f2_2
            # f2_3
            # exp
            (lambda x: np.exp(x), {"x": (0, 1)}, (np.exp(1) - 1) / (1 - 0)),
            (exp.__func__, {"x": (0, 1)}, (np.exp(1) - 1) / (1 - 0)),
            (exp.__func__, {"x": (1, 1)}, np.exp(1)),
            # sin
            (lambda x: np.sin(x), {"x": (0, np.pi)}, (1 - np.cos(np.pi)) / (np.pi - 0)),
            (sin.__func__, {"x": (0, np.pi)}, (1 - np.cos(np.pi)) / (np.pi - 0)),
            (sin.__func__, {"x": (np.pi, np.pi)}, np.sin(np.pi)),
            # discontinuous
            (lambda x: 1.0 if x < 0.5 else 2.0, {"x": (0, 1)}, 1.5),
            (discontinuous.__func__, {"x": (0, 1)}, 1.5),
        ],
    )
    def test_function(self, function, kwargs, expected):
        result, error = integral_average(function, **kwargs)
        assert result == pytest.approx(expected, rel=self.RELERR)
        assert error < self.ABSERR

    @pytest.mark.benchmark
    def test_integral_average(self, benchmark):
        """Бенчмарк для функции integral_average()"""

        def f(x):
            return (np.sin(x) / np.exp(x)) ** np.log(abs(x) + 1)

        def benchfunc(function, kwargs):
            integral_average(function, **kwargs)

        benchmark(benchfunc, f, {"x": (-100, 100)})

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


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s", "-x"])
