import numpy as np
import pytest

from src.utils import call_with_kwargs, enthalpy, integral_average


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


class TestEnthalpy:
    """Тесты для функции enthalpy"""

    @staticmethod
    def cp_constant(T):
        """Постоянная теплоемкость"""
        return 1000.0

    @staticmethod
    def cp_linear(T):
        """Линейная теплоемкость"""
        return 800.0 + 0.5 * T

    @staticmethod
    def cp_quadratic(T):
        """Квадратичная теплоемкость"""
        return 900.0 + 0.1 * T + 0.001 * T**2

    @staticmethod
    def cp_multivariate(T, P):
        """Теплоемкость от температуры и давления"""
        return 1000.0 + 0.2 * T + 0.0001 * P

    @staticmethod
    def cp_three_variables(T, P, V):
        """Теплоемкость от трех переменных"""
        return 950.0 + 0.15 * T + 0.0002 * P + 0.001 * V

    def test_constant_heat_capacity(self):
        """Тест с постоянной теплоемкостью"""
        result = enthalpy(self.cp_constant, T=(300, 500))
        expected = 1000.0 * (500 - 300)  # ∫1000 dT от 300 до 500
        assert abs(result - expected) < 1e-10

    def test_linear_heat_capacity(self):
        """Тест с линейной теплоемкостью"""
        result = enthalpy(self.cp_linear, T=(300, 500))
        # ∫(800 + 0.5T)dT от 300 до 500
        expected = 800 * (500 - 300) + 0.5 * (500**2 - 300**2) / 2
        assert abs(result - expected) < 1e-10

    def test_quadratic_heat_capacity(self):
        """Тест с квадратичной теплоемкостью"""
        result = enthalpy(self.cp_quadratic, T=(300, 500))
        # ∫(900 + 0.1T + 0.001T²)dT от 300 до 500
        expected = 900 * (500 - 300) + 0.1 * (500**2 - 300**2) / 2 + 0.001 * (500**3 - 300**3) / 3
        assert abs(result - expected) < 1e-8

    def test_multivariate_heat_capacity(self):
        """Тест с многомерной теплоемкостью"""
        result = enthalpy(self.cp_multivariate, T=(300, 500), P=(100_000, 200_000))
        # ∫∫(1000 + 0.2T + 0.0001P)dTdP
        expected = 1000 * (500 - 300) * (200000 - 100000) + 0.2 * (500**2 - 300**2) / 2 * (200000 - 100000) + 0.0001 * (500 - 300) * (200000**2 - 100000**2) / 2
        assert abs(result - expected) < 1e-5

    def test_three_variables(self):
        """Тест с тремя переменными"""
        result = enthalpy(self.cp_three_variables, T=(300, 400), P=(100000, 150000), V=(0.1, 0.2))
        # Проверяем что функция выполняется без ошибок
        assert isinstance(result, float)
        assert result > 0

    def test_reverse_integration_limits(self):
        """Тест с обратными пределами интегрирования"""
        result_forward = enthalpy(self.cp_linear, T=(300, 500))
        result_reverse = enthalpy(self.cp_linear, T=(500, 300))

        # Интеграл в обратном направлении должен быть отрицательным
        assert abs(result_forward + result_reverse) < 1e-10

    def test_zero_range(self):
        """Тест с нулевым диапазоном"""
        result = enthalpy(self.cp_constant, T=(400, 400))
        assert abs(result) < 1e-10

    def test_small_range(self):
        """Тест с малым диапазоном"""
        result = enthalpy(self.cp_linear, T=(400, 400.001))
        # Для линейной функции на малом интервале ≈ Cp_avg * ΔT
        cp_avg = (self.cp_linear(400) + self.cp_linear(400.001)) / 2
        expected = cp_avg * 0.001
        assert abs(result - expected) < 1e-8

    def test_error_cases(self):
        """Тесты обработки ошибок"""

        # Не callable объект
        with pytest.raises((AssertionError, TypeError)):
            enthalpy("not a function", T=(300, 500))

        # Отсутствующий аргумент
        with pytest.raises(Exception, match="require argument"):
            enthalpy(self.cp_linear)  # T отсутствует

        # Неправильный тип диапазона
        with pytest.raises((AssertionError, TypeError)):
            enthalpy(self.cp_linear, T="not a range")

        # Диапазон неправильной длины
        with pytest.raises((AssertionError, ValueError)):
            enthalpy(self.cp_linear, T=(300, 400, 500))

        # Лишние аргументы (должны игнорироваться)
        result = enthalpy(self.cp_linear, T=(300, 500), extra_arg=(100, 200))
        expected = enthalpy(self.cp_linear, T=(300, 500))
        assert abs(result - expected) < 1e-10

    def test_lambda_functions(self):
        """Тест с lambda-функциями"""
        # Lambda с одной переменной
        result1 = enthalpy(lambda T: 1000.0, T=(300, 500))
        expected1 = 1000.0 * 200
        assert abs(result1 - expected1) < 1e-10

        # Lambda с двумя переменными
        result2 = enthalpy(lambda T, P: 1000.0 + 0.1 * T, T=(300, 500), P=(100000, 200000))
        expected2 = (1000.0 + 0.1 * 400) * 200 * 100000  # Среднее значение
        assert abs(result2 - expected2) < 1e6  # Большая погрешность из-за приближения

    def test_edge_cases(self):
        """Тесты граничных случаев"""

        # Очень большой диапазон
        result = enthalpy(self.cp_constant, T=(0, 10000))
        expected = 1000.0 * 10000
        assert abs(result - expected) < 1e-6

        # Отрицательный диапазон
        result = enthalpy(self.cp_linear, T=(-100, 100))
        expected = 800.0 * 200 + 0.5 * (10000 - 10000) / 2  # T² симметрично
        assert abs(result - expected) < 1e-10

        # Очень маленький диапазон
        result = enthalpy(self.cp_constant, T=(400, 400 + 1e-10))
        assert abs(result) < 1e-6

    @pytest.mark.parametrize(
        "function, kwargs, expected",
        [
            (cp_constant.__func__, {"T": (0, 100)}, 100000.0),  # 1000 * 100
            (cp_linear.__func__, {"T": (0, 100)}, 80000.0 + 0.5 * 10000 / 2),  # ∫(800+0.5T)dT
            (cp_constant.__func__, {"T": (100, 200)}, 100000.0),  # 1000 * 100
        ],
    )
    def test_parametrized(self, function, kwargs, expected):
        """Параметризованные тесты"""
        result = enthalpy(function, **kwargs)
        assert abs(result - expected) < 1e-10


class TestIntegralAverage:
    """Тесты для функции integral_average"""

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
            (f0_1.__func__, {"x": (0, 1)}, 5.0),
            (f0_1.__func__, {"x": (1, 1)}, 5.0),
            # f0_2
            (f0_2.__func__, {"x": (0, 1), "y": (0, 2)}, 5.0),
            (f0_2.__func__, {"x": (0, 1), "y": (2, 2)}, 5.0),
            (f0_2.__func__, {"x": (1, 1), "y": (0, 2)}, 5.0),
            (f0_2.__func__, {"x": (1, 1), "y": (2, 2)}, 5.0),
            # f0_3
            (f0_3.__func__, {"x": (0, 1), "y": (0, 2), "z": (0, 3)}, 5.0),
            (f0_3.__func__, {"x": (1, 1), "y": (0, 2), "z": (0, 3)}, 5.0),
            (f0_3.__func__, {"x": (0, 1), "y": (2, 2), "z": (0, 3)}, 5.0),
            (f0_3.__func__, {"x": (0, 1), "y": (0, 2), "z": (3, 3)}, 5.0),
            (f0_3.__func__, {"x": (1, 1), "y": (2, 2), "z": (0, 3)}, 5.0),
            (f0_3.__func__, {"x": (1, 1), "y": (0, 2), "z": (3, 3)}, 5.0),
            (f0_3.__func__, {"x": (0, 1), "y": (2, 2), "z": (3, 3)}, 5.0),
            (f0_3.__func__, {"x": (1, 1), "y": (2, 2), "z": (3, 3)}, 5.0),
            # f1_1
            (f1_1.__func__, {"x": (0, 2)}, (2**2) / 2 / (2 - 0)),
            (f1_1.__func__, {"x": (2, 2)}, 2.0),
            # f1_2
            (f1_2.__func__, {"x": (0, 1), "y": (0, 2)}, ((1**2 - 0**2) / 2 * (2 - 0) + (2**2 - 0**2) / 2 * (1 - 0)) / (1 - 0) / (2 - 0)),
            (f1_2.__func__, {"x": (1, 1), "y": (0, 2)}, (1 * (2 - 0) + (2**2 - 0**2) / 2) / (2 - 0)),
            (f1_2.__func__, {"x": (0, 1), "y": (2, 2)}, (2 * (1 - 0) + (1**2 - 0**2) / 2) / (1 - 0)),
            (f1_2.__func__, {"x": (1, 1), "y": (2, 2)}, 1 + 2),
            # f1_3
            (f1_3.__func__, {"x": (0, 1), "y": (0, 2), "z": (0, 3)}, ((1**2 / 2 * (2 - 0) * (3 - 0)) + (2**2 / 2 * (1 - 0) * (3 - 0)) + (3**2 / 2 * (1 - 0) * (2 - 0))) / (1 - 0) / (2 - 0) / (3 - 0)),
            (f1_3.__func__, {"x": (1, 1), "y": (0, 2), "z": (0, 3)}, ((1 * (2 - 0) * (3 - 0)) + (2**2 / 2 * (1 - 0) * (3 - 0)) + (3**2 / 2 * (1 - 0) * (2 - 0))) / (2 - 0) / (3 - 0)),
            (f1_3.__func__, {"x": (0, 1), "y": (2, 2), "z": (3, 3)}, ((1**2 - 0**2) / 2 + 2 + 3) / (1 - 0)),
            (f1_3.__func__, {"x": (1, 1), "y": (2, 2), "z": (3, 3)}, 1 + 2 + 3),
            # f2_1
            (f2_1.__func__, {"x": (0, 1)}, 1**3 / 3),
            (f2_1.__func__, {"x": (1, 1)}, 1),
            # f2_2
            # f2_3
            # exp
            (exp.__func__, {"x": (0, 1)}, (np.exp(1) - 1) / (1 - 0)),
            (exp.__func__, {"x": (1, 1)}, np.exp(1)),
            # sin
            (sin.__func__, {"x": (0, np.pi)}, (1 - np.cos(np.pi)) / (np.pi - 0)),
            (sin.__func__, {"x": (np.pi, np.pi)}, np.sin(np.pi)),
            # discontinuous
            (discontinuous.__func__, {"x": (0, 1)}, 1.5),
        ],
    )
    def test_functions(self, function, kwargs, expected):
        result, error = integral_average(function, **kwargs)
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


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
