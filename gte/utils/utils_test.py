from typing import Dict

import numpy as np
import pytest

try:
    from .utils import Function, Interpolator, integral_average, integrate
except ImportError:
    from gte.utils.utils import Function, Interpolator, integral_average, integrate


def complex_function(x):
    """Сложная непрерывная функция"""
    return np.sin(x) + np.log(abs(x) + 1)


class TestFunction:
    """Тесты для класса Function"""

    @pytest.mark.parametrize(
        "function, name, args, want_name, want_args",
        [
            (complex_function, "complex_function", ("x",), "complex_function", ("x",)),
            (complex_function, "", tuple(), "complex_function", ("x",)),
            (lambda x, y: x + y, "lambda", ("x", "y"), "lambda", ("x", "y")),
            (lambda x, y: x + y, "", tuple(), "<lambda>", ("x", "y")),
        ],
    )
    def test_init(self, function, name, args, want_name, want_args):
        """Тест создания Function"""

        if name == "" and len(args) == 0:
            func = Function(function)
        else:
            func = Function(function, name=name, args=args)

        assert func.function == function
        assert func.name == want_name
        assert func.args == want_args

    def test_len(self):
        """Тест метода __len__"""

        def func1(a, b):
            return a + b

        def func2(x, y, z, w):
            return x * y * z * w

        assert len(Function(func1)) == 2
        assert len(Function(func2)) == 4

    @pytest.mark.parametrize(
        "kwargs,expected",
        [
            ({"a": 1, "b": 2}, 3),
            ({"a": 10, "b": 20}, 30),
            ({"a": -5, "b": 5}, 0),
            ({"a": 0, "b": 0}, 0),
        ],
    )
    def test_call(self, kwargs: Dict[str, float], expected: float):
        """Тест вызова функции"""

        def add(a, b):
            return a + b

        func = Function(add)
        assert func(kwargs) == expected

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"a": 1, "b": 2},
            {"a": 10, "b": 20},
            {"a": -5, "b": 5},
            {"a": 0, "b": 0},
        ],
    )
    @pytest.mark.benchmark
    def test_function_call(self, benchmark, kwargs: Dict[str, float]):
        """Бенчмарк вызова функции"""

        def add(a, b):
            return a + b

        func = Function(add)

        def benchfunc():
            return func(kwargs)

        benchmark(benchfunc)

    def test_call_with_missing_args(self):
        """Тест вызова функции с отсутствующими аргументами"""

        def multiply(x, y, z):
            return x * y * z

        func = Function(multiply)

        with pytest.raises(KeyError):
            func({"x": 5, "y": 3})  # missing 'z'

    def test_call_with_extra_args(self):
        """Тест вызова функции с избыточными параметрами"""

        def subtract(a, b):
            return a - b

        func = Function(subtract)
        result = func({"a": 10, "b": 4, "c": 100, "d": 200})

        assert result == 6

    def test_call_with_wrong_arg_types(self):
        """Тест вызова функции с неправильными типами аргументов"""

        def divide(a, b):
            return a / b

        func = Function(divide)

        with pytest.raises(TypeError):
            func({"a": 10, "b": "not a number"})

    def test_function_with_default_args(self):
        """Тест функции с аргументами по умолчанию"""

        def greet(name, greeting="Hello"):
            return f"{greeting}, {name}!"

        func = Function(greet)
        result = func({"name": "World", "greeting": "Hi"})

        assert result == "Hi, World!"

    @pytest.mark.skip
    def test_function_with_variable_args_detection(self):
        """Тест определения аргументов для функции с *args и **kwargs"""

        def variable_func(*args, **kwargs):
            return len(args) + len(kwargs)

        # Для *args и **kwargs code.co_varnames будет содержать 'args' и 'kwargs'
        func = Function(variable_func)

        # Должны быть определены 'args' и 'kwargs'
        assert "args" in func.args or "kwargs" in func.args

    def test_repr(self):
        """Тест строкового представления"""

        def test_func(data):
            return data

        func = Function(test_func)

        assert repr(func) == "test_func"
        assert str(func) == "test_func"

    def test_with_nested_function(self):
        """Тест с вложенной функцией"""

        def outer():
            def inner(a, b):
                return a * b

            return inner

        inner_func = outer()
        func = Function(inner_func)

        assert func({"a": 3, "b": 4}) == 12

    def test_function_with_non_dict_call(self):
        """Тест вызова с не-словарём (должен вызвать ошибку)"""

        def test(x):
            return x

        func = Function(test)

        with pytest.raises(TypeError):
            func([1, 2, 3])  # Передаём список вместо словаря

    def test_call_order_independent(self):
        """Тест что порядок аргументов не важен"""

        def divide(dividend, divisor):
            return dividend / divisor

        func = Function(divide)

        result1 = func({"dividend": 10, "divisor": 2})
        result2 = func({"divisor": 2, "dividend": 10})

        assert result1 == result2 == 5


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
            (Function(f0_1.__func__), {"T": (300, 500)}, 1000.0 * (500 - 300)),
            (Function(f1_1.__func__), {"T": (300, 500)}, 800 * (500 - 300) + 0.5 * (500**2 - 300**2) / 2),
            (Function(f2_1.__func__), {"T": (300, 500)}, 900 * (500 - 300) + 0.1 * (500**2 - 300**2) / 2 + 0.001 * (500**3 - 300**3) / 3),
            (Function(f1_2.__func__), {"T": (300, 500), "P": (100_000, 200_000)}, 1000 * (500 - 300) * (200000 - 100000) + 0.2 * (500**2 - 300**2) / 2 * (200000 - 100000) + 0.0001 * (500 - 300) * (200000**2 - 100000**2) / 2),
            (Function(f1_3.__func__), {"T": (300, 500), "P": (100_000, 200_000), "V": (0.3, 0.4)}, 2_080_000_700),
        ],
    )
    def test_function(self, function, kwargs, expected):
        result, error = integrate(function, **kwargs)
        assert result == pytest.approx(expected, rel=self.RELERR)
        assert error < self.ABSERR

    @pytest.mark.benchmark
    def test_integrate(self, benchmark):
        """Бенчмарк для функции integrate()"""

        def benchfunc(ranges):
            integrate(Function(complex_function), **ranges)

        benchmark(benchfunc, {"x": (-100, 100)})

    def test_reverse_integration_limits(self):
        """Тест с обратными пределами интегрирования"""
        result_forward, error_forward = integrate(Function(self.f1_1), T=(300, 500))
        result_reverse, error_reverse = integrate(Function(self.f1_1), T=(500, 300))
        assert abs(result_forward + result_reverse) < 1e-10
        assert error_forward < self.ABSERR
        assert error_reverse < self.ABSERR

    def test_zero_range(self):
        """Тест с нулевым диапазоном"""
        result, error = integrate(Function(self.f0_1), T=(400, 400))
        assert result == pytest.approx(1000.0, rel=self.RELERR)
        assert error < self.ABSERR

    def test_error_cases(self):
        """Тесты обработки ошибок"""

        # Не callable объект
        with pytest.raises(Exception):
            integrate("not a function", T=(300, 500))

        # Отсутствующий аргумент
        with pytest.raises((AssertionError, Exception)):
            integrate(Function(self.f1_1, args=("no_T",)))  # T отсутствует

        # Неправильный тип диапазона
        with pytest.raises((AssertionError, TypeError)):
            integrate(Function(self.f1_1, args=("T",)), T="not a range")

        # Диапазон неправильной длины
        with pytest.raises(Exception):
            integrate(Function(self.f1_1, args="T"), T=(300, 400, 500))

        # Лишние аргументы (должны игнорироваться)
        result, error1 = integrate(Function(self.f1_1, args=("T",)), T=(300, 500), extra_arg=(100, 200))
        expected, error2 = integrate(Function(self.f1_1, args=("T",)), T=(300, 500))
        assert result == pytest.approx(expected, rel=self.RELERR)
        assert error1 == error2

    def test_lambda_functions(self):
        """Тест с lambda-функциями"""
        result, error = integrate(
            Function(lambda T: 1000.0, args=("T",)),
            T=(300, 500),
        )
        assert result == pytest.approx(1000.0 * 200, rel=self.RELERR)
        assert error < self.ABSERR

        result, error = integrate(Function(lambda T, P: 1000.0 + 0.1 * T, args=("T", "P")), T=(300, 500), P=(100_000, 200_000))
        assert result == pytest.approx((1000.0 + 0.1 * 400) * 200 * 100000, rel=self.RELERR)
        assert error < self.ABSERR

    def test_edge_cases(self):
        """Тесты граничных случаев"""

        # Очень большой диапазон
        result, error = integrate(Function(self.f0_1, args=("T",)), T=(0, 10_000))
        assert result == pytest.approx(1000.0 * 10_000, rel=1e-6)
        assert error < self.ABSERR

        # Отрицательный диапазон
        result, error = integrate(Function(self.f1_1, args=("T",)), T=(-100, 100))
        expected = 800.0 * 200 + 0.5 * (10_000 - 10_000) / 2  # T² симметрично
        assert result == pytest.approx(expected, rel=1e-6)
        assert error < self.ABSERR

        # Очень маленький диапазон
        result, error = integrate(Function(self.f0_1, args=("T",)), T=(400, 400 + 1e-10))
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
        result, error = integral_average(Function(function), **kwargs)
        assert result == pytest.approx(expected, rel=self.RELERR)
        assert error < self.ABSERR

    @pytest.mark.benchmark
    def test_integral_average(self, benchmark):
        """Бенчмарк для функции integral_average()"""

        def benchfunc(function, kwargs):
            integral_average(Function(function), **kwargs)

        benchmark(benchfunc, complex_function, {"x": (-100, 100)})

    def test_error_cases(self):
        """Тесты обработки ошибок"""

        def test_func(x):
            return x

        # Не callable объект
        with pytest.raises(Exception):
            integral_average("not a function", x=(0, 1))

        # Отсутствующий аргумент
        with pytest.raises((AssertionError, Exception)):
            integral_average(Function(test_func), y=(0, 1))  # x отсутствует

        # Неправильный тип диапазона
        with pytest.raises((AssertionError, TypeError)):
            integral_average(Function(test_func), x="not a range")

        # Диапазон неправильной длины
        with pytest.raises((AssertionError, ValueError)):
            integral_average(Function(test_func), x=(0, 1, 2))


def f1(x):
    return float(np.sin(x))


f1_data = [{"x": x, "z": f1(x)} for x in np.linspace(-1, 1, 21)]


def f2(x, y):
    return np.sin(x) * np.exp(y)


f2_data = [{"x": x, "y": y, "z": f2(x, y)} for x in np.linspace(-1, 1, 21) for y in np.linspace(-1, 1, 21)]


class TestInterpolator:
    """Тесты для класса Interpolator"""

    @pytest.mark.parametrize(
        "datas, target, features, fill_value",
        [
            # 1D
            (f1_data, "z", ["x"], np.nan),
            (f1_data, "z", None, np.nan),
            (f1_data, "z", ["x"], 0),
            (f1_data, "z", None, 0),
            # 2D
            (f2_data, "z", ["x", "y"], np.nan),
            (f2_data, "z", None, np.nan),
            (f2_data, "z", ["x", "y"], 0),
            (f2_data, "z", None, 0),
        ],
    )
    def test_init(self, datas, target, features, fill_value):
        """Тест инициализации"""

        interpolator = Interpolator(datas, target, features, fill_value)

        assert interpolator.target == target
        if features is not None:
            assert set(interpolator.features) == set(features)
        assert interpolator.fill_value == fill_value or (np.isnan(interpolator.fill_value) and np.isnan(fill_value))
        assert callable(interpolator.function)
        assert callable(interpolator)

        for data in datas:
            result = interpolator(**data)
            assert isinstance(result, float)
            result == pytest.approx(data[target], rel=1e-6)

    @pytest.mark.parametrize(
        "datas, target, features, fill_value",
        [
            # 1D
            (f1_data, "z", ["x"], np.nan),
            (f1_data, "z", None, np.nan),
            # 2D
            (f2_data, "z", ["x", "y"], np.nan),
            (f2_data, "z", None, np.nan),
        ],
    )
    @pytest.mark.benchmark
    def test_interpolator_init(self, benchmark, datas, target, features, fill_value):
        """Бенчмарк метода init"""

        def benchfunc(datas, target, features, fill_value):
            Interpolator(datas, target, features, fill_value)

        benchmark(benchfunc, datas, target, features, fill_value)

    @pytest.mark.parametrize(
        "datas,target,features,fill_value,expected_error",
        [
            ([{"x": 1, "y": 1, "z": 1}], "z", None, np.nan, ValueError),  # недостаточно данных
            ("not a list", "z", None, np.nan, TypeError),  # неверный тип datas
            (None, "z", None, np.nan, TypeError),  # datas = None
            ([{"x": 1, "y": 1, "z": 1} for _ in range(5)], 123, None, np.nan, TypeError),  # неверный тип target
            ([{"x": 1, "y": 1, "z": 1} for _ in range(5)], "z", [], np.nan, ValueError),  # пустой список features
            ([{"x": 1, "y": 1, "z": 1} for _ in range(5)], "z", [1, 2, 3], np.nan, TypeError),  # features не из строк
            ([{"x": 1, "y": 1, "z": 1} for _ in range(5)], "nonexistent", None, np.nan, KeyError),  # отсутствие target
            ([{"x": 1, "y": 1, "z": 1} for _ in range(5)], "z", None, "invalid", TypeError),  # неверный fill_value
            ([{"x": 1, "y": 1, "z": 1} for _ in range(5)], "z", ["x", "nonexistent"], np.nan, KeyError),  # отсутствия указанных features
        ],
    )
    def test_init_error(self, datas, target, features, fill_value, expected_error):
        """Тест ошибок при инициализации"""
        with pytest.raises(expected_error):
            Interpolator(datas, target, features=features, fill_value=fill_value)

    def test_call_error(self):
        """Тест пропущенных аргументов при вызове"""
        interpolator = Interpolator(f2_data, "z")

        with pytest.raises(KeyError):
            interpolator(x=5)  # пропущен y

        with pytest.raises(Exception):
            interpolator(x="invalid", y=5)

    @pytest.mark.parametrize(
        "datas, target, features, fill_value",
        [
            # 1D
            (f1_data, "z", None, np.nan),
            # 2D
            (f2_data, "z", None, np.nan),
        ],
    )
    def test_extrapolation(self, datas, target, features, fill_value):
        """Тест экстраполяции (возврат fill_value)"""
        interpolator = Interpolator(datas, target, features=features, fill_value=fill_value)

        for data in datas:
            result = interpolator(**{k: 999_999_999 for k in data})
            assert np.isnan(result)

    @pytest.mark.parametrize(
        "datas, target",
        [
            (f1_data, "z"),
            (f2_data, "z"),
        ],
    )
    def test_repr(self, datas, target):
        """Тест строкового представления"""
        interpolator = Interpolator(datas, target)
        assert repr(interpolator) == target
        assert str(interpolator) == target

    @pytest.mark.parametrize(
        "datas, target",
        [
            (f1_data, "z"),
            (f2_data, "z"),
        ],
    )
    def test_duplicates(self, datas, target):
        """Тест обработки дублирующихся точек"""
        # Создаем данные с дубликатами
        data_with_duplicates = datas.copy()
        data_with_duplicates.append(datas[0])

        interpolator = Interpolator(data_with_duplicates, target)
        result = interpolator(**datas[0])
        assert isinstance(result, float)
        assert result == pytest.approx(datas[0][target], rel=1e-6)

    def test_different_feature_order(self):
        """Тест разного порядка features при вызове"""
        interpolator = Interpolator(f2_data, "z")

        # Порядок не должен влиять на результат
        result1 = interpolator(x=0.25, y=0.35)
        result2 = interpolator(y=0.35, x=0.25)
        assert result1 == result2


if __name__ == "__main__":
    pytest.main(
        [
            __file__,
            "-v",
            "-s",
            "-x",
            "--benchmark-columns=mean,min,max,stddev,median,rounds,outliers",
            "--benchmark-sort=name",
            "--benchmark-min-rounds=10",
        ]
    )
