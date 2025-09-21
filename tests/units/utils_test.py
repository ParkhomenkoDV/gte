import pytest

from src.utils import call_with_kwargs


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
