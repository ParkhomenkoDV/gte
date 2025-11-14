from abc import ABC, abstractmethod
from typing import Any, Dict, Tuple

from numpy import isnan
from substance import Substance
from thermodynamics import parameters as tdp

try:  # Попытка относительного импорта
    from .config import parameters as gtep
    from .errors import SUBSTANCE_ATTRIBUTE_ERROR
except ImportError:  # Резервный абсолютный импорт
    from config import parameters as gtep
    from errors import SUBSTANCE_ATTRIBUTE_ERROR

"""
Порядок расчета ТД параметров:
G -> excess_oxidizing -> gas_const -> T* -> P* -> D* -> Cp -> k -> a* -> c
"""


class GTENode(ABC):
    """Абстрактный базовый класс узла ГТД"""

    models: Dict[str, Any] = {}  # ML модели

    __slots__ = ["name", "inlet", "outlet", "leak"]  # list

    def __init__(self, name: str = "node") -> None:
        assert isinstance(name, str), TypeError(f"type name must be str, but has {type(name)}")
        self.name: str = name
        self.leak: float = 0

    def __str__(self) -> str:
        return self.name

    def __delattr__(self, name: str) -> None:
        if name == "name":
            self.name = self.__class__.__name__
        elif name == "inlet":
            self.inlet = Substance("inlet")
        elif name == "outlet":
            self.outlet = Substance("outlet")
        else:
            return super().__delattr__(name)

    @property
    @abstractmethod
    def variables(self) -> Dict[str, float]:
        """Переменные параметры"""
        return {}

    def validate_substance(self, substance: Substance) -> None:
        """Проверка параметров рабочего тела на входе"""
        assert isinstance(substance, Substance), TypeError("type substance must be Substance")
        assert substance.parameters.get(gtep.mf) is not None, AttributeError(SUBSTANCE_ATTRIBUTE_ERROR.format(substance.name, gtep.m))
        assert substance.parameters.get(gtep.TT) is not None, AttributeError(SUBSTANCE_ATTRIBUTE_ERROR.format(substance.name, gtep.TT))
        assert substance.parameters.get(gtep.PP) is not None, AttributeError(SUBSTANCE_ATTRIBUTE_ERROR.format(substance.name, gtep.PP))
        # validate functions
        tdp_keys = tdp.values()  # разрешенный список термодинамических параметров
        for name, function in substance.functions.items():
            assert name in tdp_keys, KeyError(f"function '{name}' not in {tdp_keys}")
            for func_arg in function.__code__.co_varnames:
                assert func_arg in tdp_keys, NameError(f"function '{name}' has arg '{func_arg}' not in {tdp_keys}")

    def is_solvable(self, *args) -> str:
        """Проверка резрешимости математической модели"""
        try:
            prediction = self.predict(*args, use_ml=False)  # для быстроты
        except Exception as e:
            return f"{e}"

        count_outlet_variables = sum(1 if k.startswith("outlet") else 0 for k in prediction)  # количество неизвестных выходных параметров
        args = {k: v for k, v in self.variables.items() if not isnan(v)}  # известные параметры узла
        count_variables = len(self.variables) - len(args)  # количество неизвестных параметров узла
        count_equations = len(self._equations(tuple(), args))  # количество уравнений

        if count_variables + count_outlet_variables < count_equations:
            return f"{count_variables=} + {count_outlet_variables=} < {count_equations=}"
        elif count_variables + count_outlet_variables > count_equations:
            return f"{count_variables=} + {count_outlet_variables=} > {count_equations=}"
        else:
            return ""

    @abstractmethod
    def predict(self, inlet: Substance, use_ml: bool = True) -> Dict[str, float]:
        """Начальные приближения"""
        self.validate_substance(inlet)
        assert isinstance(use_ml, bool), TypeError(f"{type(use_ml)=} must be bool")
        return {}

    @property
    def summary(self) -> Dict[str, float]:
        result = {}
        for k in self.__slots__:
            if not isinstance(getattr(self, k), Substance):
                result[k] = getattr(self, k)
            else:
                for parameter, value in getattr(self, k).parameters.items():
                    result[f"{k}_{parameter}"] = value
        return result

    def _equations(self, x: Tuple[float], args: Dict[str, Any]) -> Tuple:
        """Уравнения"""
        return tuple()

    @abstractmethod
    def solve(self, x0: Dict = None) -> Substance:
        """Расчет узла"""
        # валидация входных параметров
        # расчет входных параметров
        # расчет параметров узла
        # расчет выходных параметров
        # вывод выходных параметров
        pass

    @abstractmethod
    def is_real(self) -> bool:
        """Проверка физичной реальности"""
        return False
