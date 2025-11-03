from abc import ABC, abstractmethod
from typing import Any, Dict, Tuple

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

    __slots__ = ("name", "inlet", "outlet", "leak")

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

    @abstractmethod
    def predict(self) -> Dict[str, float]:
        """Начальные приближения"""
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

    def equations(self, x: Tuple, args: Dict) -> Tuple:
        """Уравнения"""
        pass

    @abstractmethod
    def calculate(self, x0: Dict = None) -> Substance:
        """Расчет узла"""
        # расчет входных параметров
        # расчет параметров узла
        # расчет выходных параметров
        # вывод выходных параметров
        pass

    @abstractmethod
    def is_real(self) -> bool:
        """Проверка физичной реальности"""
