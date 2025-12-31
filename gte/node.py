from abc import ABC, abstractmethod
from typing import Any, Dict, Tuple

from substance import Substance
from thermodynamics import adiabatic_index, сritical_sonic_velocity
from thermodynamics import parameters as tdp

try:  # Попытка относительного импорта
    from .config import parameters as gtep
    from .errors import SUBSTANCE_ATTRIBUTE_ERROR
    from .utils import call_with_kwargs
except ImportError:  # Резервный абсолютный импорт
    from config import parameters as gtep
    from errors import SUBSTANCE_ATTRIBUTE_ERROR
    from utils import call_with_kwargs

"""
Порядок расчета ТД параметров:
G -> excess_oxidizing -> gas_const -> T* -> P* -> D* -> Cp -> k -> a* -> c
"""


class GTENode(ABC):
    """Абстрактный базовый класс узла ГТД"""

    models: Dict[str, Any] = {}  # ML модели

    __slots__ = ["name", "characteristic"]  # list to add

    def __init__(self, name: str = "node") -> None:
        assert isinstance(name, str), TypeError(f"{type(name)=} must be str")
        self.name: str = name

    def __str__(self) -> str:
        return self.name

    def __delattr__(self, name: str) -> None:
        if name == "name":
            self.name = self.__class__.__name__
        elif name == "characteristic":
            raise Exception("deleting characteristic is prohibited!")
        else:
            return super().__delattr__(name)

    @classmethod
    def validate_substance(cls, substance: Substance) -> None:
        """Проверка обязательных параметров рабочего тела"""
        assert isinstance(substance, Substance), TypeError(f"type substance must be {Substance}")
        # validate parameters
        assert gtep.mf in substance.parameters, AttributeError(SUBSTANCE_ATTRIBUTE_ERROR.format(substance.name, gtep.m))
        assert gtep.TT in substance.parameters, AttributeError(SUBSTANCE_ATTRIBUTE_ERROR.format(substance.name, gtep.TT))
        assert gtep.PP in substance.parameters, AttributeError(SUBSTANCE_ATTRIBUTE_ERROR.format(substance.name, gtep.PP))
        # validate functions
        tdp_keys = tdp.values()  # разрешенный список термодинамических параметров
        for name, function in substance.functions.items():
            assert name in tdp_keys, KeyError(f"function '{name}' not in {tdp_keys}")
            for func_arg in function.__code__.co_varnames:
                assert func_arg in tdp_keys, NameError(f"function '{name}' has arg '{func_arg}' not in {tdp_keys}")

    @classmethod
    def calculate_substance(cls, substance: Substance) -> Substance:
        """Расчет термодинамических параметров вещества по массе, температуре, давлению"""
        cls.validate_substance(substance)

        parameters = {tdp.t: substance.parameters.get(gtep.TT), tdp.p: substance.parameters.get(gtep.PP), tdp.eo: substance.parameters.get(gtep.eo)}
        substance.parameters[gtep.gc] = call_with_kwargs(substance.functions[gtep.gc], parameters)
        substance.parameters[gtep.Cp] = call_with_kwargs(substance.functions[gtep.Cp], parameters)
        substance.parameters[gtep.DD] = substance.parameters[gtep.PP] / (substance.parameters[gtep.gc] * substance.parameters[gtep.TT])
        substance.parameters[gtep.k] = adiabatic_index(substance.parameters[gtep.gc], substance.parameters[gtep.Cp])
        substance.parameters[gtep.a_critical] = сritical_sonic_velocity(substance.parameters[gtep.k], substance.parameters[gtep.gc], substance.parameters[gtep.TT])

        return substance

    @classmethod
    @abstractmethod
    def predict(cls, inlet: Substance, use_ml: bool = True) -> Dict[str, float]:
        """Начальные приближения"""
        cls.validate_substance(inlet)

        assert isinstance(use_ml, bool), TypeError(f"{type(use_ml)=} must be bool")

        return {}

    @classmethod
    def _equations(cls, x: Tuple[float], args: Dict[str, Any]) -> Tuple:
        """Система уравнений"""
        return tuple()

    @classmethod
    @abstractmethod
    def calculate(cls, x0: Dict[str, float] = None) -> Substance:
        """Термодинамический расчет узла по СНЛАУ _equations"""
        # валидация входных параметров
        # расчет входных параметров
        # расчет параметров узла
        # расчет выходных параметров
        # вывод выходных параметров
        pass


class TurboCompressor(ABC):
    """Абстрактный класс турбокомрпессора"""

    @abstractmethod
    def get_power(self, inlet: Substance, rotation_frequency: float) -> float:
        """Расчет мощности по частоте вращения и характеристике узла"""
        return 0
