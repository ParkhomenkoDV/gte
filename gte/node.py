from abc import ABC, abstractmethod

from numpy import array, prod
from substance import Substance

try:  # Попытка относительного импорта
    from .config import parameters as gtep
    from .errors import SUBSTANCE_ATTRIBUTE_ERROR
except ImportError:  # Резервный абсолютный импорт
    from config import parameters as gtep
    from errors import SUBSTANCE_ATTRIBUTE_ERROR

"""
Порядок расчета ТД параметров:
G -> excess_oxidizing -> gas_const -> T -> P -> D -> Cp -> k -> a* -> c
"""


class GTENode(ABC):
    """Абстрактный базовый класс узла ГТД"""

    __slots__ = ("name", "inlet", "outlet", "mass_flow_leak")

    def __init__(self, name: str = "node") -> None:
        assert isinstance(name, str), TypeError(f"type name must be str, but has {type(name)}")
        self.name: str = name

        self.inlet = Substance("inlet")
        self.outlet = Substance("outlet")

        self.mass_flow_leak: float = 0

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
    def variables(self) -> dict[str:float]:
        """Переменные параметры"""
        return {}

    @property
    @abstractmethod
    def _x0(self) -> dict:
        """Начальные приближения"""
        return {}

    @property
    def summary(self) -> dict[str:float]:
        result = {
            **{k: getattr(self, k) for k in self.__slots__ if not isinstance(getattr(self, k), Substance)},
            **{f"{k}_inlet": v for k, v in self.inlet.parameters.items()},
            **{f"{k}_outlet": v for k, v in self.outlet.parameters.items()},
        }

        n = 20
        print("-" * n)
        for k in self.__slots__:
            if not k.endswith("let"):
                print(f"{k}: {getattr(self, k)}")
        for k, v in self.inlet.parameters.items():
            print(f"{k}_inlet: {v}")
        for k, v in self.outlet.parameters.items():
            print(f"{k}_outlet: {v}")
        print("-" * n)

        return result

    def get_variability(self) -> int:
        """Максимальное количество комбинаций варьируемых параметров"""
        return prod([len(value) for key, value in self.__dict__.items() if isinstance(value, (list, tuple, array)) and len(value) and not key.startswith("_")])

    def set_combination(self, combination: int, main_node) -> None:
        """Установка комбинации"""
        varible_params = [key for key, value in main_node.__dict__.items() if type(value) is list and len(value) and not key.startswith("_")]
        positions = [0] * len(varible_params)

        for i in range(combination):
            for j, param in enumerate(varible_params):
                if positions[j] == len(getattr(main_node, varible_params[j])) - 1:
                    positions[j] = 0
                else:
                    positions[j] += 1
                    continue
        for j, param in enumerate(varible_params):
            setattr(
                self,
                varible_params[j],
                getattr(main_node, varible_params[j])[positions[j]],
            )

    def validate_substance(self, substance: Substance) -> None:
        """Проверка параметров рабочего тела на входе"""
        assert isinstance(substance, Substance), TypeError("type substance must be Substance")
        assert substance.parameters.get(gtep.mf) is not None, AttributeError(SUBSTANCE_ATTRIBUTE_ERROR.format(substance.name, gtep.mf))
        assert substance.parameters.get(gtep.TT) is not None, AttributeError(SUBSTANCE_ATTRIBUTE_ERROR.format(substance.name, gtep.TT))
        assert substance.parameters.get(gtep.PP) is not None, AttributeError(SUBSTANCE_ATTRIBUTE_ERROR.format(substance.name, gtep.PP))
        # validate functions
        td_keys = gtep.values()  # разрешенный список термодинамических параметров
        for name, function in substance.functions.items():
            assert name in gtep.values(), KeyError(f"function '{name}' not in {td_keys}")
            for func_arg in function.__code__.co_varnames:
                assert func_arg in td_keys, NameError(f"function '{name}' has arg '{func_arg}' not in {td_keys}")

    def equations(self, x: tuple, args: dict) -> tuple:
        """Уравнения"""
        pass

    @abstractmethod
    def calculate(self, x0: dict = None) -> Substance:
        """Расчет узла"""
        # расчет входных параметров
        # расчет параметров узла
        # расчет выходных параметров
        # вывод выходных параметров
        pass

    @abstractmethod
    def is_real(self) -> bool:
        """Проверка физичной реальности"""
