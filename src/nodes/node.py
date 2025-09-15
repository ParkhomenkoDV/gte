from abc import ABC, abstractmethod

from numpy import array, isnan, prod
from substance import Substance

from src.config import parameters as gtep
from src.errors import SUBSTANCE_ATTRIBUTE_ERROR

"""
Порядок расчета ТД параметров:
excess_oxidizing -> gas_const -> G -> T -> P -> D -> Cp -> k ->
"""


class GTENode(ABC):
    """Абстрактный базовый класс узла ГТД"""

    def __init__(self, name: str = "node") -> None:
        self.name: str = name

        self.inlet = Substance("inlet")
        self.outlet = Substance("outlet")

        self.mass_flow_leak: float = 0

    def __str__(self) -> str:
        return self.name

    def __delattr__(self, name):
        if name == "name":
            self.name = self.__class__.__name__
        elif name == "inlet":
            self.inlet = Substance("inlet")
        elif name == "outlet":
            self.outlet = Substance("outlet")
        else:
            return super().__delattr__(name)

    @property
    # @abstractmethod
    def variables(self) -> dict[str:float]:
        return {}

    @property
    def summary(self) -> dict[str:float]:
        result = {
            **{f"{self.name}_{k}": v for k, v in self.__dict__.items()},
            **{f"{k}_inlet": v for k, v in self.inlet.parameters.items()},
            **{f"{k}_outlet": v for k, v in self.outlet.parameters.items()},
        }

        n = 20
        print("-" * n)
        for k, v in self.__dict__.items():
            if k not in ("inlet", "outlet"):
                print(f"{k}: {v}")
        for k, v in self.inlet.parameters.items():
            print(f"{k}_inlet: {v}")
        for k, v in self.outlet.parameters.items():
            print(f"{k}_outlet: {v}")
        print("-" * n)

        return result

    def get_variability(self) -> int:
        """Максимальное количество комбинаций варьируемых параметров"""
        return prod(
            [
                len(value)
                for key, value in self.__dict__.items()
                if isinstance(value, (list, tuple, array))
                and len(value)
                and not key.startswith("_")
            ]
        )

    def set_combination(self, combination: int, main_node) -> None:
        """Установка комбинации"""
        varible_params = [
            key
            for key, value in main_node.__dict__.items()
            if type(value) is list and len(value) and not key.startswith("_")
        ]
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
        assert isinstance(substance, Substance), TypeError(
            "type substance must be Substance"
        )
        assert substance.parameters.get(gtep.TT), AttributeError(
            SUBSTANCE_ATTRIBUTE_ERROR.format(substance.name, gtep.TT)
        )
        assert substance.parameters.get(gtep.PP), AttributeError(
            SUBSTANCE_ATTRIBUTE_ERROR.format(substance.name, gtep.PP)
        )
        assert substance.parameters.get(gtep.mf), AttributeError(
            SUBSTANCE_ATTRIBUTE_ERROR.format(substance.name, gtep.mf)
        )

    def equations(self, x, *args: tuple) -> tuple:
        """Уравнения"""
        pass

    @abstractmethod
    def calculate(self, Niter: int = 10) -> Substance:
        """Расчет узла"""
        # расчет входных параметров
        # расчет параметров узла
        # расчет выходных параметров
        # вывод выходных параметров
        pass
