from abc import ABC, abstractmethod

from numpy import array, nan, prod
from substance import Substance

from src.config import parameters as gtep
from src.errors import SUBSTANCE_ATTRIBUTE_ERROR


class GTENode(ABC):
    """Абстрактный базовый класс узла ГТД"""

    def __init__(self, name: str = "node") -> None:
        self.name: str = name

        self.inlet = Substance("inlet")
        self.outlet = Substance("outlet")

        self.mass_flow_leak = nan

    def __str__(self) -> str:
        return self.name

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

    @abstractmethod
    def calculate(self, Niter: int = 10) -> Substance:
        """Расчет узла"""
        # расчет входных параметров
        # расчет параметров узла
        # расчет выходных параметров
        # вывод выходных параметров
        pass

    @property
    def summary(self) -> dict:
        return {
            **{f"{self.name}_{k}": v for k, v in self.__dict__.items()},
            **{f"{k}_inlet": v for k, v in self.inlet.parameters.items()},
            **{f"{k}_outlet": v for k, v in self.outlet.parameters.items()},
        }


if __name__ == "__main__":
    n = GTENode()
    print(n.summary)
