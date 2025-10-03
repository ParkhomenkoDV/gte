from substance import Substance

try:
    from .config import parameters as params
except ImportError:
    from config import parameters as params


class Mixing:
    """Смешение"""

    __slots__ = "outlet"

    def __init__(self, *inlet_substances) -> None:
        """Расчет выходных параметров смешения"""
        self.outlet = Substance("outlet")

        for inlet_substance in inlet_substances:
            assert isinstance(inlet_substance, Substance), TypeError(f"{type(inlet_substance)=} not {type(Substance)}")

        self.outlet.parameters[params.mf] = sum(inlet_substance.parameters[params.mf] for inlet_substance in inlet_substances)


if __name__ == "__main__":
    substance_inlet_1 = Substance(
        "air",
        parameters={
            params.mf: 100,
        },
    )
    substance_inlet_2 = Substance(
        "exhaust",
        parameters={
            params.mf: 50,
        },
    )
    substance_inlet_3 = Substance(
        "air",
        parameters={
            params.mf: 3,
        },
    )

    m = Mixing(substance_inlet_1, substance_inlet_2, substance_inlet_3)

    for k, v in m.outlet.parameters.items():
        print(f"{k:<10}: {v}")
