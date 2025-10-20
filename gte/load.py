from numpy import nan

try:
    from .checks import check_efficiency
    from .config import parameters as gtep


except ImportError:
    from checks import check_efficiency
    from config import parameters as gtep


class Load:
    """Нагрузка"""

    __slots__ = (gtep.efficiency, gtep.power)

    def __init__(self):
        setattr(self, gtep.efficiency, nan)
        setattr(self, gtep.power, nan)

    @property
    def is_real(self):
        checks = (check_efficiency(getattr(self, gtep.efficiency)),)
        return all(checks)


if __name__ == "__main__":
    from colorama import Fore

    load = Load()
    setattr(load, gtep.efficiency, 0.98)
    setattr(load, gtep.power, 32 * 10**6)

    print(Fore.GREEN + f"{load.is_real = }" + Fore.RESET)
