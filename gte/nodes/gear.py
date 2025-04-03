from numpy import nan
#from thermodynamics import T0, Cp, R_gas, η_polytropic
#from tools import isnum, eps


class Gear:
    """Коробка приводов"""

    def __init__(self, name='Gear'):
        self.name = name
        self.type = ''
        self.η = nan  # КПД

        self.warnings = set()  # предупреждения

    def get_variability(self):
        return len(self.η) if type(self.η) is list and len(self.η) else 1

    def set_combination(self, combination, gear_main):
        positions = [0]
        for i in range(combination):
            if positions[0] == len(gear_main.η) - 1:
                positions[0] = 0
            else:
                positions[0] += 1
                continue

        self.η = gear_main.η[positions[0]]

    def solve(self, **kwargs):
        pass
