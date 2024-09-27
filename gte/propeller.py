from numpy import nan
import matplotlib as plt


class Propeller:
    """Пропеллер"""

    def __init__(self):
        self.N = nan  # мощность

        self.warnings = set()  # предупреждения

    def get_variability(self):
        return len(self.N) if type(self.N) is list and len(self.N) else 1

    def set_combination(self, combination, propeller_main):
        positions = [0] * 1

        for i in range(combination):
            if positions[0] == len(propeller_main.N) - 1:
                positions[0] = 0
            else:
                positions[0] += 1
                continue

        if type(propeller_main.N) is list and len(propeller_main.N): self.N = propeller_main.N[positions[0]]

    def solve(self, how='all', error=0.01, Niter=100, **kwargs):
        pass
