from numpy import nan


class Load:
    """Нагрузка"""

    def __init__(self):
        self.N = nan  # мощность

        self.warnings = set()  # предупреждения

    def get_variability(self):
        return len(self.N) if type(self.N) is list and len(self.N) else 1

    def set_combination(self, combination, load_main):
        positions = [0]
        for i in range(combination):
            if positions[0] == len(load_main.N) - 1:
                positions[0] = 0
            else:
                positions[0] += 1
                continue

        self.N = load_main.N[positions[0]]

    def solve(self, **kwargs):
        pass
