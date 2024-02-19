from numpy import linspace, nan, isnan

class Blade:
    """Лопатка или лопасть"""

    def __init__(self):
        self.material = ''
        self.height = 0
        self.amount = 1

    def solve(self, omega, *args, **kwargs):
        pass


if __name__ == '__main__':
    blade = Blade()
