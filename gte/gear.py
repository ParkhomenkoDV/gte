from numpy import nan


class Gear:
    """Коробка приводов"""

    def __init__(self, name: str = "Gear"):
        self.name = name
        self.type = ""
        self.η = nan  # КПД
