from numpy import prod
from thermodynamics import η_polytropic, Cp, R_gas
from tools import isnum, eps
from colorama import Fore

def find_node_in_scheme(scheme, node2find) -> tuple:
    """Поиск положения узла в схеме"""
    for contour in scheme:
        for i, node in enumerate(scheme[contour]):
            if scheme[contour][i] is node2find: return contour, i

class Bearing:
    """Подшипник"""

    def __init__(self, bearing_type, name='Bearing'):
        self.name = name
        self.type = bearing_type  # тип подшипника
        self.warnings = set()  # предупреждения