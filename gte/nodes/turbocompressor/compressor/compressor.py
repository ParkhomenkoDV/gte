try:
    from ....config.config import parameters as gtep
    from ...turbocompressor.rotor import Rotor
    from ...turbocompressor.stator import Stator
except ImportError:
    import os
    import sys

    sys.path.insert(0, os.getcwd())

    from gte.config.config import parameters as gtep
    from gte.nodes.turbocompressor.rotor import Rotor
    from gte.nodes.turbocompressor.stator import Stator


class Compressor:
    """Компрессор"""

    slots = ("blade_rows", gtep.pipi, gtep.titi, gtep.effeff)

    def __init__(self, *blade_rows):
        for blade_row in blade_rows:
            if not isinstance(blade_row, (Rotor, Stator)):
                raise TypeError()
