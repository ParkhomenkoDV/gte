try:
    from ...turbocompressor.rotor import Rotor
    from ...turbocompressor.stator import Stator
except ImportError:
    import os
    import sys

    sys.path.insert(0, os.getcwd())

    from gte.nodes.turbocompressor.rotor import Rotor
    from gte.nodes.turbocompressor.stator import Stator


class Compressor(Rotor):
    """Компрессор"""

    def __init__(self, *stages):
        for stage in stages:
            if not isinstance(stage, (Rotor, Stator)):
                raise TypeError
