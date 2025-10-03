try:
    from .compressor import Compressor
    from .turbine import Turbine
except ImportError:
    from compressor import Compressor
    from turbine import Turbine


class Shaft:
    """Вал"""

    __slots__ = ("name", "nodes")

    def __init__(self, name: str = "Shaft", nodes: tuple = tuple()):
        assert isinstance(name, str), TypeError(f"type name must be str not {type(name)}")
        self.name: str = name

        assert isinstance(nodes, (tuple, list)), TypeError("type nodes must be tuple")
        for node in nodes:
            assert isinstance(node, (Compressor, Turbine)), TypeError(f"type node must be {type(Compressor), type(Turbine)}")
        self.nodes: tuple = nodes


if __name__ == "__main__":
    s = Shaft("1", nodes=(Compressor(), Turbine()))
