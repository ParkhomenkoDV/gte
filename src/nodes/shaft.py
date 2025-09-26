from src.nodes.compressor import Compressor
from src.nodes.node import GTENode
from src.nodes.turbine import Turbine


class Shaft:
    """Вал"""

    __slots__ = ("name", "nodes")

    def __init__(self, name: str = "Shaft", nodes: tuple = tuple()):
        assert isinstance(name, str), TypeError(f"type name must be str not {type(name)}")
        self.name: str = name

        assert isinstance(nodes, (tuple, list)), TypeError("type nodes must be tuple")
        for node in nodes:
            assert isinstance(node, GTENode), TypeError(f"type node must be {type(GTENode)}")
        self.nodes: tuple = nodes


if __name__ == "__main__":
    s = Shaft("1", nodes=(Compressor(), Turbine()))
