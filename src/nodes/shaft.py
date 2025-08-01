class Shaft:
    """Вал"""

    __slots__ = (
        "name",
        "nodes",
    )

    def __init__(self, name: str = "Shaft", nodes: tuple = None):
        self.name: str = name
        assert isinstance(nodes, tuple), TypeError("type nodes must be tuple")
        self.nodes: tuple = nodes or tuple()


if __name__ == "__main__":
    s = Shaft("1", nodes=(1, 2, 3))
