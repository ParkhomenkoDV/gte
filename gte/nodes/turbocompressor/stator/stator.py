try:
    from ...node import Node
except ImportError:
    import os
    import sys

    sys.path.insert(0, os.getcwd())

    from gte.nodes.node import Node


class Stator:
    """Статор"""

    def __init__(self, parameters: dict[str, float], name: str = "Stator"):
        Node.__init__(self, parameters, name)
