import pytest

from src.nodes import Compressor


class TestNode:
    @pytest.mark.parametrize("Node", [Compressor])
    def test_init(self, Node):
        node = Node()
        assert isinstance(node, Node)

    @pytest.mark.parametrize("Node", [Compressor])
    def test_name(self, Node):
        node = Node()
        assert node.name == Node.__name__
        node.name = "C"
        assert node.name == "C"
        node = Node("123")
        assert node.name == "123"
