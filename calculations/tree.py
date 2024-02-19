class Tree:
    """Древо"""

    def __init__(self, data):
        self.data = data
        self.child = []
        self.parent = None

    def add_child(self, child):
        child.parent = self
        self.child.append(child)

    def get_level(self):
        level = 0
        p = self.parent
        while p:
            level += 1
            p = p.parent
        return level

    def print_tree(self):
        spaces = ' ' * self.get_level() * 3
        prefix = spaces + '|__' if self.parent else ''
        print(prefix, self.data)
        if self.child:
            for child in self.child: child.print_tree()


if __name__ == '__main__':
    root = Tree('GTE')

    H = Tree('H')
    H.add_child(Tree({'temp': 1}))
    H.add_child(Tree({'temp': 2}))
    H.add_child(Tree({'temp': 3}))
    root.add_child(H)

    TT = Tree('TT')
    TT.add_child(Tree({'temp': 4}))
    TT.add_child(Tree({'temp': 5}))
    TT.add_child(Tree({'temp': 6}))
    root.add_child(TT)

    root.print_tree()
