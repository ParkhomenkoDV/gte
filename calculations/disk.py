import numpy as np
import matplotlib.pyplot as plt


class Dick:
    def __init__(self, radius, thickness):
        assert type(radius) in (tuple, list, np.array)
        assert type(thickness) in (tuple, list, np.array)
        assert len(radius) == len(thickness)
        assert all(map(lambda r: type(r) in (int, float), radius))
        assert all(map(lambda t: type(t) in (int, float), thickness))
        assert all(map(lambda r: r >= 0, radius))
        assert all(map(lambda t: t >= 0, thickness))
        assert all(radius[i] < radius[i + 1] for i in range(len(radius) - 1))

        self.radius, self.thickness = np.array(radius), np.array(thickness)

    def calculation1(self):
        return

    def calculation2(self):
        return

    def sigma(self):
        self.calculation1()
        self.calculation2()

    def show(self, **kwargs) -> None:
        plt.figure(figsize=kwargs.pop('figsize', (8, 8)))
        plt.title("Disk", fontsize=14, fontweight='bold')
        plt.plot([-self.thickness[0] / 1.5, self.thickness[0] / 1.5], [0, 0],
                 color='orange', linestyle='dashdot', linewidth=1.5)
        plt.plot(self.thickness / 2, self.radius, color='black', linestyle='solid', linewidth=3)
        plt.plot(-self.thickness / 2, self.radius, color='black', linestyle='solid', linewidth=3)
        plt.plot([-self.thickness[-1] / 2, self.thickness[-1] / 2], [self.radius[-1], self.radius[-1]],
                 color='black', linestyle='solid', linewidth=3)
        if self.radius[0] > 0:
            plt.plot([-self.thickness[0] / 2, self.thickness[0] / 2], [self.radius[0], self.radius[0]],
                     color='black', linestyle='solid', linewidth=3)
        plt.grid(True)
        plt.axis('equal')
        plt.xlabel("Thickness", fontsize=12)
        plt.ylabel("Radius", fontsize=12)

        plt.show()

    def show_sigma(self, **kwargs) -> None:
        plt.figure(figsize=kwargs.pop('figsize', (16, 8)))

        plt.title('Sigma', fontsize=14, fontweight='bold')
        plt.grid(True)
        plt.axis('equal')
        plt.xlabel("$Sigma", fontsize=12)
        plt.ylabel("Radius", fontsize=12)

        plt.show()


if __name__ == "__main__":
    disk = Dick(radius=[0.2, 0.3, 0.7, 0.8],
                thickness=[0.6, 0.6, 0.2, 0.3])

    disk.show()
    disk.show_sigma()
