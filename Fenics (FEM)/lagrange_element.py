import matplotlib.pyplot as plt
import numpy as np


def f(x):
    return np.sin(np.pi * x)


x = np.linspace(0, 1, 100)
grid = np.linspace(0, 1, 5)

for i in range(0, len(grid) - 1):
    plt.plot([grid[i], grid[i + 1]], [0, 1], color="red")
    plt.plot([grid[i], grid[i + 1]], [1, 0], color="red")

y = f(x)

plt.figure(1)
plt.plot(x, y)
plt.scatter(grid, np.zeros_like(grid), color="green")
plt.scatter(grid, np.ones_like(grid), color="green")

plt.xlabel("x")
plt.ylabel("u")
plt.xticks(np.arange(0, 1.1, 0.1))
plt.yticks(np.arange(0, 1.1, 0.1))
plt.grid()
plt.savefig("lagrange elements.png")

plt.figure(2)
plt.plot(x, y)
for i in range(0, len(grid) - 1):
    plt.scatter(grid[i], f(grid[i]), color="green")
    plt.plot(
        [grid[i], grid[i + 1]], [0, f(grid[i + 1])], color="orange", linestyle="--"
    )
    plt.scatter(grid[i + 1], 0, color="green")
    plt.plot([grid[i + 1], grid[i]], [0, f(grid[i])], color="purple", linestyle="--")

    plt.plot([grid[i], grid[i + 1]], [f(grid[i]), f(grid[i + 1])], color="red")

plt.xlabel("x")
plt.ylabel("u")
plt.xticks(np.arange(0, 1.1, 0.1))
plt.yticks(np.arange(0, 1.1, 0.1))
plt.grid()
plt.savefig("FEM expected result.png")

plt.show()
