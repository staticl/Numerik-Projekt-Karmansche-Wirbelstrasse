import numpy as np


class LogPolarGrid:
    def __init__(self, n: int, m: int) -> None:
        self.n = n
        self.m = m
        self.h = np.pi / (m - 1)

    def grid(self) -> np.ndarray:
        xi_max = (self.n - 1) * self.h
        xi = np.linspace(0, xi_max, self.n)

        theta = np.linspace(0, np.pi, self.m)

        grid_xi, grid_theta = np.meshgrid(xi, theta, indexing='ij', sparse=False)

        self.grid_xi = grid_xi
        return grid_xi, grid_theta

    def factors(self) -> tuple[float, np.ndarray, np.ndarray]:
        exp_neg_2_xi = np.exp(-2 * self.grid_xi)
        exp_neg_xi = np.exp(-self.grid_xi)

        return self.h, exp_neg_2_xi, exp_neg_xi
