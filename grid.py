import numpy as np


class LogPolarGrid:
    def __init__(self, m: int, n: int) -> None:
        self.n = n
        self.m = m
        self.h = np.pi / (m - 1)

    def grid(self) -> np.ndarray:
        theta = np.linspace(0, np.pi, self.m)

        xi_max = (self.n - 1) * self.h
        xi = np.linspace(0, xi_max, self.n)

        grid = np.array(np.meshgrid(xi, theta, indexing='ij', sparse=False))

        self.xi = xi
        return grid

    def factors(self) -> tuple[float, np.ndarray, np.ndarray]:
        exp_neg_2_xi = np.exp(-2 * self.xi)
        exp_neg_xi = np.exp(-self.xi)

        return self.h, exp_neg_2_xi, exp_neg_xi
