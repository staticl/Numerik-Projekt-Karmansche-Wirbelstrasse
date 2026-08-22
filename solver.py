import numpy as np

import matplotlib.pyplot as plt

from scipy.sparse import eye, kron, diags
from scipy.sparse.linalg import cg, spsolve


class PoissonSolver:
    def __init__(self, h: float, exp_neg_2_xi: np.ndarray) -> None:
        self.h = h
        self.exp_pos_2_xi = 1 / exp_neg_2_xi

    def sor(self, psi: np.ndarray, omega: np.ndarray, toleranz: float, max_iterations: int, r_psi: float) -> None:
        # errors = []

        if psi.shape != omega.shape:
            raise ValueError(f'Die Matrizen des Wirbelstroms und der Wirbelstärke müssen die gleiche Shape haben. Shape-Psi: {psi.shape},\
                                Shape-Omega: {omega.shape}')

        nx, ny = psi.shape

        f = self.h**2 * self.exp_pos_2_xi[1:-1, 1:-1] * omega[1:-1, 1:-1]

        grid_y, grid_x = np.ogrid[1:nx-1, 1:ny-1]
        w_mask = ((grid_y + grid_x) % 2 == 0)
        b_mask = ~w_mask

        inner_psi = psi[1:-1, 1:-1]

        for iter in range(max_iterations):
            psi_iter = psi.copy()

            w = psi[2:, 1:-1] + psi[:-2, 1:-1] + psi[1:-1, 2:] + psi[1:-1, :-2]
            term = 1/4 * (w + f)
            inner_psi[w_mask] = (1 - r_psi) * inner_psi[w_mask] + r_psi * term[w_mask]

            b = psi[2:, 1:-1] + psi[:-2, 1:-1] + psi[1:-1, 2:] + psi[1:-1, :-2]
            term = 1/4 * (b + f)
            inner_psi[b_mask] = (1 - r_psi) * inner_psi[b_mask] + r_psi * term[b_mask]

            error_psi = np.max(np.abs(psi - psi_iter))
            # errors.append(error_psi)
            if error_psi < toleranz:
                # self.plot_error(errors)
                return
        # self.plot_error(errors)
        return

    def plot_error(self, errors: list[float]) -> None:
        fig, axis = plt.subplots(1, 1, figsize=(10, 5))
        x = np.arange(len(errors)) + 1
        axis.plot(x, errors, color='blue')
        axis.set_xlabel('Number of Iteration')
        axis.set_ylabel('Error |psi - old_psi|')
        axis.set_title('Change of the Error')

        plt.show()
