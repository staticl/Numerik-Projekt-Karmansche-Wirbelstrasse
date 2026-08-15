import numpy as np


class PoissonSolver:
    def __init__(self, Re: float, h: float, exp_neg_2_xi: np.ndarray, exp_neg_xi: np.ndarray) -> None:
        self.Re = Re

        self.h = h

        self.f1 = exp_neg_xi
        self.f2 = exp_neg_2_xi

    def sor(self, psi: np.ndarray, omega: np.ndarray, toleranz: float, max_iterations: int, r_psi: float, r_omega: float) -> int:
        n, m = psi.shape
        if psi.shape != omega.shape:
            raise ValueError(
                f'Die Matrizen des Wirbelstroms und der Wirbelstärke müssen die gleiche Shape haben.')
        for iter in max_iterations:
            psi_iter = psi.copy()
            omega_iter = omega.copy()
            for i in range(1, n - 1):
                for j in range(1, m - 1):
                    term_psi = r_psi/4 * (psi[i + 1, j] + psi[i - 1, j] + psi[i, j + 1] +
                                          psi[i, j - 1] + self.h**2 * self.f2 * omega[i, j])
                    psi[i, j] = (1 - r_psi) * psi[i, j] + term_psi

                    f = (psi[i + 1, j] - psi[j - 1, j]) * (omega[i, j + 1] - omega[i, j - 1]) + \
                        (psi[i, j + 1] - psi[i, j - 1]) * (omega[i + 1, j] - omega[j - 1, j])
                    term_omega = r_omega / 4 * \
                        (omega[i + 1, j] + omega[i - 1, j] + omega[i, j + 1] + omega[i, j - 1] + self.Re / 8 * f)
                    omega[i, j] = (1 - r_omega) * omega[i, j] + term_omega
            error_psi = np.max(np.abs(psi - psi_iter))
            error_omega = np.max(np.abs(omega - omega_iter))
            if error_psi < toleranz and error_omega < toleranz:
                return iter

        return iter
