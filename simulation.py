from solver import PoissonSolver
from grid import LogPolarGrid

import numpy as np


class VortexSimulation:
    def __init__(self, Re: float, m: int, n: int) -> None:
        self.Re = Re

        lp_grid = LogPolarGrid(n, m)
        self.grid = lp_grid.grid()
        self.h, self.exp_neg_2_xi, self.exp_neg_xi = lp_grid.factors()

        self.psi = np.zeros((n, m))
        self.omega = np.zeros((n, m))

        self.u_theta = np.zeros((n, m))
        self.u_xi = np.zeros((n, m))

    def _update_velocities(self) -> None:
        self.u_xi[:, 1:-1] = self.exp_neg_xi[:, 1:-1] * (self.psi[:, 2:] - self.psi[:, :-2]) / (2 * self.h)
        self.u_theta[1:-1, :] = -self.exp_neg_xi[1:-1, :] * (self.psi[2:, :] - self.psi[:-2, :]) / (2 * self.h)

    def _compute_vortex_boundaries(self):
        pass

    def _transport_vortex(self) -> np.ndarray:
        s = self.psi.shape
        d2omega_dtheta2 = np.zeros(s)
        d2omega_dtheta2[:, 1:-1] = (self.omega[:, 2:] - 2 * self.omega[:, 1:-1] + self.omega[:, :-2]) / (self.h**2)
        d2omega_dxi2 = np.zeros(s)
        d2omega_dxi2[1:-1, :] = (self.omega[2:, :] - 2 * self.omega[1:-1, :] + self.omega[:-2, :]) / (self.h**2)

        domega_dtheta = np.zeros(s)
        domega_dtheta[:, 1:-1] = (self.omega[:, 2:] - self.omega[:, :-2]) / (2 * self.h)
        domega_dxi = np.zeros(s)
        domega_dxi[1:-1, :] = (self.omega[2:, :] - self.omega[:-2, :]) / (2 * self.h)

        omega_dot = 2/self.Re * self.exp_neg_2_xi * \
            (d2omega_dtheta2 + d2omega_dxi2) - self.exp_neg_xi * (self.u_theta * domega_dxi + self.u_xi * domega_dtheta)

        return omega_dot

    def _step(self, r_psi: float, r_omega: float) -> np.ndarray:
        solver = PoissonSolver(self.Re, self.h, self.exp_neg_2_xi, self.exp_neg_xi)
        solver.sor(self.psi, self.omega, 1e-8, r_psi, r_omega)

        self._update_velocities()

        self._compute_vortex_boundaries()

        omega_dot = self._transport_vortex()

        return omega_dot

    def _plot(self) -> None:
        pass

    def run_simulation(self, simulation_time: float, dt: float, r_psi: float, r_omega: float) -> None:
        n_steps = int(simulation_time // dt)

        for _ in range(n_steps):
            omega_dot = self._step(r_psi, r_omega)
            self.omega += omega_dot * dt
