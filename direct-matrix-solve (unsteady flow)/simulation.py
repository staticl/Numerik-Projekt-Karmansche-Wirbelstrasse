from solver import PoissonSolver
from grid import LogPolarGrid

import matplotlib.pyplot as plt
import matplotlib.animation as animation

import numpy as np


class VortexSimulation:
    def __init__(self, Re: float, n: int, m: int, dt: float) -> None:
        self.Re = Re

        lp_grid = LogPolarGrid(n, m)
        self.grid_xi, self.grid_theta = lp_grid.grid()
        self.h, self.exp_neg_2_xi, self.exp_neg_xi = lp_grid.factors()

        self.dt = dt

        self.solver = PoissonSolver(self.h, self.exp_neg_2_xi, self.grid_xi, self.grid_theta)

        self.psi = np.exp(self.grid_xi) * np.sin(self.grid_theta)
        self.psi[0, :] = 0.0
        self.omega = np.zeros((n, m))
        self.omega += 0.05 * np.cos(self.grid_theta) * np.exp(-2 * self.grid_xi)

        self.u_theta = np.zeros((n, m))
        self.u_xi = np.zeros((n, m))

        self._compute_vortex_boundaries()

    def _update_velocities(self) -> None:
        self.u_xi[:, 1:-1] = self.psi[:, 2:] - self.psi[:, :-2]
        self.u_xi[:, 0] = self.psi[:, 1] - self.psi[:, -1]
        self.u_xi[:, -1] = self.psi[:, 0] - self.psi[:, -2]

        self.u_xi *= self.exp_neg_xi / (2 * self.h)

        self.u_theta[1:-1, :] = self.psi[2:, :] - self.psi[:-2, :]
        # self.u_theta[-1, :] = self.psi[0, :] - self.psi[-2, :]
        self.u_theta *= -self.exp_neg_xi / (2 * self.h)

        # print(f'Max-Xi: {np.max(self.u_xi)}, Min-Xi: {np.min(self.u_xi)}')
        # print(f'Max-Theta: {np.max(self.u_theta)}, Min-Theta: {np.min(self.u_theta)}')

    def _compute_vortex_boundaries(self, omega: np.ndarray | None = None):

        self.psi[-1, :] = np.exp(self.grid_xi[-1, :]) * np.sin(self.grid_theta[-1, :])
        self.psi[0, :] = 0.0

        self.omega[-1, :] = 0.0
        self.omega[0, :] = 1 / (2 * self.h**2) * (self.psi[2, :] - 8 * self.psi[1, :])

        self.u_theta[0, :] = 0.0

    def _transport_vortex(self, omega: np.ndarray) -> np.ndarray:
        s = self.psi.shape
        d2omega_dtheta2 = np.zeros(s)
        d2omega_dtheta2[:, 1:-1] = omega[:, 2:] - 2 * omega[:, 1:-1] + omega[:, :-2]
        d2omega_dtheta2[:, 0] = omega[:, 1] - 2 * omega[:, 0] + omega[:, -1]
        d2omega_dtheta2[:, -1] = omega[:, 0] - 2 * omega[:, -1] + omega[:, -2]
        d2omega_dtheta2 /= self.h**2

        d2omega_dxi2 = np.zeros(s)
        d2omega_dxi2[1:-1, :] = (omega[2:, :] - 2 * omega[1:-1, :] + omega[:-2, :]) / (self.h**2)
        d2omega_dxi2[1, :] = omega[0, :]

        domega_dtheta = np.zeros(s)
        domega_dxi = np.zeros(s)

        forward_theta = np.zeros(s)
        forward_theta[:, 1:-2] = -1 * (2 * omega[:, :-3] + 3 * omega[:, 1:-2] - 6 * omega[:, 2:-1] + omega[:, 3:])
        forward_theta[:, 0] = -1 * (2 * omega[:, -1] + 3 * omega[:, 0] - 6 * omega[:, 1] + omega[:, 2])
        forward_theta[:, -2] = -1 * (2 * omega[:, -3] + 3 * omega[:, -2] - 6 * omega[:, -1] + omega[:, 0])
        forward_theta[:, -1] = -1 * (2 * omega[:, -2] + 3 * omega[:, -1] - 6 * omega[:, 0] + omega[:, 1])
        forward_theta /= 6 * self.h

        backward_theta = np.zeros(s)
        backward_theta[:, 2:-1] = 2 * omega[:, 3:] + 3 * omega[:, 2:-1] - 6 * omega[:, 1:-2] + omega[:, :-3]
        backward_theta[:, 0] = 2 * omega[:, 1] + 3 * omega[:, 0] - 6 * omega[:, -1] + omega[:, -2]
        backward_theta[:, 1] = 2 * omega[:, 2] + 3 * omega[:, 1] - 6 * omega[:, 0] + omega[:, -1]
        backward_theta[:, -1] = 2 * omega[:, 0] + 3 * omega[:, -1] - 6 * omega[:, -2] + omega[:, -3]
        backward_theta /= 6 * self.h

        forward_xi = (omega[2:, :] - omega[1:-1, :]) / self.h
        backward_xi = (omega[1:-1, :] - omega[:-2, :]) / self.h

        domega_dtheta = np.where(self.u_theta > 0, backward_theta, forward_theta)
        domega_dxi[1:-1, :] = np.where(self.u_xi[1:-1, :] > 0, backward_xi, forward_xi)

        # print(f'd2omega_dxi2: {np.max(np.abs(d2omega_dxi2))}')
        # print(f'd2omega_dtheta2: {np.max(np.abs(d2omega_dtheta2))}')
        # print(f'domega_dxi: {np.max(np.abs(domega_dxi))}')
        # print(f'domega_dtheta: {np.max(np.abs(domega_dtheta))}')

        # print(f'omega: {omega.shape}')

        # print(f'd2omega_dxi2: {d2omega_dxi2.shape}')
        # print(f'd2omega_dtheta2: {d2omega_dtheta2.shape}')
        # print(f'domega_dxi: {domega_dxi.shape}')
        # print(f'domega_dtheta: {domega_dtheta.shape}')

        # print("exp_neg_xi shape:", self.exp_neg_xi.shape)

        diffusion = (1.0 / self.Re) * self.exp_neg_2_xi * (d2omega_dtheta2 + d2omega_dxi2)
        advection = self.exp_neg_xi * (self.u_xi * domega_dxi + self.u_theta * domega_dtheta)

        # print(f'DIFFUSION: {np.max(np.abs(diffusion))}')
        # print(f'ADVEKTION: {np.max(np.abs(advection))}')

        # print(f'diffusion: {diffusion.shape}')
        # print(f'advection: {advection.shape}')

        omega_dot = diffusion - advection

        return omega_dot

    def _runge_kutta(self) -> np.ndarray:
        k1 = self._transport_vortex(self.omega)
        # print(f'k1 - OMEGA: {np.max(np.abs(self.omega))}')
        omega_k2 = self.omega + 0.5 * self.dt * k1
        k2 = self._transport_vortex(omega_k2)
        # print(f'k2 - OMEGA: {np.max(np.abs(omega_k2))}')
        omega_k3 = self.omega + 0.5 * self.dt * k2
        k3 = self._transport_vortex(omega_k3)
        # print(f'k3 - OMEGA: {np.max(np.abs(omega_k3))}')
        omega_k4 = self.omega + self.dt * k3
        k4 = self._transport_vortex(omega_k4)
        # print(f'k4 - OMEGA: {np.max(np.abs(omega_k4))}')
        omega = self.omega + 1/6 * self.dt * (k1 + 2 * k2 + 2 * k3 + k4)

        return omega

    def _step(self) -> np.ndarray:
        self.psi = self.solver.solve(self.omega)
        # print(f'PSI: {np.abs(np.max(self.psi))}')
        self._update_velocities()

        self._compute_vortex_boundaries()

        omega = self._runge_kutta()
        
        return omega

    def init_plot(self, ax, simulation_time: float) -> None:
        r = np.exp(self.grid_xi)
        self.x_grid = r * np.cos(self.grid_theta)
        self.y_grid = r * np.sin(self.grid_theta)

        # Draw the initial mesh plot and save it to an instance variable
        # Fixing vmin and vmax is crucial for performance!
        self.mesh = ax.pcolormesh(
            self.x_grid,
            self.y_grid,
            self.omega,
            cmap="RdBu",
            shading="auto",
            vmin=-1,
            vmax=1,
        )

        ax.set_aspect("equal")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")

        ax.set_xlim(-5, 30)
        ax.set_ylim(-7, 7)

        self.sim_title = ax.set_title(
            f"Total Simulation Time: {round(simulation_time, 2)}s | t = 0.00s"
        )

    def _plot(self, frame: int, steps_per_frame: int, simulation_time: float):
        for _ in range(steps_per_frame):
            self.omega = self._step()

        current_time = (frame + 1) * steps_per_frame * self.dt

        self.sim_title.set_text(
            f"Total Simulation Time: {round(simulation_time, 2)}s | t = {current_time:.2f}s"
        )

        self.mesh.set_array(self.omega.ravel())
        return self.mesh, self.sim_title

    def run_simulation(self, simulation_time: float, steps_per_frame: int, to_plot=False) -> None:
        n_steps = int(np.ceil(simulation_time // self.dt))

        if to_plot:
            fig, ax = plt.subplots(figsize=(8, 6))
            self.init_plot(ax, simulation_time)
            anim = animation.FuncAnimation(
                fig,
                self._plot,
                fargs=(steps_per_frame, simulation_time),
                frames=n_steps // steps_per_frame + 1,
                interval=20,  # Delay between frames in ms (~50 FPS target)
                blit=False,  # Redraws only modified pixels for speed
                repeat=False
            )
            self._anim = anim
            plt.show()
        else:
            for _ in range(n_steps):
                self.omega = self._step()
                # print(f'Max: {np.max(self.psi)}, Min: {np.min(self.psi)}')
