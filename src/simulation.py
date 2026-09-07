from .solver import PoissonSolver
from .grid import LogPolarGrid

import matplotlib.pyplot as plt
import matplotlib.animation as animation

import numpy as np


class VortexSimulation:
    def __init__(self, Re: float, n: int, m: int, dt: float) -> None:
        """
        This class manages the simulation process of the Kármán vortex street. The class uses the 'run_simulation'-method to start the
        simulation.

        By initializing the class all args get saved as class variables, the 'LogPolarGrid' class gets initialized and the grid constructed,
        and the 'PoissonSolver' gets initialized.
        The initial conditions of stream function, vorticity, and speeds get set.

        Args:
            Re: Reynolds number
            n: number of grid points in ξ-direction
            m: number of grid points in θ-direction
            dt: step size in seconds
        """
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
        self.omega += 0.05 * np.sin(self.grid_theta) * np.exp(-2 * self.grid_xi)

        self.u_theta = np.zeros((n, m))
        self.u_xi = np.zeros((n, m))

        self._compute_vortex_boundaries()

    def run_simulation(self, simulation_time: float, steps_per_frame: int, to_plot=False) -> None | tuple[list[float], list[float]]:
        """
        This method is the callable method used to run the simulation. 
        When to_plot is True it starts the live simulation using matplotlib.animation.animation.FuncAnimation otherwise it calls
        the '_steps'-method without animation.

        Args:
            simulation_time: total simulation time in seconds
            steps_per_frame: number of steps taken per frame update
            to_plot: decides if the animation is shown or not

        Returns:
            A list of the stream function and voriticty near the cylinder to calculate the Strouhal-Number. Used with the 'notebooks/strouhal_number.ipynb' notebook
        """
        n_steps = int(np.ceil(simulation_time // self.dt))

        if to_plot:
            fig, ax = plt.subplots(figsize=(8, 6))
            self._init_plot(ax, simulation_time)
            anim = animation.FuncAnimation(
                fig,
                self._plot,
                fargs=(steps_per_frame, simulation_time),
                frames=n_steps // steps_per_frame + 1,
                interval=20,
                blit=False,
                repeat=False
            )
            self._anim = anim
            plt.show()
        else:
            # 'psi_for_st' and 'omega_for_st' are used to calculate the Strouhal number, both values get taken at θ = 0 and ξ = ln(r = 3)
            psi_for_st = []
            omega_for_st = []
            xi_st = int(np.log(3) / self.h)

            max_dt = 0
            for _ in range(n_steps):
                self.omega = self._step()

                this_dt = self.h / max(np.max(self.u_xi), np.max(self.u_theta))
                max_dt = this_dt if this_dt > max_dt else max_dt

                psi_for_st.append(self.psi[xi_st, 0])
                omega_for_st.append(self.omega[xi_st, 0])
            print(f"The maximum dt according to the CFL criterium for this grid size: {max_dt:.8f}")
            return omega_for_st, psi_for_st

    def _step(self) -> np.ndarray:
        """
        This method calls the function for each step inside the loop:

        time-step: solve poisson equation -> update velocities -> compute vortex boundaries -> solve vorticity equation

        Returns:
            The new calculated omega of this time step.
        """
        self.psi = self.solver.solve(self.omega)

        self._update_velocities()

        self._compute_vortex_boundaries()

        omega = self._runge_kutta()

        return omega

    def _update_velocities(self) -> None:
        """
        This methods updates the vortex velocities u_ξ = exp(-ξ) ⋅ ∂_θ ⋅ ψ and u_θ = exp(-ξ) ⋅ ∂_ξ ⋅ ψ by using a first-order 
        central derivative. 
        """
        self.u_xi[:, 1:-1] = ((self.psi[:, 2:] - self.psi[:, :-2])
                              * self.exp_neg_xi[:, 1:-1] / (2 * self.h))
        self.u_theta[1:-1, :] = (-(self.psi[2:, :] - self.psi[:-2, :])
                                 * self.exp_neg_xi[1:-1, :] / (2 * self.h))

    def _compute_vortex_boundaries(self, omega: np.ndarray | None = None) -> None:
        """
        Sets the boundary conditions for the stream function, voricity and vortex speed.

        stream function ψ: ψ_(n - 1),j = exp(ξ_(n - 1)) ⋅ sin(θ_j), ψ_0,j = 0.0
        vorticity ω: ω_(n - 1),j = 0.0, ω_0,j = 1 / (2 ⋅ h²) ⋅ (ψ_2,j - 8 ⋅ ψ_1,j)
        tangential vortex speed u_ξ: u_ξ_0,j = 0.0

        Args:
            omega: vorticity of shape n x m
        """
        self.psi[-1, :] = np.exp(self.grid_xi[-1, :]) * np.sin(self.grid_theta[-1, :])
        self.psi[0, :] = 0.0

        self.omega[-1, :] = 0.0
        self.omega[0, :] = 1 / (2 * self.h**2) * (self.psi[2, :] - 8 * self.psi[1, :])
        self.omega[:, 0] = self.omega[:, -2]
        self.omega[:, -1] = self.omega[:, 1]

        self.u_theta[0, :] = 0.0

    def _transport_vortex(self, omega: np.ndarray) -> np.ndarray:
        """
        This method solves the vorticity equation for ω_dot by calculating the convection and diffusion term.
        The diffusion term gets calculated by a second order central derivative on the vorticity in ξ and θ-direction.
        The first derivative in the convection term gets calculated by a third degree upwind scheme in θ and first degree upwind scheme 
        in ξ.

        Args:
            omega: vorticity of shape n x m

        Returns:
            The time derivative of the vorticity calculated by the vorticity equation.
        """
        omega = omega.copy()
        omega[0, :] = 1 / (2 * self.h**2) * (self.psi[2, :] - 8 * self.psi[1, :])
        omega[-1, :] = 0.0
        omega[:, 0] = omega[:, -2]
        omega[:, -1] = omega[:, 1]

        e2 = self.exp_neg_2_xi[1:-1, 1:-1]

        diffusion = (2 * e2 / (self.Re * self.h**2)) * (
            omega[2:, 1:-1] + omega[:-2, 1:-1]
            + omega[1:-1, 2:] + omega[1:-1, :-2] - 4 * omega[1:-1, 1:-1]
        )
        convection = (e2 / (4 * self.h**2)) * (
            (self.psi[2:, 1:-1] - self.psi[:-2, 1:-1])
            * (omega[1:-1, 2:] - omega[1:-1, :-2])
            - (self.psi[1:-1, 2:] - self.psi[1:-1, :-2])
            * (omega[2:, 1:-1] - omega[:-2, 1:-1])
        )

        omega_dot = np.zeros(omega.shape)
        omega_dot[1:-1, 1:-1] = diffusion + convection

        return omega_dot

    def _runge_kutta(self) -> np.ndarray:
        """
        This method applies a time step on the vorticity by using a Runge-Kutta scheme of 4th order.

        Runge-Kutta scheme of 4th order: y_k+1 = y_k + h / 6 ⋅ (k_1 + 2 ⋅ k_2 + 2 ⋅ k_3 + k_4) with:

        k_1 = f(t_n, y_n); 
        k_2 = f(t_n + h / 2, y_n + k_1 ⋅ h / 2); 
        k_3 = f(t_n + h / 2, y_n + k_2 ⋅ h / 2); 
        k_3 = f(t_n + h, y_n + k_3 ⋅ h)

        Notation: https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods

        Returns:
            The updated vorticty of this time step.
        """
        k1 = self._transport_vortex(self.omega)

        omega_k2 = self.omega + 0.5 * self.dt * k1
        k2 = self._transport_vortex(omega_k2)

        omega_k3 = self.omega + 0.5 * self.dt * k2
        k3 = self._transport_vortex(omega_k3)

        omega_k4 = self.omega + self.dt * k3
        k4 = self._transport_vortex(omega_k4)

        omega = self.omega + 1/6 * self.dt * (k1 + 2 * k2 + 2 * k3 + k4)

        return omega

    def _plot(self, frame: int, steps_per_frame: int, simulation_time: float) -> None:
        """
        This method updates the animation. It updates the animation frame every few time steps.

        Args:
            frame: the current frame of simulation
            steps_per_frame: steps taken per frame update
            simulation_time: total simulation time in seconds
        """
        for _ in range(steps_per_frame):
            self.omega = self._step()

        current_time = (frame + 1) * steps_per_frame * self.dt

        self.sim_title.set_text(
            f"Total Simulation Time: {round(simulation_time, 2)}s | t = {current_time:.2f}s"
        )

        self.mesh.set_array(self.omega.ravel())
        return self.mesh, self.sim_title

    def _init_plot(self, ax, simulation_time: float) -> None:
        """
        Initializes the plot for the simulation.
        Args:
            ax: The axis of the figure to plot on.
            simulation_time: The total simulation time in seconds.
        """
        r = np.exp(self.grid_xi)
        self.x_grid = r * np.cos(self.grid_theta)
        self.y_grid = r * np.sin(self.grid_theta)

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
