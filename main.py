"""
Simulation zur Kármanschen Wirbelstraße - Projekt Numerik
von Louis B., Linus D., Gustav H und Maximilian M.

Informationen zur Simulation, den Algorithmen und den verwendenten externen Bibliotheken sind in der 'README.md' zu finden.
"""


from src.simulation import VortexSimulation


def main():
    Re = 100  # Reynolds number, working point: 100 - 150

    n = 256  # grid points in ξ-direction
    m = 256  # grid points in θ-direction

    simulation_time = 180  # total simulation time in seconds

    # Max dt for: m=n=160, dt=0.038; m=n=256, dt=0.012
    simulation = VortexSimulation(Re, n, m, dt=0.012)
    simulation.run_simulation(simulation_time, steps_per_frame=5, to_plot=True)


if __name__ == "__main__":
    main()
