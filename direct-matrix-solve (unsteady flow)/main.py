from simulation import VortexSimulation


"""
Simulation zur Kármanschen Wirbelstraße - Projekt Numerik
von Louis B., Linus D., Gustav H und Maximilian M.

Informationen zur Simulation und den Algorithmen sind in der 'readme.md' zu finden.

Informationen zu den verwendenten externen Bibliotheken sind in der 'requirements.txt' zu finden.
"""
def main():
    Re = 1000 # Reynolds number, working point: 100 - 150

    n = 160 # grid points in ξ-direction
    m = 160 # grid points in θ-direction

    simulation_time = 60 # total simulation time in seconds

    simulation = VortexSimulation(Re, n, m, dt=0.01)
    simulation.run_simulation(simulation_time, steps_per_frame=5, to_plot=True)


if __name__ == "__main__":
    main()