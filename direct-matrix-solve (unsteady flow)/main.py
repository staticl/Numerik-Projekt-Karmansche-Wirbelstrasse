from simulation import VortexSimulation


def main():
    Re = 1000

    n = 160
    m = 160

    simulation = VortexSimulation(Re, n, m, dt=0.01)
    simulation.run_simulation(simulation_time=60, steps_per_frame=5, to_plot=True)


if __name__ == "__main__":
    main()