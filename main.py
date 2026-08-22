from simulation import VortexSimulation


def main():
    Re = 10

    n = 64
    m = 32

    # r_psi = 1.924 # for n = 128, m = 64
    r_psi = 1.852  # for n = 64, m = 32

    simulation = VortexSimulation(Re, n, m, dt=0.02)
    simulation.run_simulation(simulation_time=3, r_psi=r_psi, steps_per_frame=10, to_plot=False)


if __name__ == "__main__":
    main()
