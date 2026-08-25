from simulation import VortexSimulation


def main():
    Re = 10

    n = 128
    m = 64

    r_psi = 1.924  # for n = 128, m = 64

    simulation = VortexSimulation(Re, n, m, dt=0.005)
    simulation.run_simulation(simulation_time=30, r_psi=r_psi, steps_per_frame=2, to_plot=True)


if __name__ == "__main__":
    main()
