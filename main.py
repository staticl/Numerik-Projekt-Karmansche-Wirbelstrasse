from simulation import VortexSimulation


def main():
    Re = 100

    m = 100
    n = 100

    r_psi = 0.5
    r_omega = 0.5

    s = VortexSimulation(Re, m, n)
    s.run_simulation(simulation_time=10, dt=0.001, r_psi=r_psi, r_omega=r_omega)


if __name__ == "__main__":
    main()
