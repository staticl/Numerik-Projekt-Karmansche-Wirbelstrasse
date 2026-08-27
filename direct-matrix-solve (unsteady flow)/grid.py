import numpy as np


class LogPolarGrid:
    
    def __init__(self, n: int, m: int) -> None:
        """
        This class is used to build the math grid in logarithmic polar coordinates.

        Initializing a class object computes the grid point distance using the maximum θ value of 2π. The class methods are 'grid' which
        assembles the grid and 'factors' which calculates factors only dependent on the grid.

        Args:
            n: number of grid points in ξ-direction
            m: number of grid points in θ-direction
        """
        self.n = n
        self.m = m
        self.h = 2 * np.pi / self.m

    def grid(self) -> tuple[np.ndarray, np.ndarray]:
        """
        Assembles the logarithmic polar grid with constant point width in both ξ and θ-direction. 
        θ = [0,2π), ξ=[0,ξ_max] with ξ_max = 2π(n - 1) / m

        Returns:
            A tuple of the ξ and θ-grid both with shape n x m. 
        """
        xi_max = (self.n - 1) * self.h
        xi = np.linspace(0, xi_max, self.n)
        print(len(xi))
        theta = np.arange(self.m) * self.h
        print(len(theta))
        grid_xi, grid_theta = np.meshgrid(xi, theta, indexing='ij', sparse=False)

        self.grid_xi = grid_xi
        return grid_xi, grid_theta

    def factors(self) -> tuple[float, np.ndarray, np.ndarray]:
        """
        Calculates the factors 'exp(-2 ⋅ ξ_i,j)' and 'exp(-ξ_i,j)' using the ξ-grid.

        Returns:
            A tuple with grid point width 'h' and both factors with shape n x m.
        """
        exp_neg_2_xi = np.exp(-2 * self.grid_xi)
        exp_neg_xi = np.exp(-self.grid_xi)

        return self.h, exp_neg_2_xi, exp_neg_xi
