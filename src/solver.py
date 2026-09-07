import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


class PoissonSolver:
    def __init__(self, h: float, exp_neg_2_xi: np.ndarray, grid_xi: np.ndarray, grid_theta: np.ndarray) -> None:
        """
        This class solves the poisson equation ∇²ψ = -ω by solving the matrix equation 'Ax = b'. It uses the class method 'solve' for that.

        By initializing the class all arguments get saved as class varibales and the '_build'-method gets called, which assembles the 
        matrix A.

        Args:
            h: grid point width
            exp_neg_2_xi: factor 'exp(-2 ⋅ ξ_i,j)'
            grid_xi: meshgrid of ξ
            grid_theta: meshgrid of θ
        """
        self.h = h
        self.exp_pos_2_xi = 1 / exp_neg_2_xi

        self.grid_xi = grid_xi
        self.grid_theta = grid_theta

        self.n, self.m = self.grid_xi.shape

        self.lu = self._build()

    def _build(self) -> np.ndarray:
        """
        This method assembles the matrix A of shape m⋅n x m⋅n, used to solve for the stream function, based on the boundary 
        conditions and the FDE.

        The Finite-Difference-Equation in this case is:

        4 ⋅ ψ_i,j - ψ_i+1,j - ψ_i-1,j - ψ_i,j+1 - ψ_i,j-1 = h² ⋅ exp(2⋅ξ_i) ⋅ ω_i,j

        Returns:
            The LU decomposition of the matrix saved as a csc-sparse-array.
        """
        T_x = sp.diags([-1, 2, -1], [-1, 0, 1], shape=(self.m, self.m), format='csc')
        T_y = sp.diags([-1, 2, -1], [-1, 0, 1], shape=(self.n, self.n), format='csc')

        I_n = sp.eye(self.n, format='csc')
        I_m = sp.eye(self.m, format='csc')

        A = (sp.kron(T_x, I_n, format='csc') + sp.kron(I_m, T_y, format='csc')).tolil()

        k_L = np.arange(0, (self.m - 1) * self.n, self.n)
        A[k_L, :] = 0.0
        A[k_L, k_L] = 1.0

        k_R = np.arange(self.n - 1, self.m * self.n, self.n)
        A[k_R, :] = 0.0
        A[k_R, k_R] = 1.0

        k_B = np.arange(self.n)
        k_be = np.arange((self.m - 2) * self.n, (self.m - 1) * self.n)
        A[k_B, :] = 0.0
        A[k_B, k_B] = 1.0
        A[k_B, k_be] = -1.0

        k_T = np.arange((self.m - 1) * self.n, self.m * self.n)
        k_te = np.arange(self.n, 2 * self.n)
        A[k_T, :] = 0.0
        A[k_T, k_T] = 1.0
        A[k_T, k_te] = -1.0

        self.k_L, self.k_R, self.k_B, self.k_T = k_L, k_R, k_B, k_T

        return spla.splu(A.tocsc())

    def solve(self, omega: np.ndarray) -> np.ndarray:
        """
        This method first assembles the the right side of the equation b and applies the periodic boundary conditions on it. After it
        solves the equation for the stream function.

        Args:
            omega: vorticity of shape n x m

        Returns:
            The stream function ψ with shape n x m.
        """
        b = (self.h**2 * self.exp_pos_2_xi * omega).T.flatten()

        b[self.k_L] = 0.0
        b[self.k_R] = np.exp(self.grid_xi[-1, :]) * np.sin(self.grid_theta[-1, :])
        b[self.k_B] = 0.0
        b[self.k_T] = 0.0

        psi = self.lu.solve(b)

        return psi.reshape(self.m, self.n).T
