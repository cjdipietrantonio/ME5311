# Christian DiPietrantonio
# ME 5311: Computational Methods to Viscous Flows
# Computer Assignment 01
# 02/19/2026

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import lu_factor, lu_solve

# ------------------------------------------------------------------
# 1. Parameters and Mesh Sizing
# ------------------------------------------------------------------

def set_parameters():
    params = {}

    # Domain in y
    params["y_min"] = 0
    params["y_max"] = 1
    params["N"]     = 200  # Number of cells (finite volume method)
    params["dy"]    = (params["y_max"] - params["y_min"]) / params["N"]  # Cell size in y

    # Stepping in x
    params["x_max"] = 0.5 # Final x
    params["Nx"]    = 1000 # Number of steps in x
    params["dx"]    = params["x_max"] / params["Nx"] # Step size in x

    # Source term
    params["S"] = 2.0

    # y-array for plotting and analytical solution evaluation
    params["y"] = np.linspace(params["y_min"], params["y_max"], params["N"])

    # Initial condition u(y,0) = 0
    u0 = np.zeros(params["N"])

    return params, u0

# ------------------------------------------------------------------
# 2. Build Matrices
# ------------------------------------------------------------------
   
def build_matrices(params):
    N = params["N"]
    dx = params["dx"]
    dy = params["dy"]

    r = dx / dy**2 # Defined to simplify matrix coefficients

    A = np.zeros((N, N))
    B = np.zeros((N, N))

    # ----- Bottom Boundary: i = 0, u(x,0) = 0 ---------------------
    #LHS (coefficients from Crank-Nicholson discretization)
    A[0, 0] = 1 + 3*r
    A[0, 1] = -r

    #RHS (coefficients from Crank-Nicholson discretization)
    B[0, 0] = 1 - 3*r
    B[0, 1] = r

    # ----- Interior Cells: i = 1,..., N-2 -------------------------
    for i in range(1, N-1):
        #LHS
        A[i, i-1] = -r
        A[i, i]   = 1 + 2*r
        A[i, i+1] = -r

        #RHS
        B[i, i-1] = r
        B[i, i]   = 1 - 2*r
        B[i, i+1] = r

    # ----- Top Boundary: i = N-1, u(x,1) = 0 ----------------------
    #LHS
    A[N-1, N-2] = -r
    A[N-1, N-1] = 1 + 3*r

    #RHS
    B[N-1, N-2] = r
    B[N-1, N-1] = 1 - 3*r

    return A, B

# ------------------------------------------------------------------
# 3. Build RHS vector for a given u^n
# ------------------------------------------------------------------

def build_rhs(u_n, params, B):
    dx = params["dx"]
    S  = params["S"]
    N  = params["N"]

    rhs = B @ u_n

    # Add source term contribution
    rhs += S * dx * np.ones(N)

    return rhs

# ------------------------------------------------------------------
# 4. Step in x using LU decomposition
# ------------------------------------------------------------------

def step_in_x(params, u0):
    A, B = build_matrices(params)

    # LU factorization of A
    lu, piv = lu_factor(A)

    u = u0.copy()
    N = params["N"]
    Nx = params["Nx"]

    # Initialize storage for x-location and solution profiles
    x_vals = np.zeros(Nx + 1)
    u_store = np.zeros((Nx + 1, N))
    
    for n in range(Nx + 1):
        x_now = n * params["dx"]
        x_vals[n] = x_now
        u_store[n, :] = u.copy()

        # Break after storing solution at final x-location to avoid an extra solve
        if n == Nx:
            break

        rhs = build_rhs(u, params, B)
        u = lu_solve((lu, piv), rhs)

    return x_vals, u_store

# ------------------------------------------------------------------
# 5. Post processing
# ------------------------------------------------------------------

def analytical_solution(x, y, n_terms=100):
    """Compute the analytical solution (Eq. 19) of the parabolic PDE"""

    u_exact = 0.5 * y * (1 - y)

    for k in range(n_terms):
        n = 2*k + 1 # sum over odd integers only
        A_n = -4 / (n * np.pi)**3
        decay = np.exp(-2 * x * (n * np.pi)**2)
        u_exact += A_n * decay * np.sin(n * np.pi * y)
    
    return u_exact

def post_process(params, x_vals, u_store):
    y = params["y"]

    x_targets = [0.0, 0.01, 0.05, 0.1, 0.25, 1.0]

    errors       =[]
    x_error_vals = []

    # Plot 1: Numerical vs Analytical Solutions at selected x locations
    plt.figure(figsize=(10, 6))

    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown']

    for i, x_target in enumerate(x_targets):
        idx = np.argmin(np.abs(x_vals - x_target))
        x_now = x_vals[idx]
        u_num = u_store[idx, :]
        u_exact = analytical_solution(x_now, y)

        error = np.linalg.norm(u_num - u_exact, ord=2)
        errors.append(error)
        x_error_vals.append(x_now)

        plt.plot(y, u_num, '-d', color=colors[i], alpha=0.5, label=f'Numerical at x={x_now:.2f}', linewidth=2, markersize=4)
        plt.plot(y, u_exact, '--', color=colors[i], label=f'Analytical at x={x_now:.2f}', linewidth=1.5)

    plt.xlabel('y', fontsize=18)
    plt.ylabel('u', fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid(True, linestyle='--',alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig('profiles_comparison.png', dpi=300, bbox_inches='tight')

    # Plot 2/3: Residual, Error vs x
    all_errors = []
    residuals = []

    for n in range(len(x_vals)):
        # Error
        u_num = u_store[n, :]
        u_exact = analytical_solution(x_vals[n], y)
        all_errors.append(np.linalg.norm(u_num - u_exact, ord=2))

        # Residual
        if n < len(x_vals) - 1: # skip last point since we don't have u at x_{n+1}
            u_next = u_store[n+1, :]
            residuals.append(np.linalg.norm(u_next - u_num, ord=2))

    # Plot 2: residual vs x
    plt.figure(figsize=(5, 5))
    plt.plot(x_vals[1:], residuals, '-o', color = 'tab:blue', linewidth=2, markersize=4)
    plt.xlabel('x', fontsize=18)
    plt.ylabel(r'L2 Residual: $\|u^{n+1} - u^n\|_2$', fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig('residual_vs_x.png', dpi=300, bbox_inches='tight')

    # Plot 3: error vs x
    plt.figure(figsize=(5, 5))
    plt.plot(x_vals, all_errors, '-o', color='tab:red', linewidth=2, markersize=4)
    plt.xlabel('x', fontsize=18)
    plt.ylabel(r'L2 Error: $\|u_{num} - u_{exact}\|_2$', fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig('error_vs_x_all.png', dpi=300, bbox_inches='tight')
    
# ------------------------------------------------------------------
# 5. Main Driver
# ------------------------------------------------------------------

def main():
    params, u0 = set_parameters()
    x_vals, u_store = step_in_x(params, u0)
    post_process(params, x_vals, u_store)

if __name__ == "__main__":
    main()


