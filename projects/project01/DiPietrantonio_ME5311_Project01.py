#Christian DiPietrantonio
#ME M311: Computational Methods to Viscous Flows
#Computer Assignment 01
#DATE

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import lu_factor, lu_solve

# -----------------------------------------------------------------------------
# 1. Parameters and Mesh Sizing
# -----------------------------------------------------------------------------

def set_parameters():
    params = {}

    # Domain in y
    params["y_min"] = 0
    params["y_max"] = 1
    params["N"]     = 200  # Number of cells (finite volume method)
    params["dy"]    = (params["y_max"] - params["y_min"]) / params["N"]  #cell size in y

    # Stepping in x
    params["x_max"] = 1.0 # Final x
    params["Nx"]    = 200 # Number of steps in x
    params["dx"]    = params["x_max"] / params["Nx"] # Step size in x

    # Source term
    params["S"] = 2.0

    # y-mesh
    params["y"] = np.linspace(params["y_min"], params["y_max"], params["N"])

    #Initial condition u(y,0) = 0
    u0 = np.zeros(params["N"])

    return params, u0

# -----------------------------------------------------------------------------
# 2. Build Matrices
# -----------------------------------------------------------------------------
   
def build_matrices(params):
    N = params["N"]
    dx = params["dx"]
    dy = params["dy"]

    r = dx / dy**2

    A = np.zeros((N, N))
    B = np.zeros((N, N))

    # ----- Bottom Boundary: i = 0, u(x,0) = 0 ---------------------------------
    #LHS
    A[0, 0] = 1 + 3*r
    A[0, 1] = -r

    #RHS
    B[0, 0] = 1 - 3*r
    B[0, 1] = r

    # ----- Interior Cells: i = 1,..., N-2 -------------------------------------
    for i in range(1, N-1):
        #LHS
        A[i, i-1] = -r
        A[i, i]   = 1 + 2*r
        A[i, i+1] = -r

        #RHS
        B[i, i-1] = r
        B[i, i]   = 1 - 2*r
        B[i, i+1] = r

    # ----- Top Boundary: i = N-1, u(x,1) = 0 ----------------------------------
    #LHS
    A[N-1, N-2] = -r
    A[N-1, N-1] = 1 + 3*r

    #RHS
    B[N-1, N-2] = r
    B[N-1, N-1] = 1 - 3*r

    return A, B

# -----------------------------------------------------------------------------
# 3. Build RHS vector for a given u^n
# -----------------------------------------------------------------------------

def build_rhs(u_n, params, B):
    dx = params["dx"]
    S  = params["S"]
    N  = params["N"]

    rhs = B @ u_n

    #Add source term contribution
    rhs += S * dx * np.ones(N)

    return rhs

# -----------------------------------------------------------------------------
# 4. Step in x using LU decomposition
# -----------------------------------------------------------------------------

def step_in_x(params, u0):
    A, B = build_matrices(params)

    # LU factorization of A
    lu, piv = lu_factor(A)

    u = u0.copy()
    N = params["N"]
    Nx = params["Nx"]

    #Initialize storage lists for x location and solution profiles
    #x_vals  = []
    x_vals = np.zeros(Nx + 1)
    #u_store = []
    u_store = np.zeros((Nx + 1, N))
    
    for n in range(Nx + 1):
        x_now = n * params["dx"]
        #x_vals.append(x_now)
        #u_store.append(u.copy())
        x_vals[n] = x_now
        u_store[n, :] = u.copy()

        if n == Nx:
            break

        rhs = build_rhs(u, params, B)
        u = lu_solve((lu, piv), rhs)
        
    #return np.array(x_vals), np.array(u_store)
    return x_vals, u_store

# -----------------------------------------------------------------------------
# 5. Post processing
# -----------------------------------------------------------------------------

def analytical_solution(x, y, n_terms=100):

    u_exact = 0.5 * y * (1 - y)

    for k in range(n_terms):
        n = 2*k + 1
        A_n = -4 / (n * np.pi)**3
        decay = np.exp(-2 * x * (n * np.pi)**2)
        u_exact += A_n * decay * np.sin(n * np.pi * y)
    
    return u_exact

def post_process(params, x_vals, u_store):
    y = params["y"]

    x_targets = [0.0, 0.25, 0.5, 1.0]
    #indices = range(len(x_vals))
    #indices = [0, len(x_vals)//4, len(x_vals)//2, -1]

    errors       =[]
    x_error_vals = []

    #Plot 1: Numerical vs Analytical Solutions at selected x locations
    plt.figure(figsize=(10, 6))

    for x_target in x_targets:
        idx = np.argmin(np.abs(x_vals - x_target))
        x_now = x_vals[idx]
        u_num = u_store[idx, :]
        u_exact = analytical_solution(x_now, y)

        error = np.linalg.norm(u_num - u_exact, ord=2)
        errors.append(error)
        x_error_vals.append(x_now)

        plt.plot(y, u_num, '-o', label=f'Numerical at x={x_now:.2f}')
        plt.plot(y, u_exact, '--', label=f'Analytical at x={x_now:.2f}')

    # for idx in indices:
    #     x_now = x_vals[idx]
    #     u_num = u_store[idx, :]
    #     u_exact = analytical_solution(x_now, y)

    #     plt.plot(y, u_num, '-o', label=f'Numerical at x={x_now:.2f}', markersize=4)
    #     plt.plot(y, u_exact, '--', label=f'Analytical at x={x_now:.2f}', markersize=4)

    plt.xlabel('y')
    plt.ylabel('u')
    #plt.title()
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    #Plot 2: Error vs x
    plt.figure(figsize=(10, 6))
    plt.plot(x_error_vals, errors, '-o')
    plt.xlabel('x')
    plt.ylabel('L2 Error')
    #plt.title('L2 Error vs x')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    #Plot 3: Error vs x for every step (to demonstrate stability)
    all_errors = []

    for n in range(len(x_vals)):
        u_num = u_store[n, :]
        u_exact = analytical_solution(x_vals[n], y)
        all_errors.append(np.linalg.norm(u_num - u_exact, ord=2))

    plt.figure(figsize=(10, 6))
    plt.plot(x_vals, all_errors, '-o')
    plt.xlabel('x')
    plt.ylabel('L2 Error')
    #plt.title('L2 Error vs x (Stability Check)')
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    
# -----------------------------------------------------------------------------
# 5. Main Driver
# -----------------------------------------------------------------------------

def main():
    params, u0 = set_parameters()
    x_vals, u_store = step_in_x(params, u0)
    post_process(params, x_vals, u_store)

if __name__ == "__main__":
    main()
