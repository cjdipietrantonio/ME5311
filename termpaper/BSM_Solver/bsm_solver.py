"""
Black-Scholes-Merton Solver
----------------------------------------------------------------------
ME 5311 Term Paper - Christian DiPietrantonio

This script demonstrates three approaches to pricing a European call option:
    1. Analytical Black-Scholes formula (closed-form solution)
    2. Monte Carlo simulation (Lagrangian / stochastic particle method)
    3. Finite Difference (Eulerian / grid-based PDE solver)

Methods 2 and 3 are the Lagrangian and Eulerian numerical approaches, mirroring the two ways to solve the PDF transport equations in CFD.
"""

import numpy as np
from scipy.stats import norm
from scipy.linalg import solve_banded
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import os

# ----------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------
K          = 100.0 # Strike Price ($)
r          = 0.05  # Risk-free interest rate (5%)
sigma      = 0.2   # Volatility 
T          = 1.0   # Time to expiry (1 year)

# Monte Carlo Settings
N_paths    = 200000      # Number of simulated paths
N_steps    = 252         # Time steps (trading days in a year)
N_S_pts    = 80          # Number of starting prices to evaluate


# Stock price range
S_min = 50
S_max = 200

# Output directory
out_dir = "figures"
os.makedirs(out_dir, exist_ok=True)

# ----------------------------------------------------------------------
# Analytical Black-Scholes Formula
# ----------------------------------------------------------------------
def bs_call(S, K, r, sigma, T):
    """
    Analytical Black-Scholes price for a European call option. 
    Page 644, Black-Scholes
    Page 91,  Baxter-Rennie
    Page 336, Hull

    Parameters
    ----------
    S     : current stock price
    K     : strike price
    r     : risk free rate
    sigma : volatility
    T     : time to expiration

    Returns
    -------
    V     : option prices
    """

    d1 = (np.log(S / K) + (r + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    V = S * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)

    return V

# ----------------------------------------------------------------------
# Monte Carlo Solver - Terminal Value Method
# ----------------------------------------------------------------------
def mc_option_price(S0, r, sigma, T, N_paths):
    """
    Price a European call via Monte Carlo using the exact GBM solution. Instead of stepping through time, we jump diretly to S(T). This works because we only need the terminal stock price to compute the payoff. For the PDF evolution, we must use the full time-steppind version below.

    Page 217, Hull (Long position European Call option payoff)
    Page 336, Hull (Option expected value discounted at risk-free rate)
    Page 471, Hull (Stock price at expiration time)
    Page 473, Hull (Standard Error of estimate)
    

    Parameters
    ----------
    S0      : starting stock price
    K       : strike price
    r       : risk free rate
    sigma   : volatility
    T       : time to expiration
    N_paths : number of simulated paths

    Returns
    -------
    V_mean  : estimated option price/value
    V_std   : standard error of the estimate

    """

    Z  = np.random.standard_normal(N_paths)
    ST = S0 * np.exp((r - 0.5 * sigma**2) * T + sigma * np.sqrt(T) * Z)

    # payoff for a European call
    payoffs = np.maximum(ST - K, 0.0)

    # discount back to present value
    V_mean = np.exp(-r * T) * np.mean(payoffs)
    V_std  = np.exp(-r * T) * np.std(payoffs) / np.sqrt(N_paths)

    return V_mean, V_std

# ----------------------------------------------------------------------
# Monte Carlo Solver - Full Path Simulation
# ----------------------------------------------------------------------
def mc_full_paths(S0, r, sigma, T, N_paths, N_steps):
    """
    Simulate full GBM paths.

    Page 471, Hull (Stock price at time step)
    
    Parameters
    ----------
    S0      : starting stock price
    r       : risk free rate
    sigma   : volatility
    T       : time to expiration
    N_paths : number of simulated paths
    N_steps : number of time steps

    Returns
    -------
    S       : all price paths
    t       : time array

    """

    dt = T / N_steps
    t = np.linspace(0, T, N_steps + 1)

    S = np.zeros(N_steps + 1, N_paths)
    S[0, :] = S0

    # generate random increments 
    Z = np.random.standard_normal((N_steps, N_paths))

    # step through time
    for i in range(N_steps):
        S[i+1, :] = S[i, :] * np.exp((r - 0.5 * sigma**2) * dt + sigma * np.sqrt(dt) * Z[i, :])

    return S, t

# ----------------------------------------------------------------------
# Finite Difference Solver (Eulerian / Grid Based PDE Solver)
# ----------------------------------------------------------------------
def fd_bsm(K, r, sigma, T, S_max=500.0, M=500, N_t=5000):
    """
    Solver the BSM PDE using implicit finite difference method.
    
    This is the Eulerian approach: discritize stock price into a grid and solve the PDE directly at each point. Important to note that we march backwards in time from t = T (expiry) to t=0 (today).

    Parameters
    ----------
    K       : strike price
    r       : risk free rate
    sigma   : volatility
    T       : time to expiration
    S_max   : upper bound of stock price grid
    M       : number of spacial intervals 
    N_t     : number of time steps

    Returns
    -------
    S_grid  : stock price at grid points
    V       : option values at t = 0 (today)
    """

    # ----- Grid Setup -----
    dS      = S_max / M          # spatial step size
    dt      = T / N_t            # time step size
    S_grid  = np.linspace(0, S_max, M + 1)

    # ----- Terminal Condition -----
    # At expiry, option value = payoff = max(S - K, 0)
    V = np.maximum(S_grid - K, 0.0)

    # ----- Build Tridiagonal Coefficients for Interior Nodes j = 1, ..., M - 1 -----
    j = np.arrange(1, M)

    a_j = 0.5 * dt * (r * j - sigma**2 * j**2)
    b_j = 1.0 * dt * (sigma**2 * j**2 + r)
    c_j = -0.5 * dt * (sigma**2 * j**2 + r * j)

    # ----- Format into Banded System -----
    n_interior = M - 1 # number of unknowns
    # Banded matrix (l + u + 1, n) where l = 1 lower diagonal, u = 1 upper diagonal
    # element[row, col] maps to ab[u + row - col, col]
    ab = np.zeros((3, n_interior))
    ab[0, 1:]  = c_j[:-1] 
    ab[1,:]    = b_j
    ab[2, :-1] = a_j[1:]

    # ----- Time Marching -----
    for n in range(N_t):
        rhs = V[1:M].copy

        time_to_expiry = (n + 1) * dt
        V_top = S_max - K * np.exp(-r * time_to_expiry)

        rhs[-1] -= c_j[-1] * V_top

        # Solve the Trigiagonal System
        V[1:M] = solve_banded((1, 1), ab, rhs)

        # Apply boundary condtions to full solution array
        V[0] = 0.0
        V[M] = V_top

        return S_grid, V

# ----------------------------------------------------------------------
# Figure 1: V(S) Comparison
# ----------------------------------------------------------------------
def plot_option_value_curve():
    """
    Compare Monte Carlo, Finite Difference, and Analytical Option Prices.
    """

    print("Computing V(S) curve ...")
    S_arr = np.linspace(S_min, S_max, N_S_pts)

    #Analytical Solution
    V_analytical = bs_call(S_arr, K, r, sigma, T)

    # Monte Carlo (loop over starting prices)
    V_mc         = np.zeros(N_S_pts)
    V_err        = np.zeros(N_S_pts)

    for i, S0 in enumerate(S_arr):
        V_mc[i], V_err[i] = mc_option_price(S0, K, r, sigma, T, N_paths)

    # Finite Difference
    print("     Running finite difference solver...")
    S_fd, V_fd = fd_bsm(K, r, sigma, T, S_max=500.0, M=500, N_t=5000)

    # Interpolate FD results onto same S_arr for comparison
    V_fd_interp = np.interp(S_arr, S_fd, V_fd)

    # ----- Plot -----
    fig, ax = plt.subplots(figsize=(9, 6))

    # Analytical Solution
    ax.plot(S_arr, V_analytical, 'k-', linewidth = 2, label='Analytical (Black-Scholes)')

    # Finite Difference Solution
    ax.plot(S_arr, V_fd_interp, 's', color='tab:red', markersize=4, markerfacecolor='none', markeredgewidth=1.0, label='Finite Difference (Eulerian)', alpha=0.8)

    # Monte Carlo (with error bars)
    ax.errorbar(S_arr, V_mc, yerr=4*V_err, marker='o', color='tab:blue', markersize=4, linewidth =1.0, capsize=2, label=f'Monte Carlo / Lagrangian (N = {N_paths:,})', alpha=0.8)

    # Intrinsic Value for Reference
    intrinsic = np.maximum(S_arr - K, 0)
    ax.plot(S_arr, intrinsic, '--', color='tab:gray', linewidth=1.0, alpha=0.8, label='Intrinsic Value: max(S - K, 0)')

    ax.set_xlabel('Current Stock Price, S ($)', fontsize=12)
    ax.set_ylabel('European Call Option Value, V ($)', fontsize=12)
    ax.set_title('European Call Option Pricing: Comparison of Methods', fontsize=15)
    ax.set_xlim(S_min, S_max)
    ax.set_ylim(bottom=-2)
    ax.grid(True, alpha=0.3)

    param_text = (f'K = ${K:.0f},  r = {r:.0%},  σ = {sigma:.0%},  T = {T:.1f} yr')
    ax.text(0.02, 0.97, param_text, transform=ax.transAxes, fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', edgecolor='black', pad=0.5, alpha=0.8))

    fig.tight_layout()
    fig.savefig(f'{out_dir}/figure01_option_value_curve.png', dpi=200)
    plt.close(fig)
    print(f'     Saved:  {out_dir}/figure01_option_value_curve.png')

    return S_arr, V_analytical, V_mc, V_err, V_fd_interp

# ----------------------------------------------------------------------
# Figure 2: Sample GBM Paths
# ----------------------------------------------------------------------


