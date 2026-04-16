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

    



# 15.20 Hull pp335 (Ch. 15)