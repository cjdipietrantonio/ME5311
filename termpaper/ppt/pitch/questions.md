# Questions for Peers
## ME 5311

# Kegan Flory - Black–Scholes–Merton equation 
- Have you considered representing volatility as a function using previous market data? (i.e, have you considered non-constant volatility (which is a real limitation of the basic Black-Scholes model))
- Have you thought about the process of modeling the price/value of a portfolio comprised of many stocks? (eludes to the curse of dimensionality and the need for particle methods)


# Dwaritha Ramesh - Numerical Diffusion

- What are some downsides or limitations of nonlinear schemes like the KNP scheme? Have you encountered any of these in practice, and if so, how did you address them?

- On your last slide you stated that choosing the appropriate scheme is a balance between the physics at hand, computational cost, and implementation complexity. Could you walk me through how you navigated that process in your lab experiments with Dr. Zhao.

- How do you decide where to make the trade-off between shrinking the grid and exploring other treatment options for numerical diffusion?

- In the pressure contours from your detonation simulation, you showed the effects of numerical diffusion very clearly. Do you also encounter dispersive oscillations in practice? If so, between the two which is typically harder to mitigate?

- In your detonation simulations, have you run into cases where the flux limiter itself introduces problems — like clipping physical extrema or adding excessive dissipation near contact discontinuities? If so, how did you treat them.
