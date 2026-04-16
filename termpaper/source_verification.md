# Source Verification & Citation Tracking
## ME 5311 Term Paper — Christian DiPietrantonio

This document maps every equation, concept, and claim used in the term paper
and code to a verifiable source. Sources from the interim paper bibliography
are used wherever possible.

---

## Bibliography Key

| Cite Key           | Full Reference                                                                                      |
|--------------------|------------------------------------------------------------------------------------------------------|
| BlackScholes1973   | Black, F. & Scholes, M. (1973). "The pricing of options and corporate liabilities." *J. Political Economy*, 81(3), 637–654. |
| Merton1973         | Merton, R.C. (1973). "Theory of rational option pricing." *Bell J. Econ. & Mgmt. Sci.*, 4(1), 141–183. |
| Hull2022           | Hull, J.C. (2022). *Options, Futures, and Other Derivatives*, 10th ed. Pearson.                      |
| Shreve2004         | Shreve, S.E. (2004). *Stochastic Calculus for Finance II: Continuous-Time Models*. Springer.         |
| Pope2000           | Pope, S.B. (2000). *Turbulent Flows*. Cambridge University Press.                                   |
| BaxterRennie1996   | Baxter, M. & Rennie, A. (1996). *Financial Calculus: An Introduction to Derivative Pricing*. CUP.   |

---

## Equations & Where to Verify Them

### 1. Black-Scholes-Merton PDE

$$\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} + rS\frac{\partial V}{\partial S} - rV = 0$$

| Source             | Where to Find It                                                        |
|--------------------|-------------------------------------------------------------------------|
| BlackScholes1973   | Eq. (7) in the original paper, p. 643                                   |
| Merton1973         | Derived independently via stochastic calculus; Eq. (9), p. 164          |
| Hull2022           | Chapter 15, Eq. (15.16); also derived in Chapter 14                     |
| Shreve2004         | Chapter 4, Theorem 4.5.2 (derived from Itô's Lemma + hedging argument) |

### 2. Geometric Brownian Motion SDE

$$dS(t) = \alpha S(t)\,dt + \sigma S(t)\,dW(t)$$

(In the risk-neutral measure, α is replaced by r)

| Source             | Where to Find It                                                        |
|--------------------|-------------------------------------------------------------------------|
| Shreve2004         | Chapter 4, Eq. (4.4.1), p. 153 — used as starting point for BSM derivation |
| Hull2022           | Chapter 14, Eq. (14.1) — "model of stock price behavior"               |
| BlackScholes1973   | p. 640, described in the assumptions (stated in words, not SDE notation)|
| Pope2000           | Compare with Eq. (12.1) — same SDE structure as the Langevin equation  |

### 3. Exact GBM Solution (used in Monte Carlo code)

$$S(T) = S_0 \exp\!\left[\left(r - \tfrac{1}{2}\sigma^2\right)T + \sigma\sqrt{T}\,Z\right], \quad Z \sim \mathcal{N}(0,1)$$

| Source             | Where to Find It                                                        |
|--------------------|-------------------------------------------------------------------------|
| Shreve2004         | Chapter 4, Eq. (4.4.9) — derived by applying Itô's Lemma to ln(S)     |
| Hull2022           | Chapter 14, Eq. (14.14) — "From equation (14.3), the distribution of ln S(T)" |

**Note on the Itô correction (−½σ² term):** This comes from Itô's Lemma applied
to f(S) = ln(S). Because d(ln S) ≠ dS/S in stochastic calculus (unlike ordinary
calculus), a correction term appears. See Shreve2004 Chapter 4, Theorem 4.4.1,
or Hull2022 Chapter 14, Section 14.5.

### 4. Analytical Black-Scholes Formula for European Call

$$V = S\,N(d_1) - Ke^{-rT}\,N(d_2)$$

$$d_1 = \frac{\ln(S/K) + (r + \tfrac{1}{2}\sigma^2)T}{\sigma\sqrt{T}}, \qquad d_2 = d_1 - \sigma\sqrt{T}$$

| Source             | Where to Find It                                                        |
|--------------------|-------------------------------------------------------------------------|
| BlackScholes1973   | Eq. (13), p. 644 — the original formula                                |
| Merton1973         | Eq. (19), p. 172                                                        |
| Hull2022           | Chapter 15, Eqs. (15.20)–(15.22)                                       |
| Shreve2004         | Chapter 4, Theorem 4.5.4, p. 170                                       |

### 5. European Call Payoff Function

$$\text{Payoff} = \max(S_T - K, 0)$$

| Source             | Where to Find It                                                        |
|--------------------|-------------------------------------------------------------------------|
| BlackScholes1973   | p. 637 — defines the call option contract                              |
| Hull2022           | Chapter 10, Fig. 10.3 — payoff diagrams                                |
| Shreve2004         | Chapter 4, Section 4.5 — used as terminal condition                    |

### 6. Monte Carlo Pricing via Risk-Neutral Expectation

$$V = e^{-rT}\,\mathbb{E}^{\mathbb{Q}}\!\left[\max(S_T - K, 0)\right]$$

This is the **Feynman-Kac** representation: the PDE solution equals a
discounted expectation over stochastic paths.

| Source             | Where to Find It                                                        |
|--------------------|-------------------------------------------------------------------------|
| Shreve2004         | Chapter 4, Theorem 4.5.2 (Feynman-Kac); Ch. 4.5 risk-neutral pricing  |
| Hull2022           | Chapter 21, Section 21.6 — "Monte Carlo Simulation"                    |
| BaxterRennie1996   | Chapter 3 — risk-neutral pricing and equivalent martingale measure      |

### 7. Monte Carlo Convergence Rate

Standard error = σ / √N (error decreases as 1/√N)

| Source             | Where to Find It                                                        |
|--------------------|-------------------------------------------------------------------------|
| Hull2022           | Chapter 21, Section 21.6 — discusses MC accuracy and path count        |
| Shreve2004         | Central Limit Theorem application — standard probability result        |

### 8. Log-Normal Distribution of Stock Price

Under GBM: $\ln S(T) \sim \mathcal{N}\!\left(\ln S_0 + (r - \tfrac{1}{2}\sigma^2)T,\; \sigma^2 T\right)$

| Source             | Where to Find It                                                        |
|--------------------|-------------------------------------------------------------------------|
| Hull2022           | Chapter 14, Eq. (14.14) — explicitly states the log-normal result      |
| Shreve2004         | Chapter 4, follows from Eq. (4.4.9)                                    |

---

## CFD / PDF Transport Equations & Where to Verify Them

### 9. Langevin SDE (Generalized Diffusion Process)

$$dU^*(t) = a(U^*, X^*, t)\,dt + b(X^*, t)\,dW(t)$$

| Source             | Where to Find It                                                        |
|--------------------|-------------------------------------------------------------------------|
| Pope2000           | Chapter 12, Eq. (12.1) — general form of the Langevin model            |

### 10. Eulerian PDF Transport Equation

$$\frac{\partial f}{\partial t} + V_i\frac{\partial f}{\partial x_i} = -\frac{\partial}{\partial V_i}\left[f\left\langle\frac{DU_i}{Dt}\middle|\mathbf{V}\right\rangle\right]$$

| Source             | Where to Find It                                                        |
|--------------------|-------------------------------------------------------------------------|
| Pope2000           | Chapter 12, Eq. (12.30) — derived from Navier-Stokes                  |

### 11. Lagrangian / Fokker-Planck PDF Transport Equation

$$\frac{\partial f^*}{\partial t} + V_i\frac{\partial f^*}{\partial x_i} = -\frac{\partial}{\partial V_i}\left[f^* a_i\right] + \frac{1}{2}b^2\frac{\partial^2 f^*}{\partial V_i \partial V_i}$$

| Source             | Where to Find It                                                        |
|--------------------|-------------------------------------------------------------------------|
| Pope2000           | Chapter 12, Eq. (12.52) — Fokker-Planck for the Langevin model        |

### 12. Connection: Backward vs. Forward Kolmogorov Equations

- BSM PDE = **backward** Kolmogorov equation for GBM
- Fokker-Planck = **forward** Kolmogorov equation for the same SDE

| Source             | Where to Find It                                                        |
|--------------------|-------------------------------------------------------------------------|
| Shreve2004         | Chapter 6, Section 6.4 — Kolmogorov equations and Feynman-Kac         |
| Pope2000           | Chapter 12, Section 12.2 — derives Fokker-Planck from the Langevin SDE|

---

## Numerical Method Sources

### 13. Finite Difference Method for BSM (our Python FD solver)

The implicit Euler scheme with central differences applied to the BSM PDE,
producing a tridiagonal system solved at each time step.

| Source             | Where to Find It                                                        |
|--------------------|-------------------------------------------------------------------------|
| Hull2022           | Chapter 21, Sections 21.7–21.8 — "Finite Difference Methods" for BSM; derives the explicit, implicit, and Crank-Nicolson schemes |

**Specific items in Hull Ch. 21:**
- The grid setup with S_j = j·ΔS: Section 21.7
- Central difference stencils for ∂V/∂S and ∂²V/∂S²: Eqs. (21.8)–(21.9)
- Tridiagonal coefficient formulas (a_j, b_j, c_j): Eqs. (21.18)–(21.20)
- Boundary conditions V(0,t)=0 and V(S_max,t)≈S_max−Ke^{-r(T-t)}: Section 21.7
- Implicit vs explicit stability discussion: Section 21.8

### 14. OpenFOAM financialFoam (Finite Volume Method for BSM)

OpenFOAM solves the BSM PDE using the finite volume method with the same
advection-diffusion-reaction structure.

| Source             | Where to Find It                                                        |
|--------------------|-------------------------------------------------------------------------|
| OpenFOAM source    | github.com/OpenFOAM/OpenFOAM-3.0.x — `applications/solvers/financial/financialFoam/financialFoam.C` (link provided by Prof. Zhao) |

---

## Concept / Claim Sources

### 15. "Options give the right but not the obligation to buy/sell"

| Source             | Where to Find It                                                        |
|--------------------|-------------------------------------------------------------------------|
| BlackScholes1973   | p. 637 — opening definition                                            |
| Hull2022           | Chapter 10, Section 10.1 — "Types of Options"                          |

### 16. Assumptions of BSM (no-arbitrage, continuous trading, constant r and σ, no dividends, European exercise)

| Source             | Where to Find It                                                        |
|--------------------|-------------------------------------------------------------------------|
| BlackScholes1973   | pp. 640–641 — "ideal conditions" list                                  |
| Hull2022           | Chapter 15, Section 15.8                                                |
| BaxterRennie1996   | Chapter 3 — arbitrage-free pricing framework                           |

### 17. Volatility (σ) as a measure of stock price uncertainty

| Source             | Where to Find It                                                        |
|--------------------|-------------------------------------------------------------------------|
| Hull2022           | Chapter 14, Section 14.4 — defines and interprets volatility           |

### 18. Higher volatility → higher option value

| Source             | Where to Find It                                                        |
|--------------------|-------------------------------------------------------------------------|
| Hull2022           | Chapter 11, Table 11.1 — factors affecting option prices; σ has positive effect on calls |

### 19. Time value of an option (V > intrinsic value before expiry)

| Source             | Where to Find It                                                        |
|--------------------|-------------------------------------------------------------------------|
| Hull2022           | Chapter 11, Section 11.1 — intrinsic value vs. time value              |

### 20. Discounting (e^{-rT} factor to convert future to present value)

| Source             | Where to Find It                                                        |
|--------------------|-------------------------------------------------------------------------|
| Hull2022           | Chapter 4 — present value and continuous compounding                   |
| Shreve2004         | Chapter 4, Section 4.3 — discount factor in risk-neutral pricing       |

### 21. Curse of dimensionality (MC scales better than grids in high-D)

| Source             | Where to Find It                                                        |
|--------------------|-------------------------------------------------------------------------|
| Hull2022           | Chapter 21, Section 21.6 — discusses when MC is preferred over FD      |

---

## Potential Additional Sources to Consider

If you want to strengthen specific sections of the final paper, these could
be added to the bibliography:

- **Glasserman, P. (2003). *Monte Carlo Methods in Financial Engineering*. Springer.**
  — The definitive reference for MC methods in finance. Good for convergence
  rates, variance reduction, and path simulation details.

- **Wilmott, P., Howison, S., & Dewynne, J. (1995). *The Mathematics of Financial
  Derivatives*. Cambridge University Press.**
  — Excellent for the PDE perspective on BSM and finite difference methods.

- **Gardiner, C.W. (2009). *Stochastic Methods*, 4th ed. Springer.**
  — Bridges the Fokker-Planck / Langevin gap between physics and finance.

---

*This document will be updated as new equations or claims are added to the
paper or presentation.*
