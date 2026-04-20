# Speaker Notes — 7-Minute Presentation
## ME 5311 Term Paper: BSM ↔ Transported PDF Method

**Target: ~45-50 seconds per content slide, 8 content slides + title + references = ~7 minutes**

---

## Slide 1: Title (~30 sec)

Good morning everyone. My name is Christian DiPietrantonio and today I'll be presenting my term paper on the connection between the Black-Scholes-Merton equation for pricing financial derivatives and the transported probability density function method used in computational fluid dynamics.

The central argument is that these two equations — from completely different fields — are derived from the same type of stochastic differential equation and can be solved by the same numerical methods. The BSM PDE is the backward Kolmogorov equation, while the transported PDF equation is the forward Kolmogorov equation — both arising from SDEs of the drift-plus-diffusion form.

---

## Slide 2: The Transported PDF Method in CFD (~60 sec)

Let me start with the CFD side — the transported PDF method.

In turbulent flows, the velocity field is random. As Pope discusses in Turbulent Flows, if you repeat the same experiment under identical conditions, you get different results each time. So we use a statistical approach and work with the probability density function of velocity.

Pope derives the Eulerian PDF transport equation directly from the Navier-Stokes equations. This describes how the PDF of velocity evolves in space and time from a fixed-frame perspective.

To model the evolution of fluid particle velocities, Pope introduces the generalized Langevin model — a stochastic differential equation with a deterministic drift coefficient "a" and a diffusion coefficient "b", driven by a Wiener process. This is the Lagrangian description — we're tracking individual particles.

From this Langevin SDE, using stochastic calculus, Pope derives the Fokker-Planck equation — also called the forward Kolmogorov equation. This is the PDE that governs how the conditional PDF of particle velocity evolves. The Langevin SDE and the Fokker-Planck equation describe the same physics from two different viewpoints — Lagrangian and Eulerian.

---

## Slide 3: The BSM Equation (~50 sec)

Now the finance side. Stock prices are random — nobody knows what a stock will be worth tomorrow. Black and Scholes in 1973, and independently Merton, derived a PDE for pricing European options.

The BSM PDE is shown here. It was derived from the Geometric Brownian Motion SDE for stock prices using Ito's Lemma and a hedging argument.

For our CFD class, look at the structure. The half-sigma-squared-S-squared term is diffusion. The rS term is advection. And negative rV is a reaction or decay term. This is an advection-diffusion-reaction equation — the same type of PDE we've been solving all semester.

Here's the critical point: the BSM PDE is the BACKWARD Kolmogorov equation associated with the GBM SDE. The Fokker-Planck equation on the previous slide was the FORWARD Kolmogorov equation. Both arise from SDEs of the same form — one looks forward in probability, the other looks backward in conditional expectation.

---

## Slide 4: The Connection (~50 sec)

This is the central slide. Both equations are derived from stochastic differential equations of the same drift-plus-diffusion form.

Finance uses geometric Brownian motion: dS = rS dt + sigma S dW. CFD uses the Langevin equation: dU-star = a dt + b dW. Same mathematical structure.

The BSM PDE is the backward Kolmogorov equation — it computes the conditional expectation of a future payoff given the current state. The Fokker-Planck is the forward Kolmogorov equation — it describes how the probability density evolves forward in time.

The variable mapping is direct: stock price maps to particle velocity, risk-free rate to drift coefficient, volatility to diffusion coefficient. And the solution methods map too: Monte Carlo simulation in finance corresponds exactly to the Lagrangian particle method in CFD.

---

## Slide 5: Two Solution Approaches (~50 sec)

Both equations can be solved two ways.

The Eulerian approach solves the PDE directly on a fixed grid. For BSM, I implemented an implicit finite difference scheme following Hull Chapter 21 — using the same tridiagonal banded solver from our course projects. For CFD, this is the finite volume method. The advantage is no statistical noise. The limitation is the curse of dimensionality.

The Lagrangian approach simulates random particle trajectories and averages. For BSM, this is Monte Carlo simulation. For CFD, this is Pope's Lagrangian particle method from Chapter 12. The advantage is natural scaling to high dimensions. The limitation is slow convergence — error decreases as one over root N.

A critical distinction: for BSM with constant coefficients, Ito's Lemma gives an exact SDE step formula. For the Langevin equation with variable coefficients, only Euler-Maruyama discretization is available. And because no closed-form PDF exists for variable-coefficient SDEs, particle methods become the only practical approach in real turbulent flows.

---

## Slide 6: Results (~45 sec)

Here are the results. On the left, option value versus stock price. The black curve is the exact analytical Black-Scholes formula. Red squares are the finite difference solver — Eulerian. Blue circles are Monte Carlo — Lagrangian. All three agree, confirming the Feynman-Kac theorem.

On the right, Monte Carlo convergence. The error decreases as one over root N — the same slow convergence that characterizes Lagrangian particle methods in CFD. To cut the error in half, you need four times the paths. This is a fundamental tradeoff between Eulerian and Lagrangian approaches.

---

## Slide 7: From Particles to PDF (~50 sec)

This slide demonstrates the transported PDF method using our BSM example.

On the left: 30 individual Monte Carlo paths — particle trajectories in the Lagrangian sense. They fan out over time, showing drift and diffusion.

On the right: instead of tracking individual particles, we histogram their positions at each time snapshot. The shaded bars are Monte Carlo histograms from 100,000 particles. The solid curves are the exact analytical log-normal PDF. They match.

This IS the transported PDF method in action. Lagrangian particle statistics converge to the Eulerian PDF. Notice the drift — the peak shifts right, like advection — and the diffusion — the distribution spreads, just like a scalar field in a turbulent flow.

In real turbulent flows, you can't compute the analytical PDF because the Fokker-Planck equation has variable coefficients. The particle method is the only way to obtain the PDF.

---

## Slide 8: 3D PDF Surface (~40 sec)

This 3D surface was built entirely from Monte Carlo particle data — 100,000 simulated paths, histogrammed at each time snapshot. No analytical PDF formula was used.

This is exactly how the transported PDF method works in CFD. Simulate particles, histogram their positions, and the PDF surface emerges from the statistics.

For BSM with constant coefficients, we have the luxury of verifying against the exact log-normal PDF. For real turbulent flows with variable coefficients, the particle-based PDF is all you have — which is why the Lagrangian particle method is essential.

---

## Slide 9: Conclusions (~40 sec)

Four takeaways:

One — the BSM PDE is the backward Kolmogorov equation and the Fokker-Planck is the forward Kolmogorov equation. Both are derived from SDEs of the same drift-diffusion form, directly connecting financial derivative pricing to the transported PDF method in CFD.

Two — both Eulerian and Lagrangian methods solve these equations, producing identical results validated against the exact analytical formula.

Three — Monte Carlo particle simulations reproduce the analytical PDF, demonstrating the same Lagrangian-to-Eulerian convergence that Pope's transported PDF method relies on.

Four — BSM serves as an ideal validation problem because closed-form solutions exist for both the SDE and the PDF — a luxury unavailable in general turbulent flows, where particle methods become essential.

Thank you. I'm happy to take questions.

---

## Slide 10: References

[No speaking needed — just displayed]

---

## Slide 11: Backup — Heatmap

Available if questions arise about the PDF evolution visualization. Shows drift path (white dashed) and ±1σ bounds (white dotted) overlaid on the probability density.

---

## Potential Q&A Questions (for elevator pitch on 4/29)

**Q: What's the Wiener process / dW?**
A: It's the mathematical model for pure randomness — a continuous-time random walk. At each instant, dW is a random number from a normal distribution with mean zero and variance dt.

**Q: What is Itô's Lemma?**
A: The stochastic calculus version of the chain rule. Because the random term has a non-zero square (dW² = dt), you get an extra correction term. This is where the -½σ² correction comes from.

**Q: What is the Feynman-Kac theorem?**
A: It proves that the solution to certain PDEs equals the expected value over stochastic paths. It guarantees the FD solver and MC simulation give the same answer.

**Q: Why use Monte Carlo if the analytical solution exists?**
A: For this simple European call, you wouldn't. But for complex options or multiple correlated stocks, no analytical formula exists. Same in CFD — you need particle methods for complex turbulent flows because the Fokker-Planck equation can't be solved in closed form.

**Q: How does volatility relate to turbulence?**
A: They're analogous. Volatility sigma measures how much the stock fluctuates. Turbulence intensity measures velocity fluctuations. Both quantify randomness in their respective systems, and both appear as the diffusion coefficient in the SDE.

**Q: What's the difference between backward and forward Kolmogorov?**
A: Both come from the same SDE. The forward Kolmogorov (Fokker-Planck) describes how the probability density evolves forward in time. The backward Kolmogorov (BSM) computes the expected value of a function of the future state, given the current state. Different questions about the same stochastic process.
