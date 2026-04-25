# Speaker Notes — 2-Minute Presentation
## ME 5311 Term Paper: BSM ↔ Transported PDF Method


---

## Slide 1: The Intro

Good morning everyone. My name is Christian DiPietrantonio and today I'll be presenting my term paper on the connection between the Black-Scholes-Merton equation for pricing financial derivatives and the transported probability density function method used in computational fluid dynamics.

---

## Slide 2: The Connection

- On the left, we have the BSM PDE and on the right we have the PDF transport equation. (We can note that the BSM PDE is of the same form as many transport equations we have worked with this term.)

- The key point is that both equations are derived from SDEs of the same drift + difusion form (deterministic drift + random fluctuations driven by a Wiener process). Both equations can also be derived using Ito's Lemma and stochastic calculus. 

- One key distinction that should be made is that the BSM equation computes the conditional expectation of a future payoff given the current state, then discounts it back to give the current option price. In other words, for BSM, you use the PDF at expiry to find the expected payoff, then discount it back to today (this is done implicitly with the SDE hedging described previously, but later we show how this is done explicitly with a Monte Carlo simulation). 

- One thing to note is that the drift and difusion coefficients are constant for the Brownian Motion SDE, but depend on the particles current position and velocity for the Langevin equation. We will return to this later.

- We should also note the mapping of solution methods in addition to other variables.

---

## Slide 3: Results: Option Pricing

- To explore the different solution methods, I decided to write a python script to compare the two mentioned above for the BSM equation. The plot shows the option value for a European call versus stock prices. 
- The Eulerian approach solves the PDE directly on a fixed grid using implicit finite difference scheme, and second order central difference stencils for the first and second derivative terms. Doing so resulted in a tridiagonal system that was extremely similar to what we saw in Project 01 and Flip Classroom Exercise 01. This idea was inspired by Professor Zhao provided me with a repository for OpenFoam (more specifically finacialFoam), which executes solves the BSM equation using FVM.
- The Lagrangian approach simulated 200,000 random stock price trajectories (GBM paths) for 80 stock prices between $50 and $200 using a Monte Carlo simulation, averaged the payoffs, and discounted them, back to present day value.


---

## Slide 4: Monte Carlo Visualization

- Here we can see an example of 30 paths sampled from 100 simulated paths, along with the mean and expected price. 

- We can now extract data from the full Monte Carlo analysis to construct the transport of the PDF of the stock price

---

## Slide 5: PDF Transport

-  Here is a 3D surface resembling the transport of the PDF of the stock price (constructed entirely from Monte Carlo particle data).
- This is exactly how one can imagine the transported PDF method for a turbulent flow.

---

## Slide 6: PDF Transport

- This plot was created to help us visualize those concepts of drift and difusion.

- The shaded bars are pdfs constructed from Monte Carlo data plotted as histograms, and the solid curves are the exact analytical log-normal PDFs. Over time, the volatility spreads the distribution of possible stock prices, and the drift shifts it to the right.

- We can conceptualize this for a turbulent flow, where drift can be thought of as where the average velocity tends to go, and diffusion as how strongly turbulence randomizes the velocity around the mean

- We should note that Hull uses Ito's Lemma to show that, for our GBM model, the price of a stock at a given time follows a lognormal distribution. Therefore, we can use the standard equation for the probability density function of a lognormal variable. For the Langevin SDE, because the coefficients vary with respect to particle velocity/position and time, we don’t have this luxury, so we must construct the PDFs for given time steps using the data from the Monte Carlo simulation.

---

## Slide 7: Conclusions

Finally, to recap what we discussed:

- The BSM the Fokker-Planck (forward Kolmogorov equation) are both derived from SDEs of the same form. This directly connects financial derivative pricing to the transported PDF method in CFD.

- Second, the same results (validated against the exact analytical formula) can be obtained from both Eulerian and Lagrangian methods as demonstrated with the BSM equation.

- Third, Monte Carlo particle simulations reproduce the analytical PDF

- Fourth, BSM is actually a great validation problem because closed-form solutions exist for both the SDE and the PDF. In turbulent flows, these may not exist, which is why particle methods become essential.

- Finally, Advancements in financial modeling and computational fluid dynamics are mutually beneficial. Numerical methods/variance reduction can transfer directly between CFD and quantitative finance, so techniques/methods developed in either discipline can be utilized by/applied to problems in the other.

- One other thing we should note is although I did use a Eulerian solution method for demonstration here, Professor Zhao (and Pope) mention that this isnt really feasible for high dimensional Joint PDFs.

Thank you!

---

## Backup/Definitions:
The following assumptions were made by Black and Scholes:
- Idealized, frictionless market with continuous trading (no transaction cost)
- No dividends
- Unrestricted short selling
- Ability to borrow/lend at the risk-free rate
- European Option
- Arbitrage free market 

Other notes/Definitions:
- Volatility can be estimated from historical data
- Ito process: A stochastic process where the change in a variable during each short period of time of length t has a normal distribution. The mean and variance of the distribution are proportional to t and are not necessarily constant.
- Ito's Lemma: A result that enables the stochastic process for a function of a variable to be calculated from the stochastic process for the variable itself.
- Ito/Stochastic Calculus: Built to specifically handle random processes (define derivatives, chain rules, integrals for random processes)
- Stochastic Process: An equation describing the probabilistic behavior of a stochastic variable
Stochastic Variable: A variable whose future value is uncertain
- SDE: Defines the infinitesimal increment of a stochastic process
- Markov Process: A stochastic process where the behavior of a variable over a short period of time depends Soley on the value of the variable at the beginning of the period, not on its past history.
- Weiner Process: A stochastic process where the change in a variable during each short period of time of length ∆T has a normal distribution with a mean equal to zero and a variance equal to ∆T (Sometimes referred to as Brownian motion)
- Geometric Brownian Motion: A stochastic process often assumed for asset prices where the logarithm of the underlying variable follows a generalized Wiener process
- Risk-free rate: The rate of interest that can be earned without assuming any risks (treasuries are common, London Interbank Offered Rate, LIBOR (less common after 2008), Overnight Indexed Swap rate (OIS)).
- Vector Valued Weiner Process: a vector of Weiner Processes (random motion in multiple directions)

Mertons Derivation:
- Another note on Merton’s derivation: Merton started with the process for a stock price, which is modeled using the Stochastic Differential Equation above. Ito’s lemma is then used to express the rate of change of the value of the option, a hedged portfolio is created and its change in value over a small time interval is expressed, then the result of Ito’s lemma and the SDE are substituted to eliminate the Weiner process, making the portfolio riskless. Because of the assumptions, the change in value of the portfolio must be expressed in terms of the risk-free interest rate. Finally, the original expressions for the value of the portfolio and its change with respect to time can be substituted to obtain the BSM equation! (Hull walks through this entire derivation).

Langevin History:
- The Langevin equation was originally proposed in 1908 as a stochastic model for the velocity of a microscopic particle undergoing Brownian motion (sound familiar?). Pope claims the equation provides a good model for the velocity of a fluid particle in turbulence. The stochastic process U*(t) is called the Ornstein-Uhlenbeck (OU) process, and its PDF evolves by the Fokker-Planck equation. U*(t) is known as a diffusion process, and the Langevin equation shown above is a stochastic differential equation, where a(V, x, t) and b(x, t) are the drift and diffusion coefficients respectively. Note, the Langevin equation it is the simplest stochastic Lagrangian model according to Pope. (Pope)

Solution Methods:
Eulerian:
- Scales poorly to high dimensions, let's say we want to simulate two stocks, each axis with 500 grid points, you would now need 250,000 grid locations to cover all combinations

Lagrangian:
- If portfolio grows, we can just simulate random samples of the full vector (no giant mesh needed)
- Compatible with CPUs/accelerators/clusters because each simulation is independent so they are easy to parallelize
- Need 4x more simulations to cut error in half, and 100x more simulations to cut error by 10x

- We should note, that according to Pope: “It is not feasible, computationally, to represent the joint PDF accurately through a discretization of the seven-dimensional V-θ-x space (for the velocity-frequency joint PDF), as is required in finite-difference, finite-volume, and finite-element methods. Instead, the model PDF equations are solved by particle methods.” (Pope). Professor Zhao also mentioned this in our final lecture. 


Heat Map/Expected Stock Price:
- We clearly see the PDF shifting in line with the expected stock price, which can also be thought of as the deterministic growth of the stock in the absence of volatility. 

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
