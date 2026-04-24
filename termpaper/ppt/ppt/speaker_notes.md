# Speaker Notes — 7-Minute Presentation
## ME 5311 Term Paper: BSM ↔ Transported PDF Method

**Target: ~45-50 seconds per content slide, 8 content slides + title + references = ~7 minutes**

---

## Slide 1: Title

Good morning everyone. My name is Christian DiPietrantonio and today I'll be discussing the transported probability density function method used in computational fluid dynamics and its relation to the Black-Scholes-Merton equation for pricing financial derivatives.

The central argument is that these two equations, from completely different fields, are derived from the same type of stochastic differential equation and can be solved by the same numerical methods.

---

## Slide 2: Financial Derivatives: Options

In order to understand the content of this presentation we must define a few financial terms/concepts. A derivative as defined by Hull is listed above. Very often the variables underlying derivatives are the prices of traded assets. A stock option, for example, is a derivative whose value is dependent on the price of a stock. However, derivatives can be dependent on almost any variable, from the price of hogs to the amount of snow falling at a certain ski resort. There is now active trading in credit derivatives, electricity derivatives, weather derivatives, and insurance derivatives. (HULL)

Above we can see several definitions for different types of options. For the remainer of this presentation, we will be focusing on European Call Options. You can also long (buy) options, or short (sell) options.

As for payout rules, for a short position on an option, you can see that as the price of the stock exceeds the strike price, you lose money because you are obligated to sell the stock at the strike price. For the long position, the profit increases as the price of the underlying asset (stock) exceeds the strike price. Thus, for the remained of our analysis we can express the value of the function as shown above.


---

## Slide 3: The BSM Equation

Now, since we have shown that the value of an option depends on the price of the underlying asset/stock, we must find a way to model said asset. However, this is difficult because In financial markets, stock prices are random. Black and Scholes in 1973, and independently Merton in the same year, derived a second order linear parabolic PDE for pricing European options.

This was done by creating a riskless portfolio by combining a position in a derivative with a position in an underlying stock. Without arbitrage opportunities, the return from the portfolio must be the risk-free interest rate r. The reason this works is because the stock and derivative prices are both effected by the same underlying source of uncertainty (the stock price). When the appropriate portfolio of the stock and derivative are established, the gain or loss from the stock position always offsets the gain or loss from the derivative position, so the overall value of the portfolio at the end of the time period is always known. (Hull)

The following assumptions were made by Black and Scholes:
- Idealized, frictionless market with continuous trading (no transaction cost)
- No dividends
- Unrestricted short selling
- Ability to borrow/lend at the risk-free rate
- European Option
- Arbitrage free market 
        
As previously mentioned, Merton derived the same equation, but he did so through the use of Stochastic Calculus and Ito’s Lemma (which shows that for a variable following a stochastic process, S, the function G(s, t) follows the process shown above) . This type of calculus (stochastic, also referred to as Ito Calculus), was created specifically to handle derivatives and integrals of stochastic processes.

Now we should note that this equation looks awfully familiar to transport equations we have worked with this semester, we have a temporal term, diffusion term, advection term, and a source (decay) term.

Another note on Merton’s derivation: Merton started with the process for a stock price, which is modeled using the Stochastic Differential Equation above. Ito’s lemma is then used to express the rate of change of the value of the option, a hedged portfolio is created and its change in value over a small time interval is expressed, then the result of Ito’s lemma and the SDE are substituted to eliminate the Weiner process, making the portfolio riskless. Because of the assumptions, the change in value of the portfolio must be expressed in terms of the risk-free interest rate. Finally, the original expressions for the value of the portfolio and its change with respect to time can be substituted to obtain the BSM equation! (Hull walks through this entire derivation).

Other notes/Definitions:
- Volatility can be estimated from historical data
- Ito process: A stochastic process where the change in a variable during each short period of time of length t has a normal distribution. The mean and variance of the distribution are proportional to t and are not necessarily constant.
- Ito's Lemma: A result that enables the stochastic process for a function of a variable to be calculated from the stochastic process for the variable itself.
- Ito/Stochastic Calculus: Built to specifically handle random processes (define derivatives, chain rules, integrals for random processes)
- Stochastic Process: An equation describing the probabilistic behavior of a stochastic variable
- Stochastic Variable: A variable whose future value is uncertain
- SDE: Defines the infinitesimal increment of a stochastic process
- Markov Process: A stochastic process where the behavior of a variable over a short period of time depends Soley on the value of the variable at the beginning of the period, not on its past history.
- Weiner Process: A stochastic process where the change in a variable during each short period of time of length ∆T has a normal distribution with a mean equal to zero and a variance equal to ∆T (Sometimes referred to as Brownian motion)
Geometric Brownian Motion: A stochastic process often assumed for asset prices where the logarithm of the underlying variable follows a generalized Wiener process
- Risk-free rate: The rate of interest that can be earned without assuming any risks (treasuries are common, London Interbank Offered Rate, LIBOR (less common after 2008), Overnight Indexed Swap rate (OIS)).

---

## Slide 4: The Transported PDF Method
Pope remarks that the velocity field in a turbulent flow is random. Which means that if one were to repeat the same experiment multiple times under the same set of conditions, a given event may, but need not occur, making the event inherently unpredictable. Therefore, since U is random, its value is inherently unpredictable. Nevertheless, the probability of events can be used to create a Probability Density Function (PDF), and according to Pope, a random variable such as U is completely characterized by its PDF, which can be thought of as the probability per unit distance in the sample space.

The evolution of the Eulerian PDF of the velocity defined by Pope is shown above, and underneath we have the equation for the evolution of the conditional PDF of the particle velocity. We can think of this as the Lagrangian equivalent of the first equation (the conditional PDF, f* is also the corresponding representation of f for a particle system). 

Now, one might ask, where did this Lagrangian PDF transport equation come from? Pope began his derivation by expressing the velocity of each fluid particle as a generalized Langevin model for a diffusion process, as shown above.

The Langevin equation was originally proposed in 1908 as a stochastic model for the velocity of a microscopic particle undergoing Brownian motion (sound familiar?). Pope claims the equation provides a good model for the velocity of a fluid particle in turbulence. The stochastic process U*(t) is called the Ornstein-Uhlenbeck (OU) process, and its PDF evolves by the Fokker-Planck equation. U*(t) is known as a diffusion process, and the Langevin equation shown above is a stochastic differential equation, where a(V, x, t) and b(x, t) are the drift and diffusion coefficients respectively. Note, the Langevin equation it is the simplest stochastic Lagrangian model according to Pope. (Pope)

From this SDE, Ito calculus can be used to derive the Fokker-Planck Equation (also known as the Forward Kolmogorov Equation), which is the equation for the transport/evolution of the joint Lagrangian PDF of the particle velocity and position. This can then be divided by the marginal PDF for the particles position to obtain the equation for the evolution of the conditional PDF of the particle velocity shown above!

Definitions:
- Vector Valued Weiner Process: a vector of Weiner Processes (random motion in multiple directions)

---

## Slide 5: The Connection

Let's take a second to compare the two stochastic differential equations shown thus far. Both the BSM equation and the transported PDF equation are derived from stochastic differential equations of the same drift plus diffusion form. The BSM equation is derived from the SDE describing the price of a stock (often referred to as geometric Brownian motion), whereas the PDF transport equation is derived from the Langevin equation. Regardless, both equations have the same mathematic structure: deterministic drift plus random fluctuations driven by a Wiener process. Both equations are derived using Ito’s Lemma/Stochastic calculus.

The key distinction to be made however is that the BSM equation computes the conditional expectation of a future payoff given the current state, then discounts it back to give the current option price. On the other hand, the transported PDF equation that Pope derived from the Fokker-Planck equation describes how the PDF of velocity evolves over time. In other words, for BSM, you use the PDF at expiry to find the expected payoff, then discount it back to today (this is done implicitly with the SDE hedging described previously, but later we show how this is done explicitly with a Monte Carlo simulation). 

Regardless, both equations arise from the same type of SDE, and we can observe the direct variable mapping above, including the mapping of solution methods.

Some key differences to notice are that the drift and diffusion coefficients are constant for the Brownian Motion SDE, but depend on the particles current position and velocity for the Langevin Equation.

Note that the SDE for the velocity of a particle uses a vector-valued Wiener process because the stochastic perturbations act in multiple spatial directions (typically three-dimensional space).
---

## Slide 6: Solution Methods

It's important to take a moment to talk about both solution methods available. 

The Eulerian approach solves the PDE directly on a fixed grid.

Eulerian:
- Scales poorly to high dimensions, let's say we want to simulate two stocks, each axis with 500 grid points, you would now need 250,000 grid locations to cover all combinations

The Lagrangian approach simulates random particle trajectories.

Lagrangian:
- If portfolio grows, we can just simulate random samples of the full vector (no giant mesh needed)
- Compatible with CPUs/accelerators/clusters because each simulation is independent so they are easy to parallelize
- Need 4x more simulations to cut error in half, and 100x more simulations to cut error by 10x

We should note, that according to Pope: “It is not feasible, computationally, to represent the joint PDF accurately through a discretization of the seven-dimensional V-θ-x space (for the velocity-frequency joint PDF), as is required in finite-difference, finite-volume, and finite-element methods. Instead, the model PDF equations are solved by particle methods.” (Pope). Professor Zhao also mentioned this in our final lecture. 

---

## Slide 7: Results: Option Pricing

To explore the different solution methods, I decided to write a python script to compare the two mentioned above for the BSM equation. The plot shows the option value for a European call versus stock prices. 

The Eulerian approach solves the PDE directly on a fixed grid using implicit finite difference scheme, and second order central difference stencils for the first and second derivative terms. Doing so resulted in a tridiagonal system that was extremely similar to what we saw in Project 01 and Flip Classroom Exercise 01. This idea was inspired by Professor Zhao provided me with a repository for OpenFoam (more specifically finacialFoam), which executes solves the BSM equation using FVM.

The Lagrangian approach simulated 200,000 random stock price trajectories (GBM paths) for 80 stock prices between $50 and $200 using a Monte Carlo simulation, averaged the payoffs, and discounted them, back to present day value.

The results of both approaches are plotted here on top of the analytical solution to the BSM PDE. I have also included the intrinsic value (if the option expired today, what would it be worth at each price) to show how even though some stock prices currently sit below the strike price, they still hold value due to the possibility of their price increasing above the strike price before expiration. 

For this case, the strike price, risk free interest rate, volatility, and contract time period are listed above. 

We can observe excellent agreement between all three solutions!

---

## Slide 8: Results: Monte Carlo Visualization

To help visualize the Monte Carlo simulation, I created this plot to demonstrate the motion of different stock price paths. The plot shows 30 individual Monte Carlo paths (of 100 simulated paths). The bold blue line represents the mean of these paths; and the black dashed line is the theoretical expected return of the stock. As the number of simulations increase, the mean converges to this value.

We can now extract data from time steps in our Monte Carlo simulation to construct the transport of the PDF of the stock price! From this point on, we will focus on the case of a stock with an initial price of $100.

---

## Slide 9: Results: PDF Transport

The 3d surface show above was constructed entirely from Monte Carlo particle data (200,000 simulated GBM paths, histogram med at each time snapshot)

This is exactly how one can imagine the transported PDF method for a turbulent flow. Thousands of particles are simulated using the Langevin SDE, their velocities (conditioned on position and time) are recorded at each time step, and the PDF surface is constructed from the statistics. 

We can observe the PDF starting as a sharp spike near the initial stock price, then spreading and shifting over time due to drift and diffusion. 

---

## Slide 10: Results: PDF Transport 

The following plot was created to better illustrate this concept of drift and diffusion. The shaded bars are pdfs constructed from Monte Carlo data plotted as histograms, and the solid curves are the exact analytical log-normal PDFs. Over time, the volatility spreads the distribution of possible stock prices, and the drift shifts it to the right.

We can conceptualize this for a turbulent flow, where drift can be thought of as where the average velocity tends to go, and diffusion as how strongly turbulence randomizes the velocity around the mean

Hull uses Ito's Lemma to show that, for our GBM model, the price of a stock at a given time follows a lognormal distribution. Therefore, we can use the standard equation for the probability density function of a lognormal variable. For the Langevin SDE, because the coefficients vary with respect to particle velocity/position and time, we don’t have this luxury, so we must construct the PDFs for given time steps using the data from the Monte Carlo simulation. Regardless, both methods are used above, and we can observe excellent agreement between the two.

---

## Slide 11: Results: PDF Transport 

Finally, a top-down heat map view of the PDF evolution can help us visualize the drift of the average stock price at each time step. We clearly see the PDF shifting in line with the expected stock price, which can also be though of as the deterministic growth of the stock in the absence of volatility. 

---

## Slide 12: Monte Carlo Convergence

A quick note on convergence, as mentioned earlier, the standard error of the estimate is expressed as the standard deviation over the square root of the number of trials (to cut error in half, you need 4x more trials). With 100 trials, the estimate is noisy (each path has a larger influence so estimate can jump around) by 1,000,000 paths it approaches the exact value. Again, this is a key tradeoff between Eulerian and Lagrangian approaches.

---

## Slide 13: Conclusions

Finally, to recap what we discussed:

The BSM the Fokker-Planck (forward Kolmogorov equation) are both derived from SDEs of the same form. This directly connects financial derivative pricing to the transported PDF method in CFD.

Second, the same results (validated against the exact analytical formula) can be obtained from both Eulerian and Lagrangian methods as demonstrated with the BSM equation.

Third, Monte Carlo particle simulations reproduce the analytical PDF

Fourth, BSM is actually a great validation problem because closed-form solutions exist for both the SDE and the PDF. In turbulent flows, these may not exist, which is why particle methods become essential.

Finally, Advancements in financial modeling and computational fluid dynamics are mutually beneficial. Numerical methods/variance reduction can transfer directly between CFD and quantitative finance, so techniques/methods developed in either discipline can be utilized by/applied to problems in the other.

Thank you!

---

## Slide 14: References

No speaking necesary for this slide.

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
