# BSM Term Paper — Running Explanation & Learning Notes
## Christian DiPietrantonio | ME 5311

---

# Part 1: What Problem Are We Even Solving?

## 1.1 What Is a Stock?

You probably have some idea, but let's nail it down. A stock is a tiny piece of
ownership in a company. If a company has 1,000 shares of stock and you own 1 share,
you own 0.1% of that company. The **stock price** is what the market says one share
is worth right now. It changes constantly — every time someone buys or sells.

The key thing for us: **stock prices are random.** Nobody knows for sure what
Apple's stock will be worth tomorrow. It might go up, it might go down. This
randomness is the whole reason the math gets interesting.

## 1.2 What Is an Option?

An option is a **contract** that gives you the *right* (but not the obligation) to
buy or sell a stock at a pre-agreed price, on or before a specific date.

**Concrete example — a European Call Option:**

Say Apple stock is trading at $100 today. You buy a contract that says:

> "I have the right to buy 1 share of Apple for $100, but ONLY on the date
> exactly 1 year from now."

The key terms:
- **Strike price (K = $100):** The pre-agreed price you'd pay for the stock
- **Expiry (T = 1 year):** When the contract expires
- **"European":** You can ONLY exercise it at expiry (not before)
- **"Call":** It's the right to BUY (as opposed to a "put," which is the right to sell)

### Why would you want this?

Fast-forward 1 year. Two scenarios:

**Scenario A:** Apple stock is now $150.
- Your contract lets you buy at $100. You immediately sell at $150. Profit = $50.
- The **payoff** = $150 − $100 = $50.

**Scenario B:** Apple stock is now $80.
- Your contract lets you buy at $100... but why would you? It's cheaper on the
  open market. You just don't exercise the option. It expires worthless.
- The **payoff** = $0.

So the payoff at expiry is always:

$$\text{Payoff} = \max(S_T - K, \; 0)$$

where $S_T$ is the stock price at expiry. You either make money or you don't — you
never *lose* money on the option itself (beyond what you paid for the contract).

### The Big Question

Here's the million-dollar question (literally): **How much should this contract
cost TODAY?**

Think about it — the contract clearly has value. If there's a chance the stock
goes above $100, the contract could pay off. But how much is that chance worth?

**This is what the Black-Scholes-Merton equation answers.**

## 1.3 The Inputs We Need

To price the option, we need to know:

| Symbol   | Meaning                        | Our Value  |
|----------|--------------------------------|------------|
| $S$      | Current stock price            | varies     |
| $K$      | Strike price                   | $100       |
| $r$      | Risk-free interest rate        | 5% / year  |
| $\sigma$ | Volatility (how "jumpy" the stock is) | 20% / year |
| $T$      | Time to expiry                 | 1 year     |

**What is r?** The risk-free rate is essentially what you'd earn by putting money
in a savings account or government bonds — a "guaranteed" return with no risk.
Think of it as the baseline growth rate of money. We use r = 0.05 (5%).

**What is σ (sigma)?** Volatility measures how much the stock price bounces around.
A stock with σ = 0.20 (20%) fluctuates about ±20% per year. High σ means the stock
is wild and unpredictable; low σ means it's calm and steady.

Here's the intuition: **higher volatility makes options MORE valuable.** Why?
Because the payoff is max(S − K, 0) — you benefit from the upside but you're
protected from the downside (payoff can't go negative). More randomness = more
chances of a big upside = option is worth more.

---

# Part 2: How Does the Stock Price Move? (Geometric Brownian Motion)

## 2.1 The Analogy: A Drunk Person on an Escalator

Imagine a person standing on an escalator (moving walkway that goes upward slowly).
The escalator carries them upward at a steady rate — that's the **drift** (the
risk-free rate $r$ pulling the expected price upward over time).

But this person is also drunk, so they randomly stumble left and right. Each stumble
is random and unpredictable — that's the **diffusion** (volatility $\sigma$
causing random fluctuations).

The combination of steady upward drift + random stumbles is **exactly** what
Geometric Brownian Motion (GBM) describes for stock prices:

$$dS = r \, S \, dt + \sigma \, S \, dW$$

Let's break this apart piece by piece:

- $dS$: The tiny change in stock price over a tiny time interval
- $r \, S \, dt$: The **drift** term — the stock tends to grow at rate $r$.
  This is deterministic (predictable). Like the escalator.
- $\sigma \, S \, dW$: The **diffusion** term — random fluctuations.
  $dW$ is a Wiener process increment (a fancy way of saying "a random normal
  number scaled by $\sqrt{dt}$"). Like the drunk stumbles.

### Why "Geometric"?

Notice both terms are multiplied by $S$ — the current price. This means:
- A $100 stock might move ±$2 in a day
- A $200 stock might move ±$4 in a day

The fluctuations scale with the price level. This is called **multiplicative noise**,
which is why the resulting distribution is **log-normal** (not normal/Gaussian).

### The CFD Connection (This Is Key for Your Paper!)

This SDE is structurally identical to the Langevin equation that Pope uses for
fluid particle velocities in turbulent flows:

$$dU^*(t) = a(U^*, X^*, t) \, dt + b(X^*, t) \, dW$$

Both have the same form: **drift + diffusion driven by a Wiener process.**
The stock price $S$ plays the role of the particle velocity $U^*$. The risk-free
rate $r$ plays the role of the drift coefficient $a$. The volatility $\sigma$
plays the role of the diffusion coefficient $b$.

## 2.2 How We Simulate GBM in Code (The Exact Solution)

The SDE above has an exact solution. If you know $S$ at time $t$, then at time
$t + \Delta t$:

$$S(t + \Delta t) = S(t) \cdot \exp\!\left[\left(r - \tfrac{1}{2}\sigma^2\right)\Delta t + \sigma \sqrt{\Delta t} \; Z\right]$$

where $Z \sim \mathcal{N}(0,1)$ is a standard normal random number.

**Why $r - \frac{1}{2}\sigma^2$ and not just $r$?** This is a subtlety from
stochastic calculus (Itô's Lemma). When you take the exponential of a random
process, the average of the exponential is NOT the exponential of the average.
The $-\frac{1}{2}\sigma^2$ is a correction term (called the Itô correction) that
ensures the *expected* stock price still grows at rate $r$.

In the code, this is the line:
```python
ST = S0 * np.exp((r - 0.5 * sigma**2) * T + sigma * np.sqrt(T) * Z)
```

---

# Part 3: The Monte Carlo Method (Lagrangian / Stochastic Approach)

## 3.1 The Core Idea

Monte Carlo is beautifully simple in concept:

1. **Simulate many random futures.** Generate thousands of possible stock price
   paths, each one a different random realization of GBM.

2. **Compute the payoff for each path.** At expiry, each path gives a final stock
   price $S_T$. The payoff is $\max(S_T - K, 0)$.

3. **Average the payoffs.** The option price is the average payoff, discounted
   back to today.

That's it. You're literally asking: "If I played out the future 200,000 times,
what would I earn on average?"

### The CFD Connection

This is the **Lagrangian approach** — we track individual "particles" (each stock
price path is a particle trajectory) through the stochastic process, then extract
statistics by averaging over many particles. This is exactly how the **Lagrangian
particle method** works for solving the PDF transport equation in turbulent flows.
Pope's stochastic particle method simulates thousands of fluid particles, each
following their own random trajectory governed by the Langevin SDE, then builds
statistics from the ensemble.

## 3.2 Walking Through the Code

### The Fast Version (Terminal Value Method)

For just computing the option price, we don't need to simulate every time step.
Since we only care about where the stock ends up at expiry, we can jump directly
to $S(T)$ using the exact GBM formula:

```python
def mc_option_price(S0, K, r, sigma, T, N_paths):
    # Step 1: Generate 200,000 random standard normal numbers
    Z = np.random.standard_normal(N_paths)

    # Step 2: Compute the terminal stock price for each path
    # Each Z gives a different random outcome
    ST = S0 * np.exp((r - 0.5 * sigma**2) * T + sigma * np.sqrt(T) * Z)

    # Step 3: Compute payoff for each path
    payoffs = np.maximum(ST - K, 0.0)

    # Step 4: Average and discount back to today
    V_mean = np.exp(-r * T) * np.mean(payoffs)

    # Also compute standard error (how precise is our estimate?)
    V_std = np.exp(-r * T) * np.std(payoffs) / np.sqrt(N_paths)

    return V_mean, V_std
```

**Line by line:**

- `Z = np.random.standard_normal(N_paths)` — Draw 200,000 random numbers from
  $\mathcal{N}(0,1)$. Each one represents a different "future" for the stock.

- `ST = S0 * np.exp(...)` — Apply the exact GBM formula to jump from today's price
  $S_0$ to the expiry price $S_T$. Each random $Z$ gives a different $S_T$.

- `np.maximum(ST - K, 0.0)` — Apply the payoff rule. If $S_T > K$, you make money.
  If $S_T \leq K$, payoff is zero.

- `np.exp(-r * T) * np.mean(payoffs)` — **Discounting.** A dollar one year from now
  is worth less than a dollar today (because you could have invested it and earned
  interest). Multiplying by $e^{-rT}$ converts the future payoff to its present
  value. Then we average over all 200,000 outcomes.

- `np.std(payoffs) / np.sqrt(N_paths)` — Standard error of the mean. This tells us
  how precise our MC estimate is. It decreases as $1/\sqrt{N}$ — to cut the error
  in half, you need 4× more paths.

### The Full Path Version (for visualizations)

For the PDF evolution figures, we need to know where the stock is at every time step,
not just at the end:

```python
def mc_full_paths(S0, r, sigma, T, N_paths, N_steps):
    dt = T / N_steps
    S = np.zeros((N_steps + 1, N_paths))
    S[0, :] = S0  # All paths start at S0

    Z = np.random.standard_normal((N_steps, N_paths))

    for i in range(N_steps):
        S[i+1, :] = S[i, :] * np.exp((r - 0.5*sigma**2)*dt
                                       + sigma*np.sqrt(dt)*Z[i, :])
    return S, t
```

This is the same formula, but applied step-by-step (252 steps = 252 trading days
in a year). At each step, we multiply the current price by a random growth factor.
The result is a matrix of shape (253 × N_paths) — every row is a time step, every
column is one path.

---

# Part 4: The Analytical Black-Scholes Formula

## 4.1 The Closed-Form Solution

Black, Scholes, and Merton showed that for a European call option, there's an exact
formula (no simulation needed):

$$V(S, T) = S \cdot N(d_1) - K e^{-rT} \cdot N(d_2)$$

where:

$$d_1 = \frac{\ln(S/K) + (r + \tfrac{1}{2}\sigma^2) T}{\sigma \sqrt{T}}, \qquad d_2 = d_1 - \sigma\sqrt{T}$$

and $N(\cdot)$ is the standard normal CDF (the probability that a standard normal
random variable is less than or equal to the argument).

### Intuition for the Formula

Don't let the formula intimidate you. Here's what it's actually saying:

- $N(d_2)$ is approximately the **probability that the option expires in the money**
  (i.e., $S_T > K$). If this probability is high, the option is valuable.

- $K e^{-rT} N(d_2)$ is the present value of what you'd pay (the strike price),
  weighted by the probability you'll actually have to pay it.

- $S \cdot N(d_1)$ is the present value of what you'd receive (the stock),
  weighted by a slightly adjusted probability.

- The option value = (what you expect to receive) − (what you expect to pay).

### In Code

```python
def bs_call(S, K, r, sigma, T):
    d1 = (np.log(S / K) + (r + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    V  = S * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)
    return V
```

This is our **ground truth**. Both the Monte Carlo simulation and the OpenFOAM
finite-volume solver should produce the same answer.

### The CFD Connection

The analytical formula exists because the BSM PDE:

$$\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} + rS\frac{\partial V}{\partial S} - rV = 0$$

happens to be solvable in closed form for a European call with the payoff boundary
condition. This is analogous to how the heat equation has analytical solutions for
simple geometries and boundary conditions. The moment you change the geometry (e.g.,
American options, exotic payoffs), you need numerical methods — just like in CFD.

---

# Part 5: Understanding the Figures

## Figure 1: Option Value vs. Stock Price (V(S) Curve)

**What you're looking at:** The x-axis is the current stock price $S$. The y-axis
is how much the option contract is worth today, $V$.

**The black curve** is the analytical Black-Scholes formula — the exact answer.

**The red squares** are the Finite Difference (Eulerian) solver — our Python
implementation that solves the BSM PDE directly on a grid of stock prices.

**The blue dots** are Monte Carlo estimates, each computed by simulating 200,000
random futures for that starting price. The tiny error bars show the statistical
uncertainty (±2 standard errors ≈ 95% confidence).

**The gray dashed line** is the "intrinsic value" = max(S − K, 0). This is what
the option would be worth if it expired RIGHT NOW. The actual option value is always
above this line because there's still time for the stock to move favorably —
this difference is called **time value**.

**Key observations:**
- When S ≪ K (stock price way below strike): The option is nearly worthless. If the
  stock is at $60 and you need it above $100 to make money, that's unlikely. V ≈ $0.05.
- When S ≈ K (at-the-money): The option has significant value due to the possibility
  of ending in the money. V ≈ $10.45 when S = $100.
- When S ≫ K (deep in-the-money): The option approaches intrinsic value. If the
  stock is already at $160, you're almost certainly going to exercise it.
- **MC matches analytical perfectly** — the stochastic method reproduces the
  deterministic solution. This is the core validation.

## Figure 2: Sample GBM Paths

**What you're looking at:** 30 individual stock price trajectories, all starting
at S₀ = $100, each evolving randomly over one year.

**The black dashed line** is the expected (average) path: $E[S(t)] = S_0 e^{rt}$.
If there were no randomness (σ = 0), every path would follow this line. It's
the escalator.

**The colored spaghetti lines** are individual random paths — the drunk stumbles.
Notice how they fan out over time — early on they're all bunched near $100, but
by t = 1 year, they're spread from about $70 to $180.

**The red dotted line** is the strike price K = $100. Paths that end above this
line at T = 1 have positive payoff; paths below it expire worthless.

**Why this matters:** Each of these paths is one "particle trajectory" in the
Lagrangian sense. The Monte Carlo method runs 200,000 of these and averages. You
can see why more paths = more accurate: with only 30 paths, your average would be
noisy. With 200,000, the noise averages out.

## Figure 3: PDF Evolution (3D Surface)

**What you're looking at:** Instead of plotting individual paths, we're plotting the
*probability density* of where the stock price is at each moment in time.

- **x-axis:** Stock price S
- **y-axis:** Time t
- **z-axis (height):** Probability density — how likely the stock is to be at
  that price at that time

**At very early times (near t = 0):** The PDF is a sharp, tall spike at S = $100.
We know exactly where the stock started, so the probability is concentrated there.
This is essentially a delta function.

**As time progresses:** The spike spreads out and shifts. This is **exactly** what
you see when you release a puff of smoke in a flow — it drifts with the mean
velocity (advection) and spreads due to molecular/turbulent diffusion.

- **Drift (advection):** The peak shifts slightly to the right because the expected
  growth rate is positive ($r = 5\%$).
- **Diffusion:** The distribution gets wider and flatter because uncertainty
  accumulates over time (σ effect).
- **Skewness:** The distribution develops a right-skewed shape (long right tail).
  This is because S follows a **log-normal** distribution — it can grow to infinity
  but can never go below zero.

## Figure 4: PDF Time Slices (2D Overlay) ⭐ Best Figure for Presentation

**What you're looking at:** The same PDF evolution, but "sliced" at specific times
and overlaid on a single 2D plot. Each colored curve is the probability distribution
of the stock price at that time.

**At t = 0.02 yr (≈ 1 week):** Tall, narrow, almost symmetric peak centered near
$100. Very little uncertainty has accumulated.

**At t = 0.25 yr (3 months):** Shorter and wider. Starting to develop the
right-skew.

**At t = 1.0 yr:** Short, wide, clearly right-skewed. The stock could be anywhere
from $50 to $250, but most likely near $100–$110.

**The annotations "Drift →" and "Diffusion (spreading)"** are the punchline for
your CFD audience. You can point at this figure and say: *"This is advection and
diffusion of a probability density function — the same physics we see in scalar
transport in turbulent flows."*

## Figure 5: Monte Carlo Convergence

**What you're looking at:** How the Monte Carlo estimate improves as we use more
simulation paths.

- **x-axis:** Number of paths N (log scale)
- **y-axis:** Estimated option value
- **Red dashed line:** Exact Black-Scholes answer ($10.4506)
- **Blue dots:** MC estimate for each N
- **Blue shaded region:** ±2σ confidence interval

**Key observation:** With N = 100 paths, the estimate is all over the place
(anywhere from $6 to $15). By N = 10,000, it's close. By N = 1,000,000, it's
essentially exact.

**The convergence rate is $1/\sqrt{N}$** — this is a fundamental property of Monte
Carlo methods. To cut your error in half, you need 4× the number of paths. This
is one of the *drawbacks* of the Lagrangian/stochastic approach, and it's the same
issue in CFD: particle methods converge slowly and need many samples.

The Eulerian/grid-based approach (finite difference, finite volume — what OpenFOAM
does) doesn't have this statistical noise. It solves the PDE directly on a grid.
The tradeoff is that grid-based methods scale poorly to high dimensions (the
"curse of dimensionality"), while Monte Carlo handles high-dimensional problems
naturally. For our 1D BSM case, the grid method is more efficient. But for pricing
options on 100 correlated stocks? Monte Carlo wins.

## Figure 6: PDF Heatmap

**What you're looking at:** The same data as Figure 3, but viewed from above as a
2D color map. Think of it as a bird's-eye view of the 3D surface.

- **x-axis:** Stock price S
- **y-axis:** Time t (increasing upward)
- **Color:** Probability density (bright = high probability, dark = low)

**The white dashed line** is the expected path $E[S(t)] = S_0 e^{rt}$ — this is
the "drift" or advection velocity of the PDF.

**The white dotted lines** are the ±1σ bounds of the log-normal distribution —
they show the "diffusion cone" widening over time.

**CFD analogy:** This looks exactly like what you'd see if you simulated a puff of
scalar (dye, temperature, concentration) released into a flow with a mean velocity
and diffusivity. The center of the puff drifts downstream while it simultaneously
spreads due to diffusion. That's because *it is* the same math — both are governed
by advection-diffusion transport equations.

---

# Part 6: Connecting Everything to CFD (The Big Picture)

Here is the central argument of your paper, summarized:

## The Two Solution Approaches

| Approach       | CFD Name        | Finance Name    | What It Does                            |
|----------------|-----------------|------------------|-----------------------------------------|
| **Eulerian**   | FVM / FDM on grid | Finite Difference / OpenFOAM | Solve the PDE directly on a grid of stock prices |
| **Lagrangian** | Particle method  | Monte Carlo      | Simulate random trajectories, average the results |

## Why They Give the Same Answer

Both approaches solve the **same underlying problem**, just from different angles:

- The **Eulerian** approach solves the BSM PDE (a backward Kolmogorov equation)
  directly on a grid. The unknown is the option value $V(S,t)$ at every grid point.

- The **Lagrangian** approach simulates the SDE (Geometric Brownian Motion) that
  generates the stock price paths, then computes expectations. The unknown is
  reconstructed from statistics over many particles.

The mathematical guarantee that they agree comes from the **Feynman-Kac theorem**,
which proves that the solution to certain PDEs can be represented as an expected
value over stochastic paths. This is the same connection that links the Fokker-Planck
equation (Eulerian PDF transport) to the Langevin equation (Lagrangian particle model)
in Pope's framework.

## The Variable Mapping

| GBM / BSM (Finance)       | Langevin / PDF Transport (CFD)           |
|---------------------------|------------------------------------------|
| Stock price $S$           | Particle velocity $U^*$                  |
| Risk-free rate $r$        | Drift coefficient $a$                    |
| Volatility $\sigma$       | Diffusion coefficient $b$                |
| Option value $V(S,t)$     | Conditional expected value               |
| BSM PDE                   | Backward Kolmogorov equation             |
| Log-normal PDF of $S$     | Fokker-Planck (forward Kolmogorov) equation |
| Monte Carlo simulation    | Lagrangian particle method               |
| Finite difference solver  | Finite volume solver (OpenFOAM)          |

---

# Part 7: The Finite Difference Method (Eulerian / Grid-Based PDE Solver)

## 7.1 The Core Idea

The Monte Carlo method asked: "What happens if I simulate thousands of random
futures?" The finite difference method takes the opposite approach: "Let me solve
the PDE directly on a grid."

### The Analogy: Weather Forecasting

Think of two ways to predict tomorrow's temperature:

**Lagrangian (Monte Carlo):** Release thousands of little weather balloons. Each
one drifts and bounces randomly through the atmosphere. Collect them all tomorrow,
see what temperature each one measured, and average the results.

**Eulerian (Finite Difference):** Set up a grid of thermometers across the entire
region. Use the heat equation to compute how temperature flows from one grid point
to the next. March forward in time on the grid.

Both give you the same answer. The Eulerian approach solves the governing PDE
directly on a fixed grid — no randomness involved. This is what OpenFOAM does for
Navier-Stokes, and what our FD solver does for Black-Scholes.

## 7.2 The PDE We're Solving

The BSM PDE is:

$$\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} + rS\frac{\partial V}{\partial S} - rV = 0$$

This is a **terminal value problem** — we know the answer at the END (the payoff
at expiry) and we need to work backward to find the answer NOW. This might feel
weird, but think of it this way: at expiry, you know exactly what the option is
worth (just the payoff). The question is: what is that future payoff worth today?

### The Time Reversal Trick

To make the PDE look like a normal initial value problem (which is what we're
used to solving in CFD), we introduce:

$$\tau = T - t \quad \text{(time to expiry)}$$

When $\tau = 0$, we're at expiry (the "initial" condition we know).
When $\tau = T$, we're at the present (the answer we want).

The PDE becomes (with $\tau$ as our time-like variable):

$$\frac{\partial V}{\partial \tau} = \underbrace{\frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2}}_{\text{diffusion}} + \underbrace{rS\frac{\partial V}{\partial S}}_{\text{advection}} - \underbrace{rV}_{\text{reaction/decay}}$$

Now this is a **forward-in-time** advection-diffusion-reaction equation — the
same type you've been solving all semester! The only quirk is that the diffusion
coefficient $D = \frac{1}{2}\sigma^2 S^2$ and the advection velocity $u = rS$
both depend on the spatial variable $S$.

### Initial and Boundary Conditions

**Initial condition** (at $\tau = 0$, i.e., at expiry):
$$V(S, 0) = \max(S - K, 0)$$
This is just the payoff — the hockey-stick function you see in Figure 1 (the
gray dashed line).

**Left boundary** (at $S = 0$):
$$V(0, \tau) = 0$$
If the stock is worthless, the option to buy it is worthless.

**Right boundary** (at $S = S_{\max}$, very large $S$):
$$V(S_{\max}, \tau) \approx S_{\max} - Ke^{-r\tau}$$
If the stock price is huge, the option is almost certainly going to be exercised,
so its value is approximately the stock price minus the discounted strike price.

## 7.3 Discretization — Turning the PDE into Algebra

### The Grid

We lay down a uniform grid on the stock price axis:

$$S_j = j \cdot \Delta S, \quad j = 0, 1, 2, \ldots, M$$

where $\Delta S = S_{\max} / M$. Each grid point $j$ represents a specific stock
price. This is exactly like the $\eta$ grid in your boundary layer solver from
Project 2 — instead of height in the boundary layer, it's stock price.

### Spatial Derivatives (Central Differences — Same as Project 2!)

At interior node $j$, we approximate:

$$\frac{\partial V}{\partial S}\bigg|_j \approx \frac{V_{j+1} - V_{j-1}}{2\Delta S}$$

$$\frac{\partial^2 V}{\partial S^2}\bigg|_j \approx \frac{V_{j+1} - 2V_j + V_{j-1}}{\Delta S^2}$$

These are the same central difference stencils you used in P1 and P2.

### Substitution and Simplification

Now here's a nice trick. Since $S_j = j \cdot \Delta S$, when we plug into the
PDE:

$$\frac{1}{2}\sigma^2 S_j^2 \cdot \frac{V_{j+1} - 2V_j + V_{j-1}}{\Delta S^2} = \frac{1}{2}\sigma^2 (j\Delta S)^2 \cdot \frac{V_{j+1} - 2V_j + V_{j-1}}{\Delta S^2}$$

The $\Delta S^2$ cancels:

$$= \frac{1}{2}\sigma^2 j^2 (V_{j+1} - 2V_j + V_{j-1})$$

Similarly for the advection term:

$$rS_j \cdot \frac{V_{j+1} - V_{j-1}}{2\Delta S} = r(j\Delta S) \cdot \frac{V_{j+1} - V_{j-1}}{2\Delta S} = \frac{1}{2}rj(V_{j+1} - V_{j-1})$$

This is nice — the grid spacing $\Delta S$ drops out entirely!

### The Semi-Discrete Equation

Combining everything, at each interior node $j$:

$$\frac{dV_j}{d\tau} = \frac{1}{2}\sigma^2 j^2 (V_{j+1} - 2V_j + V_{j-1}) + \frac{1}{2}rj(V_{j+1} - V_{j-1}) - rV_j$$

### Implicit Euler Time Stepping

We use implicit Euler (backward Euler) — the same philosophy as the implicit
methods in your CFD projects. Evaluating the right-hand side at the NEW time level
$n+1$:

$$\frac{V_j^{n+1} - V_j^n}{\Delta\tau} = \frac{1}{2}\sigma^2 j^2 (V_{j+1}^{n+1} - 2V_j^{n+1} + V_{j-1}^{n+1}) + \frac{1}{2}rj(V_{j+1}^{n+1} - V_{j-1}^{n+1}) - rV_j^{n+1}$$

Why implicit? Same reason as in CFD: **unconditional stability.** An explicit
scheme would require $\Delta\tau < \Delta S^2 / (\sigma^2 S_{\max}^2)$, which for
our parameters would mean thousands upon thousands of tiny time steps. Implicit
lets us take bigger steps.

### Collecting into a Tridiagonal System

Rearranging so all $n+1$ terms are on the left and the known $n$ term is on the right:

$$a_j \, V_{j-1}^{n+1} \;+\; b_j \, V_j^{n+1} \;+\; c_j \, V_{j+1}^{n+1} \;=\; V_j^n$$

where:

$$a_j = -\Delta\tau \left(\frac{1}{2}\sigma^2 j^2 - \frac{1}{2}rj\right) \quad \text{(sub-diagonal)}$$

$$b_j = 1 + \Delta\tau \left(\sigma^2 j^2 + r\right) \quad \text{(main diagonal)}$$

$$c_j = -\Delta\tau \left(\frac{1}{2}\sigma^2 j^2 + \frac{1}{2}rj\right) \quad \text{(super-diagonal)}$$

**This is a tridiagonal system** — exactly like what you solved in Projects 1 and 2
with the banded solver! The matrix has non-zero entries only on the main diagonal
and the two adjacent diagonals.

At each time step, we solve:

$$\begin{bmatrix} b_1 & c_1 & & \\ a_2 & b_2 & c_2 & \\ & \ddots & \ddots & \ddots \\ & & a_{M-1} & b_{M-1} \end{bmatrix} \begin{bmatrix} V_1^{n+1} \\ V_2^{n+1} \\ \vdots \\ V_{M-1}^{n+1} \end{bmatrix} = \begin{bmatrix} V_1^n \\ V_2^n \\ \vdots \\ V_{M-1}^n \end{bmatrix}$$

(with the boundary conditions folded into the RHS for $j=1$ and $j=M-1$, exactly
the same technique from your projects — known boundary values get moved to the
right-hand side.)

## 7.4 Walking Through the Code

```python
def fd_bsm(K, r, sigma, T, S_max=300.0, M=200, N_t=2000):
```

**Parameters:**
- `M = 200`: 201 grid points in stock price (like your η grid in P2)
- `N_t = 2000`: 2000 time steps marching from τ=0 to τ=T
- `S_max = 300`: We truncate the infinite S domain at $300 (3× the strike price)

```python
    dS  = S_max / M           # spatial step size
    dtau = T / N_t            # time step size
    S_grid = np.linspace(0, S_max, M + 1)  # S_0, S_1, ..., S_M

    # Initial condition: the payoff at expiry
    V = np.maximum(S_grid - K, 0.0)
```

We start with V being the hockey-stick payoff function. This is our "initial
condition" in τ (= terminal condition in real time t).

```python
    # Interior node indices: j = 1, 2, ..., M-1
    j = np.arange(1, M)

    # Tridiagonal coefficients
    a_j = -dtau * (0.5 * sigma**2 * j**2 - 0.5 * r * j)   # sub-diagonal
    b_j = 1.0 + dtau * (sigma**2 * j**2 + r)               # main diagonal
    c_j = -dtau * (0.5 * sigma**2 * j**2 + 0.5 * r * j)   # super-diagonal
```

These are built once and reused at every time step (the coefficients don't change
because the PDE has constant σ and r). Notice how j appears — remember, j encodes
the stock price $S_j = j \Delta S$, so the spatial dependence of the diffusion and
advection coefficients is captured through $j^2$ and $j$.

```python
    # Pack into banded storage (same format as your P2 solver!)
    ab = np.zeros((3, n_interior))
    ab[0, 1:] = c_j[:-1]    # super-diagonal (shifted by 1)
    ab[1, :]  = b_j          # main diagonal
    ab[2, :-1] = a_j[1:]     # sub-diagonal (shifted by 1)
```

This is the banded storage format for `scipy.linalg.solve_banded` — you've done
this before. Row 0 is the super-diagonal, row 1 is the main diagonal, row 2 is
the sub-diagonal.

```python
    # Time-march forward in tau (backward in real time)
    for n in range(N_t):
        rhs = V[1:M].copy()

        # Right boundary contribution
        tau_n = (n + 1) * dtau
        V_right = S_max - K * np.exp(-r * tau_n)
        rhs[-1] -= c_j[-1] * V_right  # Move known BC to RHS

        # Solve tridiagonal system
        V[1:M] = solve_banded((1, 1), ab, rhs)

        # Apply boundary conditions
        V[0]  = 0.0
        V[M]  = V_right
```

**The time loop:** At each time step:
1. Copy the current interior values as the RHS
2. Compute the right boundary value (which changes with τ)
3. Move the known boundary value to the RHS (substitution — just like P2!)
4. Solve the tridiagonal system
5. Apply boundary conditions

After all N_t steps, V contains the option value at τ = T (which is t = 0 — today).

### What We Get

The FD solver returns V(S) at every grid point — a smooth curve across all stock
prices, computed in one pass. No randomness, no statistical noise, no need to run
it multiple times. This is the Eulerian advantage.

The spot-check table confirms it matches both the analytical formula and the
Monte Carlo to excellent precision:

```
   S ($)    Analytical    Monte Carlo    FD (Euler)
      60        0.0544        0.0528        0.0549
     100       10.4506       10.4653       10.4478
     160       64.9130       64.8843       64.9131
```

---

# Part 8: What OpenFOAM / financialFoam Does (Preview)

The `financialFoam` solver in OpenFOAM solves the BSM PDE using the **finite volume
method** — the same method used for Navier-Stokes in CFD. It treats:

- The stock price axis $S$ as the spatial dimension (like $x$ in a 1D flow)
- Time $t$ as... time (but reversed — BSM runs backward from expiry)
- The option value $V(S,t)$ as the transported scalar

Looking at the source code your professor shared:

```cpp
solve
(
    fvm::ddt(V)                    // ∂V/∂t
  + fvm::div(phi, V)              // advection: rS·∂V/∂S
  - fvm::Sp(fvc::div(phi), V)     // source correction
  - fvm::laplacian(DV, V)         // diffusion: ½σ²S²·∂²V/∂S²
 ==
  - fvm::Sp(r, V)                 // reaction: -rV (discounting)
);
```

Each line maps to a term in the BSM PDE. The advection velocity is set by the
risk-free rate $r$ and the diffusion coefficient is $D_V = \frac{1}{2}\sigma^2 S^2$.

financialFoam does essentially the same thing as our Python FD solver, but using
the **finite volume method** instead of finite differences, and with all of
OpenFOAM's industrial-strength machinery (arbitrary meshes, parallel computing,
flexible time-stepping, etc.). For our 1D problem, the results should be
identical.

If we set up financialFoam with the same parameters (K=$100, r=5%, σ=20%, T=1yr),
its V(S) curve would overlay directly on our Figure 1 — adding a fourth set of
points right on top of the other three.

---

## Updated Figure 1

Figure 1 now shows all three methods:
- **Black curve:** Analytical Black-Scholes (exact)
- **Red squares:** Finite Difference / Eulerian (our Python PDE solver)
- **Blue dots:** Monte Carlo / Lagrangian (stochastic simulation)

All three land on top of each other — confirming that the Eulerian and Lagrangian
approaches to solving the BSM equation give identical results, just as they do
for the PDF transport equation in CFD.

---

*This document will continue to grow as we build the presentation and final paper.*
