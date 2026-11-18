# Molecular Dynamics — Lecture 11: Detailed Study Notes

*Concepts, full derivations, and my questions answered, built from the handwritten lecture notes (Molecular dynamics, 8 Feb 2022).*

Equations are written in LaTeX (renders in Obsidian, Typora, VS Code preview, GitHub, Jupyter). Matrices are also given in plain text in case your viewer doesn't render them. Figures are in the `figures/` folder next to this file.

---

## Contents

1. [What an MD simulation does](#1-what-an-md-simulation-does)
2. [Time average vs ensemble average (ergodicity)](#2-time-average-vs-ensemble-average-ergodicity)
3. [MD compared with Monte Carlo](#3-md-compared-with-monte-carlo)
4. [Basic structure of an MD code](#4-basic-structure-of-an-md-code)
5. [Equations of motion and phase space](#5-equations-of-motion-and-phase-space)
6. [What makes a "good" integrator](#6-what-makes-a-good-integrator)
7. [Area preservation explained from scratch](#7-area-preservation-explained-from-scratch)
8. [Taylor series: the tool behind every integrator](#8-taylor-series-the-tool-behind-every-integrator)
9. [Velocity-Verlet](#9-velocity-verlet)
10. [Verlet](#10-verlet)
11. [Leap-frog](#11-leap-frog)
12. [All three integrators from Taylor series, side by side](#12-all-three-integrators-from-taylor-series-side-by-side)
13. [Proof: velocity-Verlet is time-reversible](#13-proof-velocity-verlet-is-time-reversible)
14. [The Jacobian determinant from scratch](#14-the-jacobian-determinant-from-scratch)
15. [Proof: velocity-Verlet is area-preserving](#15-proof-velocity-verlet-is-area-preserving)
16. [Worked example: Lennard-Jones Jacobian with real numbers](#16-worked-example-lennard-jones-jacobian-with-real-numbers)
17. [The Jacobian for a real 3D system (6N × 6N)](#17-the-jacobian-for-a-real-3d-system-6n--6n)
18. [Verlet and leap-frog are also reversible and area-preserving](#18-verlet-and-leap-frog-are-also-reversible-and-area-preserving)
19. [Computation of forces](#19-computation-of-forces)
20. [Generating initial velocities](#20-generating-initial-velocities)
21. [Why Maxwell–Boltzmann? (from counting)](#21-why-maxwellboltzmann-from-counting)
22. [How LAMMPS generates velocities](#22-how-lammps-generates-velocities)
23. [Corrections to the handwritten notes](#23-corrections-to-the-handwritten-notes)
24. [Self-check questions with answers](#24-self-check-questions-with-answers)
25. [Quick-reference summary](#25-quick-reference-summary)

---

## 1. What an MD simulation does

> *"In a molecular dynamics simulation, the system configuration, i.e. positions and momenta of all atoms, are evolved according to Newton's 2nd law of motion."*

### The configuration needs positions *and* momenta

Knowing where every atom is isn't enough to predict what happens next; you also need how fast and in which direction each is moving. For $N$ atoms in 3D:

| Quantity | Count |
|---|---|
| Positions $(x, y, z)$ per atom | $3N$ |
| Momenta $(p_x, p_y, p_z)$ per atom, $p = mv$ | $3N$ |
| **Total** | $6N$ |

A complete list of these $6N$ numbers is a **microstate**: a full snapshot. Given it, Newton's laws fix the whole future and past. (Notes say "momenta" here and "velocities" later; since $p = mv$ with fixed $m$ they carry the same information.)

### Evolving the microstate

From one microstate: compute forces from positions → use $F = ma$ to change velocities → use velocities to change positions → advance by $\Delta t$ → repeat millions of times. The sequence of microstates is the **trajectory**.

### Time averages

Measurable quantities are averages over the trajectory:

$$
\langle M \rangle = \frac{1}{\tau}\int_0^{\tau} M(t)\,dt
$$

- $M(t)$: a property evaluated on the snapshot at time $t$ (e.g. potential energy, pressure).
- $\tau$: total run length.

In a simulation time is discrete, so this is an ordinary average over saved snapshots:

$$
\langle M \rangle \approx \frac{1}{n}\sum_{k=1}^{n} M(t_k)
$$

(This is what you do when you average a column of LAMMPS thermo output.)

We average because instantaneous values fluctuate strongly (pressure in a small box can swing by hundreds of bar step to step), while experiments measure over long times and many molecules.

---

## 2. Time average vs ensemble average (ergodicity)

> *"According to postulates of statistical thermodynamics, Time average = Ensemble average."*

**Ensemble average:** imagine a huge collection of copies of the system, all at the same macroscopic conditions (e.g. same $N, V, T$) but each in a different microstate. Each microstate has a probability (canonical: $\propto e^{-E/k_BT}$). The ensemble average is the probability-weighted average of $M$ over all of them.

| | Time average (MD) | Ensemble average (theory) |
|---|---|---|
| Averaged over | one system followed through time | many copies, each in a different microstate |
| States chosen by | whatever the dynamics visits | probability weighting |

**Ergodic hypothesis:** followed long enough, one system visits every accessible microstate, spending time in each in proportion to its probability. Then the two averages are equal.

*Analogy:* average mood of a city — survey 10,000 people at noon (ensemble) or follow one person for ten years (time). Equal only if that one person eventually experiences everything the city offers, in the right proportions.

### When it fails

- **Run too short:** system trapped in one region (e.g. supercooled liquid that hasn't crystallised).
- **High barriers:** crossing takes µs, you simulate ns — the nucleation problem.
- **Non-ergodic systems:** glasses on lab timescales; a perfectly harmonic crystal with uncoupled modes.

### ❓ My question: "If the system is stuck in one basin, isn't the ensemble average also dominated by high-probability basins, so the two averages agree?"

Only if (1) one basin holds essentially all the probability **and** (2) the system happens to be stuck in *that* basin. Which basin MD gets stuck in depends on **where it started**, not which basin is most probable.

**Numerical counter-example:**

| Basin | Probability | $\langle M\rangle$ inside basin |
|---|---|---|
| A | 0.6 | 1 |
| B | 0.4 | 3 |

True ensemble average: $0.6\times1 + 0.4\times3 = 1.8$.
Stuck in A → time average 1. Stuck in B → 3. Neither equals 1.8 — even being stuck in the *more probable* basin gives the wrong answer.

**Real case:** a supercooled melt. The crystal is the high-probability basin, but MD starting from the liquid stays liquid (nucleation is rare). The time average gives liquid properties — the system is stuck in the *low*-probability basin.

**What a stuck run does give:** a **restricted ensemble average** over that basin only. That is exactly how *metastable* properties are defined ("density of the supercooled liquid"). Useful — but not the true equilibrium value. Note: basin probability is set by **free energy** (including entropy), not energy alone.

### ❓ Self-check discussed: one 1 ns run vs ten independent 100 ps runs

- The single 1 ns run is the **time average**.
- Ten runs from different random velocities, averaged, is *closer in spirit* to an ensemble average (each is itself a short time average).
- They agree if: (1) the system is ergodic on these timescales, and (2) each run is equilibrated and longer than the relaxation time. More copies only reduce scatter; they don't fix a trapped system — if all ten start in the same basin, all ten stay there and agree with each other while all missing the true answer.

---

## 3. MD compared with Monte Carlo

### What Metropolis MC does

No time, no velocities — only positions:

1. Pick one atom at random; move it by a small random displacement.
2. Compute the energy change $\Delta U$.
3. Accept if $\Delta U < 0$; if $\Delta U > 0$ accept with probability $e^{-\Delta U/k_BT}$.
4. Repeat.

The acceptance rule makes configurations appear in proportion to $e^{-U/k_BT}$, so MC computes **ensemble averages directly**. MD gets there indirectly (dynamics + ergodicity).

### The four points in the notes

**1. MD contains dynamics.** MC sequences aren't physical time. Only MD gives:
- (a) **Time correlation functions**, e.g. velocity autocorrelation $\langle \mathbf v(0)\cdot\mathbf v(t)\rangle$.
- (b) **Transport properties**: diffusion coefficient, viscosity, thermal conductivity — via Green–Kubo integrals of correlation functions, or for diffusion $D = \lim_{t\to\infty}\mathrm{MSD}/6t$.

**2. MD parallelises easily.** In one step all forces depend only on current positions → split the box among processors (domain decomposition); they only exchange boundary atoms. Metropolis moves are sequential (move 2 depends on acceptance of move 1).

**3. MD cannot change $N$.** Newton's laws don't create/destroy atoms, so plain MD can't do:
- grand canonical ($\mu VT$) — adsorption
- semi-grand — alloys, species swaps
- reaction ensemble
- osmotic, Gibbs (particle transfer between two boxes)

MC does these with insertion/deletion/swap moves. Hybrid: LAMMPS `fix gcmc`, `fix atom/swap` inside an MD run.

**4. MC is simpler for discontinuous potentials** (hard sphere, square well). MD needs $F = -dU/dr$, infinite at a jump. MC needs only $\Delta U$. (Hard-sphere MD exists but is *event-driven*, jumping collision to collision.)

| | MD | MC |
|---|---|---|
| Moves follow | Newton's laws | random trials |
| Dynamics/transport | ✅ | ❌ |
| Parallelises easily | ✅ | harder |
| $N$ can change | ❌ | ✅ |
| Discontinuous potentials | awkward | easy |

### ❓ My answer: diffusion of Cr³⁺ in a melt → MD; equilibrium composition at fixed $\mu$ → GCMC

✅ Correct. Practical addition: in a solid, atoms must relax around inserted atoms, so use **hybrid GCMC/MD** (some MC insertion/deletion attempts, then MD steps, repeated).

---

## 4. Basic structure of an MD code

1. **Specify conditions** — density, volume, temperature; this chooses the ensemble (NVE, NVT, NPT) and units. $V = N/\rho$.
2. **Read or generate initial positions** — lattice, random placement, or file. Overlapping atoms → huge repulsion → blow-up. Usually `minimize` first.
3. **Generate initial velocities** — after positions; velocities only set kinetic energy ([§20](#20-generating-initial-velocities)).
4. **Compute forces** — *most time-consuming step*. $N(N-1)/2$ pairs (10,000 atoms → ~50 million pairs per step). Reduced to ~$O(N)$ with cutoffs + neighbour lists. For ML potentials (MACE), this step is essentially the entire cost — which is why integrators are judged by **force calls per step**.
5. **Update positions** — integrate the 2nd-order ODE. $\Delta t$ ≈ 1/10–1/20 of the fastest vibration period (≈1 fs, less with H).
6. **Repeat 4 & 5** — $10^6$ steps × 1 fs = **1 ns** only.
7. **Save configurations** periodically — for RDF, MSD, correlations. Save more often than the timescale you want to resolve.

Not in the notes but essential: **equilibration** (discard) then **production** (collect).

| Step | LAMMPS |
|---|---|
| 1 | `units`, `boundary`, `fix nve/nvt/npt` |
| 2 | `read_data` or `lattice` + `create_atoms`, then `minimize` |
| 3 | `velocity all create T seed` |
| 4 | `pair_style`, `pair_coeff`, `neighbor` |
| 5 | `timestep` + integrator inside the `fix` |
| 6 | `run N` |
| 7 | `dump`, `thermo` |

### ❓ "My simulation explodes in the first steps (energy → 10⁸, lost atoms). Which step?"

| Cause | Step |
|---|---|
| Overlapping atoms / bad or skipped minimisation | 2 (most common) |
| Time step too large | 5 |
| Wrong potential: units, `pair_coeff`, species-to-type mapping | 4 |
| Velocities far too high (typo 30000 vs 300) | 3 |

**Debug trick:** `thermo 1` for the first ~100 steps. Huge energy **at step 0** → structure or potential (nothing has moved yet). Fine at step 0, blows up later → time step.

---

## 5. Equations of motion and phase space

Newton's second law for every atom:

$$
\mathbf f_i = m_i\,\mathbf a_i, \qquad \forall\, i = 1, \dots, N
$$

With $\mathbf a_i = d^2\mathbf r_i/dt^2$:

$$
\frac{d^2\mathbf r_i}{dt^2} = \frac{\mathbf f_i}{m_i} \quad\longleftarrow\quad 3N \text{ second-order ODEs}
$$

- **3N**: $N$ atoms × 3 directions.
- **Second-order**: highest derivative is $d^2/dt^2$.
- **Coupled**: $\mathbf f_i$ depends on *all* positions.
- A second-order ODE needs **two** initial conditions per coordinate (position and velocity) — that's why steps 2 *and* 3 are both needed.

Introduce velocity as a separate variable:

$$
\frac{d\mathbf r_i}{dt} = \mathbf v_i, \qquad \frac{d\mathbf v_i}{dt} = \frac{\mathbf f_i}{m_i} \quad\longleftarrow\quad 6N \text{ first-order ODEs}
$$

Same physics, but now in the standard form $d\mathbf y/dt = \mathbf g(\mathbf y)$ that every numerical method (Euler, RK4, Verlet) steps forward.

### Phase space

The $6N$ numbers form one point in a $6N$-dimensional space — **phase space**. As time passes the point traces the trajectory. For one particle in 1D: a 2D plane, $x$ across, $v$ up.

Properties of the exact motion:
1. **Deterministic** — trajectories never cross.
2. **Constant energy** — the point stays on a constant-energy surface (an ellipse for a 1D spring).
3. **Reversible** — flip all velocities and it retraces its path.

> *"One would require a 'good' algorithm to integrate the laws of motion."*

### ❓ Count for 1000 water molecules (3 atoms each, fully flexible)

- Atoms: 3000
- Second-order ODEs: $3 \times 3000 = 9000$
- First-order ODEs: 18000
- Phase-space dimensions: $6N = 18000$
- Positions alone can't start MD: each second-order ODE needs two initial conditions.
- Coupled ODEs can't be solved **analytically**, but can be solved **numerically**, step by step — that's MD.

---

## 6. What makes a "good" integrator

> *(i) time-reversible; (ii) area-preserving (volume of an element in phase space does not change with time); (iii) no long-term drift in total energy.*

### (i) Time-reversible

Exact dynamics: run forward for time $\tau$, flip $\mathbf v \to -\mathbf v$, run for $\tau$ again → back at the start (with reversed velocities). Reason: $F = m\,d^2r/dt^2$ contains $t$ twice, so $t \to -t$ leaves it unchanged.

Algorithm: $n$ steps forward, flip velocities, $n$ steps with the same formula → land exactly at the start (to rounding). Non-reversible algorithms have a built-in time direction; their errors don't undo, so they tend to push one way → drift.

### (ii) Area-preserving (Liouville's theorem)

A small cloud of nearby states in phase space may move, stretch, shear — but its **volume never changes** under exact Hamiltonian dynamics. Shrinking cloud = hidden friction; growing cloud = hidden anti-friction. An integrator that distorts volume re-weights microstates, so time averages no longer equal the correct ensemble average. Explained in detail in [§7](#7-area-preservation-explained-from-scratch). Volume-preserving integrators (with the deeper structure) are called **symplectic**: Verlet and leap-frog are; Euler and RK4 are not.

### (iii) No long-term energy drift

- **Bounded fluctuation** around the right value: harmless.
- **Drift** (steady trend): fatal over millions of steps → wrong temperature, biased averages.

Check: run **NVE**, plot total energy. A slope means $\Delta t$ too large, potential not smooth at cutoff, etc. A thermostat (NVT) **hides** drift — test $\Delta t$ in NVE first.

### How they fit together

(i) + (ii) are the cause; (iii) is the result. A reversible, area-preserving algorithm can't spiral the phase-space point in or out, so energy can only wobble. (Precisely: a symplectic integrator exactly conserves a slightly modified "shadow Hamiltonian".)

**Accuracy per step is not on the list.** Many-atom systems are chaotic (trajectories starting $10^{-10}$ apart diverge exponentially), so no integrator follows the exact path for long. The goal is staying on the right energy surface and sampling it correctly. That's why RK4 (4th order, but not symplectic, and 4 force calls per step) loses to Verlet.

### ❓ Integrator A: ±0.5% fluctuation, flat. Integrator B: +0.01% per 1000 steps, always upward. For 10 million steps?

Choose **A**. B drifts $0.01\% \times 10^4 = 100\%$ — the energy roughly doubles. A's error stays bounded forever.

---

## 7. Area preservation explained from scratch

### ❓ My question: "Which area are we talking about?"

**Not** an area in real space. It's an area on the **phase-space map** (position across, velocity up).

**Step 1 — one point = one complete state.** 1D spring, $m = k = 1$. State $(x, v)$ is one dot. Exact motion moves it round a circle (constant energy $E = \tfrac12 v^2 + \tfrac12 x^2$), clockwise.

**Step 2 — a patch of dots.** Many copies of the spring with slightly different starts, e.g. $x \in [1.0, 1.3]$, $v \in [0, 0.3]$: a small square, area $0.3\times0.3 = 0.09$. This patch *is* the "area". Why a patch?
- It's a small piece of an **ensemble** (many copies in different microstates).
- Area **counts states**: more area = more microstates. How area changes tells you how the algorithm redistributes probability.

**Step 3 — exact dynamics.** For the spring every dot turns by the same angle, so the patch **rotates rigidly**; area stays 0.09 forever. For anharmonic systems (pendulum, LJ) outer dots lag inner ones, so the patch **shears** into a long thin streak — but area is still constant (stretched dough: longer and thinner in matching proportion). That's **Liouville's theorem**.

**Step 4 — Euler.** $x' = x + v\Delta t$, $v' = v - x\Delta t$ multiplies every patch's area by $1 + \Delta t^2$ (proved in [§14](#14-the-jacobian-determinant-from-scratch)). With $\Delta t = 0.2$: ×1.04 per step → ×1.37 after 8 steps, ×3.5 after one loop (32 steps), ×$10^{17}$ after 1000. And the patch moves outwards to higher energy.

**Step 5 — why it matters.** Dot density = probability. Exact dynamics: density preserved → ensemble correctly weighted. Euler: area grows, dots spread to higher energies → energy creeps up. An area-*shrinking* algorithm squeezes states towards the centre like friction → energy drains.

**Step 6 — 6N dimensions.** Same statement with "volume" instead of "area".

> **One line:** take any small bundle of nearby states on the position–velocity map and evolve each one; the bundle may move, turn or stretch, but the room it occupies must not grow or shrink.

### ❓ An algorithm multiplies area by 0.98 per step — what happens?

Area: ×0.13 after 100 steps, ×$10^{-9}$ after 1000. The patch shrinks and spirals in towards $(0, 0)$, the zero-energy state. Energy decays like a damped spring — hidden friction (the opposite of Euler).

---

## 8. Taylor series: the tool behind every integrator

### ❓ My question: how does this connect to $f(x)=\sum_n \frac{f^{(n)}(a)}{n!}(x-a)^n$?

It's the same formula with different letters:

| General formula | MD | Meaning |
|---|---|---|
| $f$ | $\mathbf r$ | function to predict |
| $a$ | $t$ | where everything is known (now) |
| $x$ | $t + \Delta t$ | where we want to predict (next step) |
| $x - a$ | $\Delta t$ | step size |
| $f'(a)$ | $\dot{\mathbf r}(t) = \mathbf v(t)$ | velocity |
| $f''(a)$ | $\ddot{\mathbf r}(t) = \mathbf a(t) = \mathbf f(t)/m$ | acceleration (Newton) |
| $f'''(a)$ | $\dddot{\mathbf r}(t)$ | rate of change of acceleration (unknown) |

*(Letter clash: in the formula $a$ is the expansion point and $f$ the function; in physics $\mathbf a$ is acceleration and $\mathbf f$ is force.)*

Term by term:

| $n$ | $f^{(n)}(a)(x-a)^n/n!$ | MD term |
|---|---|---|
| 0 | $f(a)$ | $\mathbf r(t)$ |
| 1 | $f'(a)(x-a)$ | $\mathbf v(t)\Delta t$ |
| 2 | $f''(a)(x-a)^2/2!$ | $\dfrac{\mathbf f(t)}{2m}\Delta t^2$ |
| 3 | $f'''(a)(x-a)^3/3!$ | $\tfrac16\dddot{\mathbf r}\Delta t^3$ ← dropped |

Sum of $n = 0, 1, 2$ = equation (1). The $2m$ is the $2!$ times the $m$ from Newton.

**Why stop at $n=2$:** each extra term is smaller by a factor ~$\Delta t$ ($0.1 \to 0.01 \to 0.001$). Up to $n=2$ everything is known (and $f''$ comes free from Newton). The $n=3$ term needs the time derivative of the force — not available, expensive.

| Method | Position series kept to | Velocity series kept to |
|---|---|---|
| Euler | $n=1$ | $n=1$ |
| Velocity-Verlet | $n=2$ | $n=2$ with $\dot{\mathbf a}$ approximated |
| Higher order (Gear, RK4) | $n=3,4,\dots$ | need $\dot{\mathbf a}, \ddot{\mathbf a},\dots$ |

### Forward and backward expansions

$$
g(t+h) = g(t) + g'(t)\,h + \frac{g''(t)}{2!}h^2 + \frac{g'''(t)}{3!}h^3 + \frac{g''''(t)}{4!}h^4 + \cdots
$$

$$
g(t-h) = g(t) - g'(t)\,h + \frac{g''(t)}{2!}h^2 - \frac{g'''(t)}{3!}h^3 + \frac{g''''(t)}{4!}h^4 - \cdots
$$

**Cancellation rule:** adding cancels odd powers of $h$; subtracting cancels even powers.

Notation: $\dot{\mathbf r} = \mathbf v$, $\ddot{\mathbf r} = \mathbf a = \mathbf f/m$, $\dddot{\mathbf r} = \mathbf b$, $\ddddot{\mathbf r} = \mathbf c$.

---

## 9. Velocity-Verlet

> *Positions and velocities of all atoms are known at every step.*

$$
\mathbf r(t+\Delta t) = \mathbf r(t) + \mathbf v(t)\,\Delta t + \frac{\mathbf f(t)}{2m}\,\Delta t^2 \qquad (1)
$$

$$
\mathbf v(t+\Delta t) = \mathbf v(t) + \frac{\mathbf f(t) + \mathbf f(t+\Delta t)}{2m}\,\Delta t \qquad (2)
$$

### Derivation of (1) — forward Taylor of $\mathbf r$

$$
\mathbf r(t+\Delta t) = \mathbf r(t) + \mathbf v(t)\Delta t + \tfrac12\mathbf a(t)\Delta t^2 + \tfrac16\mathbf b(t)\Delta t^3 + \cdots
$$

Truncate after $\Delta t^2$, use $\mathbf a = \mathbf f/m$ → (1). Error per step $O(\Delta t^3)$. (It's $s = ut + \tfrac12 at^2$ over one short step.)

### Derivation of (2) — forward Taylor of $\mathbf v$

$$
\mathbf v(t+\Delta t) = \mathbf v(t) + \mathbf a(t)\Delta t + \tfrac12\dot{\mathbf a}(t)\Delta t^2 + O(\Delta t^3)
$$

$\dot{\mathbf a}$ is the problem: $\dot{\mathbf a} = \frac1m\sum_j \frac{\partial \mathbf f}{\partial \mathbf r_j}\cdot\mathbf v_j$ needs the Hessian — far more expensive than the force. **Trick:** estimate it by finite difference of forces we compute anyway:

$$
\dot{\mathbf a}(t) \approx \frac{\mathbf a(t+\Delta t) - \mathbf a(t)}{\Delta t} \qquad (\text{error } O(\Delta t))
$$

Substitute:

$$
\mathbf v(t+\Delta t) = \mathbf v(t) + \mathbf a(t)\Delta t + \tfrac12\big[\mathbf a(t+\Delta t) - \mathbf a(t)\big]\Delta t
$$

Collect $\mathbf a(t)$ terms ($1 - \tfrac12 = \tfrac12$):

$$
\mathbf v(t+\Delta t) = \mathbf v(t) + \tfrac12\big[\mathbf a(t) + \mathbf a(t+\Delta t)\big]\Delta t \;\Rightarrow\; (2)
$$

Error check: $O(\Delta t)$ error in $\dot{\mathbf a}$ × $\tfrac12\Delta t^2$ = $O(\Delta t^3)$ — same as the terms already dropped. The trick costs no accuracy.

**Second view (symmetric):** expand forward from $t$ and backward from $t+\Delta t$, average:

$$
\mathbf v(t+\Delta t) = \mathbf v(t) + \tfrac12\big[\mathbf a(t) + \mathbf a(t+\Delta t)\big]\Delta t + \tfrac14\big[\dot{\mathbf a}(t) - \dot{\mathbf a}(t+\Delta t)\big]\Delta t^2 + O(\Delta t^3)
$$

The last bracket is itself $O(\Delta t)$, so that term is $O(\Delta t^3)$. Equation (2) is the average of a prediction from each end of the step — neither end favoured — which is the origin of time-reversibility.

### The red note: "forces at current and new positions"

Order within one step:
1. New position from (1) using $\mathbf r(t), \mathbf v(t), \mathbf f(t)$.
2. New force $\mathbf f(t+\Delta t)$ from new positions (expensive).
3. New velocity from (2).

Only **one force call per step**: $\mathbf f(t+\Delta t)$ is reused as $\mathbf f(t)$ next step. At the start, $\mathbf r(0)$ (step 2), $\mathbf v(0)$ (step 3), and $\mathbf f(0)$ computed once before the loop.

### Kick–drift–kick (same algorithm, as LAMMPS writes it)

$$
\begin{aligned}
\mathbf v_{1/2} &= \mathbf v(t) + \frac{\mathbf f(t)}{2m}\Delta t &&\text{half kick}\\
\mathbf r(t+\Delta t) &= \mathbf r(t) + \mathbf v_{1/2}\Delta t &&\text{drift}\\
&\text{compute } \mathbf f(t+\Delta t)\\
\mathbf v(t+\Delta t) &= \mathbf v_{1/2} + \frac{\mathbf f(t+\Delta t)}{2m}\Delta t &&\text{half kick}
\end{aligned}
$$

Substitute line 1 into lines 2 and 4 → (1) and (2).

### ❓ "I want velocity-Verlet — why do you need half-step updates?"

You don't. Equations (1) and (2) give **identical numbers**. Kick–drift–kick is only a re-arrangement that makes it obvious *why* $\det J = 1$ (each piece is a shear). The direct route proves it by algebra ([§16](#16-worked-example-lennard-jones-jacobian-with-real-numbers)).

### Worked hand calculation (spring, $m=k=1$, $f=-x$, $x_0=1$, $v_0=0$, $\Delta t=0.1$)

| | Step 1 | Step 2 |
|---|---|---|
| $x$ | $1 + 0 + (-1/2)(0.01) = 0.995$ | $0.995 - 0.009975 - 0.004975 = 0.98005$ |
| $f$ | $-0.995$ | $-0.98005$ |
| $v$ | $0 + \frac{-1-0.995}{2}(0.1) = -0.09975$ | $-0.09975 + \frac{-0.995-0.98005}{2}(0.1) = -0.19850$ |
| exact $x=\cos t$ | 0.995004 | 0.980067 |
| exact $v=-\sin t$ | −0.099833 | −0.198669 |

Energy after step 1: $\tfrac12(0.09975)^2 + \tfrac12(0.995)^2 = 0.499988$ (started 0.5) — wobbles, doesn't drift.

**Accuracy:** local $O(\Delta t^3)$, global $O(\Delta t^2)$ — **second order** (halve $\Delta t$ → error ÷4).

---

## 10. Verlet

> *Requires positions at current and previous configurations; trajectory equivalent to velocity-Verlet.*

$$
\mathbf r(t+\Delta t) = 2\mathbf r(t) - \mathbf r(t-\Delta t) + \frac{\mathbf f(t)}{m}\Delta t^2 \qquad (5)
$$

### Derivation in the notes (from velocity-Verlet)

Write (1) one step earlier:

$$
\mathbf r(t) = \mathbf r(t-\Delta t) + \mathbf v(t-\Delta t)\Delta t + \frac{\mathbf f(t-\Delta t)}{2m}\Delta t^2
$$

$$
\Rightarrow\quad \mathbf r(t-\Delta t) = \mathbf r(t) - \mathbf v(t-\Delta t)\Delta t - \frac{\mathbf f(t-\Delta t)}{2m}\Delta t^2 \qquad (3)
$$

Add (1) and (3):

$$
\mathbf r(t+\Delta t) + \mathbf r(t-\Delta t) = 2\mathbf r(t) + \big[\mathbf v(t) - \mathbf v(t-\Delta t)\big]\Delta t + \frac{\mathbf f(t) - \mathbf f(t-\Delta t)}{2m}\Delta t^2 \qquad (4)
$$

From (2) one step earlier:

$$
\mathbf v(t) - \mathbf v(t-\Delta t) = \frac{\mathbf f(t) + \mathbf f(t-\Delta t)}{2m}\Delta t
$$

Substitute into (4); force terms combine:

$$
\frac{\mathbf f(t) + \mathbf f(t-\Delta t)}{2m} + \frac{\mathbf f(t) - \mathbf f(t-\Delta t)}{2m} = \frac{\mathbf f(t)}{m}
$$

→ equation (5). No approximation was made, so Verlet and velocity-Verlet give **identical positions**.

### Derivation directly from Taylor series

$$
\mathbf r(t\pm\Delta t) = \mathbf r \pm \mathbf v\Delta t + \tfrac12\mathbf a\Delta t^2 \pm \tfrac16\mathbf b\Delta t^3 + \tfrac1{24}\mathbf c\Delta t^4 \pm \cdots
$$

Add: $\mathbf v$ and $\mathbf b$ cancel:

$$
\mathbf r(t+\Delta t) + \mathbf r(t-\Delta t) = 2\mathbf r(t) + \mathbf a(t)\Delta t^2 + \tfrac1{12}\mathbf c(t)\Delta t^4 + \cdots
$$

→ (5), local error $O(\Delta t^4)$ — the unknown $\mathbf b$ cancelled for free.

### Velocity

Subtract instead (even terms cancel):

$$
\mathbf r(t+\Delta t) - \mathbf r(t-\Delta t) = 2\mathbf v(t)\Delta t + \tfrac13\mathbf b\Delta t^3 + \cdots
\;\Rightarrow\;
\mathbf v(t) = \frac{\mathbf r(t+\Delta t) - \mathbf r(t-\Delta t)}{2\Delta t} + O(\Delta t^2)
$$

Available only once $\mathbf r(t+\Delta t)$ is known → velocities lag one step.

### Notes

- Rearranged: $\frac{\mathbf r(t+\Delta t) - 2\mathbf r(t) + \mathbf r(t-\Delta t)}{\Delta t^2} = \frac{\mathbf f}{m}$ — the central finite-difference second derivative; Newton's law discretised directly. Error $O(\Delta t^2)$.
- Drawbacks: need $\mathbf r(-\Delta t)$ to start (from (3)); adding a small $\Delta t^2$ term to the large $2\mathbf r - \mathbf r_{\text{prev}}$ loses floating-point precision.
- Symmetric in $t\pm\Delta t$ → time-reversible by construction.

---

## 11. Leap-frog

> *Uses positions and velocities, but they differ by half a time step.*

Positions on whole steps, velocities on half steps; they take turns jumping over each other.

$$
\mathbf v\!\left(t+\tfrac{\Delta t}{2}\right) = \mathbf v\!\left(t-\tfrac{\Delta t}{2}\right) + \frac{\mathbf f(t)}{m}\Delta t \qquad (6)
$$

$$
\mathbf r(t+\Delta t) = \mathbf r(t) + \mathbf v\!\left(t+\tfrac{\Delta t}{2}\right)\Delta t \qquad (7)
$$

### Derivation in the notes (from Verlet)

Define $\mathbf v(t-\tfrac{\Delta t}{2}) = \dfrac{\mathbf r(t) - \mathbf r(t-\Delta t)}{\Delta t}$, so $\mathbf r(t-\Delta t) = \mathbf r(t) - \mathbf v(t-\tfrac{\Delta t}{2})\Delta t$. Substitute into (5):

$$
\mathbf r(t+\Delta t) = 2\mathbf r(t) - \big[\mathbf r(t) - \mathbf v(t-\tfrac{\Delta t}{2})\Delta t\big] + \frac{\mathbf f}{m}\Delta t^2
= \mathbf r(t) + \mathbf v(t-\tfrac{\Delta t}{2})\Delta t + \frac{\mathbf f}{m}\Delta t^2
$$

With the same definition one step later, $\mathbf v(t+\tfrac{\Delta t}{2}) = [\mathbf r(t+\Delta t) - \mathbf r(t)]/\Delta t$:
- Subtract $\mathbf r(t)$ from the line above and divide by $\Delta t$ → (6).
- Rearrange the definition → (7).

**Reverse direction (leap-frog → Verlet):** write (7) for two consecutive steps, subtract, use (6):

$$
\mathbf r(t+\Delta t) - 2\mathbf r(t) + \mathbf r(t-\Delta t) = \big[\mathbf v(t+\tfrac{\Delta t}{2}) - \mathbf v(t-\tfrac{\Delta t}{2})\big]\Delta t = \frac{\mathbf f}{m}\Delta t^2
$$

= Verlet. So the three are **one algorithm in three forms**.

### Derivation directly from Taylor series

**Position** — expand $\mathbf r$ about $t_{1/2} = t + \tfrac{\Delta t}{2}$ with $h = \pm\tfrac{\Delta t}{2}$:

$$
\mathbf r(t+\Delta t) = \mathbf r(t_{1/2}) + \mathbf v(t_{1/2})\tfrac{\Delta t}{2} + \tfrac12\mathbf a(t_{1/2})\left(\tfrac{\Delta t}{2}\right)^2 + \tfrac16\mathbf b(t_{1/2})\left(\tfrac{\Delta t}{2}\right)^3 + \cdots
$$

$$
\mathbf r(t) = \mathbf r(t_{1/2}) - \mathbf v(t_{1/2})\tfrac{\Delta t}{2} + \tfrac12\mathbf a(t_{1/2})\left(\tfrac{\Delta t}{2}\right)^2 - \tfrac16\mathbf b(t_{1/2})\left(\tfrac{\Delta t}{2}\right)^3 + \cdots
$$

Subtract ($\mathbf r(t_{1/2})$ and $\mathbf a$ cancel):

$$
\mathbf r(t+\Delta t) - \mathbf r(t) = \mathbf v(t_{1/2})\Delta t + \tfrac1{24}\mathbf b(t_{1/2})\Delta t^3 + \cdots \;\Rightarrow\; (7)
$$

The velocity is sampled at exactly the middle of the step — that's why a "distance = speed × time" formula is accurate.

**Velocity** — expand $\mathbf v$ about $t$ with $h = \pm\tfrac{\Delta t}{2}$, subtract ($\mathbf v(t)$ and $\dot{\mathbf a}$ cancel):

$$
\mathbf v(t+\tfrac{\Delta t}{2}) - \mathbf v(t-\tfrac{\Delta t}{2}) = \mathbf a(t)\Delta t + \tfrac1{24}\ddot{\mathbf a}(t)\Delta t^3 + \cdots \;\Rightarrow\; (6)
$$

The troublesome $\dot{\mathbf a}$ of velocity-Verlet simply disappears.

**Full-step velocity** (for kinetic energy) — add instead ($\mathbf a$ cancels):

$$
\mathbf v(t) = \tfrac12\big[\mathbf v(t-\tfrac{\Delta t}{2}) + \mathbf v(t+\tfrac{\Delta t}{2})\big] + O(\Delta t^2)
$$

**Cost:** $\mathbf r$ and $\mathbf v$ are never at the same instant, so KE at time $t$ needs this average. (GROMACS uses leap-frog; LAMMPS uses velocity-Verlet — their reported KE differs slightly.)

---

## 12. All three integrators from Taylor series, side by side

| Formula | Expanded about | Step $h$ | Operation | Cancelled | Local error |
|---|---|---|---|---|---|
| (1) VV position | $t$ | $+\Delta t$ | truncate | — | $O(\Delta t^3)$ |
| (2) VV velocity | $t$ | $+\Delta t$ | approximate $\dot{\mathbf a}$ | — | $O(\Delta t^3)$ |
| (5) Verlet | $t$ | $\pm\Delta t$ | add | $\mathbf v, \mathbf b$ | $O(\Delta t^4)$ |
| (7) LF position | $t+\tfrac{\Delta t}{2}$ | $\pm\tfrac{\Delta t}{2}$ | subtract | $\mathbf r(t_{1/2}), \mathbf a$ | $O(\Delta t^3)$ |
| (6) LF velocity | $t$ | $\pm\tfrac{\Delta t}{2}$ | subtract | $\mathbf v(t), \dot{\mathbf a}$ | $O(\Delta t^3)$ |

**Global error** over a fixed simulated time: $O(\Delta t^2)$ for all three — all second order. (Verlet's $O(\Delta t^4)$ local error is not a contradiction: its two-step recurrence accumulates errors faster; per-step error counts differently in each scheme, and the totals match — consistent with them being exact rearrangements.)

**Key idea:** forward-only expansion leaves an unknown derivative to approximate; symmetric expansion cancels it for free — and symmetry is also what makes all three time-reversible.

---

## 13. Proof: velocity-Verlet is time-reversible

**Goal:** step $(\mathbf r_0, \mathbf v_0) \to (\mathbf r_1, \mathbf v_1)$; then from $(\mathbf r_1, -\mathbf v_1)$ one step of the *same* formulas must land on $(\mathbf r_0, -\mathbf v_0)$. Write $\mathbf f_0 = \mathbf f(\mathbf r_0)$, $\mathbf f_1 = \mathbf f(\mathbf r_1)$.

Forward:

$$
\mathbf r_1 = \mathbf r_0 + \mathbf v_0\Delta t + \frac{\mathbf f_0}{2m}\Delta t^2 \quad (i)
\qquad
\mathbf v_1 = \mathbf v_0 + \frac{\mathbf f_0 + \mathbf f_1}{2m}\Delta t \quad (ii)
$$

**Identity:** multiply (ii) by $\Delta t$, substitute $\mathbf v_0\Delta t$ from (i):

$$
\mathbf v_1\Delta t = \mathbf r_1 - \mathbf r_0 + \frac{\mathbf f_1}{2m}\Delta t^2 \quad (iii)
$$

**Backward position** (velocity $-\mathbf v_1$, force $\mathbf f_1$):

$$
\mathbf r' = \mathbf r_1 - \mathbf v_1\Delta t + \frac{\mathbf f_1}{2m}\Delta t^2
= \mathbf r_1 - \Big[\mathbf r_1 - \mathbf r_0 + \frac{\mathbf f_1}{2m}\Delta t^2\Big] + \frac{\mathbf f_1}{2m}\Delta t^2 = \mathbf r_0 \;\checkmark
$$

**Backward velocity** (force at $\mathbf r'$ is $\mathbf f_0$):

$$
\mathbf v' = -\mathbf v_1 + \frac{\mathbf f_1 + \mathbf f_0}{2m}\Delta t
= -\mathbf v_0 - \frac{\mathbf f_0 + \mathbf f_1}{2m}\Delta t + \frac{\mathbf f_0 + \mathbf f_1}{2m}\Delta t = -\mathbf v_0 \;\checkmark
$$

**Why it works:** (2) uses $\mathbf f_0 + \mathbf f_1$ symmetrically; going backwards the forces swap roles but the sum is unchanged. **Euler** ($\mathbf v_1 = \mathbf v_0 + \mathbf f_0\Delta t/m$) uses only the start force; going back it would use $\mathbf f_1 \ne \mathbf f_0$, so kicks don't cancel.

---

## 14. The Jacobian determinant from scratch

### ❓ My question: what does "one step is a map $(x,v)\mapsto(x',v')$" mean?

**Step 1 — a step is a machine.** Put in two numbers $(x, v)$, get out two numbers $(x', v')$. Euler for the spring with $\Delta t = 0.5$: $x' = x + 0.5v$, $v' = v - 0.5x$.

| Input $(x,v)$ | $x'$ | $v'$ | Output |
|---|---|---|---|
| (0, 0) | 0 | 0 | (0, 0) |
| (1, 0) | 1 | −0.5 | (1, −0.5) |
| (1, 1) | 1.5 | 0.5 | (1.5, 0.5) |
| (0, 1) | 0.5 | 1 | (0.5, 1) |

A rule that moves every point of the plane to a new point is a **map**. That's all the word means.

**Step 2 — apply it to a square.** The four inputs are corners of a unit square. The outputs form a tilted parallelogram. A velocity-Verlet **kick** ($x' = x$, $v' = v - \tfrac{\Delta t}{2}x$) on the same square gives corners (0,0), (1,−0.25), (1,0.75), (0,1): vertical edges stay vertical and keep their length — a **shear**.

**Step 3 — measure area.** A parallelogram with edge vectors $(a, c)$ and $(b, d)$ has area $|ad - bc|$.
*Why:* enclosing rectangle $(a+b)(c+d)$ minus two triangles $\tfrac12 ac$, two $\tfrac12 bd$, two rectangles $bc$ = $ad - bc$.

- Euler: edges $(1, -0.5)$, $(0.5, 1)$ → $1\cdot1 - 0.5\cdot(-0.5) = 1.25$
- Kick: edges $(1, -0.25)$, $(0, 1)$ → $1\cdot1 - 0\cdot(-0.25) = 1$

**Step 4 — the edge vectors are partial derivatives.** Nudge $x$ by 1 (keep $v$): outputs change by $(\partial x'/\partial x,\ \partial v'/\partial x)$. Nudge $v$ by 1: $(\partial x'/\partial v,\ \partial v'/\partial v)$. Put these as columns:

$$
J = \begin{pmatrix}
\dfrac{\partial x'}{\partial x} & \dfrac{\partial x'}{\partial v}\\[3mm]
\dfrac{\partial v'}{\partial x} & \dfrac{\partial v'}{\partial v}
\end{pmatrix}
= \begin{pmatrix} a & b\\ c & d\end{pmatrix},
\qquad
\det J = ad - bc
$$

```
                 nudge x           nudge v
              ┌                              ┐
   x' row     │  ∂x'/∂x            ∂x'/∂v    │
   v' row     │  ∂v'/∂x            ∂v'/∂v    │
              └                              ┘
```

$$
\boxed{\det J = \frac{\text{area after one step}}{\text{area before}}}
\qquad
\text{area-preserving} \iff \det J = 1
$$

| Map | $J$ | $\det J$ |
|---|---|---|
| Euler | $\begin{pmatrix}1 & \Delta t\\ -\Delta t & 1\end{pmatrix}$ | $1 + \Delta t^2$ |
| Kick | $\begin{pmatrix}1 & 0\\ -\tfrac{\Delta t}{2} & 1\end{pmatrix}$ | $1$ |
| Drift | $\begin{pmatrix}1 & \Delta t\\ 0 & 1\end{pmatrix}$ | $1$ |

**Step 5 — real (curved) forces.** For LJ the map is non-linear; a big square comes out bent. But zoom in far enough and any smooth map looks linear, so a tiny square $dx\times dv$ still becomes a parallelogram:

$$
dx'\,dv' = \lvert\det J\rvert\,dx\,dv
$$

$J$ may vary from point to point; area-preserving means $\det J = 1$ **everywhere**.

---

## 15. Proof: velocity-Verlet is area-preserving

### Route A — kick, drift, kick (shears)

**Kick** ($x' = x$, $v' = v + \tfrac{\Delta t}{2m}f(x)$):

$$
J_{\text{kick}} = \begin{pmatrix} 1 & 0\\[1mm] \dfrac{\Delta t}{2m}f'(x) & 1\end{pmatrix},
\qquad \det = 1\cdot 1 - 0 = 1
$$

The $f'(x)$ entry is multiplied by 0 — the physics never enters the determinant.

**Drift** ($x' = x + v\Delta t$, $v' = v$): $\det J_{\text{drift}} = 1$.

**Full step** (determinant of a product = product of determinants):

$$
\det J_{\text{VV}} = 1\times1\times1 = 1
$$

for **any** force and **any** $\Delta t$.

*Why shears preserve area:* each sub-step changes one coordinate by an amount depending only on the other. Every vertical (or horizontal) line slides rigidly — like pushing a deck of cards sideways.

**Euler fails** because it changes $x$ and $v$ *simultaneously*, each using the other's old value: both off-diagonal entries are non-zero, and $-bc$ adds $+\Delta t^2$.

### Route B — directly from (1) and (2)

See the worked example in [§16](#16-worked-example-lennard-jones-jacobian-with-real-numbers); the general algebra gives $\det J = 1$ exactly.

### Direct check on the spring ($m=k=1$)

$$
x_1 = \left(1 - \tfrac{\Delta t^2}{2}\right)x_0 + \Delta t\,v_0,
\qquad
v_1 = -\Delta t\left(1 - \tfrac{\Delta t^2}{4}\right)x_0 + \left(1 - \tfrac{\Delta t^2}{2}\right)v_0
$$

$$
\det J = \left(1 - \tfrac{\Delta t^2}{2}\right)^2 + \Delta t^2\left(1 - \tfrac{\Delta t^2}{4}\right)
= 1 - \Delta t^2 + \tfrac{\Delta t^4}{4} + \Delta t^2 - \tfrac{\Delta t^4}{4} = 1
$$

---

## 16. Worked example: Lennard-Jones Jacobian with real numbers

**Setup:** one atom moving in 1D in the LJ field of an atom fixed at the origin; state $(x, v)$. Reduced units $\varepsilon = \sigma = m = 1$.

$$
u(x) = 4(x^{-12} - x^{-6}),\qquad
f(x) = -u'(x) = 48x^{-13} - 24x^{-7},\qquad
f'(x) = -624x^{-14} + 168x^{-8}
$$

$x_0 = 1.2$, $v_0 = -0.5$, $\Delta t = 0.01$. Inside the attractive well (minimum at $2^{1/6} \approx 1.122$), moving inward. $E_0 = -0.8910 + 0.125 = -0.7660$ (bound).

### One velocity-Verlet step with (1) and (2)

$$
f_0 = -2.21169,\qquad f'_0 = -9.52979
$$

$$
x_1 = 1.2 + (-0.5)(0.01) + \frac{-2.21169}{2}(0.0001) = 1.2 - 0.005 - 0.000110585 = 1.194889
$$

$$
f_1 = -2.15887,\qquad f'_1 = -11.16520
$$

$$
v_1 = -0.5 + \frac{-2.21169 - 2.15887}{2}(0.01) = -0.521853
$$

### The four partial derivatives

**Subtlety:** $v_1$ depends on $x_0$ (and $v_0$) **directly** through $f(x_0)$ and **indirectly** through $f(x_1)$, because $x_1$ depends on them. Chain rule:

$$
\frac{\partial f(x_1)}{\partial x_0} = f'(x_1)\frac{\partial x_1}{\partial x_0},
\qquad
\frac{\partial f(x_1)}{\partial v_0} = f'(x_1)\frac{\partial x_1}{\partial v_0}
$$

```
∂x1/∂x0 = 1 + (f'0/2)·Δt²                 = 1 + (−9.52979/2)(0.0001)       =  0.999524
∂x1/∂v0 = Δt                                                               =  0.01
∂v1/∂x0 = (Δt/2)·[ f'0 + f'1·(∂x1/∂x0) ]  = 0.005·[−9.52979 − 11.15989]    = −0.103448
∂v1/∂v0 = 1 + (Δt/2)·f'1·Δt               = 1 + 0.005·(−11.16520)·(0.01)   =  0.999442
```

### ❓ Derivation of $\partial v_1/\partial v_0 = 1 + \tfrac{\Delta t}{2}f'_1\Delta t$

Start from (2): $v_1 = v_0 + \tfrac{\Delta t}{2}[f(x_0) + f(x_1)]$.

1. Differentiate term by term: $\dfrac{\partial v_1}{\partial v_0} = \dfrac{\partial v_0}{\partial v_0} + \dfrac{\Delta t}{2}\Big[\dfrac{\partial f(x_0)}{\partial v_0} + \dfrac{\partial f(x_1)}{\partial v_0}\Big]$
2. $\partial v_0/\partial v_0 = 1$
3. $f(x_0)$ depends only on $x_0$ → $\partial f(x_0)/\partial v_0 = 0$
4. Chain $v_0 \to x_1 \to f(x_1)$: $\partial f(x_1)/\partial v_0 = f'(x_1)\,\partial x_1/\partial v_0$
5. From (1), $x_1 = x_0 + v_0\Delta t + \tfrac12 f(x_0)\Delta t^2$ → derivatives of the three terms w.r.t. $v_0$ are $0, \Delta t, 0$ → $\partial x_1/\partial v_0 = \Delta t$
6. So $\partial f(x_1)/\partial v_0 = f'_1\Delta t$
7. Assemble: $\dfrac{\partial v_1}{\partial v_0} = 1 + \dfrac{\Delta t}{2}f'_1\Delta t = 1 + \tfrac12 f'_1\Delta t^2$
8. Numbers: $1 + 0.005\times(-11.16520)\times0.01 = 0.999442$

*Physical meaning:* start faster by $\delta$ → travel $\Delta t\,\delta$ further → force differs by $f'_1\Delta t\,\delta$ → (2) adds $\tfrac{\Delta t}{2}$ of that. Here $f'_1 < 0$, so the end velocity changes by slightly less than $\delta$. Euler's velocity update ($v_1 = v_0 + \Delta t f(x_0)$) never sees where the atom ends up — that missing link is what breaks its area preservation.

### The matrix and its determinant

$$
J = \begin{pmatrix} 0.999524 & 0.01\\ -0.103448 & 0.999442 \end{pmatrix}
$$

```
              ┌                         ┐
              │  0.999524      0.01     │
       J  =   │                         │
              │ −0.103448      0.999442 │
              └                         ┘
```

$$
\det J = (0.999524)(0.999442) - (0.01)(-0.103448) = 0.998966 + 0.001034 = 1.000000
$$

| Entry | Value | Meaning |
|---|---|---|
| $\partial x_1/\partial x_0$ | 0.999524 | start $\delta$ further out → end ≈ $0.9995\delta$ further out |
| $\partial x_1/\partial v_0$ | 0.01 | extra speed → move $\Delta t$ further per unit |
| $\partial v_1/\partial x_0$ | −0.103448 | start further out → stronger pull → faster inward |
| $\partial v_1/\partial v_0$ | 0.999442 | extra speed mostly carries through |

### Why exactly 1 (general algebra)

Let $A = \tfrac12 f'_0\Delta t^2$, $B = \tfrac12 f'_1\Delta t^2$. Then

$$
J = \begin{pmatrix} 1 + A & \Delta t\\[1mm] \tfrac{\Delta t}{2}\big[f'_0 + f'_1(1+A)\big] & 1 + B\end{pmatrix}
$$

$$
\det J = (1+A)(1+B) - \Delta t\cdot\tfrac{\Delta t}{2}\big[f'_0 + f'_1(1+A)\big]
= (1 + A + B + AB) - (A + B + AB) = 1
$$

for **any** $f'_0, f'_1, \Delta t$ — any potential.

### Same route via kick–drift–kick (gives the same matrix)

```
J1 (kick)  = [[1, 0], [(Δt/2)f'0, 1]] = [[1, 0], [−0.047649, 1]]
J2 (drift) = [[1, 0.01], [0, 1]]
J3 (kick)  = [[1, 0], [(Δt/2)f'1, 1]] = [[1, 0], [−0.055826, 1]]

J2·J1   = [[0.999524, 0.01], [−0.047649, 1]]
J3·J2·J1 = [[0.999524, 0.01], [−0.055826×0.999524 − 0.047649, −0.055826×0.01 + 1]]
         = [[0.999524, 0.01], [−0.103448, 0.999442]]
```

### Euler at the same point

$$
J_{\text{Euler}} = \begin{pmatrix} 1 & 0.01\\ \Delta t f'_0 & 1\end{pmatrix} = \begin{pmatrix} 1 & 0.01\\ -0.095298 & 1\end{pmatrix},
\qquad \det = 1 - \Delta t^2 f'_0 = 1.000953
$$

### Numerical check with a real tiny square ($h = 0.001$)

| Corner | Start | After one VV step |
|---|---|---|
| A | (1.200, −0.500) | (1.1948894, −0.5218528) |
| B | (1.201, −0.500) | (1.1958889, −0.5219547) |
| D | (1.200, −0.499) | (1.1948994, −0.5208534) |

Edges $\mathbf a = (0.00099953, -0.00010186)$, $\mathbf b = (0.00001, 0.00099944)$ (= $h\times$ columns of $J$).
New area $= 9.99984\times10^{-7}$ → ratio **0.99999** (→ 1 as $h \to 0$; the residual is LJ curvature over a finite square). Euler: **1.00094**.

### Many steps: a whole patch (figure)

![Phase-space patch: velocity-Verlet vs Euler, LJ](figures/lj_patch.png)

A disc of 3000 states around $(1.2, -0.5)$, evolved 300 steps:

- **Velocity-Verlet:** stretched into a thin crescent (anharmonic: higher-energy states take longer round the well) but area **×1.00** at every snapshot. After 1000 steps mean energy still −0.765.
- **Euler:** area ×1.56 (step 40), ×1.62 (100), ×2.74 (200), ×3.82 (300); patch moves out to higher energy. After 1000 steps mean energy −0.765 → **+0.27** — positive, so the atom has **escaped the well**: the bond broke purely because of the integrator.

> Stretching/bending the patch is allowed (real dynamics does it). Changing *how much room the states occupy* is not.

---

## 17. The Jacobian for a real 3D system (6N × 6N)

### ❓ Why is $J$ a $6N\times6N$ matrix?

$J_{\alpha\beta} = \partial(\text{output }\alpha)/\partial(\text{input }\beta)$. **Rows = outputs, columns = inputs.**

- Reading along a **row**: how one output depends on every input.
- Reading down a **column**: how every output responds when one input is nudged (the edge vector).

| System | Inputs | Outputs | $J$ size | Entries |
|---|---|---|---|---|
| 1 atom, 1D | $x, v$ | 2 | $2\times2$ | 4 |
| 1 atom, 3D | $x,y,z,v_x,v_y,v_z$ | 6 | $6\times6$ | 36 |
| 2 atoms, 3D | 12 | 12 | $12\times12$ | 144 |
| $N$ atoms, 3D | $3N + 3N = 6N$ | $6N$ | $6N\times6N$ | $36N^2$ |
| 1000 atoms | | | $6000\times6000$ | 36 million |

Always square, because an MD step outputs a new value for every position and every velocity — and only a square matrix has a determinant.

**Example — drift for one atom in 3D:**

```
         x   y   z   vx  vy  vz
   x' [  1   0   0   Δt  0   0  ]
   y' [  0   1   0   0   Δt  0  ]
   z' [  0   0   1   0   0   Δt ]
  vx' [  0   0   0   1   0   0  ]
  vy' [  0   0   0   0   1   0  ]
  vz' [  0   0   0   0   0   1  ]
```

### Block form

$$
J = \begin{pmatrix}
\partial\mathbf R'/\partial\mathbf R & \partial\mathbf R'/\partial\mathbf V\\
\partial\mathbf V'/\partial\mathbf R & \partial\mathbf V'/\partial\mathbf V
\end{pmatrix}
$$

The only non-trivial piece is the **force Jacobian** $\mathbf K = \partial\mathbf F/\partial\mathbf R = -\partial^2 U/\partial\mathbf R\,\partial\mathbf R$ (minus the Hessian).

### Two atoms, built from scratch

**Step 0 — setup.** Two atoms, mass $m$, pair potential $u(r)$, $r = \sqrt{(x_1-x_2)^2 + (y_1-y_2)^2 + (z_1-z_2)^2}$. Place them on the $x$-axis, atom 1 to the right: $y_1 = y_2$, $z_1 = z_2$, $x_1 - x_2 = r$. (Other orientations fill more entries; the block pattern is the same.)

**Step 1 — numbering.**

| # | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| | $x_1$ | $y_1$ | $z_1$ | $x_2$ | $y_2$ | $z_2$ | $v_{x1}$ | $v_{y1}$ | $v_{z1}$ | $v_{x2}$ | $v_{y2}$ | $v_{z2}$ |

**Step 2 — the kick as 12 equations** ($\alpha = \Delta t/2m$): positions unchanged (rows 1–6); $v_{x1}' = v_{x1} + \alpha f_{x1}$, etc. (rows 7–12). Forces depend only on positions.

**Step 3 — rows 1–6:** $x_1' = x_1$ → a single 1 on the diagonal. Top-left block $= \mathbf I$, top-right $= \mathbf 0$.

**Step 4 — rows 7–12, velocity columns:** force doesn't depend on velocity → bottom-right $= \mathbf I$.

**Step 5 — rows 7–12, position columns (the physics).**

Forces:

$$
f_{x1} = -u'(r)\frac{x_1 - x_2}{r},\quad
f_{y1} = -u'(r)\frac{y_1 - y_2}{r},\quad
f_{z1} = -u'(r)\frac{z_1 - z_2}{r},\qquad
\mathbf f_2 = -\mathbf f_1
$$

*Type 1 — along the bond:* product rule, then evaluate on the axis ($x_1 - x_2 = r$):

$$
\frac{\partial f_{x1}}{\partial x_1} = -\Big[u''\Big(\frac{x_1-x_2}{r}\Big)^2 + u'\Big(\frac1r - \frac{(x_1-x_2)^2}{r^3}\Big)\Big] = -\big[u''\cdot1 + u'\cdot 0\big] = -u''(r)
$$

Meaning: stretch the bond → restoring push set by stiffness $u''$.

*Type 2 — sideways:* with $y_1 - y_2 = 0$:

$$
\frac{\partial f_{y1}}{\partial y_1} = -\Big[u''\Big(\frac{y_1-y_2}{r}\Big)^2 + u'\Big(\frac1r - \frac{(y_1-y_2)^2}{r^3}\Big)\Big] = -\frac{u'(r)}{r}
$$

Same for $z$. Meaning: push atom 1 sideways → bond tilts → tension $u'$ gets a restoring component (like a plucked string).

*Type 3 — cross terms* ($\partial f_{x1}/\partial y_1$ etc.): every term contains $(y_1-y_2) = 0$ → **0**.

*Type 4 — moving atom 2:* force depends on the difference, so $\partial/\partial x_2 = -\partial/\partial x_1$: $\partial f_{x1}/\partial x_2 = +u''$, $\partial f_{y1}/\partial y_2 = +u'/r$.

*Atom 2's forces:* $\mathbf f_2 = -\mathbf f_1$ flips every sign.

General pair stiffness block:

$$
\mathbf k = u''(r)\,\hat{\mathbf r}\hat{\mathbf r}^{\mathsf T} + \frac{u'(r)}{r}\big(\mathbf I - \hat{\mathbf r}\hat{\mathbf r}^{\mathsf T}\big),
\qquad
\mathbf K = \begin{pmatrix} -\mathbf k & \mathbf k\\ \mathbf k & -\mathbf k\end{pmatrix}
$$

On the $x$-axis $\mathbf k = \mathrm{diag}(u'', u'/r, u'/r)$.

Multiply by $\alpha$: $p = \alpha u''(r)$, $q = \alpha u'(r)/r$.

| Row | Output | $x_1$ | $y_1$ | $z_1$ | $x_2$ | $y_2$ | $z_2$ |
|---|---|---|---|---|---|---|---|
| 7 | $v_{x1}'$ | $-p$ | 0 | 0 | $+p$ | 0 | 0 |
| 8 | $v_{y1}'$ | 0 | $-q$ | 0 | 0 | $+q$ | 0 |
| 9 | $v_{z1}'$ | 0 | 0 | $-q$ | 0 | 0 | $+q$ |
| 10 | $v_{x2}'$ | $+p$ | 0 | 0 | $-p$ | 0 | 0 |
| 11 | $v_{y2}'$ | 0 | $+q$ | 0 | 0 | $-q$ | 0 |
| 12 | $v_{z2}'$ | 0 | 0 | $+q$ | 0 | 0 | $-q$ |

Each row sums to zero: moving both atoms equally changes no force (translational invariance).

**Full kick matrix:**

```
          x1  y1  z1  x2  y2  z2 | vx1 vy1 vz1 vx2 vy2 vz2
 x1'   [  1   0   0   0   0   0  |  0   0   0   0   0   0 ]
 y1'   [  0   1   0   0   0   0  |  0   0   0   0   0   0 ]
 z1'   [  0   0   1   0   0   0  |  0   0   0   0   0   0 ]
 x2'   [  0   0   0   1   0   0  |  0   0   0   0   0   0 ]
 y2'   [  0   0   0   0   1   0  |  0   0   0   0   0   0 ]
 z2'   [  0   0   0   0   0   1  |  0   0   0   0   0   0 ]
       [ ------------------------+------------------------ ]
 vx1'  [ −p   0   0   p   0   0  |  1   0   0   0   0   0 ]
 vy1'  [  0  −q   0   0   q   0  |  0   1   0   0   0   0 ]
 vz1'  [  0   0  −q   0   0   q  |  0   0   1   0   0   0 ]
 vx2'  [  p   0   0  −p   0   0  |  0   0   0   1   0   0 ]
 vy2'  [  0   q   0   0  −q   0  |  0   0   0   0   1   0 ]
 vz2'  [  0   0   q   0   0  −q  |  0   0   0   0   0   1 ]

 = [[ I , 0 ], [ αK , I ]]
```

**Full drift matrix:**

```
          x1  y1  z1  x2  y2  z2 | vx1 vy1 vz1 vx2 vy2 vz2
 x1'   [  1   0   0   0   0   0  |  Δt  0   0   0   0   0 ]
 y1'   [  0   1   0   0   0   0  |  0   Δt  0   0   0   0 ]
 z1'   [  0   0   1   0   0   0  |  0   0   Δt  0   0   0 ]
 x2'   [  0   0   0   1   0   0  |  0   0   0   Δt  0   0 ]
 y2'   [  0   0   0   0   1   0  |  0   0   0   0   Δt  0 ]
 z2'   [  0   0   0   0   0   1  |  0   0   0   0   0   Δt]
       [ ------------------------+------------------------ ]
 vx1'  [  0   0   0   0   0   0  |  1   0   0   0   0   0 ]
  ...  (identity in the velocity block, zeros elsewhere)
 vz2'  [  0   0   0   0   0   0  |  0   0   0   0   0   1 ]

 = [[ I , Δt·I ], [ 0 , I ]]
```

For LJ at $r = 1.2$ (reduced units, $\Delta t = 0.01$): $p = 0.005\times9.530 = 0.0476$, $q = 0.005\times(2.2117/1.2) = 0.0092$. The $-p = -0.0476$ matches the off-diagonal entry of $J_1$ in [§16](#16-worked-example-lennard-jones-jacobian-with-real-numbers).

### Determinant

Both matrices are **block-triangular** with identity diagonal blocks:

$$
\det\begin{pmatrix}\mathbf I & \mathbf 0\\ \alpha\mathbf K & \mathbf I\end{pmatrix} = 1,\qquad
\det\begin{pmatrix}\mathbf I & \Delta t\,\mathbf I\\ \mathbf 0 & \mathbf I\end{pmatrix} = 1,\qquad
\det J_{\text{VV}} = 1
$$

$p$ and $q$ never enter — only the **pattern of zeros** matters, so it holds for LJ, EAM, MACE, anything. Different masses: $\alpha\mathbf K \to \tfrac{\Delta t}{2}\mathbf M^{-1}\mathbf K$, still off-diagonal.

**Euler:**

$$
J_{\text{Euler}} = \begin{pmatrix}\mathbf I & \Delta t\,\mathbf I\\ \Delta t\,\mathbf M^{-1}\mathbf K & \mathbf I\end{pmatrix},
\qquad
\det = \det\big(\mathbf I - \Delta t^2\mathbf M^{-1}\mathbf K\big) \neq 1
$$

(1D spring: $\mathbf K = -k$ → $1 + \Delta t^2 k/m$.)

### In real simulations

No code builds $J$ — it's a proof tool. $\mathbf K$ is sparse (non-zero only for pairs within cutoff). $\mathbf K$ itself is used in normal-mode/phonon analysis, and sets the stability limit: velocity-Verlet is unstable for $\Delta t > 2/\omega_{\max}$.

---

## 18. Verlet and leap-frog are also reversible and area-preserving

*(Not explicitly in the notes; follows from §10–11.)*

- **Verlet (5)** is symmetric under $t+\Delta t \leftrightarrow t-\Delta t$, so reversible by construction. It produces the same trajectory as velocity-Verlet, so it inherits area preservation.
- **Leap-frog** is a full kick (6) followed by a drift (7): two shears → $\det J = 1$. It's an exact rearrangement of Verlet, so it is also reversible.
- **Symplectic Euler** ($v_1 = v_0 + \tfrac{\Delta t}{m}f(x_0)$, then $x_1 = x_0 + v_1\Delta t$): kick then drift → area-preserving ($\det = 1$), but **not** time-reversible (not a palindrome; it's first order). Adding a half kick at each end (kick–drift–kick) makes it symmetric — that's velocity-Verlet.

---

## 19. Computation of forces

The most expensive step of MD.

### Step 1 — total potential energy (pairwise)

$$
U = \frac12\sum_{i}\sum_{j\ne i} u(r_{ij}) = \sum_{i}\sum_{j>i} u(r_{ij})
$$

The ½ corrects for counting each pair twice, as $(i,j)$ and $(j,i)$.

### Step 2 — force is minus the gradient

$$
\mathbf f_i = -\nabla_i U
$$

Only terms containing atom $i$ survive. Atom $i$ appears in both $u(r_{ij})$ and $u(r_{ji})$ (equal), which cancels the ½:

$$
\mathbf f_i = -\sum_{j\ne i}\nabla_i u(r_{ij})
$$

### Step 3 — chain rule

$$
\nabla_i u(r_{ij}) = \frac{du}{dr_{ij}}\,\nabla_i r_{ij}
$$

From $r_{ij} = \sqrt{(x_i-x_j)^2 + (y_i-y_j)^2 + (z_i-z_j)^2}$:

$$
\frac{\partial r_{ij}}{\partial x_i} = \frac{x_i - x_j}{r_{ij}}\quad(\text{same for } y, z)
\qquad\Rightarrow\qquad
\nabla_i r_{ij} = \frac{\mathbf r_i - \mathbf r_j}{r_{ij}} = \hat{\mathbf r}_{ij}
$$

(⚠️ The notes write $x_i/r_{ij}$; it should be $(x_i - x_j)/r_{ij}$.) $\hat{\mathbf r}_{ij}$ is the unit vector pointing **from $j$ to $i$**.

### Step 4 — result

$$
\mathbf f_i = -\sum_{j\ne i}\frac{du}{dr_{ij}}\,\hat{\mathbf r}_{ij} = \sum_{j\ne i} f_{ij}\,\hat{\mathbf r}_{ij},
\qquad f_{ij} \equiv -\frac{du}{dr_{ij}}
$$

- $f_{ij} > 0$: pushes $i$ away from $j$ — **repulsion**.
- $f_{ij} < 0$: **attraction**.

**Newton's third law** comes free: $f_{ji}\hat{\mathbf r}_{ji} = -f_{ij}\hat{\mathbf r}_{ij}$. Codes compute each pair once and apply it to both atoms (LAMMPS `newton on`).

### Lennard-Jones

$$
u(r) = 4\varepsilon\left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right]
$$

$$
\frac{du}{dr} = 4\varepsilon\left[-\frac{12\sigma^{12}}{r^{13}} + \frac{6\sigma^{6}}{r^{7}}\right]
= -\frac{24\varepsilon}{r}\left[2\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right]
$$

(⚠️ The notes write $du/dr = +\tfrac{24\varepsilon}{r}[\dots]$ — a missing minus sign. What's written is actually the force magnitude $f_{ij}$.)

$$
f(r) = -\frac{du}{dr} = \frac{24\varepsilon}{r}\left[2\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right]
$$

Vector form (convenient in code — no square root needed for the direction):

$$
\mathbf f_i = \sum_{j\ne i}\frac{24\varepsilon}{r_{ij}^2}\left[2\left(\frac{\sigma}{r_{ij}}\right)^{12} - \left(\frac{\sigma}{r_{ij}}\right)^{6}\right](\mathbf r_i - \mathbf r_j)
$$

- Positive (repulsive) at short range.
- Zero at $r = 2^{1/6}\sigma \approx 1.122\sigma$ (the potential minimum).
- Negative (attractive) beyond.

Reduced units ($\varepsilon=\sigma=1$): $f(r) = 48r^{-13} - 24r^{-7}$; $f'(r) = -624r^{-14} + 168r^{-8}$ (used in §16).

### Cost and how codes reduce it

- Naive: $N(N-1)/2$ pairs → $O(N^2)$.
- **Cutoff** $r_c$ (e.g. $2.5\sigma$): ignore pairs beyond it.
- **Neighbour lists** (Verlet lists): store nearby pairs, rebuild every few steps → ~$O(N)$.
- **Cell lists / linked cells**: bin atoms by position to build neighbour lists in $O(N)$.
- **Periodic boundaries**: minimum-image convention, $\mathbf r_{ij} \to \mathbf r_{ij} - L\,\mathrm{round}(\mathbf r_{ij}/L)$.
- Truncation makes $u$ jump at $r_c$ → small energy drift; fix by shifting the potential (LAMMPS `pair_modify shift yes`) or smooth switching.

---

## 20. Generating initial velocities

### The target: Maxwell–Boltzmann

$$
f(v_x, v_y, v_z) = \left(\frac{m}{2\pi k_BT}\right)^{3/2}\exp\!\left[-\frac{m(v_x^2 + v_y^2 + v_z^2)}{2k_BT}\right]
$$

A Gaussian in each component: mean 0, variance $k_BT/m$.

### Derivation from the canonical ensemble

$$
P(\mathbf R, \mathbf V) \propto e^{-E/k_BT},\qquad E = K(\mathbf V) + U(\mathbf R),\qquad K = \sum_i \tfrac12 m_i v_i^2
$$

$$
e^{-E/k_BT} = \underbrace{e^{-K/k_BT}}_{\text{only }\mathbf V}\times\underbrace{e^{-U/k_BT}}_{\text{only }\mathbf R}
$$

So velocities are independent of positions, and $K$ is a sum over atoms and directions, so each component is an independent Gaussian:

$$
p(v_x) = \sqrt{\frac{m}{2\pi k_BT}}\,\exp\!\left(-\frac{mv_x^2}{2k_BT}\right),
\qquad \langle v_x\rangle = 0,\qquad \langle v_x^2\rangle = \frac{k_BT}{m}
$$

**Consequences:**
- (a) **Independent of the potential** — same for gas, liquid, solid, LJ, MACE. So velocities can be generated before knowing anything about structure or forces.
- (b) **Heavier atoms move slower**: spread $\sqrt{k_BT/m}$. Argon at 300 K: $\sqrt{1.38\times10^{-23}\times300/6.63\times10^{-26}} \approx 250$ m/s per component.
- (c) **Equipartition**: $\langle\tfrac12 mv_x^2\rangle = \tfrac12k_BT$ per degree of freedom. Defines the instantaneous temperature:

$$
T_{\text{inst}} = \frac{\sum_i m_i v_i^2}{N_f\,k_B},\qquad N_f = 3N - 3
$$

**Component vs speed:** components are Gaussian about 0; the speed $\lvert\mathbf v\rvert$ is not:

$$
p(\lvert\mathbf v\rvert) \propto \lvert\mathbf v\rvert^2\exp\!\left(-\frac{m\lvert\mathbf v\rvert^2}{2k_BT}\right)
$$

(the $\lvert\mathbf v\rvert^2$ is the sphere-surface factor $4\pi v^2$). Sample **components**, not speed.

### The notes' shortcut: uniform distribution

Only the first two moments must be right:
- $\langle v_\alpha\rangle = 0$ — centre of mass at rest (no drift).
- $\langle v_\alpha^2\rangle = k_BT/m$ — correct temperature.

Uniform on $[-a, a]$ has variance $a^2/3$ → $a = \sqrt{3k_BT/m}$.

### The recipe

1. Draw random numbers per component (Gaussian, variance 1; or uniform).
2. Scale: Gaussian $v = r\sqrt{k_BT/m}$; uniform on $[-1,1]$: $v = r\sqrt{3k_BT/m}$.
3. Remove centre-of-mass velocity: $\mathbf V_{cm} = \sum m_i\mathbf v_i/\sum m_i$, $\mathbf v_i \to \mathbf v_i - \mathbf V_{cm}$. (Otherwise the box drifts forever — momentum is conserved — and drift KE is wrongly counted as temperature. This removes 3 DOF → $N_f = 3N - 3$.)
4. Rescale to the exact target: $\mathbf v_i \to \mathbf v_i\sqrt{T_{\text{target}}/T_{\text{inst}}}$.

### Why a uniform start is good enough (figure)

![Velocity histogram relaxing to Maxwell–Boltzmann](figures/mb_hist.png)

500 LJ atoms (fcc, $\rho = 0.8$), uniform velocities at $T = 2.0$, plain velocity-Verlet NVE. The flat box becomes the Gaussian within **150 steps**. Each collision mixes velocities; a velocity becomes a sum of many random kicks → Gaussian (central limit theorem). The **shape** is forgotten fast; the **energy** is not (NVE conserves it). So only mean and variance must be right.

### The temperature drops at the start (figure)

![Temperature drop from a lattice start](figures/mb_T.png)

$T$ fell from 2.0 to ≈1.1 in ~0.1 time units. On a perfect lattice potential energy is at its minimum, so all energy starts kinetic; it then flows into potential energy. For a harmonic solid it splits equally → $T$ halves (LJ is anharmonic → ≈1.1). Remedies: start at ~$2T$ in NVE, or (usual) equilibrate with a thermostat and discard that stretch. From a disordered liquid start the drop is much smaller.

---

## 21. Why Maxwell–Boltzmann? (from counting)

### ❓ My question: "I don't understand why MB again?"

Two separate questions: **(A)** why do real atoms at equilibrium have MB velocities? **(B)** why do we generate MB in MD?

### A. The energy-sharing game (figure)

![Energy sharing game converges to an exponential](figures/game.png)

20,000 atoms all start with energy 1. Each random collision pools two atoms' energy and splits it randomly (energy conserved). After enough collisions the histogram settles into an **exponential** $P(\varepsilon)\propto e^{-\varepsilon/\langle\varepsilon\rangle}$ — nobody put it in.

### Why exponential: counting arrangements

100 atoms share 100 units. Probability that *one* atom holds $n$ units ∝ number of ways the other 99 can share $100 - n$:

| $n$ | Ways for the other 99 | Probability | Ratio to previous |
|---|---|---|---|
| 0 | $2.25\times10^{58}$ | 0.498 | |
| 1 | $1.14\times10^{58}$ | 0.251 | 0.505 |
| 2 | $5.72\times10^{57}$ | 0.126 | 0.503 |
| 3 | $2.86\times10^{57}$ | 0.063 | 0.500 |
| 4 | $1.42\times10^{57}$ | 0.031 | 0.497 |
| 5 | $7.04\times10^{56}$ | 0.016 | 0.495 |

Every extra unit cuts the rest's arrangements by the **same factor** (≈½). Repeated multiplication by a constant factor *is* an exponential:

$$
P(n) \propto \left(\tfrac12\right)^n = e^{-n\ln2}
$$

One atom hogging energy leaves the rest with far fewer arrangements; each extra unit costs the same fraction. That constant fraction is what **temperature** measures (hot: cheap to take one more unit; cold: expensive):

$$
P(\varepsilon) \propto e^{-\varepsilon/k_BT} \qquad\text{(Boltzmann factor — a result of counting, not an imposed law)}
$$

### From energy to velocity

$\varepsilon = \tfrac12 mv_x^2$ ⇒

$$
P(v_x) \propto \exp\!\left(-\frac{\tfrac12 mv_x^2}{k_BT}\right) = \exp\!\left(-\frac{v_x^2}{2\,(k_BT/m)}\right)
$$

Energy ∝ velocity², so exponential in energy → **bell curve in velocity**: $v_x = 0$ most likely (cheapest), $\pm v$ equally likely, large $\lvert v_x\rvert$ very rare. Width $\sqrt{k_BT/m}$: hotter → wider; heavier → narrower. That's Maxwell–Boltzmann. The LJ simulation in §20 plays exactly this game.

### B. Why generate MB in MD

MD aims to sample **equilibrium** at $T$, and equilibrium velocities *are* MB. Starting with MB → velocities already in shape, only positions relax. Starting with anything else still works (collisions reshape it in a few hundred steps) — just slightly slower.

> It's what the velocities will become anyway. Starting there saves time. What must be exactly right is the **total kinetic energy** (sets $T$) and **zero total momentum** (no drift). The shape fixes itself.

---

## 22. How LAMMPS generates velocities

```
velocity  all  create  1200.0  4928459  dist gaussian  mom yes  rot no  loop all
          │      │       │        │          │            │       │       │
        group  action  target   random    shape of      remove  remove  how random numbers
                        T       seed      the draw      drift   spin    map to atoms
```

Defaults: `dist uniform`, `mom yes`, `rot no`, `loop all`.

### Inside `velocity.cpp`

1. **Draw** per component: `uniform` → $r = U(0,1) - 0.5 \in [-0.5, 0.5]$; `gaussian` → $r \sim N(0,1)$. (2D: $v_z = 0$.)
2. **Divide by $\sqrt m$:** $v = r/\sqrt m$. *Shortcut:* no $\sqrt{k_BT/m}$ yet — only relative sizes between species (equal average KE per atom). So steps 1–3 need no $k_B$, $T$, or units → same code for `metal`, `real`, `lj`.
3. **`mom yes`:** subtract $\mathbf V_{cm}$ → total momentum exactly zero. **`rot yes`:** compute angular momentum $\mathbf L$, inertia tensor $\mathbf I$, $\boldsymbol\omega = \mathbf I^{-1}\mathbf L$, subtract $\boldsymbol\omega\times\mathbf r_i$ (for isolated clusters/droplets; leave off in periodic bulk).
4. **Rescale** with the group's temperature compute (default `compute temp`, same as thermo): $T_{\text{now}} = \sum m v^2/(N_f k_B)$, $N_f = 3N - 3$; $\mathbf v \to \mathbf v\sqrt{T_{\text{target}}/T_{\text{now}}}$. Units enter here. Rescaling **last** makes step-0 thermo $T$ exactly the target.

### `loop` keyword (parallel runs)

| Option | How | Same velocities on any core count? |
|---|---|---|
| `loop all` (default) | every processor runs the RNG over all atom IDs 1…N, keeps its own | ✅ (slow for huge $N$) |
| `loop local` | each processor draws for its own atoms only | ❌ (fastest) |
| `loop geom` | random numbers seeded from each atom's coordinates | ✅ and fast |

### Python equivalent

```python
import numpy as np
kB = 8.617333e-5            # eV/K (units metal)
mvv2e = 1.0364269e-4        # amu·(Å/ps)^2 -> eV

rng = np.random.default_rng(4928459)
r = rng.standard_normal((N, 3))         # step 1: dist gaussian
v = r / np.sqrt(mass)[:, None]          # step 2: divide by sqrt(m)

vcm = (mass[:, None] * v).sum(0) / mass.sum()
v -= vcm                                # step 3: mom yes

Nf = 3 * N - 3
T_now = mvv2e * (mass[:, None] * v**2).sum() / (Nf * kB)
v *= np.sqrt(1200.0 / T_now)            # step 4: rescale
```

### ❓ CrCl₃ at 1200 K: same speed or same kinetic energy?

Same average **kinetic energy** (equipartition, ½$k_BT$ per DOF), not speed:

$$
\frac{v_{\text{rms}}(\text{Cr})}{v_{\text{rms}}(\text{Cl})} = \sqrt{\frac{m_{\text{Cl}}}{m_{\text{Cr}}}} = \sqrt{\frac{35.45}{52.00}} = 0.83
$$

Cr atoms move ~17% slower.

### Related commands

| Command | Does |
|---|---|
| `velocity all scale 300` | step 4 only |
| `velocity all zero linear` | step 3 only |
| `velocity all create 300 seed sum yes` | add new random velocities to existing |
| `velocity all set 0 0 0` | all velocities zero |
| `fix momentum 100 linear 1 1 1` | remove drift during the run (e.g. with Langevin) |

### Practical tips

1. Step-0 `Temp` in thermo should equal the target exactly.
2. `fix shake` / `fix rigid` reduce $N_f$ below $3N-3$; the compute may learn this only at run init → do `run 0` then `velocity all scale T`.
3. Expect the temperature drop from a crystal start; equilibrate with a thermostat.
4. Change the seed for statistically independent runs.

---

## 23. Corrections to the handwritten notes

| Where | Notes say | Should be |
|---|---|---|
| Forces, $\partial r_{ij}/\partial x_i$ | $x_i/r_{ij}$ | $(x_i - x_j)/r_{ij}$ |
| LJ derivative | $\dfrac{du}{dr} = +\dfrac{24\varepsilon}{r}\Big[2\big(\tfrac{\sigma}{r}\big)^{12} - \big(\tfrac{\sigma}{r}\big)^6\Big]$ | $\dfrac{du}{dr} = -\dfrac{24\varepsilon}{r}[\dots]$; the written expression is the force $f_{ij} = -du/dr$ |

---

## 24. Self-check questions with answers

1. **Which step is likely to blame for an immediate blow-up?** Step 2 (overlaps / bad minimisation) most often; also step 5 (Δt), step 4 (wrong potential/units/mapping), step 3 (T typo). Check thermo at step 0.
2. **1000 flexible waters:** 3000 atoms; 9000 second-order ODEs; 18000 first-order; 18000-D phase space; positions alone insufficient.
3. **Integrator A (bounded ±0.5%) vs B (+0.01%/1000 steps) for $10^7$ steps:** A. B drifts ~100%.
4. **Area ×0.98 per step:** patch shrinks and spirals to (0,0); energy decays (hidden friction).
5. **Second velocity-Verlet step for the spring:** $x_2 = 0.98005$, $f_2 = -0.98005$, $v_2 = -0.19850$ (exact 0.980067, −0.198669).
6. **Derive the central-difference acceleration and its error:** add the forward/backward Taylor series of $\mathbf r$:
   $\mathbf a(t) = \dfrac{\mathbf r(t+\Delta t) - 2\mathbf r(t) + \mathbf r(t-\Delta t)}{\Delta t^2} - \tfrac1{12}\mathbf c\,\Delta t^2 + \cdots$ → error $O(\Delta t^2)$.
7. **Symplectic Euler — area-preserving? reversible?** Area-preserving (kick then drift, both shears, $\det = 1$). Not reversible (not symmetric; first order).
8. **Drift ($x' = x + 0.5v$, $v' = v$) on the unit square:** (0,0)→(0,0), (1,0)→(1,0), (1,1)→(1.5,1), (0,1)→(0.5,1). Horizontal shear; area $= 1\cdot1 - 0.5\cdot0 = 1$.
9. **6×6 kick Jacobian for one atom in 3D:** top-left $\mathbf I_3$; top-right $\mathbf 0$; bottom-left $\tfrac{\Delta t}{2m}\,\partial\mathbf f/\partial\mathbf r$ (3×3 force derivatives); bottom-right $\mathbf I_3$. $\det = 1$.
10. **Why do rows 7–12 of the 12×12 kick sum to zero across columns 1–6?** Translational invariance: moving both atoms by the same amount changes no force (forces depend only on separation; also Newton's third law).
11. **Atoms on the $y$-axis instead:** $\pm p$ in the $y$ rows/columns ($v_{y1}', v_{y2}'$ vs $y_1, y_2$); $\pm q$ in $x$ and $z$.
12. **Euler at the LJ wall, $x_0 = 1.0$, $f'(1.0) = -456$, $\Delta t = 0.01$:** $\det J_{\text{Euler}} = 1 - \Delta t^2 f' = 1 + 0.0456 = 1.0456$ — ~4.6% area growth per step. The wall is stiffest (largest $\lvert f'\rvert$), so error per step is largest there; it also sets the time-step stability limit.
13. **CrCl₃ at 1200 K:** same average KE; Cr ~17% slower.

---

## 25. Quick-reference summary

**Integrators**

| | Equations | Stores | Force calls/step | Reversible | Area-preserving | Global error |
|---|---|---|---|---|---|---|
| Velocity-Verlet | (1), (2) | $\mathbf r, \mathbf v, \mathbf f$ at $t$ | 1 | ✅ | ✅ | $O(\Delta t^2)$ |
| Verlet | (5) | $\mathbf r(t), \mathbf r(t-\Delta t)$ | 1 | ✅ | ✅ | $O(\Delta t^2)$ |
| Leap-frog | (6), (7) | $\mathbf r(t), \mathbf v(t-\tfrac{\Delta t}2)$ | 1 | ✅ | ✅ | $O(\Delta t^2)$ |
| Euler | — | $\mathbf r, \mathbf v$ | 1 | ❌ | ❌ ($\det = 1+\Delta t^2$ spring) | $O(\Delta t)$ |
| RK4 | — | — | 4 | ❌ | ❌ | $O(\Delta t^4)$ but drifts |

**Key formulas**

$$
\mathbf r(t+\Delta t) = \mathbf r + \mathbf v\Delta t + \frac{\mathbf f}{2m}\Delta t^2,\qquad
\mathbf v(t+\Delta t) = \mathbf v + \frac{\mathbf f(t) + \mathbf f(t+\Delta t)}{2m}\Delta t
$$

$$
\mathbf r(t+\Delta t) = 2\mathbf r(t) - \mathbf r(t-\Delta t) + \frac{\mathbf f}{m}\Delta t^2
$$

$$
\mathbf v(t+\tfrac{\Delta t}{2}) = \mathbf v(t-\tfrac{\Delta t}{2}) + \frac{\mathbf f}{m}\Delta t,\qquad
\mathbf r(t+\Delta t) = \mathbf r(t) + \mathbf v(t+\tfrac{\Delta t}{2})\Delta t
$$

$$
\mathbf f_i = \sum_{j\ne i}\frac{24\varepsilon}{r_{ij}^2}\Big[2\Big(\frac{\sigma}{r_{ij}}\Big)^{12} - \Big(\frac{\sigma}{r_{ij}}\Big)^{6}\Big](\mathbf r_i - \mathbf r_j)
$$

$$
p(v_x) = \sqrt{\frac{m}{2\pi k_BT}}\,e^{-mv_x^2/2k_BT},\qquad
T_{\text{inst}} = \frac{\sum m v^2}{(3N-3)k_B}
$$

**Big ideas**

- MD = Newton's laws + ergodicity → time averages stand in for ensemble averages (only if the run samples all relevant states).
- A good integrator copies the *structure* of exact dynamics (reversibility, volume preservation) rather than chasing per-step accuracy; bounded energy error follows.
- Symmetric Taylor expansions cancel unknown derivatives for free and give reversibility.
- Area preservation = Jacobian determinant 1; velocity-Verlet achieves it because it is built from shears (one variable updated at a time).
- Velocities: right mean (0) and variance ($k_BT/m$) matter; the Maxwell–Boltzmann shape emerges from collisions anyway, because it is what energy sharing produces.
