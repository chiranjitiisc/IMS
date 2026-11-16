# Lecture 10: NPT Monte Carlo — Complete Notes from Scratch

*Detailed notes based on 10.pdf ("NPT Simulation", 3 February 2022). Every result is derived from the basic postulate of statistical mechanics, step by step. The three volume-move acceptance rules were checked with an ideal-gas test, and the full algorithm with a working Lennard-Jones NPT code (Section 19). Follows on from the NVT, GCMC and block-averaging notes.*

---

## Contents

0. Notation
1. What NPT means and why we need it
2. The one postulate of statistical mechanics
3. Where temperature and pressure come from (two systems sharing energy and volume)
4. System plus reservoir: the probability of a state
5. What the labels mean: why $E_i$ but only $V$
6. The NPT partition function (general form)
7. Classical systems: from a probability per state to a probability density
8. Integrating out the momenta: where $\Lambda$ comes from
9. Scaled coordinates: why $d\mathbf r = L^3 d\mathbf s$
10. What $U(\mathbf s^N; L)$ means
11. Counting volumes: $V^N$ or $V^{N-1}$?
12. The final partition function and probability density, with checks
13. Monte Carlo recap: detailed balance and Metropolis
14. Move 1: displacement
15. Move 2: volume change — three ways, and the Jacobian
16. Cycles
17. The complete NPT algorithm, with a worked volume move
18. What to measure: density, pressure, compressibility, enthalpy
19. Working Python code and verification
20. Corrections and clarifications to the lecture notes
21. Common mistakes
22. Practice problems with answers
23. Summary

---

## 0. Notation

| Symbol | Meaning |
|---|---|
| $N$ | number of atoms (fixed in NPT) |
| $P$ | pressure imposed by the surroundings (fixed in NPT) |
| $T$ | temperature (fixed in NPT) |
| $V$ | box volume (fluctuates in NPT) |
| $L$ | edge of the cubic box, $L = V^{1/3}$ |
| $\rho$ | number density, $\rho = N/V$ |
| $k_B$ | Boltzmann's constant |
| $\beta$ | $1/(k_BT)$ |
| $E$, $E_i$ | energy; energy of microstate $i$ |
| $K$ | kinetic energy, $K = \sum_j \lvert\mathbf p_j\rvert^2/(2m)$ |
| $U$ | potential energy (from interactions between atoms) |
| $\mathbf r_j$, $\mathbf r^N$ | position of atom $j$; all $N$ positions together |
| $\mathbf p_j$, $\mathbf p^N$ | momentum of atom $j$; all $N$ momenta together |
| $\mathbf s_j$, $\mathbf s^N$ | scaled position $\mathbf s_j = \mathbf r_j/L$ (each component from 0 to 1); all of them together |
| $m$ | mass of one atom |
| $h$ | Planck's constant |
| $\Lambda$ | thermal de Broglie wavelength, $\Lambda = h/\sqrt{2\pi m k_BT}$ (derived in Section 8) |
| $\Omega$ | number of microstates of an isolated system |
| $\Omega_R$, $\Omega_{tot}$ | number of microstates of the reservoir; of system + reservoir together |
| $S$ | entropy, $S = k_B\ln\Omega$ |
| $Q(N,V,T)$ | canonical (NVT) partition function |
| $\Delta(N,P,T)$ | isothermal–isobaric (NPT) partition function |
| $A$, $G$ | Helmholtz free energy $A = -k_BT\ln Q$; Gibbs free energy $G = -k_BT\ln\Delta$ |
| $f(V)$ | volume scale used to count box sizes (Section 11) |
| $p(\cdot)$ | probability, or probability density |
| $o$, $n$ | old (current) and new (trial) state |
| $\alpha(o\to n)$ | probability (density) of **proposing** $n$ from $o$ |
| $\text{acc}(o\to n)$ | probability of **accepting** that proposal |
| $\Delta U$, $\Delta V$ | $U_n - U_o$; $V_n - V_o$ |
| $u$ | a uniform random number in $[0,1)$ (the lecture's `ran2(iseed)`) |
| $\Delta r_m$ | maximum displacement per coordinate |
| $\Delta V_{max}$, $\Delta\ln V_{max}$, $\Delta L_{max}$ | maximum step in $V$, $\ln V$, $L$ |
| $r_c$ | cutoff distance of the pair potential |
| $N_c$, $M$ | number of cycles; total number of attempted moves |

> **Notation warning.** The lecture uses $f$ for two different things: the volume scale $f(V)$ and the probabilities $f(o)$, $f(n)$. In these notes **$f(V)$ is only the volume scale**; probabilities are always $p$.

> **Symbol warning.** $\Delta$ (capital delta, alone) is the partition function. $\Delta U$, $\Delta V$ (delta in front of something) mean a change.

---

## 1. What NPT means and why we need it

In **NVT** the number of atoms $N$, the volume $V$ and the temperature $T$ are fixed.

In **NPT** the pressure $P$ is fixed instead of $V$. Picture the box closed by a frictionless piston with a constant weight on top, all sitting in a heat bath:

- the bath fixes $T$ (heat flows in and out),
- the weight fixes $P$ (the piston moves until the inside pushes back equally hard),
- the walls are closed to atoms, so $N$ is fixed.

| Fixed | Fluctuates |
|---|---|
| $N$, $P$, $T$ | $V$ and the energy |

**Why NPT is useful.**

- Experiments are usually done at fixed pressure (often 1 atm), not fixed volume.
- You don't need to know the density in advance: the simulation finds the equilibrium density at the given $P$ and $T$.
- Crystals and dense phases can relax their box size instead of being squeezed or stretched by a wrong guess of $V$.
- Volume fluctuations give the compressibility directly (Section 18).

---

## 2. The one postulate of statistical mechanics

**Postulate (equal a priori probabilities).** For an **isolated** system (fixed energy $E$, volume $V$, number of atoms $N$), every microstate consistent with those values is equally likely.

A **microstate** is one complete specification of the system: every atom's position and momentum (or, quantum mechanically, one quantum state).

Let $\Omega(E,V,N)$ be the number of such microstates. Each one has probability $1/\Omega$. Boltzmann's definition of entropy is

$$S = k_B\ln\Omega$$

Everything below — the NPT partition function, the probability density, the acceptance rules — follows from this postulate plus careful counting.

---

## 3. Where temperature and pressure come from

Before we can use a reservoir, we must know what "the reservoir's temperature and pressure" mean in terms of counting. Take two systems, 1 and 2, separated by a wall that **moves** and **conducts heat**. Together they are isolated:

$$E_1 + E_2 = E, \qquad V_1 + V_2 = V \qquad (\text{both totals fixed})$$

### 3.1 Counting joint states

For a given split $(E_1, V_1)$, each microstate of system 1 can be combined with each microstate of system 2. The number of joint microstates is the product:

$$\Omega_{joint}(E_1,V_1) = \Omega_1(E_1,V_1)\;\Omega_2(E-E_1,\;V-V_1)$$

### 3.2 The most probable split

All joint microstates are equally likely (Section 2), so the split that is actually observed is the one with the **most** joint microstates. For macroscopic systems this maximum is extraordinarily sharp. Maximise the logarithm (same maximum, easier algebra):

$$\ln\Omega_{joint} = \ln\Omega_1(E_1,V_1) + \ln\Omega_2(E-E_1,\;V-V_1)$$

Set the derivatives to zero. Because $E_2 = E - E_1$, $dE_2/dE_1 = -1$; the same for $V$:

$$\frac{\partial\ln\Omega_1}{\partial E_1} = \frac{\partial\ln\Omega_2}{\partial E_2}, \qquad \frac{\partial\ln\Omega_1}{\partial V_1} = \frac{\partial\ln\Omega_2}{\partial V_2}$$

### 3.3 What these derivatives are

- **Energy derivative.** Heat stops flowing when temperatures are equal. So $\partial\ln\Omega/\partial E$ is a quantity that is equal exactly when temperatures are equal. With $S = k_B\ln\Omega$ and the thermodynamic definition $\partial S/\partial E = 1/T$:
  $$\frac{\partial\ln\Omega}{\partial E} = \frac{1}{k_BT} = \beta$$
- **Volume derivative.** The wall stops moving when pressures balance. With $\partial S/\partial V = P/T$:
  $$\frac{\partial\ln\Omega}{\partial V} = \frac{P}{k_BT} = \beta P$$

These identifications are not borrowed from nowhere: counting produces two quantities that must be equal at equilibrium, and they behave exactly like temperature and pressure.

---

## 4. System plus reservoir: the probability of a state

### 4.1 The setup

Our system (N atoms) sits inside a **reservoir** R that is enormously larger. They can exchange energy (heat) and volume (through the piston), but not atoms. System + reservoir together are isolated, with fixed totals $E_{tot}$ and $V_{tot}$.

If the system has energy $E_i$ and volume $V$, the reservoir has $E_{tot}-E_i$ and $V_{tot}-V$.

### 4.2 What $\Omega_R$ and $\Omega_{tot}$ are

- $\Omega_R(E_R, V_R)$: the number of microstates of the **reservoir alone** when it holds energy $E_R$ and volume $V_R$.
- $\Omega_{tot}$: the number of microstates of the **whole isolated universe** (system + reservoir). One "total microstate" is a **pair**: (one microstate of the system, one microstate of the reservoir).

To count $\Omega_{tot}$, sort the total microstates by what the system is doing. If the system is in microstate $i$ with volume $V$, the reservoir can be in any of $\Omega_R(E_{tot}-E_i, V_{tot}-V)$ states. Adding over every system state counts everything once:

$$\Omega_{tot} = \sum_V\sum_i \Omega_R(E_{tot}-E_i,\;V_{tot}-V)$$

### 4.3 The probability of one system state

By the postulate, every total microstate has probability $1/\Omega_{tot}$. The system is in state $(i, V)$ whenever the reservoir is in any of its matching states, so ("favourable over total")

$$p(i,V) = \frac{\Omega_R(E_{tot}-E_i,\;V_{tot}-V)}{\Omega_{tot}}$$

This is exact.

### 4.4 Toy example with numbers

To keep it small, fix the volume and exchange only energy.

- System: one atom with three microstates, A ($E=0$), B ($E=1$), C ($E=2$).
- Reservoir: suppose $\Omega_R(E) = 2^E$.
- $E_{tot} = 4$.

| System state | Reservoir energy | $\Omega_R$ |
|---|---|---|
| A ($E=0$) | 4 | $2^4 = 16$ |
| B ($E=1$) | 3 | $2^3 = 8$ |
| C ($E=2$) | 2 | $2^2 = 4$ |
| **Total** | | $\Omega_{tot} = 28$ |

$$p(A) = \tfrac{16}{28},\quad p(B) = \tfrac{8}{28},\quad p(C) = \tfrac{4}{28}$$

Two lessons:

1. **The Boltzmann factor appears.** Each extra unit of system energy halves the probability: $p \propto 2^{-E} = e^{-E\ln 2}$, i.e. $e^{-\beta E}$ with $\beta = \ln 2$. It comes from the reservoir losing states when the system takes energy.
2. **$\Omega_{tot}$ is the partition function in disguise.** $\Omega_{tot}/\Omega_R(E_{tot}) = 28/16 = 1.75 = 1 + \tfrac12 + \tfrac14 = \sum_i e^{-\beta E_i}$.

### 4.5 Taylor expansion of $\ln\Omega_R$

The reservoir is huge, so $E_i \ll E_{tot}$ and $V \ll V_{tot}$. We expand $\ln\Omega_R$ (not $\Omega_R$ itself — see the remark below) to first order:

$$\ln\Omega_R(E_{tot}-E_i,\;V_{tot}-V) \approx \ln\Omega_R(E_{tot},V_{tot}) - \frac{\partial\ln\Omega_R}{\partial E}E_i - \frac{\partial\ln\Omega_R}{\partial V}V$$

Using Section 3.3, the derivatives are $\beta$ and $\beta P$ — the **reservoir's** values, evaluated at one fixed point $(E_{tot},V_{tot})$, so they are fixed numbers that do not depend on the system's state:

$$\ln\Omega_R(E_{tot}-E_i,\;V_{tot}-V) \approx \ln\Omega_R(E_{tot},V_{tot}) - \beta E_i - \beta PV$$

The second-order terms involve derivatives such as $\partial\beta/\partial E_R$, which scale as one over the reservoir's size and vanish for an infinite reservoir.

> **Why expand $\ln\Omega_R$ and not $\Omega_R$?** $\Omega_R$ changes by enormous factors for tiny energy changes, so a straight-line approximation of $\Omega_R$ fails immediately. Its logarithm ($S/k_B$) varies smoothly, so the straight line is accurate. Exponentiating afterwards turns the straight line into exponentials.

### 4.6 Exponentiate

Take $e^{(\dots)}$ of both sides, using $e^{a+b+c} = e^ae^be^c$ and $e^{\ln x} = x$:

$$\Omega_R(E_{tot}-E_i,\;V_{tot}-V) \approx \Omega_R(E_{tot},V_{tot})\;e^{-\beta E_i}\;e^{-\beta PV}$$

### 4.7 Normalise: what happens to the constants

Put this into Section 4.3:

$$p(i,V) = \underbrace{\frac{\Omega_R(E_{tot},V_{tot})}{\Omega_{tot}}}_{C}\;e^{-\beta E_i}\,e^{-\beta PV}$$

$C$ is the same number for every system state. Probabilities must add to 1:

$$\sum_V\sum_i C\,e^{-\beta E_i}e^{-\beta PV} = 1 \quad\Rightarrow\quad C = \frac{1}{\sum_V\sum_i e^{-\beta E_i}e^{-\beta PV}} \equiv \frac{1}{\Delta}$$

$$\boxed{\;p(i,V) = \frac{e^{-\beta E_i}\,e^{-\beta PV}}{\Delta}\;}$$

The huge reservoir numbers never need to be computed: normalisation replaces them with $\Delta$, which depends only on the system. Compared with NVT there is one new factor, $e^{-\beta PV}$: growing the box costs work $PV$ against the surroundings.

---

## 5. What the labels mean: why $E_i$ but only $V$

- **$E_i$ carries the subscript $i$** because the energy depends on the arrangement: microstate $i$ (all positions and momenta) has its own energy.
- **$V$ does not** because the volume is not a property each arrangement has; it is chosen **first**. Fix the piston (choose $V$), then list the arrangements $i$ that fit in that box. Every arrangement in that list has the same $V$, but different energies.

So the full state is the **pair** $(i, V)$: which box size, and which arrangement inside it.

**Example.** Two atoms in a 1D box. For box length 2, the arrangements are all ways to place both atoms between 0 and 2: each has its own $E_i$ (depends on their separation), all have $V = 2$. For box length 3: a new, larger list of arrangements, all with $V = 3$.

To be fully precise, $E_i = E_i(V)$: stretching the box changes the distances and so the energies. Section 10 makes this explicit.

---

## 6. The NPT partition function (general form)

$$\Delta = \sum_V\sum_i e^{-\beta E_i}\,e^{-\beta PV}$$

Inside the inner sum $V$ is fixed, so $e^{-\beta PV}$ is a common factor ($ab_1 + ab_2 + ab_3 = a(b_1+b_2+b_3)$):

$$\boxed{\;\Delta = \sum_V e^{-\beta PV}\underbrace{\sum_i e^{-\beta E_i}}_{Q(N,V,T)}\;}$$

**The NPT partition function is the NVT partition function for each box size, weighted by $e^{-\beta PV}$, summed over all box sizes.** This holds for quantum or classical systems. Sections 7–11 turn it into the classical form used in simulations.

---

## 7. Classical systems: from a probability per state to a probability density

### 7.1 Why we need a density

Positions, momenta and volume are continuous. The probability of one exact value (say $x = 1.3000\ldots$) is zero. Instead we use a **probability density** $p(x)$: the probability of finding $x$ in a small interval of width $dx$ is $p(x)\,dx$.

Example: an atom equally likely anywhere in a 1D box of length 4 has $p(x) = 1/4$; the chance of finding it between 1.0 and 1.1 is $(1/4)(0.1) = 0.025$.

### 7.2 The bridge: states in a region

If a small region contains many states, all with (nearly) the same probability, then

$$\text{prob(region)} = (\text{number of states in the region}) \times (\text{probability of one state})$$

(10 states of probability 0.01 each make a region of probability 0.1.) So we need: how many states fit in a small region of positions, momenta and volume.

### 7.3 Phase-space cells: the factor $h$

One atom in 1D has a state $(x, p_x)$: a point in a plane called **phase space**. Quantum mechanics (the uncertainty principle) says a state cannot be pinned down more finely than a cell of area $h$:

$$\text{number of states in } dx\,dp_x = \frac{dx\,dp_x}{h}$$

In 3D an atom has three such pairs ($x$ with $p_x$, $y$ with $p_y$, $z$ with $p_z$): $d\mathbf r\,d\mathbf p/h^3$. For $N$ atoms: $d\mathbf r^N d\mathbf p^N/h^{3N}$.

### 7.4 Identical atoms: the factor $N!$

Two identical atoms, one at $x = 1$ and one at $x = 2$: the integral counts "A at 1, B at 2" and "B at 1, A at 2" separately, but physically they are the same state. With $N$ atoms there are $N!$ relabellings, so every state is counted $N!$ times. Divide by $N!$:

$$\sum_i \;\longrightarrow\; \frac{1}{N!\,h^{3N}}\int d\mathbf r^N\,d\mathbf p^N$$

### 7.5 Counting volumes: the factor $f(V)$

$V$ is continuous too. $dV$ has units of volume, but a count of states must be a pure number, so divide by some volume $f(V)$:

$$\sum_V \;\longrightarrow\; \int\frac{dV}{f(V)}$$

(The choice of $f(V)$ is discussed in Section 11.)

### 7.6 The full density

A small region $dV\,d\mathbf r^N d\mathbf p^N$ contains $\dfrac{dV}{f(V)}\dfrac{d\mathbf r^N d\mathbf p^N}{N!\,h^{3N}}$ states, each with probability $e^{-\beta(K+U)}e^{-\beta PV}/\Delta$ (since $E = K + U$). Multiplying:

$$\text{prob(region)} = \underbrace{\frac{e^{-\beta K}\,e^{-\beta U}\,e^{-\beta PV}}{N!\,h^{3N}\,f(V)\,\Delta}}_{\text{probability density}}\;dV\,d\mathbf r^N\,d\mathbf p^N$$

The factor multiplying $dV\,d\mathbf r^N d\mathbf p^N$ is, by definition, the probability density. We used $e^{-\beta E} = e^{-\beta K}e^{-\beta U}$: kinetic energy $K$ depends only on momenta, potential energy $U$ only on positions.

---

## 8. Integrating out the momenta: where $\Lambda$ comes from

### 8.1 Why we may

Monte Carlo never moves momenta. We want the probability of positions and volume **whatever the momenta are** — so we add up (integrate) the density over all momenta, like finding a row total in a table.

### 8.2 Why it is easy

Only $e^{-\beta K}$ contains momenta. $K = \sum_j \lvert\mathbf p_j\rvert^2/(2m)$ is a sum of $3N$ separate terms like $p_x^2/(2m)$, so $e^{-\beta K}$ is a product of $3N$ identical factors, and the integral is a product of $3N$ identical one-dimensional integrals.

### 8.3 One Gaussian integral

First $I = \int_{-\infty}^\infty e^{-x^2}dx = \sqrt\pi$. Square it and switch to polar coordinates $(r,\theta)$:

$$I^2 = \iint e^{-(x^2+y^2)}dx\,dy = \int_0^{2\pi}\!\!\int_0^\infty e^{-r^2}\,r\,dr\,d\theta = 2\pi\cdot\tfrac12 = \pi$$

Now substitute $x = \sqrt a\,p_x$ ($dx = \sqrt a\,dp_x$):

$$\int_{-\infty}^\infty e^{-a p_x^2}\,dp_x = \frac{1}{\sqrt a}\int e^{-x^2}dx = \sqrt{\frac{\pi}{a}}$$

Our integrand is $e^{-\beta p_x^2/(2m)}$, so $a = \beta/(2m)$:

$$\int_{-\infty}^\infty e^{-\beta p_x^2/(2m)}\,dp_x = \sqrt{\frac{2\pi m}{\beta}} = \sqrt{2\pi m k_BT}$$

### 8.4 All $3N$ integrals together

$$\frac{(2\pi m k_BT)^{3N/2}}{h^{3N}} = \left(\frac{\sqrt{2\pi mk_BT}}{h}\right)^{3N} = \frac{1}{\Lambda^{3N}}, \qquad \boxed{\Lambda \equiv \frac{h}{\sqrt{2\pi mk_BT}}}$$

$\Lambda$ (a length, the thermal de Broglie wavelength) is simply what is left after integrating out the momenta. The density in positions and volume is now

$$p(V,\mathbf r^N) = \frac{e^{-\beta U(\mathbf r^N)}\,e^{-\beta PV}}{N!\,\Lambda^{3N}\,f(V)\,\Delta}$$

Equivalently, the NVT partition function becomes

$$Q(N,V,T) = \frac{1}{N!\,\Lambda^{3N}}\int_V e^{-\beta U(\mathbf r^N)}\,d\mathbf r^N$$

---

## 9. Scaled coordinates: why $d\mathbf r = L^3 d\mathbf s$

### 9.1 Why we need them

In NPT the box changes size. With real positions, the range of each position ($0$ to $L$) changes whenever $V$ changes — $V$ hides inside the integration limits. Scaled positions

$$\mathbf s_j = \frac{\mathbf r_j}{L}$$

always run from 0 to 1, whatever the box does.

### 9.2 One direction

$d\mathbf r$ is shorthand for a small **volume**, $d\mathbf r = dx\,dy\,dz$ (not a small vector). With $x = L s_x$, and $L$ held fixed while we integrate over positions (the position integral sits inside $Q$ for one particular box):

$$\frac{dx}{ds_x} = L \quad\Rightarrow\quad dx = L\,ds_x$$

Same for $y$ and $z$.

### 9.3 Multiply the three

$$d\mathbf r = dx\,dy\,dz = (L\,ds_x)(L\,ds_y)(L\,ds_z) = L^3\,d\mathbf s = V\,d\mathbf s$$

Each direction is stretched by $L$, so a small cube's volume is stretched by $L^3$.

**Numbers.** $L = 10$, scaled cube of side 0.01: scaled volume $10^{-6}$; real side $0.1$, real volume $10^{-3}$; ratio $1000 = L^3$. ✓

**Whole box.** $\int d\mathbf s = 1$ (unit cube), $\int d\mathbf r = L^3 = V$. ✓

For $N$ atoms each contributes $L^3$: $\;d\mathbf r^N = V^N\,d\mathbf s^N$.

### 9.4 How densities change with variables (the Jacobian rule)

Probability in a small region cannot depend on how you label it:

$$p_s\,d\mathbf s^N = p_r\,d\mathbf r^N \quad\Rightarrow\quad p_s = p_r\,\frac{d\mathbf r^N}{d\mathbf s^N} = V^N\,p_r$$

1D example: an atom uniform in a box of length $L$ has $p_x = 1/L$. With $s = x/L$: $p_s = p_x\,(dx/ds) = (1/L)(L) = 1$ — uniform on $[0,1]$, as it should be. **Changing variables multiplies the density by the stretch factor.** This rule is used again for the volume moves (Section 15).

Applying it:

$$p(V,\mathbf s^N) = \frac{V^N\,e^{-\beta PV}\,e^{-\beta U(\mathbf s^N;L)}}{N!\,\Lambda^{3N}\,f(V)\,\Delta}, \qquad Q = \frac{V^N}{N!\,\Lambda^{3N}}\int_0^1 e^{-\beta U(\mathbf s^N;L)}\,d\mathbf s^N$$

$V^N$ is the "a bigger box has more room" factor, and it came purely from the change of variables.

---

## 10. What $U(\mathbf s^N; L)$ means

The energy depends on **real** distances. Scaled positions alone do not give real distances; you also need the box size, since $\mathbf r_j = L\mathbf s_j$:

$$U(\mathbf s^N; L) \equiv U(\mathbf r^N = L\,\mathbf s^N)$$

The **semicolon** separates the two roles: $\mathbf s^N$ are the variables integrated over in $\int d\mathbf s^N$; $L$ is a parameter held fixed during that integral (it only changes in the outer integral over $V$).

**Example.** Lennard-Jones pair energy $u(r) = 4\varepsilon[(\sigma/r)^{12} - (\sigma/r)^6]$, where $r$ is the real distance, $\varepsilon$ the well depth, $\sigma$ the atom size. Atom 1 at scaled $(0.1,0,0)$, atom 2 at $(0.2,0,0)$: scaled separation 0.1, real separation $0.1L$.

| $L$ | $r = 0.1L$ | $U$ |
|---|---|---|
| $10\sigma$ | $1.0\sigma$ | $0$ |
| $11.22\sigma$ | $1.122\sigma$ (minimum of $u$) | $-\varepsilon$ |
| $15\sigma$ | $1.5\sigma$ | $-0.32\varepsilon$ |

Same scaled positions, different energies.

**Periodic boundaries** work the same way: minimum image in scaled units, $\Delta s \leftarrow \Delta s - \text{round}(\Delta s)$, then multiply by $L$.

**Why this matters for Monte Carlo.** A volume move keeps all $\mathbf s_j$ and changes only $L$. Every pair distance is multiplied by $L_n/L_o$, so every pair energy changes, and $\Delta U = U(\mathbf s^N;L_n) - U(\mathbf s^N;L_o)$ needs a full recalculation over all pairs.

---

## 11. Counting volumes: $V^N$ or $V^{N-1}$?

### 11.1 The two choices

After Section 9 the partition function is

$$\Delta = \frac{1}{N!\,\Lambda^{3N}}\int\frac{e^{-\beta PV}\,V^N}{f(V)}\int e^{-\beta U}\,d\mathbf s^N\,dV$$

This is **exactly the first line of the lecture.** The only open question is $f(V)$.

| Choice | Result | Used by |
|---|---|---|
| constant $f(V) = V_0$ | $V^N$ (constant prefactor) | most textbooks, e.g. Frenkel & Smit |
| $f(V) = V$ | $V^N/V = V^{N-1}$ | **this lecture**, "for periodic systems" |

Statistical mechanics does not by itself say how to count continuous box sizes, so neither choice follows from the postulate alone. The two differ by one power of $V$ out of $N$: effects of order $1/N$, which vanish for large systems. **What matters is consistency:** the acceptance rules must come from the same $\Delta$ you start with. For this course, use the lecture's $V^{N-1}$.

### 11.2 Reason 1 for $f(V) = V$: in a periodic box, sliding everything together changes nothing

**Picture a periodic 1D box as a ring.** An atom leaving the right end re-enters at the left, so the two ends are the same point: a periodic box of length $L$ is a ring of circumference $L$.

**Example: a ring with 6 spots (0–5), two atoms 2 spots apart.**

| Atom 1 | Atom 2 |
|---|---|
| 0 | 2 |
| 1 | 3 |
| 2 | 4 |
| 3 | 5 |
| 4 | 0 |
| 5 | 1 |

Six placements ($= L$). Drawn on a ring, they are **all the same picture, rotated**. The numbers 0–5 are labels we wrote on; the ring has no marked spot. Distance, energy, forces — everything measurable — is identical. There is one physical arrangement, counted $L$ times.

**Compare a box with walls** (spots 0–5, no wrap-around): atoms at (0, 2) press against the left wall; atoms at (2, 4) sit in the middle. Different energies — the walls are a landmark, so these really are different states.

**In the integral.** For two atoms in a periodic box, the energy depends only on the separation $d = x_2 - x_1$. Change variables from $(x_1, x_2)$ to $(x_1, d)$:

$$\int_0^L\!\!\int_0^L e^{-\beta u}\,dx_1\,dx_2 = \underbrace{\int_0^L dx_1}_{=\,L}\;\times\int e^{-\beta u(d)}\,dd$$

The first factor counts where the pair sits as a whole — the "same picture, rotated" freedom. Removing this overcount removes one factor of $L$ (in 3D: one factor of $V$, since the whole group can slide in $x$, $y$ and $z$). What remains are $N-1$ **relative** position vectors.

### 11.3 Two variables, two integrals — don't mix them up

Two different things share the letter $V$ (or $L$):

- $x_1$ (in 3D, $\mathbf r_1$): the **position of atom 1**, running from 0 to $L$;
- $L$ (in 3D, $V$): the **size of the box**, which fluctuates in NPT.

The degree of freedom removed is $x_1$, a position. Its integral **happens to equal** $L$ ($\int_0^L dx_1 = L$), the way "how many positions are in a box of length 5?" gives 5 — you integrated over position, not over box length. The full 1D NPT expression has three integrals:

$$\Delta \propto \underbrace{\int dL}_{\text{box size}}\;e^{-\beta PL}\;\underbrace{\int_0^L dx_1}_{\text{atom 1's position}}\;\underbrace{\int dd}_{\text{separation}}\;e^{-\beta u(d)}$$

Removing the overcount deletes only the middle one:

$$\Delta \propto \int dL\;e^{-\beta PL}\int dd\;e^{-\beta u(d)}$$

As nested loops:

```
for each box size L:              # outer integral, dL (dV in 3D) -- stays
    for each position x1:         # removed: same physical state, counted L times
        for each separation d:    # inner integral, dd -- stays
            add exp(-beta*P*L) * exp(-beta*u(d))
```

In scaled form ($d = L\sigma$, $dd = L\,d\sigma$): $\;\Delta \propto \int dL\,e^{-\beta PL}\,L^{1}\int e^{-\beta u(L\sigma)}d\sigma$, i.e. $L^{N-1}$ with $N = 2$. ✓

**Why the lecture still writes $\int d\mathbf s^N$ over all $N$ atoms.** In scaled coordinates the box has size 1, so the integral over where the whole group sits, $\int d\mathbf s_1$, equals 1. The overcount factor is invisible inside the scaled integral; it appears only as the outer power of $V$.

### 11.4 Where the $1/V$ is written is bookkeeping

$$\underbrace{V^{N-1}}_{\text{positions}}\;e^{-\beta PV}\cdots\;\underbrace{dV}_{\text{volume}} \;=\; \underbrace{V^{N}}_{\text{positions}}\;e^{-\beta PV}\cdots\;\underbrace{\frac{dV}{V}}_{\text{volume}}$$

Left: the physical origin (positions give $V^{N-1}$; volume integrated with plain $dV$). Right: the lecture's notation ($dV/f(V) = dV/V = d(\ln V)$), which is exactly the form used in move (b). Same expression.

### 11.5 Reason 2: no arbitrary volume unit

A constant $V_0$ requires picking a reference volume that nothing in the problem specifies. $dV/V = d(\ln V)$ has no units and gives equal weight to equal **percentage** changes (100 → 101 counts the same as 1000 → 1010). A periodic box has no preferred size, so a scale-free measure fits it.

### 11.6 Reason 3: the ideal gas comes out exactly

With $V^{N-1}$: $\langle V\rangle = Nk_BT/P$ exactly (Section 12.3). With $V^N$: $(N+1)k_BT/P$.

---

## 12. The final partition function and probability density, with checks

### 12.1 The results

$$\boxed{\;\Delta = \frac{1}{N!\,\Lambda^{3N}}\int_0^\infty e^{-\beta PV}\,V^{N-1}\int_0^1 e^{-\beta U(\mathbf s^N;L)}\,d\mathbf s^N\,dV\;}$$

$$\boxed{\;p(V,\mathbf s^N) = \frac{V^{N-1}\,e^{-\beta PV}\,e^{-\beta U(\mathbf s^N;L)}}{N!\,\Lambda^{3N}\,\Delta}\;}$$

$p(V,\mathbf s^N)$ is a density **per unit $dV$ and per unit $d\mathbf s^N$**. It is simply the integrand of $\Delta$ divided by the whole integral — the same rule as $p_i = e^{-\beta E_i}/Q$ in NVT: each state's probability is its own term in the sum divided by the total.

| Piece | Origin |
|---|---|
| $e^{-\beta U}$ | Boltzmann factor for the positions (Sections 4, 7) |
| $e^{-\beta PV}$ | work against the reservoir's pressure (Section 4) |
| $V^N$ | change to scaled coordinates (Section 9) |
| $1/f(V) = 1/V$ | counting volume states (Sections 7.5, 11) |
| $1/N!$ | identical atoms (Section 7.4) |
| $1/\Lambda^{3N}$ | phase-space cells $h$ + momentum integrals (Sections 7.3, 8) |
| $1/\Delta$ | normalisation |

### 12.2 Check: normalisation and units

$\int dV\int d\mathbf s^N\,p = 1$ gives back exactly the boxed $\Delta$, so density and partition function are consistent.

Units: $\mathbf s$ is dimensionless, $\Lambda^{3N}$ has units volume$^N$, and $dV\,V^{N-1}$ also volume$^N$, so $\Delta$ is dimensionless and $p$ has units 1/volume — right for a density per $dV$.

### 12.3 Check: the ideal gas

$U = 0 \Rightarrow \int d\mathbf s^N = 1$. Use $\int_0^\infty V^{N-1}e^{-aV}dV = (N-1)!/a^N$ (integrate by parts $N-1$ times), with $a = \beta P$:

$$\Delta = \frac{1}{N!\,\Lambda^{3N}}\cdot\frac{(N-1)!}{(\beta P)^N} = \frac{1}{N\,\Lambda^{3N}(\beta P)^N}$$

To get $\langle V\rangle$: differentiating $\Delta$ (in its sum form) with respect to $P$ brings down $-\beta V$ from each $e^{-\beta PV}$:

$$\frac{\partial\Delta}{\partial P} = -\beta\langle V\rangle\Delta \quad\Rightarrow\quad \langle V\rangle = -k_BT\,\frac{\partial\ln\Delta}{\partial P}$$

$\ln\Delta = -N\ln P + (\text{terms without } P)$, so

$$\langle V\rangle = \frac{Nk_BT}{P}$$

the ideal gas law, exactly.

### 12.4 Connection to thermodynamics: $G = -k_BT\ln\Delta$

Write $Q = e^{-\beta A}$. Then $\Delta = \int\frac{dV}{f(V)}e^{-\beta[A(V) + PV]}$. $A + PV$ grows with $N$, so the integrand is extremely sharply peaked at the $V^*$ minimising $A + PV$, and

$$-k_BT\ln\Delta \approx \min_V\,[A(V) + PV] = G$$

The minimum condition $\partial A/\partial V = -P$ says the box settles where the internal pressure equals the applied pressure — the piston picture. $\Delta$ plays the role at fixed $(N,P,T)$ that $Q$ plays at fixed $(N,V,T)$.

---

## 13. Monte Carlo recap: detailed balance and Metropolis

(From the Metropolis and NVT notes.) A Markov chain visits states with probability $p$ if **detailed balance** holds:

$$p(o)\,\alpha(o\to n)\,\text{acc}(o\to n) = p(n)\,\alpha(n\to o)\,\text{acc}(n\to o)$$

If the proposal is **symmetric**, $\alpha(o\to n) = \alpha(n\to o)$, the $\alpha$'s cancel:

$$\frac{\text{acc}(o\to n)}{\text{acc}(n\to o)} = \frac{p(n)}{p(o)} \equiv R$$

**Metropolis choice:** $\text{acc}(o\to n) = \min\{1, R\}$. Check:

- $R < 1$: $\text{acc}(o\to n) = R$; reverse ratio $1/R > 1$ so $\text{acc}(n\to o) = 1$; ratio $R/1 = R$. ✓
- $R \ge 1$: $\text{acc}(o\to n) = 1$, $\text{acc}(n\to o) = 1/R$; ratio $R$. ✓

**In NPT the state is $(V, \mathbf s^N)$**, sampled with the density of Section 12.1. Every constant ($N!$, $\Lambda^{3N}$, $\Delta$) cancels in $p(n)/p(o)$, which is why the lecture simply writes $f(o) = V_o^{N-1}e^{-\beta PV_o}e^{-\beta U_o}$.

Each derivation below has two jobs: (1) show the proposal is symmetric, (2) compute $p(n)/p(o)$ **in the variable being stepped in**.

---

## 14. Move 1: displacement

Exactly the NVT move: pick atom $k$, shift each coordinate by $(2u-1)\Delta r_m$, apply periodic boundaries. $V$ is unchanged, so $V^{N-1}$ and $e^{-\beta PV}$ cancel:

$$\text{acc}(o\to n) = \min\left\{1,\ e^{-\beta\Delta U}\right\}$$

$\Delta U$ needs only the pairs involving atom $k$.

---

## 15. Move 2: volume change — three ways, and the Jacobian

### 15.1 What a volume move does

Propose a new box size; **rescale every atom** with the box, $\mathbf r_j \leftarrow \mathbf r_j\,(L_n/L_o)$, so all scaled positions $\mathbf s_j$ stay fixed; recompute $U$ from scratch. The positions get no random proposal of their own (the rescaling is deterministic and identical in both directions), so only the proposal of the box size matters for symmetry.

### 15.2 The key idea: the density depends on the variable you step in

$p(V,\mathbf s^N)$ is a density **per unit $V$**. If random steps are taken in some other variable $y$ (like $\ln V$ or $L$), the symmetric proposal is uniform **in $y$**, so detailed balance needs the density **per unit $y$**. By the Jacobian rule (Section 9.4):

$$p_y\,dy = p_V\,dV \quad\Rightarrow\quad p_y = p_V\,\frac{dV}{dy}$$

**Example ($y = \ln V$, so $dV = V\,dy$).** A step of 0.01 in $\ln V$ covers $dV \approx 0.02$ near $V = 2$ but $dV \approx 0.5$ near $V = 50$ — 25 times more volume. A fixed-width interval in $\ln V$ holds more probability at large $V$, in proportion to $V$; so $p_{\ln V} = V\,p_V$.

### 15.3 (a) Random walk in $V$

**Proposal.** $V_n = V_o + (2u-1)\Delta V_{max}$. $(2u-1)$ is uniform in $[-1,1)$, so $V_n$ is uniform in a window of width $2\Delta V_{max}$ centred on $V_o$: $\alpha(o\to n) = 1/(2\Delta V_{max})$. The reverse window, centred on $V_n$, has the same width and contains $V_o$, so $\alpha(n\to o)$ is the same. Symmetric. ✓

**Ratio** (here $y = V$, no Jacobian):

| Factor | $p(n)/p(o)$ |
|---|---|
| $V^{N-1}$ | $(V_n/V_o)^{N-1}$ |
| $e^{-\beta PV}$ | $e^{-\beta P\Delta V}$, $\Delta V = V_n - V_o$ |
| $e^{-\beta U}$ | $e^{-\beta\Delta U}$, $\Delta U = U(\mathbf s^N;L_n) - U(\mathbf s^N;L_o)$ |

$$\text{acc}(o\to n) = \min\left\{1,\ \left[\frac{V_n}{V_o}\right]^{N-1}e^{-\beta P\Delta V}\,e^{-\beta\Delta U}\right\}$$

If $V_n \le 0$: reject.

### 15.4 (b) Random walk in $\ln V$

**Proposal.** $\ln V_n = \ln V_o + (2u-1)\Delta\ln V_{max}$, symmetric in $\ln V$ by the same window argument.

**Density per unit $\ln V$.** $dV/d(\ln V) = V$:

$$p_{\ln V} = V^{N-1}e^{-\beta PV}e^{-\beta U}\times V = V^{N}e^{-\beta PV}e^{-\beta U}$$

This is the lecture's line $V^N\cdots\frac{dV}{V} = V^N\cdots d(\ln V)$.

$$\text{acc}(o\to n) = \min\left\{1,\ \left[\frac{V_n}{V_o}\right]^{N}e^{-\beta P\Delta V}\,e^{-\beta\Delta U}\right\}$$

The power rises by one **only** because of the Jacobian. $V_n = V_o e^{(\dots)}$ is always positive.

### 15.5 (c) Random walk in $L$

**Proposal.** $L_n = L_o + (2u-1)\Delta L_{max}$, symmetric in $L$.

**Density per unit $L$.** $V = L^3$, $dV/dL = 3L^2$, and $V^{N-1} = L^{3(N-1)}$:

$$p_L = L^{3(N-1)}\times 3L^2\times e^{-\beta PV}e^{-\beta U} = 3\,L^{3N-1}e^{-\beta PV}e^{-\beta U}$$

The constant 3 cancels:

$$\text{acc}(o\to n) = \min\left\{1,\ \left[\frac{L_n}{L_o}\right]^{3N-1}e^{-\beta P\Delta V}\,e^{-\beta\Delta U}\right\}$$

The pressure term still uses $\Delta V = L_n^3 - L_o^3$: the work against the piston is always $P$ times the change in **volume**. If $L_n \le 0$: reject.

### 15.6 Three alternatives, not three parts

(a), (b) and (c) are **alternative ways of doing the same volume move. Pick one and use only that one.** Each, with its own acceptance rule, samples exactly the same distribution. They differ in efficiency and convenience. The key difference is a **fixed amount** versus a **fixed percentage**. Two systems from the NVT notes, both $N = 256$: liquid $V \approx 330$, gas $V \approx 26{,}500$:

| Step choice | Liquid ($V\approx 330$) | Gas ($V\approx 26{,}500$) |
|---|---|---|
| (a) $\Delta V_{max} = 5$ | up to 1.5% — sensible | up to 0.02% — far too small |
| (b) $\Delta\ln V_{max} = 0.02$ | up to 2% | up to 2% |

| Version | Strengths | Weaknesses |
|---|---|---|
| (a) $V$ | simplest | can propose $V \le 0$; step doesn't adapt to system size |
| (b) $\ln V$ | fixed-percentage steps, one setting works at any density; $V$ always positive | — (usual default) |
| (c) $L$ | natural if the code works with $L$; extends to non-cubic boxes (each edge separately) | can propose $L \le 0$ |

**Why the lecture shows all three:** to teach that the acceptance rule depends on the variable you step in, not only on the physics. They can even be mixed (choose (a) or (b) at random each time), provided each uses its own rule and the choice doesn't depend on the current state — but there is no need.

### 15.7 Same numbers, three rules

$N = 100$, $\beta = 1$, $P = 1$, $V_o = 200$, $V_n = 202$ (so $V_n/V_o = 1.01$, $\ln 1.01 = 0.009950$, $\Delta V = 2$), $\Delta U = -0.3$. The common part is $-\beta(\Delta U + P\Delta V) = -1.7$.

| Rule | Volume term | Exponent | Acceptance |
|---|---|---|---|
| (a) | $99\ln 1.01 = 0.9851$ | $-0.7149$ | 0.489 |
| (b) | $100\ln 1.01 = 0.9950$ | $-0.7050$ | 0.494 |
| (c) | $299\ln(1.01^{1/3}) = 0.9917$ | $-0.7083$ | 0.492 |

Nearly identical — but each is correct only for its own kind of step.

---

## 16. Cycles

**A cycle is a unit for counting simulation length. In NPT one cycle is $N + 1$ attempted moves.** During a cycle each move is a displacement with probability $N/(N+1)$ or a volume move with probability $1/(N+1)$. On average: $N$ displacements (about one per atom) and one volume move.

$$M = N_c \times (N+1)$$

Example: $N = 256$, $N_c = 20{,}000$: $M = 5{,}140{,}000$ moves.

**"On average" matters.** The move type is drawn at random every time; a given cycle may contain 0, 1 or 2 volume attempts. A fixed pattern ("N displacements, then one volume move") breaks detailed balance.

**Why count in cycles.**

- **Comparable run lengths:** $10^6$ moves is long for 100 atoms, short for 10,000; "10,000 cycles" means every atom had ~10,000 chances to move, whatever $N$.
- **When to record data:** consecutive moves change one atom, so data are recorded once per cycle (or every few cycles).
- **Run settings:** equilibration, production, tuning intervals and block sizes are all stated in cycles.

**Why one volume move per cycle: cost balance.**

| Per cycle | Cost (pair calculations) |
|---|---|
| $N$ displacements | $N\times N \approx N^2$ |
| 1 volume move (full recompute) | $\approx N^2/2$ |

The single volume move costs about as much as all the displacements together. More volume moves would slow the run without much gain; far fewer would make the volume relax slowly.

---

## 17. The complete NPT algorithm, with a worked volume move

### 17.1 Steps

1. **Inputs:** $N$, $T$, $P$; starting volume $V_0$ (or density $\rho_0$, $V_0 = N/\rho_0$); $\Delta r_m$; $\Delta\ln V_{max}$ (we use move (b)); PRNG seed.
2. **Start:** FCC lattice in a cube of edge $L = V^{1/3}$; compute $U$ including the tail correction.
3. **Choose the move:** draw $u$, set $k = \text{INT}(u(N+1)) + 1$ (a whole number from 1 to $N+1$). If $k \le N$: displace atom $k$. If $k = N+1$: volume move. One draw gives the right probabilities and picks the atom.
4. **Displacement:** Section 14.
5. **Volume move (b):**
   1. Save $V_o$, $L_o$, all positions, $U_o$.
   2. $\ln V_n = \ln V_o + (2u-1)\Delta\ln V_{max}$; $V_n = e^{\ln V_n}$.
   3. $L_n = V_n^{1/3}$.
   4. Rescale: $\mathbf r_j \leftarrow \mathbf r_j\,(L_n/L_o)$ for every atom.
   5. Recompute $U_n$ over all pairs, **including the tail correction at the new density** $N/V_n$; $\Delta U = U_n - U_o$.
   6. Acceptance in logarithms (avoids overflow, since $(V_n/V_o)^N$ can be huge):
      $$\text{arg} = N\ln\frac{V_n}{V_o} - \beta(\Delta U + P\,\Delta V)$$
      Draw a fresh $u$; accept if $\ln u < \text{arg}$ (always if $\text{arg} \ge 0$).
   7. If rejected, **restore everything** saved in step 1.
6. **Cutoff check:** $r_c$ must stay below $L/2$ (minimum image). Reject a volume move that would violate it, or use a bigger system.
7. **Tuning (equilibration only):** $\Delta r_m$ towards ~50% displacement acceptance; $\Delta\ln V_{max}$ towards ~30–50% volume acceptance. Freeze both before production.
8. **Sampling:** after each cycle record $V$, $\rho = N/V$, $U/N$, the virial pressure. Average over production only; error bars by block averaging (Lecture 9 notes).

### 17.2 One volume move with numbers

Reduced units, $N = 256$, $T^* = 1$ ($\beta = 1$), $P^* = 1$, $V_o = 400$, $\Delta\ln V_{max} = 0.02$.

| Step | Calculation | Result |
|---|---|---|
| draw $u = 0.75$ | $(2u-1)\Delta\ln V_{max} = 0.5\times 0.02$ | $\ln V$ changes by 0.01 |
| new volume | $V_n = 400\,e^{0.01}$ | 404.02 |
| $\Delta V$ | $404.02 - 400$ | 4.02 |
| box edges | $L_o = 400^{1/3}$, $L_n = 404.02^{1/3}$ | 7.368 → 7.393 |
| scale factor | $L_n/L_o = e^{0.01/3}$ | 1.00334 |
| energy | recompute; suppose | $\Delta U = +0.50$ |
| arg | $256\times 0.01 - (0.50 + 4.02)$ | $2.56 - 4.52 = -1.96$ |
| acceptance | $e^{-1.96}$ | 0.141 |
| draw $u = 0.10$ | $\ln 0.10 = -2.30 < -1.96$ | **accept** |

The volume factor (+2.56) favours expansion; energy and work against the pressure (−4.52) oppose it.

---

## 18. What to measure: density, pressure, compressibility, enthalpy

### 18.1 Density

$\langle\rho\rangle = \langle N/V\rangle$. This is not exactly $N/\langle V\rangle$ (the average of $1/V$ is not $1/\langle V\rangle$); the difference is small but real. Report $\langle N/V\rangle$, or report $\langle V\rangle$ and say so.

### 18.2 The exact pressure check (new)

Integrate a derivative over all $V$. The integrand vanishes at $V = 0$ (for $N \ge 2$) and at $V\to\infty$, so

$$\int_0^\infty dV\,\frac{\partial}{\partial V}\Big[V^{N-1}e^{-\beta PV}Z(V)\Big] = 0, \qquad Z(V) = \int d\mathbf s^N e^{-\beta U(\mathbf s^N;L)}$$

Carry out the derivative:

$$\frac{\partial}{\partial V}\Big[\cdots\Big] = \Big[\cdots\Big]\left(\frac{N-1}{V} - \beta P + \frac{Z'(V)}{Z(V)}\right)$$

and $Z'/Z = -\beta\langle\partial U/\partial V\rangle$ at fixed $\mathbf s^N$. For pair potentials $r_{ij} = L\,s_{ij}$ gives $\partial r_{ij}/\partial V = r_{ij}/(3V)$, so

$$P_{ex} \equiv -\frac{\partial U}{\partial V} = -\frac{1}{3V}\sum_{i<j} r_{ij}\,u'(r_{ij})$$

the usual virial (excess) pressure. Dividing by the normalisation and averaging:

$$\boxed{\;P = \left\langle \frac{(N-1)\,k_BT}{V} + P_{ex}\right\rangle\;}$$

With the textbook $V^N$ convention the same derivation gives $N$ instead of $N-1$. **The two conventions trade exactness:** the lecture's gives $\langle V\rangle = Nk_BT/P$ exactly for the ideal gas, the textbook's gives $\langle Nk_BT/V\rangle = P$ exactly. Both differences are of order $1/N$.

> **Cutoff caveat.** The derivation needs $U$ to change smoothly with $V$. A plainly truncated potential jumps by $u(r_c)$ at the cutoff, and a volume move pushes pairs across it; the virial misses that contribution. Use a cut-and-**shifted** potential for this check, or expect a small gap (Section 19.3).

### 18.3 Isothermal compressibility from volume fluctuations

$\kappa_T = -\frac{1}{\langle V\rangle}\frac{\partial\langle V\rangle}{\partial P}$. Differentiating $\langle V\rangle = \frac{1}{\Delta}\sum V(\dots)e^{-\beta PV}$ with respect to $P$ gives $\partial\langle V\rangle/\partial P = -\beta(\langle V^2\rangle - \langle V\rangle^2)$, so

$$\kappa_T = \frac{\langle V^2\rangle - \langle V\rangle^2}{k_BT\,\langle V\rangle}$$

This needs long runs: fluctuations converge more slowly than averages.

### 18.4 Enthalpy

$$H = \langle K\rangle + \langle U\rangle + P\langle V\rangle, \qquad \langle K\rangle = \tfrac32 Nk_BT$$

---

## 19. Working Python code and verification

### 19.1 The code (Lennard-Jones, reduced units, move (b))

```python
"""NPT Monte Carlo for Lennard-Jones atoms, reduced units (epsilon, sigma, m = k_B = 1).
Lecture-10 convention: p(V, s^N) ~ V^(N-1) exp(-beta P V) exp(-beta U).
Volume move (b): random walk in ln V, acceptance factor (V_n/V_o)^N."""
import numpy as np, time

def fcc(N, L):
    n = int(round((N/4)**(1/3))); a = L/n
    base = np.array([[0,0,0],[.5,.5,0],[.5,0,.5],[0,.5,.5]])
    pts = [(np.array([i,j,k])+b)*a for i in range(n) for j in range(n) for k in range(n) for b in base]
    return np.array(pts) - L/2

def tail(N, V, rc):                       # tail corrections U_LR and P_LR
    rho = N/V; sr3 = (1/rc)**3
    U = 8/3*np.pi*N*rho*(sr3**3/3 - sr3)
    P = 16/3*np.pi*rho**2*(2/3*sr3**3 - sr3)
    return U, P

def total(r, L, rc):                      # total pair energy U and W = -(1/3) sum r u'(r)
    d = r[:,None,:] - r[None,:,:]; d -= L*np.rint(d/L)
    r2 = (d**2).sum(-1); iu = np.triu_indices(len(r), 1); r2 = r2[iu]
    r2 = r2[r2 < rc*rc]; s6 = 1/r2**3
    U = np.sum(4*(s6*s6 - s6)); W = np.sum(24*(2*s6*s6 - s6))/3
    return U, W

def one(r, i, ri, L, rc):                 # energy of atom i at position ri with all others
    d = r - ri; d -= L*np.rint(d/L); r2 = (d**2).sum(1); r2[i] = np.inf
    r2 = r2[r2 < rc*rc]; s6 = 1/r2**3
    return np.sum(4*(s6*s6 - s6))

def run(N=108, T=2.0, P=2.0, rho0=0.6, rc=2.5, neq=2000, nprod=6000, seed=11):
    rng = np.random.default_rng(seed); beta = 1/T
    V = N/rho0; L = V**(1/3); r = fcc(N, L)
    U, _ = total(r, L, rc); dr, dlnV = 0.15, 0.02
    na = nt = va = vt = 0; samples = []
    for cyc in range(neq + nprod):
        for _ in range(N + 1):
            k = int(rng.random()*(N + 1))          # 0..N ; k == N means volume move
            if k < N:                               # ---- displacement
                rn = r[k] + (2*rng.random(3) - 1)*dr; rn -= L*np.rint(rn/L)
                dU = one(r, k, rn, L, rc) - one(r, k, r[k], L, rc); nt += 1
                if dU <= 0 or rng.random() < np.exp(-beta*dU):
                    r[k] = rn; U += dU; na += 1
            else:                                   # ---- volume move (b)
                Vn = V*np.exp((2*rng.random() - 1)*dlnV); Ln = Vn**(1/3); vt += 1
                if Ln/2 <= rc: continue             # cutoff must stay below L/2
                rn = r*(Ln/L); Un, _ = total(rn, Ln, rc)
                dU = (Un + tail(N, Vn, rc)[0]) - (U + tail(N, V, rc)[0])
                arg = N*np.log(Vn/V) - beta*(dU + P*(Vn - V))
                if arg >= 0 or rng.random() < np.exp(arg):
                    r, V, L, U = rn, Vn, Ln, Un; va += 1
        if cyc < neq and (cyc + 1) % 50 == 0:      # tune during equilibration only
            dr *= 1.05 if na/nt > 0.5 else 0.95
            if vt: dlnV *= 1.05 if va/vt > 0.4 else 0.95
            na = nt = va = vt = 0
        if cyc == neq: na = nt = va = vt = 0
        if cyc >= neq:
            Uc, W = total(r, L, rc); Ut, Pt = tail(N, V, rc)
            samples.append([V, N/V, (Uc + Ut)/N, W/V + Pt])
    s = np.array(samples)
    return s, na/nt, va/max(vt, 1), dr, dlnV, U

def blockerr(x, nb=20):
    m = np.array([b.mean() for b in np.array_split(x, nb)]); return m.std(ddof=1)/np.sqrt(nb)

if __name__ == "__main__":
    t = time.time(); N, T, P = 108, 2.0, 2.0
    s, accd, accv, dr, dlnV, Urun = run(N, T, P)
    V, rho, u, Pex = s.T
    print(f"time {time.time()-t:.0f}s  acc disp {accd:.2f}  acc vol {accv:.2f}")
    print(f"<V> = {V.mean():.2f} +- {blockerr(V):.2f}   <rho> = {rho.mean():.4f} +- {blockerr(rho):.4f}")
    print(f"<U/N> = {u.mean():.4f} +- {blockerr(u):.4f}")
    p1 = (N-1)*T/V + Pex
    print(f"<(N-1)kT/V + P_ex> = {p1.mean():.4f} +- {blockerr(p1):.4f}   (should equal P = {P})")
    print(f"kappa_T = {V.var()/(T*V.mean()):.4f}")
```

Notes on the code: atoms are numbered from 0, so `k == N` is the volume move; the running $U$ stores the pair energy only, and the tail correction is added explicitly when computing $\Delta U$ for a volume move (it is constant during displacements). A first attempt with $r_c = 3.0$ and $N = 108$ at $\rho = 0.6$ gave **0%** volume acceptance: $L/2 = 2.82 < r_c$, so every move was rejected by the cutoff check. Watch for this.

### 19.2 Test 1: ideal gas, all three move types

$U = 0$, $N = 10$, $\beta = 1$, $P = 1$; exact answer $\langle V\rangle = N k_BT/P = 10$. Small $N$ makes an exponent error of one easy to see. $2\times10^6$ volume moves each, error bars by block averaging:

| Move | Exponent used | $\langle V\rangle$ | Verdict |
|---|---|---|---|
| (a) $V$ | $N-1 = 9$ | $9.98 \pm 0.03$ | correct |
| (b) $\ln V$ | $N = 10$ | $9.98 \pm 0.03$ | correct |
| (c) $L$ | $3N-1 = 29$ | $9.97 \pm 0.03$ | correct |
| (b) $\ln V$ | $N-1 = 9$ (**wrong**) | $8.97 \pm 0.03$ | shifted by $-1$ |
| (a) $V$ | $N = 10$ (textbook $V^N$) | $10.97 \pm 0.04$ | gives $(N+1)k_BT/P$ |

The last row is not wrong — it is the textbook convention, consistent with its own $\Delta$ — but mixing it with the lecture's partition function would be. Also from these runs: $\langle (N-1)k_BT/V\rangle = 0.999$, $0.997$, $0.994$ for the three correct rules, versus $\langle Nk_BT/V\rangle \approx 1.11$ — confirming Section 18.2's identity for the lecture's convention.

### 19.3 Test 2: Lennard-Jones, $N = 108$, $T^* = 2.0$, $P^* = 2.0$, $r_c = 2.5$

2,000 equilibration + 6,000 production cycles (about 50 s in NumPy). Tuned: $\Delta r_m = 0.17$–$0.19$, $\Delta\ln V_{max} = 0.077$.

| Quantity | Truncated LJ + tail corrections | Cut-and-shifted LJ, no tail |
|---|---|---|
| displacement / volume acceptance | 0.51 / 0.42 | 0.50 / 0.44 |
| $\langle\rho\rangle$ | $0.6207 \pm 0.0018$ | $0.5870 \pm 0.0017$ |
| $\langle U/N\rangle$ | $-3.862 \pm 0.014$ | $-3.041 \pm 0.010$ |
| $\langle(N-1)k_BT/V + P_{ex}\rangle$ | $1.950 \pm 0.011$ | $\mathbf{1.984 \pm 0.015}$ |
| imposed $P$ | 2.0 | 2.0 |
| $\kappa_T$ (short run, rough) | 0.13 | 0.17 |

- **Shifted potential:** the pressure identity holds within about one error bar. ✓ (The two models are different potentials, so their densities and energies differ — that is expected.)
- **Truncated + tail:** a gap of about 0.05. This is the cutoff caveat of Section 18.2: the energy jumps at $r_c$, which the virial does not see during volume moves. The sampling itself is fine; it is the pressure *estimator* that needs an extra (impulsive) correction for truncated potentials.

---

## 20. Corrections and clarifications to the lecture notes

1. **$f$ is used twice:** $f(V)$ (volume scale) and $f(o)$, $f(n)$ (probabilities). Keep them apart.
2. **The boxed ratios are $\text{acc}(o\to n)/\text{acc}(n\to o)$, not the acceptance itself.** The acceptance is $\min\{1, \text{ratio}\}$.
3. **$V^{N-1}$ is a convention.** Most textbooks use $V^N$; the acceptance exponents then all shift up by one. Never mix the two.
4. **Unstated but essential for volume moves:** rescale all positions with the box; recompute $U$ over all pairs; include the volume dependence of the tail correction; reject $V_n \le 0$ (a) or $L_n \le 0$ (c); keep $r_c < L/2$.
5. **Each proposal must be symmetric in its own variable** — that is why (a), (b), (c) have different exponents. Using a $\ln V$ step with the (a) rule is wrong (Section 19.2, row 4).
6. **The move type must be chosen at random** with probabilities $N/(N+1)$ and $1/(N+1)$, not in a fixed sequence.
7. **Tuning $\Delta V_{max}$ (or $\Delta\ln V_{max}$, $\Delta L_{max}$)** is not mentioned; do it during equilibration only.
8. **Pressure check:** with the lecture's convention the exact identity is $P = \langle(N-1)k_BT/V + P_{ex}\rangle$ (Section 18.2).

---

## 21. Common mistakes

| Mistake | Consequence |
|---|---|
| Not restoring positions after a rejected volume move | atoms stay rescaled while $V$ goes back: garbage |
| Stepping in $\ln V$ with the exponent from (a) | biased towards small volumes by $\sim 1/N$ — no obvious alarm |
| Forgetting the tail correction changes with $V$ | wrong density |
| Volume moves in a fixed pattern | breaks detailed balance |
| Tuning step sizes during production | breaks detailed balance |
| Computing $(V_n/V_o)^N$ directly | overflow for large $N$; use logarithms |
| Reporting $N/\langle V\rangle$ as $\langle\rho\rangle$ | small systematic difference |
| $r_c > L/2$ after the box shrinks | minimum image fails |
| Reusing $\Delta U$ from only the "moved" atoms in a volume move | every pair changes; recompute all |

---

## 22. Practice problems with answers

**22.1 Reservoir counting.** System states with $E = 0, 1, 2$; reservoir with $\Omega_R(E) = 3^E$; $E_{tot} = 5$. Find $\Omega_{tot}$, the three probabilities, $\beta$ and $\sum_i e^{-\beta E_i}$.

> $\Omega_R = 243, 81, 27$; $\Omega_{tot} = 351$; $p = 0.692, 0.231, 0.077$; $\beta = \ln 3 = 1.099$; $\sum e^{-\beta E_i} = 1 + \tfrac13 + \tfrac19 = 1.444 = 351/243$. ✓

**22.2 Scaled volume.** $L = 5$, a scaled cube of side 0.02. Real and scaled volumes, and their ratio?

> Scaled $8\times10^{-6}$; real side 0.1, real volume $10^{-3}$; ratio $125 = L^3$. ✓

**22.3 $U(\mathbf s;L)$.** Two LJ atoms with scaled separation 0.2. Energy for $L = 5$ and $L = 6$?

> $L = 5$: $r = 1.0$, $u = 0$. $L = 6$: $r = 1.2$, $(1/1.2)^6 = 0.3349$, $(1/1.2)^{12} = 0.1122$, $u = 4(0.1122 - 0.3349) = -0.891$.

**22.4 Ideal gas.** $N = 100$, $T^* = 1.5$, $P^* = 0.5$. $\langle V\rangle$ with the lecture's convention? With the textbook's?

> $Nk_BT/P = 300$; $(N+1)k_BT/P = 303$.

**22.5 Volume move (a).** $N = 100$, $\beta = 1$, $P = 1$, $V_o = 200$, $V_n = 202$, $\Delta U = -0.3$. Acceptance?

> $99\ln 1.01 - (-0.3 + 2) = 0.985 - 1.7 = -0.715$; acc $= 0.489$. (Same move proposed in $\ln V$: 0.494; in $L$: 0.492.)

**22.6 Cycles.** $N = 500$, 10,000 cycles. Total moves? Average number of volume moves?

> $10{,}000\times 501 = 5{,}010{,}000$ moves; about 10,000 volume moves.

**22.7 Debugging.** An ideal-gas NPT test with $N = 10$, $k_BT/P = 1$ and $\ln V$ steps gives $\langle V\rangle = 9.0$. What is wrong?

> The $\ln V$ rule was used with exponent $N - 1$ instead of $N$ (Jacobian forgotten).

**22.8 Cutoff.** $N = 108$, $\rho = 0.6$, $r_c = 3.0$. Can volume moves work?

> $L = (108/0.6)^{1/3} = 5.65$, $L/2 = 2.82 < 3.0$: no — every move fails the cutoff check. Use $r_c \le 2.5$ or more atoms.

---

## 23. Summary

- **Statistical mechanics:** for a system exchanging energy and volume with a huge reservoir, $p(i,V) = e^{-\beta E_i}e^{-\beta PV}/\Delta$, from counting reservoir states ($\Omega_R/\Omega_{tot}$) and expanding $\ln\Omega_R$; the derivatives $\beta$ and $\beta P$ are the reservoir's temperature and pressure.
- **Partition function:** $\Delta = \sum_V e^{-\beta PV}Q(N,V,T)$; classically, with the lecture's $f(V) = V$,
  $$\Delta = \frac{1}{N!\,\Lambda^{3N}}\int_0^\infty e^{-\beta PV}\,V^{N-1}\int_0^1 e^{-\beta U(\mathbf s^N;L)}\,d\mathbf s^N\,dV$$
  $V^N$ from scaled coordinates, $1/V$ from removing the overall slide of a periodic system (equivalently, measuring volume in $d\ln V$). $G = -k_BT\ln\Delta$.
- **Density sampled:** $p(V,\mathbf s^N) \propto V^{N-1}e^{-\beta PV}e^{-\beta U}$.
- **Moves:** displacement, $\min\{1, e^{-\beta\Delta U}\}$; volume, with the exponent set by the stepping variable (Jacobian): $V$ → $N-1$, $\ln V$ → $N$, $L$ → $3N-1$; all with $e^{-\beta P\Delta V}e^{-\beta\Delta U}$.
- **Cycle:** $N + 1$ moves, type chosen at random with probabilities $N/(N+1)$ and $1/(N+1)$; one volume move costs about as much as $N$ displacements.
- **Checks:** ideal gas gives $\langle V\rangle = Nk_BT/P$ exactly; interacting systems satisfy $P = \langle(N-1)k_BT/V + P_{ex}\rangle$ (smooth potentials).
