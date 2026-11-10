# Lecture 8: Grand Canonical Monte Carlo — Complete Notes with Every Derivation

*Based on 8.pdf ("GCMC Simulation", 1 February 2022). This file combines and extends the earlier lecture-8 notes and the partition-function derivation: every formula in the lecture is derived here from first principles, step by step, with numerical checks. Numbers marked "simulation" come from a tested Python GCMC code (Part VI).*

---

## Contents

**Part 0. Orientation**
- 0.1 Notation
- 0.2 The lecture in one page

**Part I. The grand canonical ensemble from scratch**
- 1. Open systems: the physical picture
- 2. Probability of a microstate (counting the reservoir)
- 3. The grand partition function Ξ
- 4. Classical partition function for identical atoms (Λ and N!)
- 5. The lecture's formula, scaled coordinates, microstate probability

**Part II. Thermodynamics from Ξ**
- 6. Average N, fluctuations of N, and the grand potential
- 7. The ideal gas: every result in closed form
- 8. Chemical potential of a real fluid: βμ = ln(ρΛ³) + βμᵉˣ
- 9. Activity z: concept, physical meaning, derivation
- 10. Fugacity f and z = βf

**Part III. Why GCMC: adsorption**
- 11. Adsorption in a slit pore; why NVT is impractical
- 12. Analytical example: the Langmuir isotherm from Ξ

**Part IV. The Monte Carlo moves**
- 13. The detailed-balance recipe
- 14. Displacement
- 15. Insertion (and a slip in the lecture)
- 16. Deletion
- 17. Two bookkeeping conventions, one result
- 18. Reading the acceptance rules with numbers

**Part V. Implementation**
- 19. Insertion and deletion algorithms (the array trick)
- 20. The tail correction changes with N
- 21. The complete GCMC algorithm, line by line

**Part VI. Tests and working code**
- 22. Sanity checks (ideal gas, low density, Widom)
- 23. Python code and results

**Part VII. Wrap-up**
- 24. Practical issues
- 25. Corrections and clarifications to the lecture notes
- 26. Practice questions with answers
- 27. Summary: the whole chain

---

# Part 0. Orientation

## 0.1 Notation

| Symbol | Meaning |
|---|---|
| S, R | the system (simulation box) and the reservoir around it |
| $N$ | number of molecules in the box |
| $V$, $L$ | box volume and edge length, $V = L^3$ |
| $T$ | temperature |
| $k_B$ | Boltzmann's constant, $1.380649\times10^{-23}$ J/K |
| $\beta$ | $1/(k_BT)$ |
| $E$, $E_s$ | energy; energy of system microstate $s$ |
| $E_{tot}$, $N_{tot}$ | total energy and total molecules of system + reservoir (fixed) |
| $\Omega_R(E, N)$ | number of reservoir microstates with energy $E$ and $N$ molecules |
| $S$ | entropy, $S = k_B\ln\Omega$ (in Part I, "S" also labels the system; the meaning is clear from context) |
| $P$ | pressure |
| $\mu$ | chemical potential |
| $\mu^{\text{id}}$, $\mu^{\text{ex}}$ | ideal-gas and excess parts of $\mu$ |
| $F$ | Helmholtz free energy, $F = E - TS$ |
| $\Phi$ | grand potential, $\Phi = -k_BT\ln\Xi$ |
| $Q(N, V, T)$ | canonical partition function |
| $\Xi(\mu, V, T)$ | grand partition function |
| $H$, $K$, $U$ | total energy (Hamiltonian), kinetic energy, potential energy |
| $\mathbf r^N$, $\mathbf p^N$ | positions and momenta of all $N$ molecules |
| $m$ | mass of one molecule |
| $h$ | Planck's constant |
| $\Lambda$ | thermal de Broglie wavelength, $h/\sqrt{2\pi m k_BT}$ |
| $Z_N$ | configurational integral, $\int e^{-\beta U(\mathbf r^N)}\,d\mathbf r^N$ |
| $\mathbf s$ | scaled coordinates, $\mathbf s = \mathbf r/V^{1/3}$ |
| $\rho$ | number density $N/V$ |
| $z$ | activity, $e^{\beta\mu}/\Lambda^3$ |
| $f$ | fugacity |
| $\Delta U$ | potential-energy change of a move |
| $p(o)$, $p(n)$ | probability of the old / new microstate (the lecture writes $f$; here $f$ is reserved for fugacity) |
| $\alpha(o\to n)$, $\text{acc}(o\to n)$ | proposal and acceptance probabilities |
| $\Delta r_m$ | maximum displacement |

## 0.2 The lecture in one page

1. **Grand canonical ensemble:** an open box at fixed $\mu$, $V$, $T$; $N$ fluctuates.
2. **Partition function:** $\displaystyle\Xi = \sum_{N=0}^{\infty}\frac{e^{\beta\mu N}}{N!\,\Lambda^{3N}}\int e^{-\beta U}\,d\mathbf r^N = \sum_{N}\frac{V^N e^{\beta\mu N}}{N!\,\Lambda^{3N}}\int e^{-\beta U}\,d\mathbf s^N$.
3. **Application:** adsorption in a pore; simulate only the pore at the $\mu$ of the outside gas.
4. **Moves and acceptance rules:**
   - displacement: $\min\lbrace 1, e^{-\beta\Delta U}\rbrace$
   - insertion: $\min\lbrace 1, \frac{zV}{N+1}e^{-\beta\Delta U}\rbrace$, with $z = e^{\beta\mu}/\Lambda^3 = \rho e^{\beta\mu^{ex}} = \beta f$
   - deletion: $\min\lbrace 1, \frac{N}{zV}e^{-\beta\Delta U}\rbrace$
5. **Implementation:** insert at a uniformly random point; delete a random molecule by copying the last one into its slot; choose each move type with probability 1/3; cycle length $\max(10, N)$.

Everything below derives these results.

---

# Part I. The grand canonical ensemble from scratch

## 1. Open systems: the physical picture

The **grand canonical ensemble** describes a box that exchanges **molecules and energy** with a very large **reservoir**.

| Fixed by the reservoir | Fluctuates |
|---|---|
| chemical potential $\mu$ | number of molecules $N$ |
| volume $V$ (of the box) | energy |
| temperature $T$ | |

At equilibrium the box and the reservoir have the same $T$ (no net heat flow) and the same $\mu$ (no net flow of molecules). $\mu$ is the free-energy change on adding one molecule: a high $\mu$ in the reservoir pushes molecules into the box, a low $\mu$ pulls them out.

| Ensemble | Fixed | Monte Carlo moves |
|---|---|---|
| Canonical (NVT) | $N, V, T$ | displacement (rotation for molecules) |
| Grand canonical ($\mu VT$) | $\mu, V, T$ | displacement **+ insertion + deletion** |

## 2. Probability of a microstate

**Goal:** the probability $P(s)$ that the box is in one particular microstate $s$, which has $N$ molecules and energy $E_s$.

### Step 2.1 Build an isolated system

Box + reservoir together are isolated, so

$$E_s + E_R = E_{tot}, \qquad N + N_R = N_{tot}$$

where $E_R$, $N_R$ belong to the reservoir. The reservoir is huge: $E_s \ll E_{tot}$, $N \ll N_{tot}$.

### Step 2.2 The one postulate

In an isolated system, **every microstate with the given $E_{tot}$ and $N_{tot}$ is equally likely** (microcanonical ensemble).

### Step 2.3 Count the reservoir's arrangements

If the box is in state $s$, the reservoir must have energy $E_{tot} - E_s$ and $N_{tot} - N$ molecules, and can be in any of its $\Omega_R$ microstates with those values. All combined microstates are equally likely, so

$$P(s) \propto \Omega_R(E_{tot} - E_s,\; N_{tot} - N)$$

**The box prefers states that leave the reservoir the most ways to arrange itself.**

### Step 2.4 Work with $\ln\Omega_R$, not $\Omega_R$

We do not know the function $\Omega_R$, but we only need how it changes when the box takes a tiny energy $E_s$ and a few molecules $N$. $\Omega_R$ itself changes far too steeply to expand directly; its logarithm changes smoothly.

**Concrete reservoir.** For an ideal monatomic gas of $N_R$ atoms at fixed volume, $\Omega_R \propto E^{3N_R/2}$, so

$$\ln\Omega_R = \tfrac{3N_R}{2}\ln E + \text{const}$$

$$\frac{\partial\ln\Omega_R}{\partial E} = \frac{3N_R}{2E} = \frac{3N_R}{2\cdot\frac32 N_R k_BT} = \frac{1}{k_BT} = \beta$$

using $E_{tot} = \frac32 N_R k_BT$. Compare the two ways of expanding when the box takes energy $E_s$:

| $E_s$ | True ratio $\Omega_{new}/\Omega_{old}\approx e^{-\beta E_s}$ | Linear in $\Omega$: $1 - \beta E_s$ | Linear in $\ln\Omega$: $e^{-\beta E_s}$ |
|---|---|---|---|
| $1\,k_BT$ | 0.368 | 0 | 0.368 ✓ |
| $5\,k_BT$ | 0.0067 | −4 | 0.0067 ✓ |
| $10\,k_BT$ | 0.000045 | −9 | 0.000045 ✓ |

Taking one molecule's worth of thermal energy from $10^{23}$ molecules already cuts $\Omega_R$ by a factor $e$. Expanding $\Omega$ predicts negative numbers of states; expanding $\ln\Omega$ is accurate.

### Step 2.5 Taylor-expand to first order

**One variable:** $g(a - h) \approx g(a) - h\,g'(a)$ for small $h$. Check: $\ln 998 = 6.905753$; $\ln 1000 - 2/1000 = 6.905755$.

**Two variables:** $g(a - h, b - k) \approx g(a,b) - h\,\partial g/\partial a - k\,\partial g/\partial b$, each partial derivative being the slope in one direction with the other variable held fixed.

With $g = \ln\Omega_R$, $a = E_{tot}$, $h = E_s$, $b = N_{tot}$, $k = N$:

$$\ln\Omega_R(E_{tot} - E_s,\, N_{tot} - N)$$

$$\approx \ln\Omega_R(E_{tot}, N_{tot}) - E_s\frac{\partial\ln\Omega_R}{\partial E} - N\frac{\partial\ln\Omega_R}{\partial N}$$

**Why first order is enough.** The second-order energy term is $\tfrac12E_s^2\,\partial^2\ln\Omega_R/\partial E^2$. For the ideal-gas reservoir this is $E_s^2/(2\cdot\frac32N_R(k_BT)^2)$: for $E_s = 10\,k_BT$ and $N_R = 6\times10^{23}$, about $6\times10^{-23}$. In general it scales as 1/(heat capacity of the reservoir) and vanishes for a huge reservoir.

### Step 2.6 The slopes are $\beta$ and $-\beta\mu$

**(a) Rearrange the fundamental relation.** Start from $dE = T\,dS - P\,dV + \mu\,dN$. Move the last two terms and divide by $T$:

$$dS = \frac{1}{T}\,dE + \frac{P}{T}\,dV - \frac{\mu}{T}\,dN \qquad (\ast)$$

**(b) Read off the partial derivatives.** $S$ is a function $S(E, V, N)$, so calculus gives

$$dS = \left(\frac{\partial S}{\partial E}\right)_{V,N}dE + \left(\frac{\partial S}{\partial V}\right)_{E,N}dV$$

$$\qquad + \left(\frac{\partial S}{\partial N}\right)_{E,V}dN \qquad (\ast\ast)$$

Both describe the same $dS$ for **any** small changes. Change only $E$ ($dV = dN = 0$): $(\ast)$ gives $dS = dE/T$ and $(\ast\ast)$ gives $dS = (\partial S/\partial E)\,dE$, so $\partial S/\partial E = 1/T$. Change only $N$: $\partial S/\partial N = -\mu/T$. (Change only $V$: $\partial S/\partial V = P/T$.)

*Check of "reading off" on a simple function:* $f(x,y) = x^2y$ has $df = 2xy\,dx + x^2\,dy$; at $x=3$, $y=2$, $df = 12\,dx + 9\,dy$, and indeed $\partial f/\partial x = 12$.

**(c) Convert to $\ln\Omega$.** Boltzmann: $S = k_B\ln\Omega$, so $\ln\Omega = S/k_B$:

$$\boxed{\;\frac{\partial\ln\Omega_R}{\partial E} = \frac{1}{k_BT} = \beta\;}$$

$$\boxed{\;\frac{\partial\ln\Omega_R}{\partial N} = -\frac{\mu}{k_BT} = -\beta\mu\;}$$

Units check: $\ln\Omega$ is a pure number, so $\partial\ln\Omega/\partial E$ has units 1/energy (as $\beta$ does) and $\partial\ln\Omega/\partial N$ is a pure number (as $\beta\mu$ is). These are the **reservoir's** $T$ and $\mu$: this is how the reservoir imposes them on the box.

### Step 2.7 Substitute, exponentiate, drop the constant

$$\ln\Omega_R \approx C - \beta E_s + \beta\mu N, \qquad C \equiv \ln\Omega_R(E_{tot}, N_{tot})$$

Exponentiate both sides, using $e^{\ln x} = x$ and $e^{a+b+c} = e^ae^be^c$:

$$\Omega_R(E_{tot} - E_s,\, N_{tot} - N) \approx e^{C}\,e^{-\beta E_s}\,e^{\beta\mu N}$$

$e^{C}$ depends only on the fixed totals, so it is **the same for every state $s$** and goes into the proportionality. Factor out $-\beta$:

$$P(s) \propto e^{-\beta E_s}e^{\beta\mu N} = e^{-\beta(E_s - \mu N)}$$

- $e^{-\beta E_s}$: low energy favoured (as in NVT).
- $e^{\beta\mu N}$: each molecule in the box multiplies the weight by $e^{\beta\mu}$.

## 3. The grand partition function Ξ

### Step 3.1 Normalize

Write $P(s) = c\,w(s)$ with weight $w(s) = e^{-\beta(E_s - \mu N)}$. Probabilities over **all** states (every possible $N$) must add to 1:

$$c\sum_{\text{all } s} w(s) = 1 \;\Rightarrow\; c = \frac{1}{\Xi}$$

$$\boxed{\;\Xi = \sum_{\text{all } s} e^{-\beta(E_s - \mu N)}\;}$$

$$\boxed{\;P(s) = \frac{e^{-\beta(E_s - \mu N)}}{\Xi}\;}$$

Keeping $e^C$ would not have mattered: it would appear in every term of $\Xi$ and cancel.

**Numerical example** (energies in $k_BT$, $\beta\mu = -0.5$):

| State | $N$ | $\beta E_s$ | Weight $e^{-\beta E_s}e^{\beta\mu N}$ | $P = w/\Xi$ |
|---|---|---|---|---|
| empty | 0 | 0 | 1.0000 | 0.348 |
| one molecule, favourable spot | 1 | −1 | $e^{1}e^{-0.5} = 1.6487$ | 0.574 |
| one molecule, unfavourable spot | 1 | +1 | $e^{-1}e^{-0.5} = 0.2231$ | 0.078 |
| | | | $\Xi = 2.8718$ | sum = 1.000 |

### Step 3.2 Group the states by N

For fixed $N$, $e^{\beta\mu N}$ is common to all states, so it comes out of the inner sum; what remains is the canonical partition function $Q(N, V, T) = \sum_{s \text{ with } N} e^{-\beta E_s}$:

$$\boxed{\;\Xi(\mu, V, T) = \sum_{N=0}^{\infty} e^{\beta\mu N}\,Q(N, V, T)\;}$$

**The grand partition function is a weighted sum of canonical partition functions, one for each $N$.**

## 4. Classical partition function for identical atoms

### Step 4.1 Phase-space form

For atoms obeying classical mechanics, the sum over states becomes an integral over positions and momenta:

$$Q(N,V,T) = \frac{1}{N!\,h^{3N}}\int d\mathbf p^N\int d\mathbf r^N\,e^{-\beta H}$$

$$H = \sum_{i=1}^{N}\frac{p_i^2}{2m} + U(\mathbf r^N)$$

| Factor | Why it is there |
|---|---|
| $h^{3N}$ | quantum mechanics assigns one state to each phase-space cell of size $h$ per coordinate–momentum pair; dividing by $h^{3N}$ counts states and makes $Q$ dimensionless |
| $N!$ | identical atoms: swapping two gives the same physical state, but the integral over labelled positions counts it $N!$ times |

### Step 4.2 Separate kinetic and potential parts

$e^{-\beta H} = e^{-\beta K}e^{-\beta U}$, so

$$Q = \frac{1}{N!\,h^{3N}}\left[\int e^{-\beta\sum_i p_i^2/2m}\,d\mathbf p^N\right]\left[\int e^{-\beta U(\mathbf r^N)}\,d\mathbf r^N\right]$$

### Step 4.3 The Gaussian integral

$$I = \int_{-\infty}^{\infty}e^{-ax^2}\,dx$$

$$I^2 = \int\!\!\int e^{-a(x^2+y^2)}\,dx\,dy$$

$$= \int_0^{2\pi}\!\!\int_0^\infty e^{-ar^2}\,r\,dr\,d\theta = 2\pi\cdot\frac{1}{2a} = \frac{\pi}{a}$$

(polar coordinates: $x^2 + y^2 = r^2$, $dx\,dy = r\,dr\,d\theta$). So $I = \sqrt{\pi/a}$.

### Step 4.4 The momentum integral and Λ

Each of the $3N$ momentum components gives, with $a = \beta/(2m)$,

$$\int_{-\infty}^{\infty}e^{-\beta p_x^2/2m}\,dp_x = \sqrt{\frac{2\pi m}{\beta}} = \sqrt{2\pi m k_BT}$$

All $3N$ together: $(2\pi m k_BT)^{3N/2}$. Combined with $h^{3N}$:

$$\frac{(2\pi m k_BT)^{3N/2}}{h^{3N}} = \frac{1}{\Lambda^{3N}}, \qquad \Lambda = \frac{h}{\sqrt{2\pi m k_BT}}$$

$\Lambda$ is the **thermal de Broglie wavelength**. Argon: 0.25 Å at 120 K, 0.16 Å at 300 K, tiny compared with the ~3.6 Å spacing between atoms in liquid argon, so classical mechanics is fine.

$$\boxed{\;Q(N,V,T) = \frac{Z_N}{N!\,\Lambda^{3N}}, \qquad Z_N = \int e^{-\beta U(\mathbf r^N)}\,d\mathbf r^N\;}$$

## 5. The lecture's formula, scaled coordinates, microstate probability

### Step 5.1 Assemble

Put Section 4 into Section 3:

$$\boxed{\;\Xi(\mu,V,T) = \sum_{N=0}^{\infty}\frac{e^{\beta\mu N}}{N!\,\Lambda^{3N}}\int\!\!\cdots\!\!\int e^{-\beta U(\mathbf r^N)}\,d\mathbf r^N\;}$$

This is the lecture's first line.

### Step 5.2 Scaled coordinates

Define $\mathbf s = \mathbf r/V^{1/3} = \mathbf r/L$, so each component of $\mathbf s$ lies in $[0, 1)$ whatever the box size. Each molecule's volume element is $d^3\mathbf r = L^3\,d^3\mathbf s = V\,d^3\mathbf s$, so $d\mathbf r^N = V^N d\mathbf s^N$:

$$\Xi = \sum_{N=0}^{\infty}\frac{V^N e^{\beta\mu N}}{N!\,\Lambda^{3N}}\int\!\!\cdots\!\!\int e^{-\beta U}\,d\mathbf s^N$$

This is the lecture's second line. It makes the $V$-dependence explicit and keeps the integration region fixed (the unit cube), which is convenient whenever $N$ or $V$ changes.

### Step 5.3 Probability of a microstate

After integrating out the momenta, the probability (density) of $N$ molecules at scaled positions $\mathbf s^N$ is the term inside the sum divided by $\Xi$:

$$p(\mathbf s^N; N) = \frac{1}{\Xi}\,\frac{V^N e^{\beta\mu N}}{N!\,\Lambda^{3N}}\,e^{-\beta U(\mathbf s^N)}$$

In real coordinates the same probability has no $V^N$ (it is absorbed into $d\mathbf r^N$). $\Xi$ is never needed in a simulation: acceptance rules use ratios, in which it cancels.

---

# Part II. Thermodynamics from Ξ

## 6. Average N, fluctuations of N, and the grand potential

### 6.1 Probability of having N molecules

Summing $P(s)$ over all states with the same $N$:

$$P(N) = \frac{e^{\beta\mu N}\,Q(N,V,T)}{\Xi}$$

### 6.2 Average number of molecules

Differentiate a single term of $\Xi$ with respect to $\mu$: $\frac{\partial}{\partial\mu}e^{\beta\mu N} = \beta N\,e^{\beta\mu N}$. Therefore

$$\frac{\partial\Xi}{\partial\mu} = \sum_N \beta N\,e^{\beta\mu N}Q(N) = \beta\,\Xi\sum_N N\,P(N) = \beta\,\Xi\,\langle N\rangle$$

$$\boxed{\;\langle N\rangle = \frac{1}{\beta}\frac{\partial\ln\Xi}{\partial\mu}\;}$$

### 6.3 Fluctuations of N

Differentiating twice brings down $\beta^2N^2$, so $\frac{1}{\Xi}\frac{\partial^2\Xi}{\partial\mu^2} = \beta^2\langle N^2\rangle$. Using the quotient rule,

$$\frac{\partial^2\ln\Xi}{\partial\mu^2} = \frac{1}{\Xi}\frac{\partial^2\Xi}{\partial\mu^2} - \left(\frac{1}{\Xi}\frac{\partial\Xi}{\partial\mu}\right)^2$$

$$= \beta^2\langle N^2\rangle - \beta^2\langle N\rangle^2$$

$$\boxed{\;\langle N^2\rangle - \langle N\rangle^2 = \frac{1}{\beta^2}\frac{\partial^2\ln\Xi}{\partial\mu^2} = k_BT\,\frac{\partial\langle N\rangle}{\partial\mu}\;}$$

Large fluctuations mean $\langle N\rangle$ responds strongly to $\mu$ (high compressibility, e.g. near a critical point).

### 6.4 The grand potential Φ = −k_BT ln Ξ = −PV

**(a) Entropy of the ensemble.** Using the Gibbs entropy $S = -k_B\sum_s P(s)\ln P(s)$ and $\ln P(s) = -\beta(E_s - \mu N_s) - \ln\Xi$:

$$S = k_B\beta\left(\langle E\rangle - \mu\langle N\rangle\right) + k_B\ln\Xi$$

$$\Rightarrow\quad -k_BT\ln\Xi = \langle E\rangle - TS - \mu\langle N\rangle \equiv \Phi$$

$\Phi$ is the **grand potential**.

**(b) Euler relation.** Energy is extensive: scaling the system by a factor $\lambda$ scales $S$, $V$, $N$ and $E$ by $\lambda$, i.e. $E(\lambda S, \lambda V, \lambda N) = \lambda E(S, V, N)$. Differentiate with respect to $\lambda$ at $\lambda = 1$ and use $\partial E/\partial S = T$, $\partial E/\partial V = -P$, $\partial E/\partial N = \mu$:

$$E = TS - PV + \mu N$$

**(c) Combine:**

$$\Phi = E - TS - \mu N = -PV \qquad\Rightarrow\qquad \boxed{\;PV = k_BT\ln\Xi\;}$$

(In the same way, the canonical ensemble gives the Helmholtz free energy $F = E - TS = -k_BT\ln Q$, used in Section 8.)

## 7. The ideal gas: every result in closed form

No interactions: $U = 0$, so $Z_N = \int d\mathbf r^N = V^N$. Define $z = e^{\beta\mu}/\Lambda^3$ (Section 9). Then

$$\Xi_{\text{id}} = \sum_{N=0}^{\infty}\frac{e^{\beta\mu N}V^N}{N!\,\Lambda^{3N}} = \sum_{N=0}^{\infty}\frac{(zV)^N}{N!} = e^{zV}$$

using the exponential series $\sum_N x^N/N! = e^x$.

| Quantity | Derivation | Result |
|---|---|---|
| average $N$ | $\ln\Xi_{\text{id}} = zV$, $\partial z/\partial\mu = \beta z$, Section 6.2 | $\langle N\rangle = zV$ |
| density | $\rho = \langle N\rangle/V$ | $\rho = z$ |
| chemical potential | $e^{\beta\mu}/\Lambda^3 = \rho$, take ln | $\beta\mu^{\text{id}} = \ln(\rho\Lambda^3)$ |
| pressure | $PV = k_BT\ln\Xi = k_BT\,zV$ | $PV = \langle N\rangle k_BT$ |
| fluctuations | Section 6.3 | variance $= zV = \langle N\rangle$ |
| distribution | $P(N) = (zV)^N/(N!\,e^{zV})$ | Poisson |

**Meaning of $\rho\Lambda^3$:** the number of molecules in a cube of side $\Lambda$. For argon at 300 K and 1 bar, $\rho = P/k_BT = 2.41\times10^{25}$ m⁻³ and $\Lambda^3 = 4.1\times10^{-33}$ m³, so $\rho\Lambda^3 = 9.8\times10^{-8}$ and $\beta\mu = -16.1$. A classical gas has $\rho\Lambda^3 \ll 1$, which is why its chemical potential is strongly negative.

## 8. Chemical potential of a real fluid: βμ = ln(ρΛ³) + βμᵉˣ

### Step 8.1 μ is the free-energy cost of one molecule

With $F = E - TS$ and $dE = T\,dS - P\,dV + \mu\,dN$:

$$dF = dE - T\,dS - S\,dT = -S\,dT - P\,dV + \mu\,dN$$

$$\Rightarrow\quad \mu = \left(\frac{\partial F}{\partial N}\right)_{V,T} \approx F(N+1) - F(N)$$

(adding one molecule to a large system is a tiny change).

### Step 8.2 Use F = −k_BT ln Q

$$\beta\mu = -\ln\frac{Q(N+1)}{Q(N)}$$

### Step 8.3 Ratio of partition functions

With $Q = Z_N/(N!\,\Lambda^{3N})$:

| Piece | $Q(N+1)/Q(N)$ |
|---|---|
| configurational integrals | $Z_{N+1}/Z_N$ |
| identical molecules | $N!/(N+1)! = 1/(N+1)$ |
| momenta | $\Lambda^{3N}/\Lambda^{3(N+1)} = 1/\Lambda^3$ |

$$\frac{Q(N+1)}{Q(N)} = \frac{Z_{N+1}}{Z_N}\cdot\frac{1}{(N+1)\Lambda^3}$$

### Step 8.4 Ideal gas gives the first term

$Z_N = V^N$, so $Z_{N+1}/Z_N = V$:

$$\beta\mu^{\text{id}} = \ln\frac{(N+1)\Lambda^3}{V} \approx \ln(\rho\Lambda^3) \qquad (N + 1 \approx N)$$

### Step 8.5 Interactions give the second term

Add molecule $N+1$ at position $\mathbf r_{new}$. Its interaction energy with the others is $\Delta U$, so $U(N+1) = U(N) + \Delta U$ and

$$\frac{Z_{N+1}}{Z_N} = \int d\mathbf r_{new}\;\frac{\int d\mathbf r^N e^{-\beta U(N)}e^{-\beta\Delta U}}{\int d\mathbf r^N e^{-\beta U(N)}}$$

$$= \int d\mathbf r_{new}\,\left\langle e^{-\beta\Delta U}\right\rangle_N = V\left\langle e^{-\beta\Delta U}\right\rangle_N$$

The fraction is the Boltzmann-weighted average, over configurations of the $N$-molecule fluid, of the Boltzmann factor felt by a "ghost" molecule at $\mathbf r_{new}$. A uniform fluid looks the same everywhere, so the average does not depend on $\mathbf r_{new}$ and the remaining integral gives $V$.

### Step 8.6 Result

$$\beta\mu = -\ln\frac{V\langle e^{-\beta\Delta U}\rangle}{(N+1)\Lambda^3}$$

$$= \ln\frac{(N+1)\Lambda^3}{V} - \ln\left\langle e^{-\beta\Delta U}\right\rangle$$

With $(N+1)/V \approx \rho$:

$$\boxed{\;\beta\mu = \ln(\rho\Lambda^3) + \beta\mu^{\text{ex}}\;}$$

$$\boxed{\;\beta\mu^{\text{ex}} = -\ln\left\langle e^{-\beta\Delta U}\right\rangle\;}$$

The second formula is **Widom's test-particle method**: $\mu^{\text{ex}}$ can be measured in an NVT simulation by trying random ghost insertions and averaging $e^{-\beta\Delta U}$.

| Term | Meaning |
|---|---|
| $\ln(\rho\Lambda^3)$ | cost of fitting one more molecule among the others as if they did not interact |
| $\beta\mu^{\text{ex}}$ | extra cost (or gain) from interactions, averaged over every possible insertion point |

| Situation | $\langle e^{-\beta\Delta U}\rangle$ | $\beta\mu^{\text{ex}}$ | Effect |
|---|---|---|---|
| ideal gas | 1 | 0 | none |
| every spot attractive, $\Delta U = -1\,k_BT$ | 2.72 | −1 | easier to insert |
| half the spots overlap, half neutral | 0.5 | +0.69 | harder to insert (dense fluid) |
| Lennard-Jones, $T^{\ast} = 2$, $\rho^{\ast} = 0.135$ (simulation) | 1.35 | −0.30 | attraction wins |

### Step 8.7 Low-density limit: βμᵉˣ ≈ 2B₂ρ

**(a) Gibbs–Duhem.** Differentiate the Euler relation $E = TS - PV + \mu N$ and subtract $dE = T\,dS - P\,dV + \mu\,dN$: $0 = S\,dT - V\,dP + N\,d\mu$. At constant $T$: $d\mu = dP/\rho$.

**(b) Virial equation of state.** At low density, $\beta P = \rho + B_2\rho^2$, where the **second virial coefficient** is

$$B_2(T) = -2\pi\int_0^\infty\left(e^{-\beta u(r)} - 1\right)r^2\,dr$$

and $u(r)$ is the pair potential. Then $\beta\,dP = (1 + 2B_2\rho)\,d\rho$.

**(c) Integrate.** $\beta\,d\mu = \beta\,dP/\rho = (1/\rho + 2B_2)\,d\rho$, so $\beta\mu = \ln\rho + 2B_2\rho + \text{const}(T)$. The $\ln\rho + \text{const}$ part is the ideal term $\ln(\rho\Lambda^3)$, leaving

$$\beta\mu^{\text{ex}} \approx 2B_2\rho$$

For Lennard-Jones at $T^{\ast} = 2$, $B_2^{\ast} = -1.314$ (negative: attraction dominates).

## 9. Activity z: concept, physical meaning and derivation

### 9.1 The problem activity solves

The reservoir's "push" on molecules is measured by its chemical potential $\mu$. But $\mu$ is awkward to work with:

- it is an energy, usually large and negative (argon gas at 300 K and 1 bar: $\beta\mu = -16.1$, Section 7);
- in every formula it comes together with the constant $\Lambda$.

The **activity** $z$ repackages $\mu$ into something that can be pictured: a **concentration**, in molecules per volume.

### 9.2 The idea: an effective concentration

$$z = \rho\times\gamma$$

- $\rho$ is the **actual** number density.
- $\gamma$ is a **correction factor for interactions**. In GCMC notation it is $\gamma = e^{\beta\mu^{\text{ex}}}$ (Section 9.6).

This is the same idea as the activity of a component in a liquid mixture in chemical-engineering thermodynamics, $a = \gamma x$ (activity coefficient times mole fraction): the concentration a substance **behaves as if** it had, once interactions are accounted for.

| Molecules … | $\gamma$ | Activity vs density | Meaning |
|---|---|---|---|
| do not interact (ideal gas) | 1 | $z = \rho$ | behaves exactly as crowded as it is |
| attract each other | below 1 | $z < \rho$ | each molecule is held by its neighbours and is less eager to leave |
| repel each other (very dense) | above 1 | $z > \rho$ | molecules are pushed out; the fluid behaves more crowded than it is |

### 9.3 Is activity a property of the system? Density versus activity

Yes. Activity is a property of the system's **state** (an intensive property, like $T$, $P$ and $\mu$). But it is a different kind of property from density:

- **density** says **how many molecules are here**;
- **activity** says **how strongly the molecules here tend to leave** (their *escaping tendency*).

The pair behaves exactly like energy content and temperature:

| | "How much is here" | "How strongly it tends to leave" |
|---|---|---|
| heat | energy content | **temperature** |
| molecules | density $\rho$ | **activity** $z$ (equivalently $\mu$) |

When a hot iron block is put in water, heat flows until the **temperatures** are equal; the energy contents are not equal (a kilogram of water stores far more heat than a kilogram of iron at the same temperature). In the same way, when two regions can exchange molecules, molecules flow from **high activity to low activity** until the **activities** are equal; the densities are not equal.

**This is what adsorption is.** A pore in contact with a gas at 300 K and 1 bar:

| | Gas outside | Inside an attractive pore |
|---|---|---|
| activity $z$ | 0.024 nm⁻³ | 0.024 nm⁻³ (equal: no net flow) |
| correction factor $\gamma$ | about 1 (nearly ideal) | 0.01 (strong wall attraction, illustrative value) |
| density $\rho = z/\gamma$ | 0.024 nm⁻³ | 2.4 nm⁻³ (100 times more) |

Same activity, a hundred times the density: the adsorbed molecules "feel" as crowded as the gas outside, but there are far more of them.

### 9.4 A precise physical definition, and how to "measure" activity

> **The activity of a system is the density that an ideal gas would need to have in order to be in equilibrium with that system.**

Reason: an ideal gas has no interactions, so its activity is its density ($z = \rho$). Bring any system into equilibrium with an ideal gas and their activities become equal, so the system's activity equals that gas's density.

This gives a way to "measure" activity, just as a thermometer measures temperature: let a dilute gas reach equilibrium with the system and read off the gas density,

$$z_{\text{system}} = \rho_{\text{gas}} = \frac{P}{k_BT}$$

(for a non-ideal gas, replace $P$ by the fugacity $f$, Section 10). Nobody can put a pressure gauge inside a nanometre-sized pore, but its activity is known: it equals the activity of the gas it is in contact with. **In GCMC you set $z$, the reservoir's property, and the simulation finds the density the box must have to match it.**

### 9.5 Molecular meaning: density divided by a "welcome factor"

From Section 8.6, $\beta\mu^{\text{ex}} = -\ln\langle e^{-\beta\Delta U}\rangle$, so $e^{\beta\mu^{\text{ex}}} = 1/\langle e^{-\beta\Delta U}\rangle$ and

$$z = \rho\,e^{\beta\mu^{\text{ex}}} = \frac{\rho}{\left\langle e^{-\beta\Delta U}\right\rangle}$$

Here $\Delta U$ is the energy change when one extra molecule is placed at a random point, and $\langle e^{-\beta\Delta U}\rangle$ (averaged over many random placements) is a **"welcome factor"**: how comfortable an extra molecule would be in this fluid. So

$$\text{escaping tendency} = \frac{\text{how many molecules there are}}{\text{how comfortable each one is}}$$

| Situation | $\rho$ | Welcome factor | Activity $z$ | Meaning |
|---|---|---|---|---|
| ideal gas | 0.1 | 1 | 0.1 | as crowded as it is |
| Lennard-Jones, $T^{\ast} = 2$ (simulation) | 0.134 | 1.35 (attraction) | 0.0998 | comfortable molecules escape less than their number suggests |
| dense fluid with overlaps (illustrative) | 0.8 | 0.5 (repulsion) | 1.6 | uncomfortable molecules escape more than their number suggests |

Many comfortable molecules can have the **same** activity as a few uncomfortable ones. That is why a dense adsorbed layer and a dilute gas can coexist.

### 9.6 The lecture's three lines, derived step by step

The lecture writes:

$$\beta\mu = \ln(\rho\Lambda^3) + \beta\mu^{\text{ex}}$$

$$\frac{e^{\beta\mu}}{\Lambda^3} = \rho\,e^{\beta\mu^{\text{ex}}} = z$$

$$\text{acc}(o\to n) = \min\left\lbrace 1,\; \frac{zV}{N+1}\,e^{-\beta\Delta U}\right\rbrace$$

**Line 1: chemical potential = ideal part + interaction part.**

*Ideal part.* For an ideal gas, Section 7 gives $\langle N\rangle = V e^{\beta\mu}/\Lambda^3$. Divide by $V$:

$$\rho = \frac{e^{\beta\mu}}{\Lambda^3}$$

Multiply by $\Lambda^3$ and take the natural logarithm:

$$\beta\mu^{\text{id}} = \ln(\rho\Lambda^3)$$

*Interaction part.* A real fluid's $\mu$ differs from this. The difference is, by definition, the **excess chemical potential**:

$$\mu^{\text{ex}} = \mu - \mu^{\text{id}}$$

Multiply by $\beta$ and rearrange:

$$\beta\mu = \ln(\rho\Lambda^3) + \beta\mu^{\text{ex}}$$

which is line 1. (Section 8 derives the same result from the free energy and shows that $\beta\mu^{\text{ex}} = -\ln\langle e^{-\beta\Delta U}\rangle$, the Widom formula.)

**Line 2: the activity.**

*Step 1.* Exponentiate both sides of line 1, using $e^{a+b} = e^a e^b$ and $e^{\ln x} = x$:

$$e^{\beta\mu} = e^{\ln(\rho\Lambda^3)}\,e^{\beta\mu^{\text{ex}}} = \rho\Lambda^3\,e^{\beta\mu^{\text{ex}}}$$

*Step 2.* Divide both sides by $\Lambda^3$:

$$\frac{e^{\beta\mu}}{\Lambda^3} = \rho\,e^{\beta\mu^{\text{ex}}}$$

*Step 3.* Name this combination the activity:

$$\boxed{\;z \equiv \frac{e^{\beta\mu}}{\Lambda^3} = \rho\,e^{\beta\mu^{\text{ex}}}\;}$$

Comparing with Section 9.2, the correction factor is $\gamma = e^{\beta\mu^{\text{ex}}}$: equal to 1 for an ideal gas ($\mu^{\text{ex}} = 0$), below 1 for attraction ($\mu^{\text{ex}} < 0$), above 1 for repulsion ($\mu^{\text{ex}} > 0$).

**Line 3: the insertion rule in terms of z.**

*Step 1.* Start from the insertion rule derived from detailed balance (Section 15):

$$\text{acc}(o\to n) = \min\left\lbrace 1,\; \frac{V\,e^{\beta\mu}}{(N+1)\,\Lambda^3}\,e^{-\beta\Delta U}\right\rbrace$$

(The lecture's line just before this one has $e^{\beta\mu N}$; that is a slip, see Section 15.1.)

*Step 2.* Regroup the fraction so that $\mu$ and $\Lambda$ sit together:

$$\frac{V\,e^{\beta\mu}}{(N+1)\,\Lambda^3} = \frac{V}{N+1}\times\frac{e^{\beta\mu}}{\Lambda^3}$$

*Step 3.* Recognise the second factor as $z$ (line 2):

$$\frac{V}{N+1}\times z = \frac{zV}{N+1}$$

*Step 4.* Substitute:

$$\boxed{\;\text{acc}(o\to n) = \min\left\lbrace 1,\; \frac{zV}{N+1}\,e^{-\beta\Delta U}\right\rbrace\;}$$

**Why the lecture makes this change.** $\mu$ and $\Lambda$ only ever appear together as $e^{\beta\mu}/\Lambda^3$, so specifying $z$ means you never need the molecular mass or Planck's constant separately. And $zV$ has a direct meaning: the number of molecules an ideal gas with this activity would put in the box, i.e. the number the reservoir "wants". In reduced units, $z^{\ast} = z\sigma^3$.

### 9.7 Numerical check (simulation)

Lennard-Jones, $T^{\ast} = 2$: GCMC at $z^{\ast} = 0.1$ gave $\rho^{\ast} = 0.1343$, about 34% denser than an ideal gas with the same activity. An independent NVT Widom run at $\rho^{\ast} = 0.1348$ gave $\beta\mu^{\text{ex}} = -0.30$ (welcome factor 1.35). Check:

$$z^{\ast} = 0.1348\times e^{-0.30} = 0.0998 \approx 0.1 \;\checkmark$$

### 9.8 Summary: density versus activity

| | Density $\rho$ | Activity $z$ |
|---|---|---|
| property of the system? | yes | yes (of its state, like $T$, $P$, $\mu$) |
| what it says | how many molecules per volume | how strongly molecules tend to leave |
| units | molecules per volume | molecules per volume (the density of the ideal gas it would be in equilibrium with) |
| equal across connected regions at equilibrium? | **no** | **yes** |
| analogy | energy content | temperature |
| relation | $\rho = z\,\langle e^{-\beta\Delta U}\rangle$ | $z = \rho\,e^{\beta\mu^{\text{ex}}} = e^{\beta\mu}/\Lambda^3$ |

**In one sentence:** density is how full a system is; activity is how hard it pushes molecules out. GCMC fixes the push ($z$, set by the reservoir) and lets the simulation find how full the box must be to push back equally.

## 10. Fugacity f and z = βf

Experiments specify the **pressure** of the outside gas, not $\mu$. Fugacity is the bridge.

### Step 10.1 Definition

$$k_BT\ln\left(\frac{f}{P}\right) = \mu(T, P) - \mu^{\text{id}}(T, P)$$

where $\mu^{\text{id}}(T, P)$ is the chemical potential an ideal gas would have at the **same $T$ and $P$**. For an ideal gas $f = P$. The ratio $\phi = f/P$ is the **fugacity coefficient**.

### Step 10.2 Derivation of z = βf

1. Ideal gas at pressure $P$: $\rho = P/(k_BT)$, so $\beta\mu^{\text{id}}(T,P) = \ln(\rho\Lambda^3) = \ln\!\left(\dfrac{P\Lambda^3}{k_BT}\right)$.
2. Divide the definition by $k_BT$: $\ln(f/P) = \beta\mu - \ln\!\left(\dfrac{P\Lambda^3}{k_BT}\right)$.
3. Exponentiate: $\dfrac{f}{P} = e^{\beta\mu}\,\dfrac{k_BT}{P\Lambda^3}$.
4. Cancel $P$ and rearrange:

$$\boxed{\;z = \frac{e^{\beta\mu}}{\Lambda^3} = \frac{f}{k_BT} = \beta f\;}$$

### Step 10.3 Relation between fugacity and excess chemical potential

Combine $z = \rho e^{\beta\mu^{\text{ex}}}$ and $z = \beta f$, and define the compressibility factor $\mathcal Z = \beta P/\rho$:

$$\frac{f}{P} = \frac{\rho\,e^{\beta\mu^{\text{ex}}}}{\beta P} \quad\Rightarrow\quad \ln\phi = \beta\mu^{\text{ex}} - \ln\mathcal Z$$

The two "corrections" differ only because $\mu^{\text{ex}}$ compares with an ideal gas at the same **density**, while $\phi$ compares with one at the same **pressure**.

### Step 10.4 Numbers

Methane at 300 K and 1 bar ($f \approx P$): $z = P/(k_BT) = 10^5/(1.381\times10^{-23}\times300) = 2.41\times10^{25}$ m⁻³ $= 0.0241$ nm⁻³. For a $3\times3\times3$ nm box, $zV = 0.65$; at 10 bar, $zV = 6.5$. In practice $f$ is obtained from an equation of state of the bulk gas (e.g. Peng–Robinson) at the experimental $T$ and $P$.

### Step 10.5 Three equivalent ways to specify the reservoir

| Specify | Insertion | Deletion |
|---|---|---|
| $\mu$ | $\min\lbrace 1, \frac{Ve^{\beta\mu}}{(N+1)\Lambda^3}e^{-\beta\Delta U}\rbrace$ | $\min\lbrace 1, \frac{N\Lambda^3}{Ve^{\beta\mu}}e^{-\beta\Delta U}\rbrace$ |
| $z$ | $\min\lbrace 1, \frac{zV}{N+1}e^{-\beta\Delta U}\rbrace$ | $\min\lbrace 1, \frac{N}{zV}e^{-\beta\Delta U}\rbrace$ |
| $f$ | $\min\lbrace 1, \frac{\beta fV}{N+1}e^{-\beta\Delta U}\rbrace$ | $\min\lbrace 1, \frac{N}{\beta fV}e^{-\beta\Delta U}\rbrace$ |

(These rules are derived in Part IV.) This is why the lecture's algorithm begins "specify $\mu$ (or $f$, $z$)".

---

# Part III. Why GCMC: adsorption

## 11. Adsorption in a slit pore; why NVT is impractical

**Adsorption:** gas molecules accumulate inside the pores of a solid (activated carbon, zeolites, MOFs). The lecture considers a **slit pore**: two parallel solid walls a distance $d$ apart. Question: how many molecules sit in the pore when the outside gas is at a given pressure and temperature?

### 11.1 NVT approach (left sketch in the lecture)

Place the adsorbent particle, with its pore, inside a large box of gas and run NVT. Two difficulties:

1. **The entire adsorbent must be simulated.** Real adsorbent particles are macroscopic: far too many atoms.
2. **The external gas must be simulated.** Most of the computer time goes into gas far from the pore, which is not of interest.

### 11.2 GCMC approach (right sketch)

Simulate **only the pore**. The pore is in equilibrium with the outside gas exactly when both have the same $\mu$ and $T$ (Section 1), so the outside gas can be replaced by a single number, its $\mu$ (or $z$, or $f$). Molecules are inserted into and deleted from the pore until $N$ fluctuates around the equilibrium value. **This is computationally much cheaper.**

### 11.3 Output: the adsorption isotherm

One GCMC run per outside pressure gives $\langle N\rangle$ in the pore. Plotting $\langle N\rangle$ against $P$ (or $f$) at fixed $T$ gives the **adsorption isotherm**, directly comparable with experiment. The walls act on the molecules through a wall–molecule potential added to $U$; attractive walls make insertions easier ($\Delta U < 0$), so the pore holds many more molecules than the same volume of bulk gas.

## 12. Analytical example: the Langmuir isotherm from Ξ

The machinery of Part I solved exactly for a simple model.

### 12.1 One site

An adsorption site is either empty or holds one molecule with binding energy $\varepsilon_b > 0$:

| State | $N$ | Energy | Weight $e^{-\beta(E-\mu N)}$ |
|---|---|---|---|
| empty | 0 | 0 | 1 |
| occupied | 1 | $-\varepsilon_b$ | $x \equiv e^{\beta(\mu + \varepsilon_b)}$ |

$$\Xi_{\text{site}} = 1 + x, \qquad \theta = \frac{x}{1 + x}$$

where $\theta$ is the fractional occupancy.

### 12.2 Many independent sites

Two sites have four states, (0,0), (1,0), (0,1), (1,1), with weights $1, x, x, x^2$: $\Xi = (1 + x)^2$. For $M$ independent sites, $\Xi = (1 + x)^M$, and the average occupancy per site is still $\theta = x/(1+x)$ (Section 6.2 applied to $\ln\Xi = M\ln(1 + x)$).

### 12.3 Connect to the gas pressure

For an ideal gas at pressure $P$: $z = \rho = \beta P$, so $e^{\beta\mu} = z\Lambda^3 = \beta P\Lambda^3$ and $x = KP$ with $K = \beta\Lambda^3e^{\beta\varepsilon_b}$:

$$\boxed{\;\theta = \frac{KP}{1 + KP}\;}$$

the **Langmuir isotherm**. The site is half full at $P = 1/K$; stronger binding means larger $K$, so filling happens at lower pressure. GCMC does the same job for a continuous, interacting pore where no closed form exists.

---

# Part IV. The Monte Carlo moves

## 13. The detailed-balance recipe

Every move is accepted or rejected so that detailed balance holds with respect to the target probability $p$:

$$p(o)\,\alpha(o\to n)\,\text{acc}(o\to n) = p(n)\,\alpha(n\to o)\,\text{acc}(n\to o)$$

so that, with the Metropolis choice,

$$\text{acc}(o\to n) = \min\left\lbrace 1,\; \frac{p(n)\,\alpha(n\to o)}{p(o)\,\alpha(o\to n)}\right\rbrace$$

Here the target is the grand canonical probability of Section 5.3:

$$p(\mathbf s^N; N) \propto \frac{V^N e^{\beta\mu N}}{N!\,\Lambda^{3N}}\,e^{-\beta U(\mathbf s^N)}$$

Recipe for each move: (1) write $\alpha$ forward and backward, (2) write $p$ for old and new, (3) take the ratio.

## 14. Displacement

**Move:** choose "displace" (probability 1/3 in the lecture's algorithm), pick one of the $N$ molecules, shift it uniformly within $\pm\Delta r_m$ in each direction.

- $\alpha(o\to n) = \frac13\cdot\frac1N\cdot\frac{1}{(2\Delta r_m)^3} = \alpha(n\to o)$ (the reverse picks the same molecule and the opposite shift).
- $N$ does not change, so every prefactor of $p$ cancels: $p(n)/p(o) = e^{-\beta\Delta U}$.

$$\boxed{\;\text{acc}_{\text{disp}} = \min\left\lbrace 1,\; e^{-\beta\Delta U}\right\rbrace \;}$$

the same as in NVT ("as derived earlier" in the lecture).

## 15. Insertion (N → N + 1)

**Probabilities of old and new states:**

$$p(o) \propto \frac{V^N e^{\beta\mu N}}{N!\,\Lambda^{3N}}\,e^{-\beta U(o)}$$

$$p(n) \propto \frac{V^{N+1} e^{\beta\mu(N+1)}}{(N+1)!\,\Lambda^{3(N+1)}}\,e^{-\beta U(n)}$$

**Proposal probabilities.** If an insertion is attempted as often as a deletion, $\alpha(o\to n) = \alpha(n\to o)$ (Section 17 explains the bookkeeping).

**Ratio, factor by factor:**

| Factor | New / old | Result |
|---|---|---|
| volume | $V^{N+1}/V^N$ | $V$ |
| chemical potential | $e^{\beta\mu(N+1)}/e^{\beta\mu N}$ | $e^{\beta\mu}$ |
| identical molecules | $N!/(N+1)!$ | $1/(N+1)$ |
| momenta | $\Lambda^{3N}/\Lambda^{3(N+1)}$ | $1/\Lambda^3$ |
| energy | $e^{-\beta U(n)}/e^{-\beta U(o)}$ | $e^{-\beta\Delta U}$, $\Delta U = U(n) - U(o)$ |

$$\frac{\text{acc}(o\to n)}{\text{acc}(n\to o)} = \frac{V e^{\beta\mu}}{(N+1)\Lambda^3}\,e^{-\beta\Delta U}$$

$$\boxed{\;\text{acc}_{\text{ins}} = \min\left\lbrace 1,\; \frac{V e^{\beta\mu}}{(N+1)\Lambda^3}\,e^{-\beta\Delta U}\right\rbrace\;}$$

Using $z = e^{\beta\mu}/\Lambda^3 = \beta f$ (Sections 9 and 10), the same rule reads

$$\text{acc}_{\text{ins}} = \min\left\lbrace 1,\; \frac{zV}{N+1}\,e^{-\beta\Delta U}\right\rbrace$$

$$\text{acc}_{\text{ins}} = \min\left\lbrace 1,\; \frac{\beta fV}{N+1}\,e^{-\beta\Delta U}\right\rbrace$$

$\Delta U$ is the interaction energy of the new molecule with all the others, plus the change in the tail correction (Section 20).

**Where every factor comes from:**

| Factor | Origin |
|---|---|
| $e^{\beta\mu}$ | the reservoir's chemical potential (Section 2.6) |
| $1/\Lambda^3$ | the new molecule's momentum integral (Section 4.4) |
| $1/(N+1)$ | $N!/(N+1)!$, identical molecules (Section 4.1) |
| $V$ | the new molecule's position, in scaled coordinates (Section 5.2) |
| $e^{-\beta\Delta U}$ | the configurational Boltzmann factor |

### 15.1 A slip in the lecture: the factor is exp(βμ), not exp(βμN)

One line of the lecture writes the insertion rule, **incorrectly**, as

$$\text{acc}(o\to n) = \min\left\lbrace 1,\; \frac{V\,e^{\beta\mu N}}{(N+1)\,\Lambda^3}\,e^{-\beta\Delta U}\right\rbrace$$

The correct factor is $e^{\beta\mu}$.

**Derivation.** The chemical-potential part of $p(n)/p(o)$, for $N+1$ versus $N$ molecules, is

$$\frac{e^{\beta\mu(N+1)}}{e^{\beta\mu N}} = e^{\beta\mu(N+1) - \beta\mu N} = e^{\beta\mu}$$

Only **one** molecule is added, so only **one** factor $e^{\beta\mu}$ survives. The lecture's preceding line (the ratio of acceptances) already has $e^{\beta\mu}$; the $N$ slipped in when copying.

**Three checks that $e^{\beta\mu N}$ cannot be right.**

1. *The next step would fail.* The lecture next replaces $e^{\beta\mu}/\Lambda^3$ by $z$ to get $zV/(N+1)$. That substitution only works if the formula contains $e^{\beta\mu}$.
2. *An empty box would ignore the reservoir.* For the first insertion, $N = 0$ and $e^{\beta\mu N} = e^0 = 1$: the chance of filling an empty box would not depend on the reservoir's chemical potential, although a high-$\mu$ reservoir must push molecules in far more readily than a low-$\mu$ one.
3. *Insertions would become effectively impossible.* For a gas $\beta\mu$ is negative, e.g. $-5$. With $N = 10$, $e^{\beta\mu N} = e^{-50} \approx 2\times10^{-22}$, so insertions would essentially never be accepted once a few molecules were present. The correct $e^{\beta\mu} = e^{-5}$ is the same for every insertion; the $N$-dependence enters only through the gentle factor $1/(N+1)$.

**Rule of thumb for catching such slips.** Each insertion or deletion changes $N$ by exactly one, so every factor in the acceptance ratio must be a **one-molecule** factor: one $e^{\beta\mu}$, one $1/\Lambda^3$, one $V$, one $1/(N+1)$. Anything raised to the power $N$ is a mistake.

## 16. Deletion (N → N − 1)

$$p(o) \propto \frac{V^N e^{\beta\mu N}}{N!\,\Lambda^{3N}}\,e^{-\beta U(o)}$$

$$p(n) \propto \frac{V^{N-1} e^{\beta\mu(N-1)}}{(N-1)!\,\Lambda^{3(N-1)}}\,e^{-\beta U(n)}$$

| Factor | New / old | Result |
|---|---|---|
| volume | $V^{N-1}/V^N$ | $1/V$ |
| chemical potential | $e^{\beta\mu(N-1)}/e^{\beta\mu N}$ | $e^{-\beta\mu}$ |
| identical molecules | $N!/(N-1)!$ | $N$ |
| momenta | $\Lambda^{3N}/\Lambda^{3(N-1)}$ | $\Lambda^3$ |
| energy | | $e^{-\beta\Delta U}$ |

$$\frac{\text{acc}(o\to n)}{\text{acc}(n\to o)} = \frac{N\Lambda^3}{Ve^{\beta\mu}}\,e^{-\beta\Delta U}$$

$$= \frac{N}{zV}\,e^{-\beta\Delta U} = \frac{N}{\beta fV}\,e^{-\beta\Delta U}$$

$$\boxed{\;\text{acc}_{\text{del}} = \min\left\lbrace 1,\; \frac{N}{zV}\,e^{-\beta\Delta U}\right\rbrace\;}$$

$$\text{acc}_{\text{del}} = \min\left\lbrace 1,\; \frac{N}{\beta fV}\,e^{-\beta\Delta U}\right\rbrace$$

$$\Delta U = U(N-1) - U(N)$$

Removing molecule $j$ removes all its interactions, so $\Delta U = -(\text{energy of } j \text{ with all others})$. Deleting a strongly bound molecule ($\Delta U > 0$) is hard; deleting one in an unfavourable spot is easy.

**Consistency check:** deleting a molecule from $N$ molecules and inserting it back into $N-1$ molecules connect the same pair of states in opposite directions. The insertion ratio (starting from $N-1$ molecules) is $\frac{zV}{N}e^{-\beta\Delta U_{\text{ins}}}$; the deletion ratio is $\frac{N}{zV}e^{-\beta\Delta U_{\text{del}}}$ with $\Delta U_{\text{del}} = -\Delta U_{\text{ins}}$. Their product is exactly 1: they are reciprocals, as detailed balance requires.

## 17. Two bookkeeping conventions, one result

The factor $V/(N+1)$ can be placed in two different spots; both give the same rule and the same code.

| | Lecture's convention | Alternative convention |
|---|---|---|
| probability $p$ | includes $1/N!$ (labelled molecules) and $V^N$ (scaled coordinates) | no $1/N!$ (identical molecules), real coordinates |
| insertion proposal | "insert" chosen as often as "delete"; position uniform in the unit cube | random point in the box: $\alpha = 1/V$ |
| deletion proposal | same | pick one of the $N+1$ molecules: $\alpha = 1/(N+1)$ |
| where $V/(N+1)$ comes from | the ratio $p(n)/p(o)$ | the ratio $\alpha(n\to o)/\alpha(o\to n)$ |
| final rule | $\min\lbrace 1, \frac{zV}{N+1}e^{-\beta\Delta U}\rbrace$ | identical |

In both, the code inserts at a uniformly random point and deletes a uniformly random molecule.

## 18. Reading the acceptance rules with numbers

### 18.1 The three factors

| Factor | Meaning |
|---|---|
| $zV$ | number of molecules an ideal gas at this $\mu$ would put in the box: what the reservoir "wants" |
| $1/(N+1)$ | the box population after the insertion |
| $e^{-\beta\Delta U}$ | energy correction: attractive spots ($\Delta U < 0$) help, overlaps ($\Delta U > 0$) are penalised |

### 18.2 Example: $zV = 50$, $N = 40$ (so $zV/(N+1) = 1.22$)

| $\Delta U$ ($k_BT$) | $e^{-\beta\Delta U}$ | $1.22\times e^{-\beta\Delta U}$ | acc (insertion) |
|---|---|---|---|
| −2 | 7.39 | 9.01 | 1 |
| 0 | 1 | 1.22 | 1 |
| +0.5 | 0.607 | 0.740 | 0.740 |
| +3 | 0.050 | 0.061 | 0.061 |

Deleting a molecule bound by $2\,k_BT$ ($\Delta U = +2$): $\frac{40}{50}e^{-2} = 0.108$.

**A second example, comparing box populations.** Let the reservoir "want" $zV = 10$ molecules:

| Box holds $N$ | $\Delta U$ | $\frac{zV}{N+1}\,e^{-\beta\Delta U}$ | Insertion acceptance |
|---|---|---|---|
| 4 (too few) | 0 | $10/5\times1 = 2$ | 1 (always insert) |
| 19 (too many) | 0 | $10/20\times1 = 0.5$ | 0.5 |
| 9 (about right) | $+1\,k_BT$ | $10/10\times0.37 = 0.37$ | 0.37 |
| 9 | $-1\,k_BT$ | $10/10\times2.72 = 2.72$ | 1 |

The deletion rule mirrors this with $N/(zV)$: when the box holds more than the reservoir wants, deletions are favoured. Together the two rules keep $N$ near the value that equilibrium with the reservoir demands.

### 18.3 Ideal gas: the rules drive N to a Poisson distribution with mean zV

With $\Delta U = 0$, consider the chain of values of $N$. Insertions ($N\to N+1$) and deletions ($N+1\to N$) are attempted equally often, so detailed balance between $N$ and $N+1$ reads

$$P(N)\,\min\left\lbrace 1, \tfrac{zV}{N+1}\right\rbrace = P(N+1)\,\min\left\lbrace 1, \tfrac{N+1}{zV}\right\rbrace$$

$$\Rightarrow\quad \frac{P(N+1)}{P(N)} = \frac{zV}{N+1}$$

(check both cases $N+1 < zV$ and $N+1 > zV$: the ratio is the same). Iterating from $P(0)$: $P(N) = P(0)\,(zV)^N/N!$, the Poisson distribution of Section 7, with mean $zV$. When $N < zV$ the box gains molecules; when $N > zV$ it loses them.

### 18.4 What if insertion and deletion are not attempted equally often?

If insertions were attempted twice as often as deletions with the rules unchanged, the left side above gains a factor 2: $P(N+1)/P(N) = 2zV/(N+1)$. The simulation then behaves as if the activity were $2z$: too many molecules. Unequal attempt probabilities must appear in the acceptance rules.

---

# Part V. Implementation

## 19. Insertion and deletion algorithms

### 19.1 Insertion

1. Choose a **random point uniformly** in the box: $\mathbf r_{new} = (\mathbf u - 0.5)L$, where $\mathbf u$ is three uniform random numbers in $[0,1)$ (box from $-L/2$ to $L/2$).
2. Compute $\Delta U$: energy of a molecule at $\mathbf r_{new}$ with all $N$ existing molecules (minimum image, cutoff) plus $\Delta U_{\text{tail}}$.
3. Accept with $\min\lbrace 1, \frac{zV}{N+1}e^{-\beta\Delta U}\rbrace$. If accepted, store the coordinates in position $N+1$ and set $N \leftarrow N+1$.

### 19.2 Deletion: the coordinate-array trick

The first $N$ entries of the coordinate arrays (`rx`, `ry`, `rz`) are in use:

```
index:   1   2   3   4   5   6   7   8   9  10 ...
        [x] [x] [x] [x] [x] [x] [x] [x] [ ] [ ] ...      N = 8
```

1. **Choose a molecule $j$ at random** (each with probability $1/N$):

```
        [x] [x] [x] [x] [X] [x] [x] [x] [ ] [ ]          j = 5
```

2. **Compute $\Delta U$** $= -(\text{energy of } j \text{ with all others}) + \Delta U_{\text{tail}}$.
3. **If accepted**, slot $j$ becomes empty. Instead of shifting every later molecule one place left, **copy the last molecule into slot $j$** and shorten the list:

```
before:  [x] [x] [x] [x] [ ] [x] [x] [x]                 hole at j = 5
         copy molecule N = 8 into slot 5
after:   [x] [x] [x] [x] [x] [x] [x]                     N = 7
```

```fortran
rx(j) = rx(N)
ry(j) = ry(N)
rz(j) = rz(N)
N = N - 1
```

The cost does not depend on $N$, the occupied part of the array stays contiguous, and the order of identical molecules does not matter. In Python (indices start at 0, last molecule is `pos[N-1]`):

```python
pos[j] = pos[N - 1]
pos = pos[:N - 1]
N -= 1
```

## 20. The tail correction changes with N

In NVT the long-range correction is constant and never enters $\Delta U$. In GCMC it depends on $N$.

**Derivation.** From the cutoff/tail-correction lecture, $U_{\text{tail}} = \frac{2\pi N^2}{V}\int_{r_c}^\infty r^2u(r)\,dr$. For Lennard-Jones in reduced units, $u = 4(r^{-12} - r^{-6})$, so $r^2u = 4(r^{-10} - r^{-4})$, $\int_{r_c}^\infty r^{-10}dr = r_c^{-9}/9$, $\int_{r_c}^\infty r^{-4}dr = r_c^{-3}/3$, and

$$U_{\text{tail}}(N) = \frac{8}{3}\pi\,\frac{N^2}{V}\left[\frac13 r_c^{-9} - r_c^{-3}\right]$$

**Change on insertion and deletion** (using $(N+1)^2 - N^2 = 2N+1$ and $(N-1)^2 - N^2 = -(2N-1)$):

$$\Delta U_{\text{tail}}^{\text{ins}} = \frac{8}{3}\pi\,\frac{2N+1}{V}\left[\frac13 r_c^{-9} - r_c^{-3}\right]$$

$$\Delta U_{\text{tail}}^{\text{del}} = -\frac{8}{3}\pi\,\frac{2N-1}{V}\left[\frac13 r_c^{-9} - r_c^{-3}\right]$$

**Example:** $N = 68$, $V^{\ast} = 512$, $r_c^{\ast} = 3$: $\Delta U_{\text{tail}}^{\text{ins}} = 8.378\times\frac{137}{512}\times(-0.03702) = -0.083\,\varepsilon$. At $T^{\ast} = 2$ this changes the insertion acceptance by about 4%, so it must be included.

## 21. The complete GCMC algorithm (lecture page 3), line by line

### 21.1 The algorithm

1. **Specify** $\mu$ (or $f$, or $z$), $V$, $T$.
2. **Read** initial coordinates or **generate** them (an empty box is fine in GCMC).
3. **Specify** the simulation length, $\Delta r_m$, and the random-number seed.
4. **Perform** `nequil + ncycles` cycles:

```fortran
do l = 1, (nequil + ncycles)
   cycle_length = max(10, N)
   do m = 1, cycle_length
      ! choose the move type
      xrand = ran2(iseed)
      if (xrand < 0.33) then
         ! translation (displacement) move
      else if (xrand < 0.67) then
         ! insertion move
      else
         ! deletion move
      endif
      ! sample properties of interest
   enddo
   ! optimise Delta r_m during the equilibration phase
enddo
```

### 21.2 Line by line

| Line | Purpose and reasoning |
|---|---|
| `cycle_length = max(10, N)` | about $N$ moves per cycle, so each molecule gets roughly one attempt; at least 10 when the box is nearly empty (e.g. $N = 0$ at the start) |
| `xrand = ran2(iseed)` | uniform random number in $[0,1)$ (Numerical Recipes generator) |
| three branches | each move type chosen with probability 1/3 |
| equal insertion and deletion probability | makes $\alpha(o\to n) = \alpha(n\to o)$ in Section 15; otherwise Section 18.4 applies |
| "sample properties" | record $N$, $U$, density profiles … after **every** move, accepted or rejected (repeats carry the Boltzmann weighting) |
| "optimise $\Delta r_m$" | adjust towards ~40–50% displacement acceptance, **only during equilibration** (a changing $\Delta r_m$ breaks the symmetry of $\alpha$) |

### 21.3 Edge cases

- **Deletion at $N = 0$:** nothing to delete; count as rejected, still sample.
- **Displacement at $N = 0$:** nothing to move; count as rejected.
- **Cap $\Delta r_m$ at $L/2$.** At low density the displacement acceptance stays high, so automatic tuning keeps enlarging $\Delta r_m$. A step longer than half the box is meaningless under periodic boundaries, and an unbounded $\Delta r_m$ destroys floating-point precision (two molecules can end up at exactly the same coordinates). This happened while testing the code in Part VI.

---

# Part VI. Tests and working code

## 22. Sanity checks

| Check | What theory says | Section |
|---|---|---|
| ideal gas ($U = 0$) | $\langle N\rangle = zV$, variance $= zV$, Poisson histogram | 7, 18.3 |
| low density | $\rho \approx z\,e^{-2B_2\rho}$ (from $z = \rho e^{\beta\mu^{ex}}$ with $\beta\mu^{ex} \approx 2B_2\rho$) | 8.7 |
| any density | $z$ used in GCMC $= \rho/\langle e^{-\beta\Delta U}\rangle$ from an NVT Widom run at the GCMC density | 8.6, 9 |

## 23. Python code and results

### 23.1 GCMC code (reduced units, Lennard-Jones, follows Section 21)

```python
import numpy as np

T        = 2.0          # temperature T*
z        = 0.05         # activity z* = z sigma^3
L        = 8.0          # box edge L*
rc       = 3.0          # cutoff r_c*
dr_max   = 0.5          # maximum displacement (tuned during equilibration)
n_equil  = 2000         # equilibration cycles
n_cycles = 5000         # production cycles

beta, V = 1.0 / T, L**3
rng = np.random.default_rng(1)

def lj(r2):
    s3 = (1.0 / r2) ** 3
    return 4.0 * (s3 * s3 - s3)

def tail(N):            # U_tail(N): depends on N in GCMC (Section 20)
    return (8.0 / 3.0) * np.pi * N * N / V * ((1.0 / 3.0) * rc**-9 - rc**-3)

def energy_with_others(r, pos, skip=None):
    """Energy of a molecule at r with every molecule in pos (except index skip)."""
    if len(pos) == 0:
        return 0.0
    d = pos - r
    d -= L * np.rint(d / L)                   # minimum image
    r2 = np.einsum("ij,ij->i", d, d)
    if skip is not None:
        r2[skip] = np.inf                     # no self-interaction
    return lj(r2[r2 < rc * rc]).sum()

pos = np.empty((0, 3))                        # start from an empty box
N = 0
N_samples = []

for cycle in range(n_equil + n_cycles):
    n_disp_try = n_disp_acc = 0
    for m in range(max(10, N)):               # cycle length = max(10, N)
        x = rng.random()

        if x < 1 / 3:                         # ---- displacement ----
            if N > 0:
                n_disp_try += 1
                j = rng.integers(N)
                r_old = pos[j].copy()
                r_new = r_old + (2 * rng.random(3) - 1) * dr_max
                r_new -= L * np.rint(r_new / L)
                dU = (energy_with_others(r_new, pos, skip=j)
                      - energy_with_others(r_old, pos, skip=j))
                if rng.random() < np.exp(-beta * dU):
                    pos[j] = r_new
                    n_disp_acc += 1

        elif x < 2 / 3:                       # ---- insertion ----
            r_new = (rng.random(3) - 0.5) * L
            dU = energy_with_others(r_new, pos) + tail(N + 1) - tail(N)
            if rng.random() < z * V / (N + 1) * np.exp(-beta * dU):
                pos = np.vstack([pos, r_new])
                N += 1

        else:                                 # ---- deletion ----
            if N > 0:
                j = rng.integers(N)
                dU = -energy_with_others(pos[j], pos, skip=j) + tail(N - 1) - tail(N)
                if rng.random() < N / (z * V) * np.exp(-beta * dU):
                    pos[j] = pos[N - 1]       # copy the last molecule into slot j
                    pos = pos[:N - 1]
                    N -= 1

        if cycle >= n_equil:
            N_samples.append(N)               # sample after EVERY move

    if cycle < n_equil and n_disp_try > 0:    # tune dr_max during equilibration only
        dr_max *= 1.05 if n_disp_acc / n_disp_try > 0.4 else 0.95
        dr_max = min(dr_max, L / 2)           # never more than half a box

N_samples = np.array(N_samples)
print("<N> =", N_samples.mean(), " var(N) =", N_samples.var(), " rho* =", N_samples.mean() / V)
```

To run the ideal-gas test, make `energy_with_others` and `tail` return 0.

### 23.2 Widom test-particle check (NVT)

```python
# inside an NVT Monte Carlo run, after each production sweep:
for _ in range(50):
    r = (rng.random(3) - 0.5) * L                  # ghost molecule, never actually added
    dU = energy_with_others(r, pos) + tail(N + 1) - tail(N)
    widom.append(np.exp(-beta * dU))
# at the end:
beta_mu_ex = -np.log(np.mean(widom))               # Section 8.6
z_check = (N / V) * np.exp(beta_mu_ex)             # should equal the GCMC z
```

### 23.3 Results ($T^{\ast} = 2.0$, $L^{\ast} = 8$, $V^{\ast} = 512$, $r_c^{\ast} = 3$; error bars from 10 blocks)

**Ideal gas** ($z^{\ast} = 0.05$, so $zV = 25.6$): $\langle N\rangle = 25.60$, variance $= 26.39$ (theory: 25.6 and 25.6).

**Lennard-Jones:**

| $z^{\ast}$ | Ideal gas $\rho^{\ast} = z^{\ast}$ | GCMC $\rho^{\ast}$ | Second-virial estimate | Displacement acc. | Insertion acc. |
|---|---|---|---|---|---|
| 0.01 | 0.0100 | 0.0101 ± 0.0001 | 0.0103 | 0.93 | 0.82 |
| 0.03 | 0.0300 | 0.0323 ± 0.0004 | 0.0327 | 0.81 | 0.83 |
| 0.10 | 0.1000 | 0.1343 ± 0.0009 | 0.147 | 0.47 | 0.60 |

- At low $z$ the fluid is almost ideal ($\rho \approx z$).
- As $z$ grows, attractions make it denser than an ideal gas ($B_2^{\ast} = -1.314 < 0$).
- The second-virial estimate is good at low density and overshoots at $\rho^{\ast} \approx 0.13$, where higher virial terms matter.
- Insertion acceptance falls as the fluid gets denser.

**Widom check at $z^{\ast} = 0.10$:** NVT at $\rho^{\ast} = 0.1348$ ($N = 69$) gave $\beta\mu^{ex} = -0.300 \pm 0.002$, so $z^{\ast} = 0.1348\times e^{-0.300} = 0.0998$, matching the GCMC input of 0.10.

---

# Part VII. Wrap-up

## 24. Practical issues

- **Insertions fail in dense phases.** In a dense liquid ($\rho^{\ast} \approx 0.8$) almost every random point overlaps a molecule, so $\Delta U$ is huge and insertions are almost never accepted; deletions of well-bound molecules are also rare. $N$ then changes very slowly. Configurational-bias or cavity-biased insertions address this.
- **Insertion and deletion must be attempted equally often**, or the ratio must enter the acceptance rules (Section 18.4).
- **Equilibration** is judged by watching $N$ (not only $U$) settle.
- **Sampling** happens after every move, including rejected moves and failed deletions at $N = 0$.
- **Box size:** $L > 2r_c$ for the minimum image, and large enough that $\langle N\rangle$ is not tiny.
- **Molecules** (not atoms) need a random orientation on insertion, and its probability must be included in $\alpha$.
- **Tail correction** must be in $\Delta U$ for insertions and deletions (Section 20).

## 25. Corrections and clarifications to the lecture notes

1. **Two meanings of $f$.** The lecture uses $f$ for the microstate probability and for the fugacity. Here the probability is $p$.
2. **Microstate probability (page 1).** The $V^N$ factor belongs with the scaled coordinates $\mathbf s$; in real coordinates $\mathbf r$ the probability has no $V^N$.
3. **Insertion rule (page 1).** One line has $e^{\beta\mu N}$ in the numerator; it should be $e^{\beta\mu}$, because only one molecule is added. Section 15.1 gives the derivation and three independent checks.

## 26. Practice questions (with answers)

1. **Why does $1/(N+1)$ appear in the insertion rule?**
   From $N!/(N+1)!$ (identical molecules). Equivalently, from the probability $1/(N+1)$ of picking the newly inserted molecule in the reverse deletion (Section 17).

2. **Why do we Taylor-expand $\ln\Omega_R$ instead of $\Omega_R$?**
   $\Omega_R$ changes by a factor $e$ when the reservoir gives away just $k_BT$ of energy; a linear expansion of $\Omega_R$ predicts negative numbers of states. $\ln\Omega_R$ changes linearly over such small steps (Section 2.4).

3. **Show that $\langle N\rangle = \frac{1}{\beta}\frac{\partial\ln\Xi}{\partial\mu}$.**
   Section 6.2.

4. **For an ideal gas with $zV = 5$, what are $\langle N\rangle$ and its variance?**
   Both 5 (Poisson).

5. **In an NVT run, ghost insertions give $\langle e^{-\beta\Delta U}\rangle = 0.5$. What is $\beta\mu^{ex}$, and is insertion easier or harder than in an ideal gas?**
   $\beta\mu^{ex} = -\ln 0.5 = +0.69$; harder (repulsion dominates).

6. **An ideal-gas GCMC run gives $\langle N\rangle = 12.3$ in 100 nm³ at 300 K. What pressure was simulated?**
   $z = 0.123$ nm⁻³ $= 1.23\times10^{26}$ m⁻³; $P = zk_BT = 1.23\times10^{26}\times1.381\times10^{-23}\times300 = 5.1\times10^5$ Pa $= 5.1$ bar.

7. **$N = 10$, $zV = 10$, a deletion has $\Delta U = 0$. Acceptance?**
   $\min\lbrace 1, 10/10\rbrace  = 1$.

8. **Why must the tail correction be in $\Delta U$ for insertions but not displacements?**
   It depends on $N$, which only insertions and deletions change.

9. **Insertions are attempted twice as often as deletions with unchanged rules. What happens?**
   The simulation behaves as if the activity were $2z$ (Section 18.4).

10. **For the Langmuir model, at what pressure is the surface half covered?**
    $P = 1/K$, with $K = \beta\Lambda^3e^{\beta\varepsilon_b}$ (Section 12).

11. **What is the difference between $\mu^{ex}$ and the fugacity coefficient $\phi$?**
    $\mu^{ex}$ compares with an ideal gas at the same density; $\phi$ compares at the same pressure. $\ln\phi = \beta\mu^{ex} - \ln\mathcal Z$ (Section 10.3).

12. **(Conceptual) Why is GCMC much cheaper than NVT for adsorption in a slit pore?**
    Only the pore is simulated. The macroscopic adsorbent and the outside gas are replaced by one number, the gas's $\mu$ (or $z$, or $f$), because pore and gas are in equilibrium when their $\mu$ and $T$ are equal (Section 11).

13. **A GCMC box holds $N = 24$ molecules and $zV = 30$. A trial insertion has $\Delta U = +1.0\,k_BT$. Acceptance?**
    $\min\lbrace 1, \frac{30}{25}e^{-1}\rbrace = 1.2\times0.368 = 0.441$. (Common slips: forgetting $zV/(N+1)$ gives 0.368; using $N$ instead of $N+1$ gives 0.460.)

14. **$N = 20$, $zV = 25$. You try to delete a molecule whose interaction energy with the others is $-2\,k_BT$. Acceptance?**
    Removing it gives $\Delta U = +2\,k_BT$, so $\min\lbrace 1, \frac{20}{25}e^{-2}\rbrace = 0.8\times0.135 = 0.108$.

15. **An ideal-gas reservoir is at 300 K and 2 bar ($f = P$). The box volume is 50 nm³. What is $\langle N\rangle$?**
    $z = \beta f = 2\times10^5/(1.381\times10^{-23}\times300) = 4.83\times10^{25}$ m⁻³ $= 0.0483$ nm⁻³, so $\langle N\rangle = zV = 2.41$ (and the variance is also 2.41).

16. **(Conceptual) Is activity a property of the system? When a pore and a gas reach equilibrium, which becomes equal, density or activity?**
    Yes, it is a property of the system's state: its escaping tendency. The **activities** become equal; the densities can differ enormously (Section 9.3), just as temperatures, not energy contents, become equal in heat exchange.

17. **An NVT Widom run gives a welcome factor $\langle e^{-\beta\Delta U}\rangle = 1.35$ at $\rho^{\ast} = 0.135$. What activity $z^{\ast}$ would a GCMC run need to reproduce this density?**
    $z^{\ast} = \rho^{\ast}/\langle e^{-\beta\Delta U}\rangle = 0.135/1.35 = 0.100$.


## 27. Summary: the whole chain

1. **Isolated box + reservoir** → all combined microstates equally likely.
2. **Count the reservoir:** $P(s) \propto \Omega_R(E_{tot} - E_s, N_{tot} - N)$.
3. **Take the log and Taylor-expand** to first order; the slopes are $\beta$ and $-\beta\mu$ (from $dS = dE/T + P\,dV/T - \mu\,dN/T$ and $S = k_B\ln\Omega$).
4. **Exponentiate:** $P(s) \propto e^{-\beta(E_s - \mu N)}$; **normalize:** $\Xi = \sum_s e^{-\beta(E_s - \mu N)} = \sum_N e^{\beta\mu N}Q(N,V,T)$.
5. **Classical $Q$:** $Q = Z_N/(N!\Lambda^{3N})$, with $\Lambda$ from the Gaussian momentum integrals and $N!$ from identical molecules.
6. **Lecture formula:** $\Xi = \sum_N \frac{e^{\beta\mu N}}{N!\Lambda^{3N}}\int e^{-\beta U}d\mathbf r^N = \sum_N \frac{V^Ne^{\beta\mu N}}{N!\Lambda^{3N}}\int e^{-\beta U}d\mathbf s^N$.
7. **Thermodynamics:** $\langle N\rangle = \beta^{-1}\partial\ln\Xi/\partial\mu$, variance $= \beta^{-2}\partial^2\ln\Xi/\partial\mu^2$, $PV = k_BT\ln\Xi$.
8. **Chemical potential:** $\beta\mu = \ln(\rho\Lambda^3) + \beta\mu^{ex}$ with $\beta\mu^{ex} = -\ln\langle e^{-\beta\Delta U}\rangle$ (Widom).
9. **Activity and fugacity:** $z = e^{\beta\mu}/\Lambda^3 = \rho e^{\beta\mu^{ex}} = \beta f$; for an ideal gas $z = \rho = \beta P$.
10. **Acceptance rules** (detailed balance on $p(\mathbf s^N;N)$): displacement $\min\lbrace 1, e^{-\beta\Delta U}\rbrace$; insertion $\min\lbrace 1, \frac{zV}{N+1}e^{-\beta\Delta U}\rbrace$; deletion $\min\lbrace 1, \frac{N}{zV}e^{-\beta\Delta U}\rbrace$.
11. **Algorithm:** each move with probability 1/3, cycle length $\max(10, N)$, array trick for deletion, $N$-dependent tail correction, sample after every move, tune $\Delta r_m$ (capped at $L/2$) only during equilibration.
12. **Application:** adsorption isotherms. Simulate only the pore at the $\mu$ of the outside gas; the Langmuir isotherm is the exactly solvable special case.
