# Derivation of the Grand Canonical Partition Function

*Companion notes to Lecture 8 (GCMC). Starting from one postulate, these notes derive the grand partition function used in the lecture, show where every factor of the GCMC acceptance rules comes from, and check the result on the ideal gas and on Langmuir adsorption.*

---

## Contents

0. Notation
1. The idea in one sentence
2. Part A: General derivation (system + reservoir)
   - Step 1: The set-up
   - Step 2: The one postulate
   - Step 3: Probability of one system microstate
   - Step 4: The Taylor expansion of ln Ω (explained in detail)
   - Step 5: The derivatives are temperature and chemical potential
   - Step 6: Undo the logarithm
   - Step 7: Normalize to get Ξ
3. Part B: Classical canonical partition function for N identical atoms
4. Part C: Assembling the lecture's formula
5. Part D: Checks on the ideal gas
6. Part E: Worked example, the Langmuir adsorption isotherm
7. The whole chain

---

## 0. Notation

| Symbol | Meaning |
|---|---|
| S, R | the system (simulation box) and the reservoir around it |
| $E_s$ | energy of system microstate $s$ |
| $N$ | number of molecules in the system |
| $E_{tot}$, $N_{tot}$ | total energy and total number of molecules of system + reservoir (fixed) |
| $\Omega_R(E, N)$ | number of reservoir microstates with energy $E$ and $N$ molecules |
| $P(s)$ | probability that the system is in microstate $s$ |
| $k_B$ | Boltzmann's constant |
| $T$ | temperature |
| $\beta$ | $1/(k_BT)$ |
| $S$ | entropy, $S = k_B\ln\Omega$ |
| $P$ | pressure (in thermodynamic formulas) |
| $\mu$ | chemical potential |
| $V$ | system volume |
| $Q(N, V, T)$ | canonical partition function for $N$ molecules |
| $\Xi(\mu, V, T)$ | grand partition function |
| $H$, $K$, $U$ | total energy (Hamiltonian), kinetic energy, potential energy |
| $\mathbf r^N$, $\mathbf p^N$ | positions and momenta of all $N$ atoms |
| $m$ | mass of one atom |
| $h$ | Planck's constant |
| $\Lambda$ | thermal de Broglie wavelength, $h/\sqrt{2\pi m k_BT}$ |
| $Z_N$ | configurational integral, $\int e^{-\beta U}\,d\mathbf r^N$ |
| $\mathbf s$ | scaled coordinates, $\mathbf r/V^{1/3}$ |
| $z$ | activity, $e^{\beta\mu}/\Lambda^3$ |
| $\rho$ | number density |
| $\Phi$ | grand potential, $-k_BT\ln\Xi$ |

---

## 1. The idea in one sentence

**Put a small box in contact with a huge reservoir, count how many ways the reservoir can arrange itself for each state of the box, and that count becomes $e^{-\beta(E - \mu N)}$.**

---

## 2. Part A: General derivation

### Step 1: The set-up

- **System S:** the simulation box, volume $V$. It can exchange **energy** and **molecules** with its surroundings.
- **Reservoir R:** a huge amount of the same substance around the box.
- **S + R together** are isolated from everything else, so the totals are fixed:

$$E_s + E_R = E_{tot}, \qquad N + N_R = N_{tot}$$

The reservoir is enormous compared with the box: $E_s \ll E_{tot}$ and $N \ll N_{tot}$.

### Step 2: The one postulate

For an isolated system, **every microstate with the given total energy and particle number is equally likely**. This is the basic postulate of statistical mechanics (the microcanonical ensemble). Everything else follows from counting.

### Step 3: Probability of one system microstate

Pick one particular microstate $s$ of the box, with $N$ molecules and energy $E_s$.

The box's state is fixed, so the reservoir must hold the rest: energy $E_{tot} - E_s$ and $N_{tot} - N$ molecules. It can be in any of its microstates with those values. Since every combined microstate is equally likely, the probability of $s$ is proportional to the number of ways the reservoir can accommodate it:

$$P(s) \propto \Omega_R(E_{tot} - E_s,\; N_{tot} - N)$$

**The box prefers states that leave the reservoir the most ways to arrange itself.**

### Step 4: The Taylor expansion of ln Ω

$$\boxed{\;\ln\Omega_R(E_{tot} - E_s,\, N_{tot} - N) \approx \ln\Omega_R(E_{tot}, N_{tot}) \;-\; E_s\,\frac{\partial \ln\Omega_R}{\partial E} \;-\; N\,\frac{\partial \ln\Omega_R}{\partial N}\;}$$

#### 4.1 What we are trying to do

We do not know the function $\Omega_R$: the reservoir has around $10^{23}$ molecules and could be any complicated fluid. But we do not need the whole function. We only need **how much $\Omega_R$ changes when the box takes away a tiny energy $E_s$ and a few molecules $N$**. That question ("how does a function change when its inputs change a little?") is what a Taylor expansion answers.

#### 4.2 Taylor expansion with one variable

For a smooth function $g$, changing its input from $a$ to $a - h$ with $h$ small:

$$g(a - h) \approx g(a) - h\,g'(a)$$

New value ≈ old value − (step) × (slope), where $g'$ is the derivative $dg/dx$.

**Numerical check:** $g(x) = \ln x$, $a = 1000$, $h = 2$. The slope is $g'(x) = 1/x$, so $g'(1000) = 0.001$.

- Exact: $\ln 998 = 6.905753$
- Approximation: $\ln 1000 - 2\times0.001 = 6.907755 - 0.002 = 6.905755$

They agree to six digits because the step (2) is tiny compared with where we are (1000).

#### 4.3 Two variables

If both inputs change a little, add one correction for each:

$$g(a - h,\; b - k) \approx g(a, b) - h\,\frac{\partial g}{\partial a} - k\,\frac{\partial g}{\partial b}$$

Each **partial derivative** is the slope in one direction while the other input is held fixed. (On a hillside: change in height ≈ steps east × slope eastward + steps north × slope northward.)

#### 4.4 Matching it to our line

| General formula | Our line |
|---|---|
| $g$ | $\ln\Omega_R$ |
| $a$ (where we are) | $E_{tot}$ |
| $h$ (small step in $a$) | $E_s$, the energy the box takes |
| $b$ | $N_{tot}$ |
| $k$ (small step in $b$) | $N$, the molecules in the box |

Substituting gives the boxed line. It is accurate because $E_s \ll E_{tot}$ and $N \ll N_{tot}$.

#### 4.5 Why expand ln Ω and not Ω itself

$\Omega$ changes **enormously** with even the smallest change in energy. Take an ideal monatomic gas of $N_R$ atoms at fixed volume as the reservoir. A standard result is $\Omega_R \propto E^{3N_R/2}$, so

$$\ln\Omega_R = \frac{3N_R}{2}\ln E + \text{const} \quad\Rightarrow\quad \frac{\partial\ln\Omega_R}{\partial E} = \frac{3N_R}{2E}$$

The energy of the gas is $E_{tot} = \frac{3}{2}N_R k_BT$, so

$$\frac{\partial\ln\Omega_R}{\partial E} = \frac{3N_R}{2\cdot\frac32 N_R k_BT} = \frac{1}{k_BT} = \beta$$

Now let the box take energy $E_s$ and compare the two ways of expanding:

| Energy taken $E_s$ | True change $\Omega_{new}/\Omega_{old} \approx e^{-\beta E_s}$ | Linear expansion of $\Omega$: $1 - \beta E_s$ | Linear expansion of $\ln\Omega$: $e^{-\beta E_s}$ |
|---|---|---|---|
| $1\,k_BT$ | 0.368 | 0 | 0.368 ✓ |
| $5\,k_BT$ | 0.0067 | −4 | 0.0067 ✓ |
| $10\,k_BT$ | 0.000045 | −9 | 0.000045 ✓ |

Taking away just one molecule's worth of thermal energy, out of $10^{23}$ molecules, already cuts the reservoir's number of states by a factor of $e$.

- Expanding $\Omega$ directly predicts zero and then **negative** numbers of states: nonsense.
- Expanding $\ln\Omega$ is accurate, because $\ln\Omega$ changes **linearly** over such tiny steps.

#### 4.6 Why the first derivative is enough

The next Taylor term is $\tfrac12 E_s^2\,\partial^2\ln\Omega_R/\partial E^2$. For the ideal-gas reservoir, $\partial^2\ln\Omega_R/\partial E^2 = -\frac{3N_R}{2E^2}$, so its size is

$$\frac12 E_s^2\,\frac{3N_R/2}{E_{tot}^2} = \frac{E_s^2}{2\cdot\frac32 N_R\,(k_BT)^2}$$

For $E_s = 10\,k_BT$ and $N_R = 6\times10^{23}$ this is about $6\times10^{-23}$. In general it is proportional to 1/(heat capacity of the reservoir), which vanishes for a huge reservoir. Only the first-order terms survive.

### Step 5: The derivatives are temperature and chemical potential

Boltzmann's entropy is $S = k_B\ln\Omega$. The fundamental thermodynamic relation $dE = T\,dS - P\,dV + \mu\,dN$, rearranged for $dS$:

$$dS = \frac{1}{T}\,dE + \frac{P}{T}\,dV - \frac{\mu}{T}\,dN$$

Divide by $k_B$ and read off the coefficients:

$$\boxed{\;\frac{\partial\ln\Omega_R}{\partial E} = \frac{1}{k_BT} = \beta, \qquad \frac{\partial\ln\Omega_R}{\partial N} = -\frac{\mu}{k_BT} = -\beta\mu\;}$$

**Meaning:**

- **$\beta$:** how fast the reservoir loses states when it gives away energy. A cold reservoir (large $\beta$) loses many states per unit of energy given up, so low-energy states of the box are favoured more strongly at low $T$.
- **$-\beta\mu$:** how fast the reservoir's number of states changes when it gives away molecules. For a gas $\mu$ is negative. For argon gas at 300 K and 1 bar, $\beta\mu = \ln(\rho\Lambda^3) \approx -16$ (this formula is derived in Part D). Each molecule moved into the box multiplies the reservoir's count by $e^{\beta\mu} \approx 10^{-7}$, so a low $\mu$ keeps the box emptier.

Both are properties of the **reservoir**: this is how the reservoir imposes its $T$ and $\mu$ on the box.

### Step 6: Undo the logarithm

Substitute the derivatives into the expansion:

$$\ln\Omega_R(E_{tot} - E_s,\, N_{tot} - N) \approx \ln\Omega_R(E_{tot}, N_{tot}) - \beta E_s + \beta\mu N$$

Exponentiate:

$$\Omega_R(E_{tot} - E_s,\, N_{tot} - N) \approx \underbrace{\Omega_R(E_{tot}, N_{tot})}_{\text{same for every } s}\; e^{-\beta E_s}\, e^{\beta\mu N}$$

The first factor is the same for every state $s$, so it goes into the normalization:

$$\boxed{\;P(s) = \frac{e^{-\beta(E_s - \mu N)}}{\Xi}\;}$$

- **$e^{-\beta E_s}$:** low energy is favoured, as in the canonical ensemble.
- **$e^{\beta\mu N}$:** each molecule in the box multiplies the weight by $e^{\beta\mu}$. High reservoir $\mu$ favours a full box; low $\mu$ an empty one.

### Step 7: Normalize to get Ξ

The probabilities must add to 1 over **all** states, including every possible $N$:

$$\Xi = \sum_{N=0}^{\infty}\;\sum_{s \text{ with } N} e^{-\beta(E_s - \mu N)}$$

For fixed $N$, $e^{\beta\mu N}$ is a common factor and comes out of the inner sum. What is left is the **canonical partition function** $Q(N,V,T) = \sum_s e^{-\beta E_s}$:

$$\boxed{\;\Xi(\mu, V, T) = \sum_{N=0}^{\infty} e^{\beta\mu N}\, Q(N, V, T)\;}$$

**The grand partition function is a weighted sum of canonical partition functions, one for each $N$.**

---

## 3. Part B: Classical canonical partition function for N identical atoms

For atoms obeying classical mechanics, the sum over states becomes an integral over positions and momenta:

$$Q(N, V, T) = \frac{1}{N!\,h^{3N}} \int\! d\mathbf p^N \!\int\! d\mathbf r^N\; e^{-\beta H(\mathbf p^N, \mathbf r^N)}, \qquad H = \sum_{i=1}^{N}\frac{p_i^2}{2m} + U(\mathbf r^N)$$

| Factor | Why it is there |
|---|---|
| $h^{3N}$ | quantum mechanics assigns one state to each phase-space cell of size $h$ per coordinate–momentum pair; dividing by $h^{3N}$ turns the integral into a count of states and makes it dimensionless |
| $N!$ | the atoms are identical; swapping two gives the same physical state, but integrating over all labelled positions counts it $N!$ times |

### Separate kinetic and potential parts

Because $e^{-\beta H} = e^{-\beta K}\,e^{-\beta U}$:

$$Q = \frac{1}{N!\,h^{3N}}\underbrace{\int e^{-\beta\sum_i p_i^2/2m}\,d\mathbf p^N}_{\text{momenta}}\;\underbrace{\int e^{-\beta U(\mathbf r^N)}\,d\mathbf r^N}_{Z_N}$$

### The momentum integral

There are $3N$ independent momentum components, each giving the same Gaussian integral. Using

$$\int_{-\infty}^{\infty} e^{-a x^2}\,dx = \sqrt{\frac{\pi}{a}}$$

(proof: square the integral, write it as a 2D integral, switch to polar coordinates, and it becomes $\pi/a$) with $a = \beta/(2m)$:

$$\int_{-\infty}^{\infty} e^{-\beta p_x^2/2m}\,dp_x = \sqrt{\frac{2\pi m}{\beta}} = \sqrt{2\pi m k_BT}$$

All $3N$ components together give $(2\pi m k_BT)^{3N/2}$. Combined with $h^{3N}$:

$$\frac{(2\pi m k_BT)^{3N/2}}{h^{3N}} = \frac{1}{\Lambda^{3N}}, \qquad \Lambda = \frac{h}{\sqrt{2\pi m k_BT}}$$

$\Lambda$ is the **thermal de Broglie wavelength**, roughly the quantum "size" of an atom at temperature $T$. For argon at 300 K, $\Lambda = 0.16$ Å, tiny compared with the roughly 3.6 Å spacing between atoms in liquid argon. That is why classical mechanics works for argon.

$$\boxed{\;Q(N, V, T) = \frac{Z_N}{N!\,\Lambda^{3N}}, \qquad Z_N = \int e^{-\beta U(\mathbf r^N)}\,d\mathbf r^N\;}$$

---

## 4. Part C: Assembling the lecture's formula

Put Part B into Step 7:

$$\Xi(\mu, V, T) = \sum_{N=0}^{\infty}\frac{e^{\beta\mu N}}{N!\,\Lambda^{3N}}\int\!\!\cdots\!\!\int e^{-\beta U(\mathbf r^N)}\,d\mathbf r^N$$

This is the first line of the lecture.

### Scaled coordinates

With $\mathbf s = \mathbf r/V^{1/3}$, every component of $\mathbf s$ runs from 0 to 1. Each atom's volume element becomes $d^3\mathbf r = V\,d^3\mathbf s$, so $d\mathbf r^N = V^N d\mathbf s^N$:

$$\Xi = \sum_{N=0}^{\infty}\frac{V^N e^{\beta\mu N}}{N!\,\Lambda^{3N}}\int\!\!\cdots\!\!\int e^{-\beta U}\,d\mathbf s^N$$

This is the lecture's second line.

### Probability of one microstate

After integrating out the momenta, the probability of $N$ atoms at scaled positions $\mathbf s^N$ is the term inside the sum divided by $\Xi$:

$$p(\mathbf s^N; N) = \frac{1}{\Xi}\,\frac{V^N e^{\beta\mu N}}{N!\,\Lambda^{3N}}\,e^{-\beta U(\mathbf s^N)}$$

### Where every factor of the insertion rule comes from

The ratio $p(N+1)/p(N)$ gives the insertion rule $\dfrac{Ve^{\beta\mu}}{(N+1)\Lambda^3}\,e^{-\beta\Delta U}$:

| Factor | Origin in the derivation |
|---|---|
| $e^{\beta\mu}$ | the reservoir's chemical potential (Step 5) |
| $1/\Lambda^3$ | one more atom's momentum integral (Part B) |
| $1/(N+1)$ | $N!/(N+1)!$, identical atoms (Part B) |
| $V$ | one more atom's position, in scaled coordinates (Part C) |
| $e^{-\beta\Delta U}$ | the configurational Boltzmann factor |

---

## 5. Part D: Checks on the ideal gas

With no interactions, $U = 0$ and $Z_N = \int d\mathbf r^N = V^N$. With the activity $z = e^{\beta\mu}/\Lambda^3$:

$$\Xi_{ig} = \sum_{N=0}^{\infty}\frac{(zV)^N}{N!} = e^{zV}$$

using the exponential series $\sum_N x^N/N! = e^x$.

### (a) Average number of molecules

The probability of $N$ molecules is $e^{\beta\mu N}Q(N)/\Xi$. Differentiate $\Xi$ with respect to $\mu$: each term gains a factor $\beta N$, so

$$\frac{\partial\Xi}{\partial\mu} = \sum_N \beta N\,e^{\beta\mu N}Q(N) = \beta\,\Xi\,\langle N\rangle \quad\Rightarrow\quad \langle N\rangle = \frac{1}{\beta}\frac{\partial\ln\Xi}{\partial\mu}$$

For the ideal gas, $\ln\Xi_{ig} = zV$ and $\partial z/\partial\mu = \beta z$, so $\langle N\rangle = zV$. The density is $\rho = \langle N\rangle/V = z$, which gives

$$\boxed{\;\beta\mu = \ln(\rho\Lambda^3)\;}$$

This is the ideal-gas chemical potential that Lecture 8 "recalls"; it comes straight out of $\Xi$. (Argon at 300 K and 1 bar: $\rho = P/k_BT = 2.41\times10^{25}$ m⁻³, $\Lambda^3 = 4.1\times10^{-33}$ m³, so $\beta\mu = \ln(9.8\times10^{-8}) = -16.1$.)

### (b) The ideal gas law

The grand potential is $\Phi = -k_BT\ln\Xi$. From the Euler relation of thermodynamics, $E = TS - PV + \mu N$, so $\Phi = E - TS - \mu N = -PV$. Therefore

$$PV = k_BT\ln\Xi_{ig} = k_BT\,zV = \langle N\rangle k_BT$$

which is the ideal gas law.

### (c) Fluctuations in N

Differentiating once more gives the variance:

$$\langle N^2\rangle - \langle N\rangle^2 = \frac{1}{\beta^2}\frac{\partial^2\ln\Xi}{\partial\mu^2} = zV = \langle N\rangle$$

The probability of $N$ molecules is $e^{-zV}(zV)^N/N!$, a **Poisson distribution**. The GCMC test code from Lecture 8 gave $\langle N\rangle = 25.60$ and variance 26.4, against 25.6 for both.

---

## 6. Part E: Worked example, the Langmuir adsorption isotherm

The same machinery applied to a discrete system.

### One site

An adsorption site can be empty or hold one molecule with binding energy $\varepsilon_b > 0$:

| State | $N$ | Energy | Weight $e^{-\beta(E - \mu N)}$ |
|---|---|---|---|
| empty | 0 | 0 | 1 |
| occupied | 1 | $-\varepsilon_b$ | $e^{\beta(\mu + \varepsilon_b)}$ |

Let $x = e^{\beta(\mu + \varepsilon_b)}$. Then

$$\Xi_{site} = 1 + x, \qquad \theta = \frac{x}{1 + x}$$

where $\theta$ is the fraction of time the site is occupied.

### Two and M sites

Two independent sites have four states: (0,0) with weight 1, (1,0) and (0,1) with weight $x$ each, (1,1) with weight $x^2$. So $\Xi = 1 + 2x + x^2 = (1 + x)^2$. For $M$ independent sites, $\Xi = (1 + x)^M$: independent sites multiply.

### Connecting to gas pressure

For an ideal gas at pressure $P$, $z = \rho = P/(k_BT)$, so $e^{\beta\mu} = z\Lambda^3 = \beta P\Lambda^3$. Then $x = KP$ with $K = \beta\Lambda^3 e^{\beta\varepsilon_b}$, and

$$\boxed{\;\theta = \frac{KP}{1 + KP}\;}$$

which is the **Langmuir isotherm**. The site is half full at $P = 1/K$, and stronger binding (larger $\varepsilon_b$) makes $K$ larger, so it fills at lower pressure. GCMC in a pore does the same thing for a continuous, interacting system where no closed-form answer exists.

---

## 7. The whole chain

1. System + reservoir are isolated, so every combined microstate is equally likely.
2. $P(s) \propto \Omega_R(E_{tot} - E_s,\; N_{tot} - N)$.
3. Taylor-expand $\ln\Omega_R$ (not $\Omega_R$, which changes far too violently); the slopes are $\beta$ and $-\beta\mu$; higher orders vanish for a huge reservoir.
4. Exponentiate: $P(s) \propto e^{-\beta(E_s - \mu N)}$.
5. Normalize: $\Xi = \sum_N e^{\beta\mu N}Q(N, V, T)$.
6. Classically, $Q = Z_N/(N!\Lambda^{3N})$, with $\Lambda$ from the momentum integrals and $N!$ from identical atoms.
7. Together: $\Xi = \sum_N \frac{e^{\beta\mu N}}{N!\Lambda^{3N}}\int e^{-\beta U}d\mathbf r^N$, the lecture's formula. Ratios of its terms give the GCMC acceptance rules.
8. Checks: the ideal gas gives $\beta\mu = \ln\rho\Lambda^3$, $PV = Nk_BT$ and Poisson fluctuations; independent sites give the Langmuir isotherm.
