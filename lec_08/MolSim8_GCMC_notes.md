# Lecture 8: Grand Canonical Monte Carlo (GCMC)

*Detailed notes based on 8.pdf ("GCMC Simulation", 1 February 2022). The acceptance rules were checked with a working Python GCMC code (Section 9), an ideal-gas test, and an independent Widom-insertion test.*

---

## Contents

0. Notation
1. Open systems and the grand canonical ensemble
2. The grand partition function, term by term
3. Why GCMC: adsorption in a slit pore
4. The three Monte Carlo moves and their acceptance rules
5. Chemical potential, activity and fugacity
6. How the insertion and deletion moves are implemented
7. The complete GCMC algorithm (lecture page 3)
8. Sanity checks
9. Working Python code and results
10. Practical issues
11. Corrections and clarifications to the lecture notes
12. Practice questions (with answers)
13. Summary

---

## 0. Notation

| Symbol | Meaning |
|---|---|
| $N$ | number of molecules in the simulation box (changes during GCMC) |
| $V$, $L$ | box volume and box edge length ($V = L^3$) |
| $T$ | temperature |
| $k_B$ | Boltzmann's constant |
| $\beta$ | $1/(k_BT)$ |
| $\mu$ | chemical potential (fixed by the reservoir) |
| $U(\mathbf r^N)$ | potential energy of the configuration $\mathbf r^N = (\mathbf r_1, \dots, \mathbf r_N)$ |
| $\Delta U$ | $U(\text{new}) - U(\text{old})$ |
| $h$ | Planck's constant |
| $m$ | mass of one molecule |
| $\Lambda$ | thermal de Broglie wavelength, $h/\sqrt{2\pi m k_BT}$ |
| $\mathbf s$ | scaled coordinates, $\mathbf s = \mathbf r/V^{1/3} = \mathbf r/L$ (each component between 0 and 1) |
| $\Xi$ | grand partition function |
| $p(o)$, $p(n)$ | probability (density) of the old and new microstate |
| $\alpha(o\to n)$ | probability of proposing the move $o \to n$ |
| $\text{acc}(o\to n)$ | probability of accepting it |
| $\rho$ | number density $N/V$ |
| $\mu^{\text{ex}}$ | excess chemical potential (the part due to interactions) |
| $z$ | activity, $e^{\beta\mu}/\Lambda^3$ |
| $f$ | fugacity (a "corrected pressure") |
| $P$ | pressure |
| $\Delta r_m$ | maximum displacement in a translation move |

> **Notation warning.** The lecture writes $f$ both for the probability of a microstate and for the fugacity. In these notes the probability is $p(o)$ and **$f$ always means fugacity**.

---

## 1. Open systems and the grand canonical ensemble

The **grand canonical ensemble** describes an **open system**: a box that exchanges **molecules** (and energy) with a very large **reservoir**. The reservoir fixes three things:

| Fixed | Fluctuates |
|---|---|
| chemical potential $\mu$ | number of molecules $N$ |
| volume $V$ | energy $U$ |
| temperature $T$ | |

**Physical picture.** The box is a small region inside a huge amount of gas. Molecules wander in and out. At equilibrium the box and the reservoir have the **same $\mu$ and $T$**, so molecules have no net tendency to flow either way. $N$ fluctuates around its average $\langle N\rangle$.

**What $\mu$ means.** $\mu$ is the free-energy change when one molecule is added. A high $\mu$ in the reservoir pushes molecules into the box; a low $\mu$ pulls them out.

Compare with the ensembles used so far:

| Ensemble | Fixed | Monte Carlo moves |
|---|---|---|
| Canonical (NVT) | $N, V, T$ | displacement (and rotation) |
| Grand canonical ($\mu VT$) | $\mu, V, T$ | displacement **+ insertion + deletion** |

---

## 2. The grand partition function, term by term

### 2.1 The formula

$$\Xi(\mu, V, T) = \sum_{N=0}^{\infty} \frac{e^{\beta\mu N}}{N!\,\Lambda^{3N}} \int\!\!\cdots\!\!\int e^{-\beta U(\mathbf r^N)}\, d\mathbf r^N$$

Read it as: **sum over every possible number of molecules $N$**, and for each $N$ take the usual canonical (NVT) weight, multiplied by $e^{\beta\mu N}$.

| Piece | Where it comes from | Meaning |
|---|---|---|
| $\sum_{N=0}^{\infty}$ | $N$ is not fixed | every number of molecules is allowed |
| $e^{\beta\mu N}$ | the reservoir | each molecule in the box contributes a factor $e^{\beta\mu}$; large $\mu$ favours large $N$ |
| $\dfrac{1}{\Lambda^{3N}}$ | integrating out the momenta | each molecule's kinetic energy gives $\int e^{-\beta p^2/2m}\,d^3p/h^3 = 1/\Lambda^3$ |
| $\dfrac{1}{N!}$ | identical molecules | swapping two identical molecules gives the same physical state, so integrating over all labelled positions overcounts by $N!$ |
| $\int e^{-\beta U}\,d\mathbf r^N$ | configurations | the Boltzmann weight of every arrangement, as in NVT |

Equivalently, $\Xi = \sum_N e^{\beta\mu N} Q(N, V, T)$, where $Q(N,V,T)$ is the ordinary canonical partition function for $N$ molecules.

### 2.2 Scaled coordinates

The lecture rewrites the integral with **scaled coordinates** $\mathbf s = \mathbf r/V^{1/3}$, so every component of $\mathbf s$ runs over the unit interval whatever the box size. Each molecule's volume element becomes $d^3\mathbf r = V\, d^3\mathbf s$, so for $N$ molecules $d\mathbf r^N = V^N d\mathbf s^N$:

$$\Xi = \sum_{N=0}^{\infty} \frac{V^N e^{\beta\mu N}}{N!\,\Lambda^{3N}} \int\!\!\cdots\!\!\int e^{-\beta U}\, d\mathbf s^N$$

Why bother: the volume now appears explicitly as $V^N$, and the integration region (the unit cube) never changes. This is convenient whenever $N$ or $V$ changes (GCMC, and constant-pressure simulations).

### 2.3 Probability of a microstate

A **microstate** is "$N$ molecules at scaled positions $\mathbf s_1, \dots, \mathbf s_N$". Its probability is the term inside the sum:

$$p(\mathbf s^N; N) \propto \frac{V^N e^{\beta\mu N}}{N!\,\Lambda^{3N}}\; e^{-\beta U(\mathbf s^N)}$$

(The lecture writes this with $\mathbf r_1, \dots, \mathbf r_N$; the $V^N$ factor belongs with the scaled coordinates $\mathbf s$. In real coordinates the same probability has no $V^N$, because the $V^N$ is absorbed into $d\mathbf r^N$.)

This is the target distribution that Metropolis must sample. **The normalization $\Xi$ is never needed**, because acceptance rules use only ratios.

---

## 3. Why GCMC: adsorption in a slit pore

### 3.1 The problem

**Adsorption** is gas molecules accumulating inside the pores of a solid (activated carbon, zeolites, MOFs). The lecture considers a **slit pore**: two parallel solid walls a distance $d$ apart, with gas between them. The question is: how many molecules sit in the pore when the outside gas is at a given pressure and temperature?

### 3.2 The NVT approach (left sketch in the lecture)

Put the whole adsorbent particle, with its pore, inside a large box of gas and run NVT. Two difficulties:

1. **The entire adsorbent has to be simulated.** Real adsorbent particles are macroscopic; that is far too many atoms.
2. **The external gas has to be simulated too.** Most of the computer time goes into gas far from the pore, which we do not care about.

### 3.3 The GCMC approach (right sketch)

Simulate **only the pore** (a small box between the two walls). Replace the external gas by its **chemical potential $\mu$**: the pore is in equilibrium with the outside gas exactly when both have the same $\mu$ and $T$. Molecules are inserted into and deleted from the pore so that $N$ settles at the value that equilibrium demands.

**GCMC is computationally much cheaper**, because the reservoir is never simulated: it is represented by a single number, $\mu$.

### 3.4 What you get: the adsorption isotherm

Run one GCMC simulation per outside pressure. Each run gives $\langle N\rangle$ in the pore. Plotting $\langle N\rangle$ against pressure at fixed $T$ gives the **adsorption isotherm**, which can be compared directly with experiments. The walls act on the gas through a wall–molecule potential added to $U$; strongly attractive walls make insertions easier (negative $\Delta U$), so the pore holds far more molecules than the same volume of bulk gas.

The same idea applies to any region in equilibrium with a gas or vapour at known conditions, for example vapour in contact with a surface or a substrate.

---

## 4. The three Monte Carlo moves and their acceptance rules

All three follow the detailed-balance recipe from the Metropolis notes:

$$p(o)\,\alpha(o\to n)\,\text{acc}(o\to n) = p(n)\,\alpha(n\to o)\,\text{acc}(n\to o)$$

$$\Rightarrow\quad \text{acc}(o\to n) = \min\left\{1,\; \frac{p(n)\,\alpha(n\to o)}{p(o)\,\alpha(o\to n)}\right\}$$

### 4.1 Move 1: displacement (as in NVT)

$N$ does not change, so every factor in $p$ except $e^{-\beta U}$ cancels. The proposal (pick a molecule, shift it uniformly within $\pm\Delta r_m$) is symmetric. As derived earlier:

$$\text{acc}(o\to n) = \min\left\{1,\; e^{-\beta\Delta U}\right\}$$

### 4.2 Move 2: insertion ($N \to N+1$)

**Old and new probabilities:**

$$p(o) \propto \frac{V^N e^{\beta\mu N}}{N!\,\Lambda^{3N}}\, e^{-\beta U(o)}, \qquad p(n) \propto \frac{V^{N+1} e^{\beta\mu(N+1)}}{(N+1)!\,\Lambda^{3(N+1)}}\, e^{-\beta U(n)}$$

**Proposal probabilities.** If an insertion is attempted as often as a deletion, then $\alpha(o\to n) = \alpha(n\to o)$ and the $\alpha$'s cancel (Section 4.4 explains the bookkeeping).

**The ratio, factor by factor:**

| Factor | New / old | Result |
|---|---|---|
| volume | $V^{N+1}/V^N$ | $V$ |
| chemical potential | $e^{\beta\mu(N+1)}/e^{\beta\mu N}$ | $e^{\beta\mu}$ |
| identical molecules | $N!/(N+1)!$ | $1/(N+1)$ |
| momenta | $\Lambda^{3N}/\Lambda^{3(N+1)}$ | $1/\Lambda^3$ |
| energy | $e^{-\beta U(n)}/e^{-\beta U(o)}$ | $e^{-\beta\Delta U}$ |

$$\frac{\text{acc}(o\to n)}{\text{acc}(n\to o)} = \frac{p(n)}{p(o)} = \frac{V e^{\beta\mu}}{(N+1)\,\Lambda^3}\, e^{-\beta\Delta U}, \qquad \Delta U = U(n) - U(o)$$

$$\boxed{\;\text{acc}_{\text{ins}}(N\to N+1) = \min\left\{1,\; \frac{V e^{\beta\mu}}{(N+1)\,\Lambda^3}\, e^{-\beta\Delta U}\right\}\;}$$

Here $\Delta U$ is the interaction energy of the new molecule with all the others (plus the change in the tail correction, Section 6.3).

### 4.3 Move 3: deletion ($N \to N-1$)

$$p(o) \propto \frac{V^N e^{\beta\mu N}}{N!\,\Lambda^{3N}}\, e^{-\beta U(o)}, \qquad p(n) \propto \frac{V^{N-1} e^{\beta\mu(N-1)}}{(N-1)!\,\Lambda^{3(N-1)}}\, e^{-\beta U(n)}$$

| Factor | New / old | Result |
|---|---|---|
| volume | $V^{N-1}/V^N$ | $1/V$ |
| chemical potential | $e^{\beta\mu(N-1)}/e^{\beta\mu N}$ | $e^{-\beta\mu}$ |
| identical molecules | $N!/(N-1)!$ | $N$ |
| momenta | $\Lambda^{3N}/\Lambda^{3(N-1)}$ | $\Lambda^3$ |
| energy | | $e^{-\beta\Delta U}$ |

$$\boxed{\;\text{acc}_{\text{del}}(N\to N-1) = \min\left\{1,\; \frac{N\Lambda^3}{V e^{\beta\mu}}\, e^{-\beta\Delta U}\right\}\;}, \qquad \Delta U = U(N-1) - U(N)$$

Removing molecule $j$ removes all its interactions, so $\Delta U = -(\text{energy of } j \text{ with all the others})$. Deleting a strongly bound molecule means $\Delta U > 0$, which is hard; deleting a molecule in an unfavourable spot is easy.

### 4.4 Two bookkeeping conventions, one result

The $(N+1)$ in the insertion rule can be placed in two different spots; both give the same acceptance rule.

| | Lecture's convention | Alternative convention (earlier Metropolis notes) |
|---|---|---|
| Probability $p$ | includes $1/N!$ (numbered molecules) and $V^N$ (scaled coordinates) | no $1/N!$ (identical molecules), real coordinates |
| Insertion proposal | "insert" chosen with the same probability as "delete" | random point: $\alpha = 1/V$ |
| Deletion proposal | same | pick one of $N+1$ molecules: $\alpha = 1/(N+1)$ |
| Where $V/(N+1)$ comes from | the ratio $p(n)/p(o)$ | the ratio $\alpha(n\to o)/\alpha(o\to n)$ |
| Final rule | $\min\{1, \frac{zV}{N+1}e^{-\beta\Delta U}\}$ | the same |

The code is identical either way: insert at a uniformly random point, delete a uniformly random molecule.

---

## 5. Chemical potential, activity and fugacity

The acceptance rules contain the combination $e^{\beta\mu}/\Lambda^3$. This section gives it two physical meanings.

### 5.1 Activity: split $\mu$ into ideal and excess parts

For an ideal gas (no interactions), $\beta\mu^{\text{ig}} = \ln(\rho\Lambda^3)$. A real fluid has an extra part from interactions, the **excess chemical potential** $\mu^{\text{ex}}$:

$$\beta\mu = \ln(\rho\Lambda^3) + \beta\mu^{\text{ex}}$$

Exponentiate and divide by $\Lambda^3$:

$$\boxed{\;z \equiv \frac{e^{\beta\mu}}{\Lambda^3} = \rho\, e^{\beta\mu^{\text{ex}}}\;}$$

$z$ is called the **activity**. It has units of number density (molecules per volume).

**Meaning:** for an ideal gas, $\mu^{\text{ex}} = 0$, so $z = \rho$. The activity is **the density an ideal gas would have at this chemical potential**, and **$zV$ is the number of molecules the reservoir "wants" in the box**.

With $z$, the acceptance rules become:

$$\text{acc}_{\text{ins}} = \min\left\{1,\; \frac{zV}{N+1}\, e^{-\beta\Delta U}\right\}, \qquad \text{acc}_{\text{del}} = \min\left\{1,\; \frac{N}{zV}\, e^{-\beta\Delta U}\right\}$$

**Reading the rules:** if the box holds fewer molecules than the reservoir wants ($N+1 < zV$), insertions are favoured; if it holds more, deletions are favoured. The energy factor then adjusts this: favourable spots (attractive walls, $\Delta U < 0$) pull extra molecules in.

**Numerical example:** $zV = 50$, $N = 40$.

- Insert where $\Delta U = +0.5\,k_BT$: acc $= \frac{50}{41}e^{-0.5} = 0.740$.
- Delete a molecule bound by $2\,k_BT$ ($\Delta U = +2\,k_BT$): acc $= \frac{40}{50}e^{-2} = 0.108$.

### 5.2 Fugacity: a "corrected pressure"

Experiments specify the **pressure** of the outside gas, not $\mu$. Fugacity connects the two.

**Definition.** The fugacity $f$ is defined by

$$k_BT\ln\left(\frac{f}{P}\right) = \mu(T, P) - \mu^{\text{ig}}(T, P)$$

where $\mu^{\text{ig}}(T, P)$ is the chemical potential an ideal gas would have at the same $T$ and $P$. For an ideal gas $f = P$; for a real gas $f/P$ (the fugacity coefficient) measures the deviation.

**Derivation of the link to activity.**

1. For an ideal gas, $\rho = P/(k_BT)$, so $\beta\mu^{\text{ig}} = \ln(\rho\Lambda^3) = \ln\left(\dfrac{P\Lambda^3}{k_BT}\right)$.
2. Divide the definition by $k_BT$:
$$\ln\left(\frac{f}{P}\right) = \beta\mu - \ln\left(\frac{P\Lambda^3}{k_BT}\right)$$
3. Exponentiate:
$$\frac{f}{P} = e^{\beta\mu}\,\frac{k_BT}{P\Lambda^3}$$
4. Cancel $P$ and rearrange:
$$\boxed{\;z = \frac{e^{\beta\mu}}{\Lambda^3} = \frac{f}{k_BT} = \beta f\;}$$

So the acceptance rules can be written with fugacity:

$$\text{acc}_{\text{ins}} = \min\left\{1,\; \frac{\beta f V}{N+1}\, e^{-\beta\Delta U}\right\}, \qquad \text{acc}_{\text{del}} = \min\left\{1,\; \frac{N}{\beta f V}\, e^{-\beta\Delta U}\right\}$$

**In practice:** take the outside pressure $P$, get $f$ from an equation of state of the bulk gas (at low pressure $f \approx P$), set $z = \beta f$, and run GCMC.

### 5.3 Numbers for a real gas

Methane at 300 K and 1 bar (low enough that $f \approx P$):

$$z = \frac{P}{k_BT} = \frac{10^5}{1.381\times10^{-23}\times300} = 2.41\times10^{25}\ \text{m}^{-3} = 0.0241\ \text{nm}^{-3}$$

For a $3\times3\times3$ nm box ($V = 27$ nm³): $zV = 0.65$. In bulk gas the box would hold less than one molecule on average. At 10 bar, $zV = 6.5$. In a strongly attractive pore of the same volume, $\langle N\rangle$ can be many times larger: that excess is the adsorption.

### 5.4 Three equivalent ways to write the input

| Specify | Insertion acceptance | Deletion acceptance |
|---|---|---|
| chemical potential $\mu$ | $\min\{1, \frac{Ve^{\beta\mu}}{(N+1)\Lambda^3}e^{-\beta\Delta U}\}$ | $\min\{1, \frac{N\Lambda^3}{Ve^{\beta\mu}}e^{-\beta\Delta U}\}$ |
| activity $z$ | $\min\{1, \frac{zV}{N+1}e^{-\beta\Delta U}\}$ | $\min\{1, \frac{N}{zV}e^{-\beta\Delta U}\}$ |
| fugacity $f$ | $\min\{1, \frac{\beta fV}{N+1}e^{-\beta\Delta U}\}$ | $\min\{1, \frac{N}{\beta fV}e^{-\beta\Delta U}\}$ |

This is why step 1 of the GCMC algorithm says "specify $\mu$ (or $f$, $z$)".

---

## 6. How the insertion and deletion moves are implemented

### 6.1 Insertion

1. Choose a **random point uniformly** inside the simulation box: $\mathbf r_{\text{new}} = (\mathbf u - 0.5)\,L$, where $\mathbf u$ is three uniform random numbers in $[0, 1)$.
2. Compute $\Delta U$: the energy of a molecule at $\mathbf r_{\text{new}}$ with all $N$ existing molecules (plus the tail-correction change).
3. Accept with $\min\{1, \frac{zV}{N+1}e^{-\beta\Delta U}\}$. If accepted, store the new coordinates in position $N+1$ of the coordinate array and update $N \leftarrow N+1$.

### 6.2 Deletion: the coordinate-array trick

The coordinates are stored in arrays (`rx`, `ry`, `rz`) whose first $N$ entries are in use:

```
index:   1   2   3   4   5   6   7   8   9  10 ...
        [x] [x] [x] [x] [x] [x] [x] [x] [ ] [ ] ...      N = 8
```

1. **Choose a molecule $j$ at random** (each of the $N$ with probability $1/N$):

```
        [x] [x] [x] [x] [X] [x] [x] [x] [ ] [ ]          j = 5
```

2. **Compute $\Delta U$** $= -(\text{energy of } j \text{ with all others})$ (plus the tail-correction change).

3. **If accepted**, the slot $j$ becomes empty. Instead of shifting every molecule after $j$ one place left, **copy the last molecule into slot $j$** and shorten the list:

```
before:  [x] [x] [x] [x] [ ] [x] [x] [x]                 hole at j = 5
         copy the last one (index N = 8) into slot j
after:   [x] [x] [x] [x] [x] [x] [x]                     N = 7
```

```fortran
rx(j) = rx(N)
ry(j) = ry(N)
rz(j) = rz(N)
N = N - 1
```

This takes the same time whatever $N$ is, and keeps the occupied part of the array contiguous. The order of molecules does not matter because they are identical.

In Python (indices start at 0, so the last molecule is `pos[N-1]`):

```python
pos[j] = pos[N - 1]
pos = pos[:N - 1]
N -= 1
```

### 6.3 The tail correction changes when $N$ changes

In NVT the long-range (tail) correction is a constant and never enters $\Delta U$. In GCMC it depends on $N$. For Lennard-Jones (reduced units):

$$U_{\text{tail}}(N) = \frac{8}{3}\pi\,\frac{N^2}{V}\left[\frac{1}{3}r_c^{-9} - r_c^{-3}\right]$$

so an insertion changes it by

$$\Delta U_{\text{tail}} = U_{\text{tail}}(N+1) - U_{\text{tail}}(N) = \frac{8}{3}\pi\,\frac{2N+1}{V}\left[\frac{1}{3}r_c^{-9} - r_c^{-3}\right]$$

Example: $N = 68$, $V = 512$, $r_c = 3$: $\Delta U_{\text{tail}} = -0.083\,\varepsilon$. At $T^* = 2$ this changes the insertion acceptance by about 4%. It must be included in $\Delta U$ for insertions and deletions.

---

## 7. The complete GCMC algorithm (lecture page 3)

### 7.1 The steps

1. **Specify** $\mu$ (or $f$, or $z$), $V$, $T$.
2. **Read** initial coordinates from a file, or **generate** them (an empty box is fine in GCMC).
3. **Specify** the length of the simulation, the maximum displacement $\Delta r_m$, and the random number seed.
4. **Perform** `nequil + ncycles` cycles of moves:

```fortran
do l = 1, (nequil + ncycles)
   cycle_length = max(10, N)
   do m = 1, cycle_length
      ! choose the move type
      xrand = ran2(iseed)
      if (xrand < 0.33) then
         ! perform translation (displacement) move
      else if (xrand < 0.67) then
         ! perform insertion move
      else
         ! perform deletion move
      endif
      ! sample properties of interest
   enddo
   ! optimise Delta r_m during the equilibration phase
enddo
```

### 7.2 Line by line

| Line | Purpose |
|---|---|
| `cycle_length = max(10, N)` | one cycle is about $N$ moves, so each molecule gets about one attempt per cycle. The `max(10, …)` keeps at least 10 attempts per cycle when the box is nearly empty (for example at the start, with $N = 0$). |
| `xrand = ran2(iseed)` | uniform random number in $[0, 1)$ (`ran2` is the Numerical Recipes generator) |
| three branches | each move type is chosen with probability $1/3$ |
| equal probability for insertion and deletion | this is what makes $\alpha(o\to n) = \alpha(n\to o)$ in Section 4.2. If they were chosen with different probabilities, the ratio would have to appear in the acceptance rules. |
| "sample properties of interest" | record $N$, $U$, density profiles, … after **every** move, accepted or rejected |
| "optimise $\Delta r_m$" | adjust the maximum displacement towards about 40–50% acceptance, **only during equilibration**, then keep it fixed |

### 7.3 Edge cases

- **Deletion when $N = 0$:** nothing to delete, so the move is rejected (the state is still sampled).
- **Displacement when $N = 0$:** nothing to move; count it as rejected.
- **Cap $\Delta r_m$ at $L/2$.** At low density the displacement acceptance stays high, so automatic tuning keeps increasing $\Delta r_m$. A step longer than half the box means nothing under periodic boundaries, and if $\Delta r_m$ grows without limit, floating-point precision is lost (this actually happened while testing the code in Section 9).

---

## 8. Sanity checks

### 8.1 Ideal gas: a Poisson distribution

With no interactions ($U = 0$), the configurational integral is $V^N$, so

$$\Xi_{\text{ig}} = \sum_{N=0}^{\infty}\frac{(zV)^N}{N!} = e^{zV}, \qquad P(N) = \frac{e^{-zV}(zV)^N}{N!}$$

where $P(N)$ is the probability of finding $N$ molecules. This is a **Poisson distribution**: $\langle N\rangle = zV$ and the variance is also $zV$.

Test (Section 9 code, interactions switched off, $z^* = 0.05$, $V^* = 512$, so $zV = 25.6$):

| | Simulation | Exact |
|---|---|---|
| $\langle N\rangle$ | 25.60 | 25.6 |
| variance of $N$ | 26.39 | 25.6 |

### 8.2 Low-density limit for a real fluid

At low density, $\beta\mu^{\text{ex}} \approx 2B_2\rho$, where $B_2$ is the second virial coefficient. Then $z = \rho e^{2B_2\rho}$. For Lennard-Jones at $T^* = 2$, $B_2^* = -1.314$ (negative: attraction wins), so the fluid is **denser than an ideal gas** at the same $z$.

### 8.3 Widom insertion (independent check)

The excess chemical potential can be measured in an **NVT** simulation by attempting ghost insertions at random points and averaging the Boltzmann factor (Widom's method):

$$\beta\mu^{\text{ex}} = -\ln\left\langle e^{-\beta\Delta U}\right\rangle$$

Then $z = \rho e^{\beta\mu^{\text{ex}}}$ must reproduce the activity that GCMC used. See Section 9.2 for the result.

---

## 9. Working Python code and results

### 9.1 The code (reduced units, Lennard-Jones, follows Section 7)

```python
import numpy as np

T        = 2.0          # temperature T*
z        = 0.05         # activity z* = z sigma^3
L        = 8.0          # box edge L*
rc       = 3.0          # cutoff
dr_max   = 0.5          # maximum displacement (tuned during equilibration)
n_equil  = 2000         # equilibration cycles
n_cycles = 5000         # production cycles

beta, V = 1.0 / T, L**3
rng = np.random.default_rng(1)

def lj(r2):
    s3 = (1.0 / r2) ** 3
    return 4.0 * (s3 * s3 - s3)

def tail(N):            # U_tail depends on N in GCMC
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
        dr_max = min(dr_max, L / 2)

N_samples = np.array(N_samples)
print("<N> =", N_samples.mean(), " rho* =", N_samples.mean() / V)
```

### 9.2 Results ($T^* = 2.0$, $L^* = 8$, $V^* = 512$, $r_c^* = 3$)

| $z^*$ | Ideal gas $\rho^* = z^*$ | GCMC $\rho^*$ | Second-virial estimate | Insertion acceptance |
|---|---|---|---|---|
| 0.01 | 0.0100 | 0.0101 ± 0.0001 | 0.0103 | 0.82 |
| 0.03 | 0.0300 | 0.0323 ± 0.0004 | 0.0327 | 0.83 |
| 0.10 | 0.1000 | 0.1343 ± 0.0009 | 0.147 | 0.60 |

- At low $z$ the fluid behaves almost ideally ($\rho \approx z$).
- As $z$ increases, attractions make the fluid denser than an ideal gas.
- The second-virial estimate works at low density and overshoots at $\rho^* \approx 0.13$, where higher virial terms matter.

**Widom check at $z^* = 0.10$:** an independent NVT simulation at $\rho^* = 0.1348$ with ghost insertions gave $\beta\mu^{\text{ex}} = -0.300 \pm 0.002$, so $z^* = \rho^* e^{\beta\mu^{\text{ex}}} = 0.0998$. That matches the $z^* = 0.10$ used in GCMC.

---

## 10. Practical issues

- **Insertions fail in dense phases.** In a dense liquid ($\rho^* \approx 0.8$), almost every random point overlaps an existing molecule, so $\Delta U$ is huge and insertions are almost never accepted; deletions of well-bound molecules are also rare. The number of molecules then changes very slowly. Special techniques (configurational-bias insertion, cavity-biased insertion) exist for this.
- **Insertion and deletion must be attempted equally often** (or the ratio must be put in the acceptance rule).
- **Equilibration** is judged by watching $N$ (not only $U$) settle down.
- **Sampling:** record after every move, including rejected moves and failed deletions at $N = 0$.
- **Finite size:** the box must be larger than $2r_c$ and large enough that $\langle N\rangle$ is not tiny.
- **Molecules** (not single atoms) need random orientations on insertion, and the orientation probability must be accounted for in $\alpha$.

---

## 11. Corrections and clarifications to the lecture notes

1. **The same letter $f$** is used for the microstate probability and for the fugacity. These notes use $p$ for the probability.
2. **Microstate probability (page 1):** the factor $V^N$ goes with the scaled coordinates $\mathbf s$; in real coordinates $\mathbf r$ the probability has no $V^N$.
3. **Insertion rule (page 1):** the line before the activity form shows $e^{\beta\mu N}$ in the numerator. It should be $e^{\beta\mu}$ (one molecule is added), as in the line above it.
4. **Final deletion line (page 2):** the lecture labels it acc$(o\to n)$; here $o$ has $N$ molecules and $n$ has $N-1$.

---

## 12. Practice questions (with answers)

1. **Why does $1/(N+1)$ appear in the insertion rule?**
   From $N!/(N+1)!$ in the ratio of probabilities: identical molecules (lecture convention). Equivalently, from the chance $1/(N+1)$ of picking the newly inserted molecule in the reverse deletion (Section 4.4).

2. **For an ideal gas with $zV = 5$, what are $\langle N\rangle$ and its variance?**
   Both are 5 (Poisson).

3. **A GCMC run of an ideal gas gives $\langle N\rangle = 12.3$ in a box of 100 nm³ at 300 K. What pressure was simulated?**
   $z = 0.123$ nm⁻³ $= 1.23\times10^{26}$ m⁻³, so $P = zk_BT = 1.23\times10^{26}\times1.381\times10^{-23}\times300 = 5.1\times10^5$ Pa $= 5.1$ bar.

4. **The box holds $N = 10$ molecules and $zV = 10$. A deletion has $\Delta U = 0$. What is the acceptance?**
   $\min\{1, 10/10\} = 1$.

5. **Why must the tail correction be included in $\Delta U$ for insertions but not for displacements?**
   It depends on $N$, which changes only in insertions and deletions.

6. **What happens if insertions are attempted twice as often as deletions but the acceptance rules are unchanged?**
   Detailed balance is broken: the $\alpha$ ratio of 2 is missing, so the box ends up with too many molecules.

---

## 13. Summary

- **GCMC** simulates an open system at fixed $\mu$, $V$, $T$. $N$ fluctuates, and the reservoir is represented by one number ($\mu$, $z$ or $f$), never simulated.
- **Grand partition function:** $\Xi = \sum_N \frac{e^{\beta\mu N}}{N!\Lambda^{3N}}\int e^{-\beta U}d\mathbf r^N$. Each term is the probability weight of states with $N$ molecules.
- **Three moves**, chosen with probability 1/3 each:
  - displacement: $\min\{1, e^{-\beta\Delta U}\}$
  - insertion: $\min\{1, \frac{zV}{N+1}e^{-\beta\Delta U}\}$
  - deletion: $\min\{1, \frac{N}{zV}e^{-\beta\Delta U}\}$
- **Activity** $z = e^{\beta\mu}/\Lambda^3 = \rho e^{\beta\mu^{\text{ex}}} = \beta f$. For an ideal gas $z = \rho = P/k_BT$, and $zV$ is the number of molecules the reservoir wants in the box.
- **Implementation:** insert at a uniformly random point; delete a random molecule by copying the last one into its slot; include the $N$-dependent tail correction; record after every move.
- **Main application:** adsorption isotherms. Simulate only the pore, at the $\mu$ of the outside gas.
