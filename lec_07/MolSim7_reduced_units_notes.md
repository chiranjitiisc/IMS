# Lecture 7: Reduced Units

*Detailed notes based on MolSim7.pdf ("Reduced Units", 27 January 2022). Every number for argon below was recomputed and checked.*

---

## Contents

0. Notation
1. Why reduced units?
2. The Lennard-Jones potential and its two parameters
3. Making the potential dimensionless
4. Why one simulation covers every Lennard-Jones material
5. Building every unit from σ, ε and m
6. Worked conversions for argon
7. The phase diagram of the Lennard-Jones fluid in reduced units
8. Reduced units in practice
9. Python converter
10. Common mistakes
11. Practice questions (with answers)
12. Summary

---

## 0. Notation

| Symbol | Meaning | SI unit |
|---|---|---|
| $r$ | distance between two atoms | m |
| $u(r)$ | pair potential energy of two atoms | J |
| $\sigma$ | Lennard-Jones size parameter (distance where $u = 0$) | m |
| $\varepsilon$ | Lennard-Jones energy parameter (depth of the well) | J |
| $m$ | mass of one atom | kg |
| $k_B$ | Boltzmann's constant, $1.380649\times10^{-23}$ J/K | J/K |
| $U$ | total potential energy of the system | J |
| $L$ | a length (for example the box edge) | m |
| $P$ | pressure | Pa |
| $\rho$ | number density, atoms per unit volume ($N/V$) | m⁻³ |
| $T$ | temperature | K |
| $t$ | time | s |
| $\gamma$ | surface tension (energy per unit area) | N/m = J/m² |
| $\tau$ | the reduced unit of time, $\sqrt{m\sigma^2/\varepsilon}$ | s |
| $X^*$ | any quantity $X$ written in reduced units (a pure number) | none |

---

## 1. Why reduced units?

The lecture gives two reasons.

### Reason 1: convenient numbers

In SI units, the numbers for a single argon atom are awkward:

| Quantity | SI value |
|---|---|
| $\varepsilon$ | $1.654\times10^{-21}$ J |
| $\sigma$ | $3.405\times10^{-10}$ m |
| $m$ | $6.634\times10^{-26}$ kg |

Every calculation then drags along large powers of ten, which makes the code hard to read and mistakes easy to make. The lecture's point is to choose units so that the numbers are **neither too large nor too small**. In reduced units, typical values are between about 0.01 and 10: density 0.8, temperature 1.0, energy per atom −5.5.

### Reason 2: one simulation for every Lennard-Jones material

The properties of a Lennard-Jones fluid depend on $\sigma$ and $\varepsilon$, which differ from substance to substance. But if every quantity is measured in units built from $\sigma$ and $\varepsilon$, the dependence on the particular substance disappears. **One simulation in reduced units gives the results for every Lennard-Jones material at once**; you only convert the answer back with that material's $\sigma$ and $\varepsilon$. Section 4 shows why.

---

## 2. The Lennard-Jones potential and its two parameters

$$u(r) = 4\varepsilon\left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right]$$

- The $(\sigma/r)^{12}$ term is **repulsion** (electron clouds overlapping at short range).
- The $(\sigma/r)^{6}$ term is **attraction** (van der Waals dispersion).

### Key points on the curve (the sketch in the lecture)

| Point | Where | Value of $u$ | What it tells you |
|---|---|---|---|
| Zero crossing | $r = \sigma$ | $0$ | $\sigma$ is the effective **size** of the atom |
| Minimum | $r = 2^{1/6}\sigma = 1.1225\,\sigma$ | $-\varepsilon$ | $\varepsilon$ is the **strength** of the attraction |
| Large $r$ | $r \to \infty$ | $\to 0$ from below | atoms far apart do not interact |
| Small $r$ | $r < \sigma$ | rises steeply | atoms cannot overlap |

**Where the minimum comes from.** At the bottom of the well the slope is zero:

$$\frac{du}{dr} = 4\varepsilon\left[-12\,\frac{\sigma^{12}}{r^{13}} + 6\,\frac{\sigma^{6}}{r^{7}}\right] = 0 \;\Rightarrow\; r^6 = 2\sigma^6 \;\Rightarrow\; r_{\min} = 2^{1/6}\sigma$$

At $r_{\min}$, $(\sigma/r)^6 = \tfrac12$ and $(\sigma/r)^{12} = \tfrac14$, so $u(r_{\min}) = 4\varepsilon\left[\tfrac14 - \tfrac12\right] = -\varepsilon$.

---

## 3. Making the potential dimensionless

Define a **reduced energy** and a **reduced distance**:

$$u^* = \frac{u}{\varepsilon}, \qquad r^* = \frac{r}{\sigma}$$

Divide the Lennard-Jones formula by $\varepsilon$ and use $\sigma/r = 1/r^*$:

$$\frac{u}{\varepsilon} = 4\left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right] \quad\Longrightarrow\quad \boxed{\;u^* = 4\left[\left(\frac{1}{r^*}\right)^{12} - \left(\frac{1}{r^*}\right)^{6}\right]\;}$$

**$\sigma$ and $\varepsilon$ have disappeared.** In reduced units there is only **one** Lennard-Jones curve, the same for argon, krypton, methane or any other Lennard-Jones material:

- zero at $r^* = 1$,
- minimum $u^* = -1$ at $r^* = 1.1225$,
- $u^*(1.5) = -0.3203$, $u^*(3) = -0.0055$.

This is exactly what writing $\sigma = \varepsilon = 1$ in a code means.

---

## 4. Why one simulation covers every Lennard-Jones material

### 4.1 Equilibrium (Monte Carlo) properties

Everything in equilibrium statistical mechanics comes from the Boltzmann factor $e^{-U/k_BT}$. For a Lennard-Jones system:

$$\frac{U}{k_BT} = \frac{1}{k_BT}\sum_{i<j} 4\varepsilon\left[\left(\frac{\sigma}{r_{ij}}\right)^{12} - \left(\frac{\sigma}{r_{ij}}\right)^{6}\right] = \frac{\varepsilon}{k_BT}\sum_{i<j} u^*(r^*_{ij}) = \frac{U^*}{T^*}$$

where $T^* = k_BT/\varepsilon$ and $U^* = U/\varepsilon$.

The positions enter only through $r^*_{ij} = r_{ij}/\sigma$, and the box only through $L^* = L/\sigma$, so the number of atoms per $\sigma^3$, $\rho^* = \rho\sigma^3$, is all that matters about the volume.

**Conclusion:** every reduced equilibrium property ($U^*/N$, $P^*$, the structure as a function of $r^*$, …) depends **only on $T^*$ and $\rho^*$** (and on $N$ for small systems). Two Lennard-Jones substances at the same $T^*$ and $\rho^*$ have identical reduced properties. This is the **law of corresponding states**.

The Metropolis acceptance rule shows this directly:

$$\text{acc} = \min\left\{1,\; e^{-\Delta U/k_BT}\right\} = \min\left\{1,\; e^{-\Delta U^*/T^*}\right\}$$

That is why the Monte Carlo codes in the earlier notes (with $\sigma = \varepsilon = k_B = 1$) were already in reduced units.

### 4.2 Dynamics (molecular dynamics), and where the time unit comes from

Newton's equation for atom $i$ is

$$m\frac{d^2 \mathbf r_i}{dt^2} = -\frac{\partial U}{\partial \mathbf r_i}$$

Substitute $\mathbf r = \sigma\,\mathbf r^*$, $U = \varepsilon\,U^*$ and $t = \tau\,t^*$, where $\tau$ is a time unit still to be chosen:

$$m\,\frac{\sigma}{\tau^2}\,\frac{d^2\mathbf r_i^*}{dt^{*2}} = -\frac{\varepsilon}{\sigma}\,\frac{\partial U^*}{\partial \mathbf r_i^*} \quad\Longrightarrow\quad \frac{d^2\mathbf r_i^*}{dt^{*2}} = -\frac{\varepsilon\,\tau^2}{m\,\sigma^2}\,\frac{\partial U^*}{\partial \mathbf r_i^*}$$

Choose $\tau$ so that the constant in front is exactly 1:

$$\boxed{\;\tau = \sqrt{\frac{m\sigma^2}{\varepsilon}}\;}\qquad\Longrightarrow\qquad \frac{d^2\mathbf r_i^*}{dt^{*2}} = -\frac{\partial U^*}{\partial \mathbf r_i^*}$$

The equation of motion is now free of $m$, $\sigma$ and $\varepsilon$, so the **dynamics** in reduced units is also the same for every Lennard-Jones material. This is the origin of the time row in the table below.

---

## 5. Building every unit from σ, ε and m

Choose three **base units**:

| Base quantity | Reduced unit |
|---|---|
| length | $\sigma$ |
| energy | $\varepsilon$ |
| mass | $m$ |

and use $k_B$ to turn temperature into an energy. Every other unit follows from dimensions.

### 5.1 The conversion table (lecture table, with derivations and argon values)

| Property | Real quantity | Reduced quantity | Unit | How the unit is built | Argon: 1 reduced unit = |
|---|---|---|---|---|---|
| Energy | $U$ | $U^* = U/\varepsilon$ | $\varepsilon$ | base unit | $1.654\times10^{-21}$ J = 0.996 kJ/mol = 0.0103 eV |
| Length | $L$ | $L^* = L/\sigma$ | $\sigma$ | base unit | 3.405 Å |
| Pressure | $P$ | $P^* = P\sigma^3/\varepsilon$ | $\varepsilon/\sigma^3$ | pressure = force/area = energy/volume | 41.9 MPa |
| Number density | $\rho$ | $\rho^* = \rho\sigma^3$ | $1/\sigma^3$ | number/volume | $2.533\times10^{28}$ m⁻³ |
| (Mass density) | $\rho_m$ | $\rho_m\sigma^3/m$ | $m/\sigma^3$ | mass/volume | 1680 kg/m³ |
| Temperature | $T$ | $T^* = \dfrac{T}{\varepsilon/k_B} = \dfrac{k_BT}{\varepsilon}$ | $\varepsilon/k_B$ | $k_BT$ is an energy | 119.8 K |
| Time | $t$ | $t^* = t\big/\sqrt{m\sigma^2/\varepsilon}$ | $\tau = \sqrt{m\sigma^2/\varepsilon}$ | energy = mass × length²/time² | 2.156 ps |
| Surface tension | $\gamma$ | $\gamma^* = \gamma\sigma^2/\varepsilon$ | $\varepsilon/\sigma^2$ | energy/area | 14.3 mN/m |
| (Force) | $F$ | $F\sigma/\varepsilon$ | $\varepsilon/\sigma$ | energy/length | $4.86\times10^{-12}$ N |
| (Velocity) | $v$ | $v\big/\sqrt{\varepsilon/m}$ | $\sqrt{\varepsilon/m} = \sigma/\tau$ | length/time | 158 m/s |

Rows in brackets are not in the lecture but follow the same logic.

### 5.2 The rule for every row

$$\text{reduced value} = \frac{\text{real value}}{\text{unit}}, \qquad \text{real value} = \text{reduced value}\times\text{unit}$$

### 5.3 Why temperature is quoted as $\varepsilon/k_B$ in kelvin

$\varepsilon$ is an energy (J), but $\varepsilon/k_B$ has units of kelvin. Quoting $\varepsilon/k_B = 119.8$ K means "the well depth equals the thermal energy $k_BT$ at 119.8 K". So $T^* = 1$ means $k_BT = \varepsilon$: thermal energy equals the well depth.

### 5.4 A consistency check: the ideal gas law

$P = \rho k_BT$. Multiply both sides by $\sigma^3/\varepsilon$:

$$\frac{P\sigma^3}{\varepsilon} = \rho\sigma^3\cdot\frac{k_BT}{\varepsilon} \quad\Longrightarrow\quad P^* = \rho^*\,T^*$$

The ideal gas law has the same form in reduced units, as it should.

---

## 6. Worked conversions for argon

**Argon parameters (from the lecture):** $\sigma = 3.405$ Å, $\varepsilon/k_B = 119.8$ K. Atomic mass 39.948 u, so $m = 39.948\times1.6605\times10^{-27} = 6.634\times10^{-26}$ kg.

### 6.1 The three examples in the lecture

**Temperature, $T^* = 1.0$:**

$$T = T^*\times\frac{\varepsilon}{k_B} = 1.0\times119.8\text{ K} = \mathbf{119.8\ K}$$

**Density, $\rho^* = 1.0$.** First the number density, then multiply by the mass of one atom:

$$\sigma^3 = (3.405\times10^{-10}\text{ m})^3 = 3.948\times10^{-29}\text{ m}^3$$

$$\rho = \frac{\rho^*}{\sigma^3} = \frac{1.0}{3.948\times10^{-29}} = 2.533\times10^{28}\text{ atoms/m}^3$$

$$\rho_m = \rho\times m = 2.533\times10^{28}\times6.634\times10^{-26} = \mathbf{1680\ kg/m^3}$$

**Pressure, $P^* = 1.0$:**

$$\varepsilon = 119.8\times1.380649\times10^{-23} = 1.654\times10^{-21}\text{ J}$$

$$P = P^*\times\frac{\varepsilon}{\sigma^3} = \frac{1.654\times10^{-21}}{3.948\times10^{-29}} = 4.19\times10^{7}\text{ Pa} = \mathbf{41.9\ MPa}$$

### 6.2 Going the other way

| Real argon state | Reduced value |
|---|---|
| Normal boiling point, 87.3 K | $T^* = 87.3/119.8 = 0.729$ |
| Liquid at the boiling point, about 1395 kg/m³ | $\rho^* = 1395/1680 = 0.830$ |
| Atmospheric pressure, 0.101 MPa | $P^* = 0.101/41.9 = 0.0024$ |
| An MD time step of 2 fs | $\Delta t^* = 2\times10^{-15}/2.156\times10^{-12} = 0.00093$ |

**What $\rho^* = 1$ means physically:** one atom per $\sigma^3$, i.e. atoms packed at about one diameter apart. Ordinary liquids sit around $\rho^* \approx 0.7$–$0.85$.

### 6.3 The Monte Carlo problem from the earlier notes

| Reduced | Argon |
|---|---|
| $T^* = 0.85$ | 101.8 K |
| $\rho^* = 0.776$ | 1304 kg/m³ |
| $\langle U\rangle/N = -5.51\,\varepsilon$ | $-5.49$ kJ/mol |
| Typical MD time step $0.005\,\tau$ | 10.8 fs |

So that simulation was liquid argon at about 102 K.

---

## 7. The phase diagram of the Lennard-Jones fluid in reduced units

### 7.1 Reading the sketch

The lecture sketches the vapour–liquid equilibrium (VLE) of the Lennard-Jones fluid as a dome, with density along the horizontal axis. The top of the dome is the critical temperature and the bottom is the triple point, so the vertical axis is temperature.

| Region | Where | Lecture's description |
|---|---|---|
| Vapour | left of the dome | $\rho^* < 0.01$ (at low temperature) |
| Liquid | right of the dome | $\rho^* > 0.7$ |
| Two-phase region | inside the dome | vapour and liquid coexist |
| Critical point | top of the dome | $T_c^* = k_BT_c/\varepsilon \approx 1.4$ in the lecture |
| Triple point | bottom of the dome | $T^* \approx 0.7$ |

Because of corresponding states, this **single diagram** applies to every Lennard-Jones material.

### 7.2 Accurate literature values

The lecture's "≈ 1.4" is a rough figure. Modern simulations of the full Lennard-Jones potential (no cutoff, or with tail corrections) give:

| Quantity | Reduced value |
|---|---|
| Critical temperature $T_c^*$ | 1.312 |
| Critical density $\rho_c^*$ | 0.316 |
| Critical pressure $P_c^*$ | about 0.13–0.14 |
| Triple-point temperature $T_t^*$ | about 0.69 |
| Liquid density at the triple point | about 0.85 |

These values depend on how the potential is cut off. Truncating it (for example at $r_c = 2.5\sigma$ without a tail correction) lowers $T_c^*$ noticeably, to roughly 1.1–1.2. Always state the cutoff when you quote a phase-diagram number.

### 7.3 Checking against real argon

| Property | Lennard-Jones prediction (× argon units) | Real argon |
|---|---|---|
| Critical temperature | $1.312\times119.8 = 157$ K | 150.7 K |
| Critical density | $0.316\times1680 = 531$ kg/m³ | about 536 kg/m³ |
| Critical pressure | $0.13$–$0.14\times41.9 \approx 5.5$–$5.9$ MPa | 4.86 MPa |
| Triple-point temperature | $0.69\times119.8 = 83$ K | 83.8 K |

The Lennard-Jones model with these two parameters reproduces argon's phase diagram to within a few percent for temperature and density. The agreement is partly by construction, since $\sigma$ and $\varepsilon$ were fitted to argon data.

### 7.4 Where the earlier Monte Carlo state point sits

$T^* = 0.85$ lies between the triple point (0.69) and the critical point (1.31). At $\rho^* = 0.776$ the system is on the **liquid side, right next to coexistence**: the NIST reference pressure at this state point is almost zero ($P^* = 0.0068$).

---

## 8. Reduced units in practice

### 8.1 In your own codes

Set $\sigma = \varepsilon = m = k_B = 1$ and write everything in reduced form:

- potential: `4*(r**-12 - r**-6)`
- Boltzmann factor: `exp(-dU / T)` with `T` meaning $T^*$
- tail correction: $\frac{8}{3}\pi N\rho^*\left[\frac13 r_c^{*-9} - r_c^{*-3}\right]$

Convert to SI only at the end, using the table in Section 5.

### 8.2 In LAMMPS

`units lj` is exactly this system: distances in $\sigma$, energies in $\varepsilon$, masses in $m$, time in $\tau$, and $k_B = 1$. A typical MD time step is `timestep 0.005` (that is, $0.005\,\tau \approx 10.8$ fs for argon). All thermo output (temperature, pressure, energy) is then reduced and must be multiplied by the units above to get SI values.

### 8.3 When reduced units do not fit naturally

- **Mixtures:** choose one species as the reference ($\sigma_A$, $\varepsilon_A$) and describe the others by ratios such as $\sigma_B/\sigma_A$ and $\varepsilon_B/\varepsilon_A$.
- **Molecules with charges or bonds:** electrostatics and bond potentials bring in extra energy and length scales, so a single $\sigma$, $\varepsilon$ no longer removes all constants.
- **Machine-learned potentials (for example MACE):** there is no single $\sigma$ or $\varepsilon$, so real units are used (in LAMMPS, `units metal`: Å, eV, ps, K, bar).

### 8.4 Limits of corresponding states

The law of corresponding states is exact only for a pure Lennard-Jones model. Real substances deviate when:

- the atoms are light, so quantum effects matter (helium, hydrogen, neon),
- the molecules are not spherical,
- many-body forces are significant.

---

## 9. Python converter

```python
kB  = 1.380649e-23        # Boltzmann constant, J/K
amu = 1.66053906660e-27   # atomic mass unit, kg

sigma = 3.405e-10         # m
eps   = 119.8 * kB        # J   (eps/kB = 119.8 K)
m     = 39.948 * amu      # kg  (argon)

units = {
    "length (m)":            sigma,
    "energy (J)":            eps,
    "temperature (K)":       eps / kB,
    "number density (1/m3)": 1 / sigma**3,
    "mass density (kg/m3)":  m / sigma**3,
    "pressure (Pa)":         eps / sigma**3,
    "time (s)":              (m * sigma**2 / eps) ** 0.5,
    "surface tension (N/m)": eps / sigma**2,
}
for name, value in units.items():
    print(f"{name:24s} 1 reduced unit = {value:.4g}")

# example: the state point of the Monte Carlo problem
T_star, rho_star = 0.85, 0.776
print("T   =", round(T_star * units["temperature (K)"], 1), "K")
print("rho =", round(rho_star * units["mass density (kg/m3)"]), "kg/m3")
```

Output:

```
length (m)               1 reduced unit = 3.405e-10
energy (J)               1 reduced unit = 1.654e-21
temperature (K)          1 reduced unit = 119.8
number density (1/m3)    1 reduced unit = 2.533e+28
mass density (kg/m3)     1 reduced unit = 1680
pressure (Pa)            1 reduced unit = 4.19e+07
time (s)                 1 reduced unit = 2.156e-12
surface tension (N/m)    1 reduced unit = 0.01427
T   = 101.8 K
rho = 1304 kg/m3
```

To use another Lennard-Jones material, change only `sigma`, `eps` and `m`.

---

## 10. Common mistakes

- **Mixing up number density and mass density.** $\rho^* = \rho\sigma^3$ uses the **number** density. To get kg/m³, multiply by the atomic mass $m$.
- **Mixing up $\varepsilon$ and $\varepsilon/k_B$.** Tables quote $\varepsilon/k_B$ in kelvin; for energies and pressures you need $\varepsilon$ in joules ($\varepsilon = 119.8\times k_B$).
- **Writing $T^* = T/\varepsilon$.** This mixes kelvin with joules. The correct form is $T^* = k_BT/\varepsilon$.
- **Forgetting the mass in the time unit.** $\tau = \sqrt{m\sigma^2/\varepsilon}$ depends on the atomic mass, so two materials with the same $\sigma$ and $\varepsilon$ but different masses have different time units.
- **Comparing results with different cutoffs.** Reduced phase-diagram numbers (like $T_c^*$) change with the cutoff and with whether tail corrections are used.

---

## 11. Practice questions (with answers)

1. **Convert $T^* = 0.69$ and $T^* = 1.312$ to kelvin for argon.**
   Answer: 82.7 K and 157.2 K.

2. **Liquid argon at its boiling point has a density of about 1395 kg/m³. What is $\rho^*$?**
   Answer: $1395/1680 = 0.830$.

3. **What is atmospheric pressure (0.101 MPa) in reduced units for argon?**
   Answer: $0.101/41.9 = 0.0024$.

4. **An MD run uses `timestep 0.005` in `units lj`. How long is one step for argon?**
   Answer: $0.005\times2.156$ ps $= 10.8$ fs.

5. **Show that the ideal gas law becomes $P^* = \rho^*T^*$.**
   Answer: Section 5.4.

6. **Two Lennard-Jones liquids, A and B, are simulated at $T^* = 0.85$, $\rho^* = 0.776$. Liquid B has twice the $\varepsilon$ of A. How do their potential energies per atom compare (in J)?**
   Answer: the reduced energy per atom is the same ($-5.51$), so in joules B's is twice A's. B's real temperature is also twice A's.

---

## 12. Summary

- Choose $\sigma$ (length), $\varepsilon$ (energy) and $m$ (mass) as units, and use $k_B$ to turn temperature into energy. Every other unit follows from dimensions.
- In these units the Lennard-Jones potential becomes $u^* = 4[r^{*-12} - r^{*-6}]$, the same curve for every material.
- The Boltzmann factor depends only on $U^*/T^*$, and Newton's equation becomes constant-free when time is measured in $\tau = \sqrt{m\sigma^2/\varepsilon}$. So **all** Lennard-Jones materials share the same reduced properties at the same $T^*$ and $\rho^*$ (law of corresponding states).
- One simulation in reduced units gives results for argon, krypton, methane, … : multiply by each material's units at the end.
- Key reduced numbers: critical point $T_c^* \approx 1.31$, $\rho_c^* \approx 0.32$; triple point $T^* \approx 0.69$; liquid $\rho^* \gtrsim 0.7$; vapour $\rho^* \lesssim 0.01$.
- Argon units: $\sigma = 3.405$ Å, $\varepsilon/k_B = 119.8$ K, 1680 kg/m³, 41.9 MPa, 2.156 ps, 14.3 mN/m.

**Sources for the literature values in Section 7:** Heyes (2015), *Computational Methods in Science and Technology* 21(4), which compiles $T_c^* = 1.3120$ (Potoff & Panagiotopoulos), $\rho_c^* = 0.316$, $P_c^* = 0.141$ and $T_t^* \approx 0.69$; NIST Lennard-Jones benchmark data, https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm
