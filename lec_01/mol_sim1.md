# Statistical Thermodynamics --- Study Notes

## Correctly formatted equations for Markdown and LaTeX

**Rendering tip:** Open this file in a Markdown editor that supports
LaTeX math, such as VS Code with a Markdown math extension, Obsidian, or
a compatible previewer. Use `$...$` for inline math and `$$...$$` for
displayed equations. Some basic Markdown viewers show the raw LaTeX
rather than rendered equations.

## 1. Statistical thermodynamics

Statistical thermodynamics connects the microscopic behavior of
particles to macroscopic properties such as temperature, pressure,
energy, and free energy.

A classical microscopic state is specified by the positions and momenta
of all particles. Different microscopic states can correspond to the
same macrostate.

## 2. Macrostate, microstate, and phase space

For a system of (N) particles,

$$
\mathbf r^N=(\mathbf r_1,\ldots,\mathbf r_N), \qquad
\mathbf p^N=(\mathbf p_1,\ldots,\mathbf p_N).
$$

The superscript (N) means "the collection for all (N) particles," not a
power applied to a vector.

In three dimensions, the phase space has (3N) position coordinates and
(3N) momentum coordinates, for a total of (6N) dimensions.

Define the phase-space volume element

$$
d\Gamma_N \equiv d\mathbf r^N\,d\mathbf p^N.
$$

Written explicitly,

$$
d\mathbf r^N=\prod_{i=1}^{N} dx_i\,dy_i\,dz_i,
\qquad
d\mathbf p^N=\prod_{i=1}^{N} dp_{ix}\,dp_{iy}\,dp_{iz}.
$$

## 3. Mechanical properties and the Hamiltonian

A mechanical property (M) depends on the microscopic state. Examples
include kinetic energy, potential energy, total energy, and
instantaneous pressure.

The classical Hamiltonian is

$$
H_N(\mathbf r^N,\mathbf p^N)=K_N(\mathbf p^N)+U_N(\mathbf r^N).
$$

The kinetic energy is

$$
K_N=\sum_{i=1}^{N}\frac{\mathbf p_i^2}{2m_i},
$$

where (m_i) is the mass of particle (i). The potential energy (U_N)
depends on particle positions and their interaction model.

## 4. Discrete and continuous systems

### Discrete states

For discrete states (i), with property (M_i) and probability (P_i),

$$
\boxed{\langle M\rangle=\sum_i M_iP_i},
\qquad
\sum_iP_i=1.
$$

Example: for values (M=(10,20,30)) and probabilities (P=(0.2,0.5,0.3)),

$$
\langle M\rangle=10(0.2)+20(0.5)+30(0.3)=21.
$$

A quantum particle in a one-dimensional infinite potential well has
energy levels

$$
E_n=\frac{n^2h^2}{8mL^2},\qquad n=1,2,3,\ldots
$$

A quantum system can also have a continuous energy spectrum; discrete
energies are not universal to all quantum systems.

### Continuous variables

For a continuous variable (x) with probability density (f(x)),

$$
\boxed{\langle x\rangle=\int x f(x)\,dx}.
$$

Normalization requires

$$
f(x)\geq 0,\qquad \int f(x)\,dx=1.
$$

The probability of finding (x) in a small interval from (x) to (x+dx) is
approximately (f(x),dx).

## 5. General ensemble average

For a normalized probability density
(f(`\mathbf`{=tex}r^N,`\mathbf`{=tex}p^N)), the ensemble average of a
property (M) is

$$
\boxed{
\langle M\rangle_{\mathrm{ens}}=
\int M(\mathbf r^N,\mathbf p^N)
f(\mathbf r^N,\mathbf p^N)\,d\Gamma_N
}.
$$

The probability density must satisfy

$$
\int f(\mathbf r^N,\mathbf p^N)\,d\Gamma_N=1.
$$

For discrete states, the equivalent expression is

$$
\langle M\rangle=\sum_i M_iP_i.
$$

## 6. Time average and ergodic hypothesis

For one trajectory observed over duration (`\tau`{=tex}),

$$
\boxed{
\overline M_{\mathrm{time}}=\frac{1}{\tau}\int_0^\tau M(t)\,dt
}.
$$

For sampled data, an estimate is

$$
\overline M_{\mathrm{time}}\approx \frac{1}{K}\sum_{k=1}^{K}M(t_k).
$$

Under suitable equilibrium and ergodicity conditions, the long-time
average equals the corresponding ensemble average:

$$
\boxed{\overline M_{\mathrm{time}}=\langle M\rangle_{\mathrm{ens}}}.
$$

Ergodicity does not mean every imaginable state is visited. It concerns
states accessible under the system's constraints and the relevant
statistical measure. A long trajectory does not guarantee good sampling
if the system remains trapped in a metastable state.

## 7. Main ensembles

  -----------------------------------------------------------------------
  Ensemble                Controlled or fixed     Quantities that can
                          quantities              fluctuate
  ----------------------- ----------------------- -----------------------
  NVE (microcanonical)    (N,V,E)                 Temperature is not
                                                  imposed; microscopic
                                                  observables fluctuate

  NVT (canonical)         (N,V,T)                 Energy

  NPT                     (N,P,T)                 Volume and energy
  (isothermal-isobaric)                           

  (`\mu`{=tex}VT) (grand  (`\mu`{=tex},V,T)       Particle number and
  canonical)                                      energy
  -----------------------------------------------------------------------

## 8. NVE --- microcanonical ensemble

The Hamiltonian is restricted to the constant-energy surface:

$$
H_N(\mathbf r^N,\mathbf p^N)=E.
$$

A common classical microcanonical density is

$$
\boxed{
f_E(\mathbf r^N,\mathbf p^N)=
\frac{\delta(E-H_N)}{\omega_N(E)}
}
$$

where

$$
\boxed{
\omega_N(E)=\int \delta(E-H_N)\,d\Gamma_N
}
$$

is the density of states on the constant-energy surface under this
convention. The Dirac delta function (`\delta`{=tex}) restricts the
integral to energy (E).

The ensemble average is

$$
\boxed{
\langle M\rangle_{NVE}=
\frac{\int M\,\delta(E-H_N)\,d\Gamma_N}
{\int \delta(E-H_N)\,d\Gamma_N}
}.
$$

Some texts use (`\Omega`{=tex}\_N(E)) for the cumulative phase-space
volume

$$
\Omega_N(E)=\int \Theta(E-H_N)\,d\Gamma_N,
$$

where (`\Theta`{=tex}) is the Heaviside step function. This differs from
the density of states (`\omega`{=tex}\_N(E)). Where differentiable,

$$
\omega_N(E)=\frac{d\Omega_N(E)}{dE}.
$$

Check the textbook's definition before using (1/`\Omega`{=tex}).

## 9. NVT --- canonical ensemble

For the canonical ensemble,

$$
\beta=\frac{1}{k_BT},
$$

where (k_B) is Boltzmann's constant.

### Discrete physical states

For discrete physical states (i) with energies (E_i),

$$
\boxed{Z_N(V,T)=\sum_i e^{-\beta E_i}}.
$$

If (E_i) denotes an energy level with degeneracy (g_i), this can be
written as

$$
Z_N=\sum_i g_i e^{-\beta E_i}.
$$

The probability of state (i) is

$$
\boxed{P_i=\frac{e^{-\beta E_i}}{Z_N}},
$$

and the ensemble average is

$$
\boxed{
\langle M\rangle_{NVT}
=\frac{\sum_i M_i e^{-\beta E_i}}{Z_N}
=\sum_i M_iP_i
}.
$$

The quantum state sum should count properly defined physical
many-particle states; do not automatically multiply it by the classical
factor (1/N!).

### Classical phase-space form

For (N) identical classical particles, a conventional canonical
partition function is

$$
\boxed{
Z_N(V,T)=\frac{1}{N!h^{3N}}
\int e^{-\beta H_N(\mathbf r^N,\mathbf p^N)}\,d\Gamma_N
}.
$$

Here (h) is Planck's constant and (1/N!) corrects for permutation
overcounting of identical particles in classical phase space.

The classical ensemble average can be written as

$$
\boxed{
\langle M\rangle_{NVT}=
\frac{\int M e^{-\beta H_N}\,d\Gamma_N}
{\int e^{-\beta H_N}\,d\Gamma_N}
}.
$$

The same constant normalization factors cancel between numerator and
denominator.

### Integrating out momenta

For (H_N=K_N+U_N), define the thermal de Broglie wavelength

$$
\Lambda=\frac{h}{\sqrt{2\pi m k_BT}}.
$$

For (N) identical particles of mass (m),

$$
\boxed{
Z_N(V,T)=\frac{1}{N!\Lambda^{3N}}
\int_V e^{-\beta U(\mathbf r^N)}\,d\mathbf r^N
}.
$$

### Binary mixture

For (N_A) identical particles of species A and (N_B) identical particles
of species B, with (N=N_A+N_B),

$$
\boxed{
Z_{N_A,N_B}(V,T)=
\frac{1}{N_A!N_B!\Lambda_A^{3N_A}\Lambda_B^{3N_B}}
\int_V e^{-\beta U(\mathbf r^{N_A},\mathbf r^{N_B})}
\,d\mathbf r^{N_A}\,d\mathbf r^{N_B}
}.
$$

For species (s=A,B),

$$
\Lambda_s=\frac{h}{\sqrt{2\pi m_s k_BT}}.
$$

The factorials are (N_A!N_B!), not ((N_A+N_B)!), because A and B are
distinct species, while particles within each species are treated as
identical.

### When does (N!) appear?

-   Classical phase-space partition function for (N) identical
    particles: use (1/N!) under the conventional indistinguishability
    correction.
-   Classical mixture: use (1/`\prod`{=tex}\_s N_s!).
-   Genuinely distinguishable classical particles: do not add (1/N!)
    just as an indistinguishability correction.
-   Quantum state sums: count physical many-particle states correctly;
    do not blindly insert the classical factorial.

## 10. NPT --- isothermal-isobaric ensemble

NPT fixes (N,P,T); volume and energy can fluctuate. It is useful when
the target state is specified by temperature and pressure.

Let (Z_N(V,T)) be the canonical partition function at volume (V). A
common isothermal-isobaric partition function convention is

$$
\boxed{
\Delta_N(P,T)=\int_0^\infty Z_N(V,T)e^{-\beta PV}\,dV
}.
$$

The factor (e\^{-`\beta`{=tex}PV}) weights the volumes. The equilibrium
distribution reflects both this factor and the volume dependence of
(Z_N(V,T)).

The volume probability density, under this convention, is

$$
\boxed{
p(V)=\frac{Z_N(V,T)e^{-\beta PV}}{\Delta_N(P,T)}
}.
$$

If (`\langle`{=tex}M`\rangle`{=tex}\_{NVT,V}) is the canonical average
at a particular volume, then

$$
\boxed{
\langle M\rangle_{NPT}=
\frac{
\int_0^\infty Z_N(V,T)e^{-\beta PV}
\langle M\rangle_{NVT,V}\,dV
}{
\int_0^\infty Z_N(V,T)e^{-\beta PV}\,dV
}
}.
$$

The expression assumes the stated volume measure and consistent
normalization. Flexible-cell algorithms may require additional
cell-shape variables and the appropriate measure.

In MD, a thermostat controls temperature and a barostat adjusts volume
or cell dimensions to sample the target pressure. The precise algorithm
affects how sampling is implemented.

## 11. Grand canonical --- (`\mu`{=tex}VT)

The grand canonical ensemble fixes chemical potential (`\mu`{=tex}),
volume (V), and temperature (T). Particle number (N) and energy can
fluctuate.

### Discrete physical states

Let (i) label a physical state with (N) particles and energy (E\_{N,i}).
The grand partition function is

$$
\boxed{
\Xi(\mu,V,T)=\sum_{N=0}^{\infty}\sum_i
e^{-\beta(E_{N,i}-\mu N)}
}.
$$

The probability of state ((N,i)) is

$$
\boxed{
P_{N,i}=\frac{e^{-\beta(E_{N,i}-\mu N)}}{\Xi}
}.
$$

The ensemble average is

$$
\boxed{
\langle M\rangle_{\mu VT}=
\frac{1}{\Xi}\sum_{N=0}^{\infty}\sum_i
M_{N,i}e^{-\beta(E_{N,i}-\mu N)}
}.
$$

### Classical phase-space form

For a classical system of identical particles,

$$
\boxed{
\Xi(\mu,V,T)=
\sum_{N=0}^{\infty}\frac{1}{N!h^{3N}}
\int e^{-\beta[H_N(\mathbf r^N,\mathbf p^N)-\mu N]}\,d\Gamma_N
}.
$$

The ensemble average is

$$
\boxed{
\langle M\rangle_{\mu VT}=
\frac{1}{\Xi}\sum_{N=0}^{\infty}\frac{1}{N!h^{3N}}
\int M_N e^{-\beta(H_N-\mu N)}\,d\Gamma_N
}.
$$

For mixtures, use a product of species factorials, such as (N_A!N_B!),
with separate particle-number sums for species that can exchange with
reservoirs.

## 12. Compare the ensemble weights

The following formulas summarize the statistical weight and
normalization for each ensemble. They are written separately rather than
inside a table so that Markdown math renderers can display the equations
reliably.

### NVE --- microcanonical ensemble

**Statistical weight / restriction**

$$
\delta(E-H_N)
$$

**Normalization**

$$
\omega_N(E)=\int \delta(E-H_N)\,d\Gamma_N
$$

### NVT --- canonical ensemble, discrete states

**Statistical weight**

$$
e^{-\beta E_i}
$$

**Partition function**

$$
Z_N(V,T)=\sum_i e^{-\beta E_i}
$$

### NVT --- canonical ensemble, classical phase space

**Statistical weight**

$$
e^{-\beta H_N}
$$

**Partition function**

$$
Z_N(V,T)=\frac{1}{N!h^{3N}}\int e^{-\beta H_N}\,d\Gamma_N
$$

### NPT --- isothermal-isobaric ensemble

**Volume-dependent statistical weight**

$$
Z_N(V,T)e^{-\beta PV}
$$

**Isothermal-isobaric partition function**

$$
\Delta_N(P,T)=\int_0^\infty Z_N(V,T)e^{-\beta PV}\,dV
$$

### Grand canonical --- discrete states

**Statistical weight**

$$
e^{-\beta(E_{N,i}-\mu N)}
$$

**Grand partition function**

$$
\Xi(\mu,V,T)=\sum_{N=0}^{\infty}\sum_i e^{-\beta(E_{N,i}-\mu N)}
$$

### General principle

For any ensemble, the average is the weighted sum or integral of the
property divided by the sum or integral of the statistical weights:

$$
\boxed{
\text{Ensemble average}=
\frac{\text{weighted sum or integral of the property}}
{\text{sum or integral of the statistical weights}}
}
$$

## 13. Thermodynamic potentials

### Helmholtz free energy (canonical ensemble)

$$
\boxed{F=-k_BT\ln Z_N}
$$

### Gibbs free energy (isothermal-isobaric ensemble)

$$
\boxed{G=-k_BT\ln\Delta_N}
$$

### Grand potential (grand canonical ensemble)

$$
\boxed{\mathcal J=-k_BT\ln\Xi}
$$

The grand potential is also denoted (`\Omega`{=tex}\_G) or
(`\Phi`{=tex}\_G) in some texts. Do not confuse it with a microcanonical
phase-space volume (`\Omega`{=tex}(E)). Additive constants can depend on
partition-function normalization; use a consistent convention for
thermodynamic derivatives and physical predictions.

## 14. Quick review

1.  **NVE:** fixed (N,V,E).
2.  **NVT:** fixed (N,V,T); energy fluctuates.
3.  **NPT:** fixed (N,P,T); volume and energy fluctuate.
4.  **Grand canonical:** fixed (`\mu`{=tex},V,T); particle number and
    energy fluctuate.
5.  **Discrete states:** use sums.
6.  **Continuous phase-space variables:** use integrals.
7.  **Ergodic hypothesis:** under suitable conditions, long-time and
    ensemble averages agree.
8.  **Classical (1/N!):** corrects permutation overcounting for
    identical particles.
