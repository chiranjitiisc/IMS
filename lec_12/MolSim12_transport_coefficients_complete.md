# Lecture 12: Transport Coefficients — Complete Notes from Scratch

*Based on 12.pdf ("Transport coefficients", 10 February 2022). Everything in the lecture is explained here from the beginning, in simple steps: the diffusion equation, its Gaussian solution, the mean squared displacement (MSD), the lecture's derivation of $d\langle r^2\rangle/dt = 6D$, ballistic and diffusive motion, how the MSD looks in a gas, a liquid, a supercooled liquid and a crystal, the Green–Kubo formula for diffusion, and shear viscosity from the pressure tensor. Every question asked while studying this lecture is answered where it belongs (look for the boxes marked **Your question**). All numbers were checked with small molecular-dynamics (MD) simulations written for these notes. The three figures are in the `figures` folder next to this file.*

---

## Contents

**Part 0. Orientation**
- 0.1 What a transport coefficient is
- 0.2 Notation
- 0.3 The lecture in one page
- 0.4 What the lecture covers, and what these notes add

**Part I. The diffusion equation**
- 1. Concentration and flux
- 2. Conservation: what comes in minus what goes out
- 3. Fick's law: flow goes downhill
- 4. The diffusion equation and what the second derivative means
- 5. The same thing in 3D: gradient, divergence, Laplacian

**Part II. Spreading from a point**
- 6. The starting condition: the Dirac delta
- 7. The 1D solution is a Gaussian (and checking it)
- 8. The 3D solution
- 9. Why the width grows as √t: the random walk
- 10. What $c$ means in an MD simulation

**Part III. Mean squared displacement**
- 11. What the MSD is, and why we square
- 12. MSD = 6Dt straight from the Gaussian
- 13. The lecture's derivation of $d\langle r^2\rangle/dt = 6D$
- 14. Computing the MSD from an MD trajectory (time origins explained)
- 15. How LAMMPS computes the MSD
- 16. Practical traps: periodic boundaries, drift, species

**Part IV. Reading the MSD curve**
- 17. Ballistic motion: what it means and why MSD ∝ t²
- 18. Diffusive motion: why MSD ∝ t
- 19. The simulations used in these notes
- 20. Log–log plots and the local slope α
- 21. Gas versus liquid
- 22. Supercooled liquid: the cage
- 23. Crystal: overshoot and plateau
- 24. Other things the MSD tells you, and warning signs

**Part V. Green–Kubo relation for diffusion**
- 25. The velocity autocorrelation function (VACF)
- 26. Deriving $D = \frac13\int\langle\vec v(t)\cdot\vec v(0)\rangle dt$
- 27. Checks with the simulations, and what D means
- 28. MSD or VACF in practice

**Part VI. Shear viscosity**
- 29. What viscosity is
- 30. The pressure tensor
- 31. The Green–Kubo formula for viscosity
- 32. Our result, its error bar, and two cross-checks
- 33. Viscosity in LAMMPS

**Part VII. Wrap-up**
- 34. Corrections and clarifications to the lecture notes
- 35. Your questions, in one place
- 36. Practice questions with answers
- 37. Python code used in these notes
- 38. Summary: the whole chain

---

# Part 0. Orientation

## 0.1 What a transport coefficient is

A **transport coefficient** is a number that says **how fast something spreads or relaxes** in a material.

| Coefficient | Symbol | What spreads | Everyday picture |
|---|---|---|---|
| diffusion coefficient | $D$ | molecules | a drop of ink spreading in still water |
| shear viscosity | $\eta$ | momentum (sideways motion) | honey resisting being stirred |
| thermal conductivity | $\kappa$ | heat | a metal spoon warming up in hot tea (not in this lecture) |

The lecture shows **two ways** to get them from an ordinary equilibrium MD run, without pushing the system:

1. **Einstein route:** watch how far molecules wander (the mean squared displacement). Its slope gives $D$.
2. **Green–Kubo route:** watch how long a fluctuation is "remembered" (an autocorrelation function). Its integral gives $D$, or $\eta$.

## 0.2 Notation

Every symbol is also defined again where it first appears.

| Symbol | Meaning |
|---|---|
| $t$ | time (in MSD and correlation functions: the **lag**, the time elapsed since the start) |
| $x, y, z$ | Cartesian coordinates |
| $\vec r = (x, y, z)$ | position vector, or (in MD) displacement vector |
| $r^2 = x^2 + y^2 + z^2$ | squared length of $\vec r$ |
| $c(x,t)$, $c(\vec r,t)$ | concentration (molecules per unit volume); later, probability density of a displacement |
| $j$, $\vec j$ | flux: net molecules crossing unit area per unit time |
| $D$ | diffusion coefficient (m²/s) |
| $\partial$ | partial derivative (vary one variable, keep the others fixed) |
| $\nabla$ | "del": the vector $(\partial/\partial x,\ \partial/\partial y,\ \partial/\partial z)$ |
| $\nabla c$, $\nabla\cdot\vec j$, $\nabla^2 c$ | gradient, divergence, Laplacian (Section 5) |
| $\delta(\vec r)$ | Dirac delta: a spike of total area 1 at the origin |
| $\sigma$ | standard deviation of a Gaussian (Part II); the Lennard-Jones length (Parts IV–VI) |
| $\langle\ \rangle$ | average |
| $N$ | number of molecules (Part III); number of steps (Part IV, random walk) |
| $i$ | molecule index |
| $\vec r_i(t)$ | position of molecule $i$ at time $t$ |
| $t_0$ | time origin (the moment the stopwatch is started) |
| $K$ | number of time origins used for one lag |
| $L$ | box length |
| $m$ | mass of a molecule |
| $k_B$ | Boltzmann's constant, $1.381\times10^{-23}$ J/K |
| $T$ | temperature |
| $\vec v$, $v_x$ | velocity, and its x-component |
| $\ell$ | step (hop) length in a random walk |
| $t_h$, $t_c$ | time per hop; time per step (about the time between collisions) |
| $\alpha$ | local slope of log MSD versus log t |
| $\vec u$ | displacement of an atom from its lattice site (crystal) |
| $C(t)$ | velocity autocorrelation function, $\langle\vec v(t)\cdot\vec v(0)\rangle$ |
| $\tau_v$ | velocity memory time |
| $\eta$ | shear viscosity |
| $P_{\alpha\beta}$ | pressure tensor; $\alpha, \beta \in \lbrace x, y, z\rbrace$ |
| $V$ | volume of the simulation box |

## 0.3 The lecture in one page

1. **Diffusion equation.** Conservation, $\partial c/\partial t = -\nabla\cdot\vec j$, plus Fick's law, $\vec j = -D\nabla c$, give $\partial c/\partial t = D\nabla^2 c$.
2. **Point source.** If everything starts at one point, $c(\vec r, 0) = \delta(\vec r)$, the solution is a spreading Gaussian, $c(\vec r,t) = (4\pi Dt)^{-3/2}e^{-r^2/4Dt}$.
3. **MSD.** $\langle r^2(t)\rangle = \int r^2 c\,d\vec r$. Using the diffusion equation and integrating by parts twice, $d\langle r^2\rangle/dt = 6D$. So **the slope of the MSD is $6D$**.
4. **MSD curve of a fluid.** At short times $\langle r^2\rangle \propto t^2$ ("ballistic", before collisions); at long times $\langle r^2\rangle \propto t$ ("diffusive"), with slope $6D$.
5. **Green–Kubo.** Transport coefficients are integrals of autocorrelation functions:
$$D = \frac13\int_0^\infty\langle\vec v(t)\cdot\vec v(0)\rangle\,dt, \qquad \eta = \frac{V}{k_BT}\int_0^\infty\langle P_{xy}(t)P_{xy}(0)\rangle\,dt$$
6. **Pressure tensor.** $P_{\alpha\beta} = \frac1V\left(\sum_i p_{i\alpha}p_{i\beta}/m + \sum_i r_{i\alpha}f_{i\beta}\right)$.

## 0.4 What the lecture covers, and what these notes add

| Topic | Lecture | Section |
|---|---|---|
| diffusion equation | stated | 1–5 (derived, with intuition) |
| Gaussian solution | stated | 6–9 (checked, random-walk meaning) |
| MSD definition | stated | 11 |
| $d\langle r^2\rangle/dt = 6D$ | derived | 12–13 (every step, 1D first) |
| what $c$ and $r$ mean in MD | not stated | 10, 14 (correction: $r$ is a displacement) |
| computing the MSD, time origins, LAMMPS | added | 14–16 |
| ballistic $t^2$ and diffusive $t$ | sketch | 17–18 (derived) |
| MSD of gas, liquid, supercooled liquid, crystal | partly (fluid only) | 19–24 (simulated) |
| Green–Kubo for $D$ | stated | 25–28 (derived and simulated) |
| viscosity and pressure tensor | stated | 29–33 (corrected and simulated) |
| your questions | — | boxes throughout, collected in 35 |

---

# Part I. The diffusion equation

**The picture.** Put a drop of ink into still water. Nobody stirs it. Each ink molecule is just knocked about at random by water molecules. Yet the ink spreads from where there is a lot of it to where there is little. The diffusion equation is the bookkeeping of this spreading. It needs two ingredients:

- (a) molecules are not created or destroyed (**conservation**);
- (b) a rule for which way they flow (**Fick's law**).

We do everything in **1D first** (a long, thin tube along $x$), then move to 3D.

## 1. Concentration and flux

- **Concentration** $c(x,t)$: the number of molecules per unit volume at position $x$ and time $t$.
- **Flux** $j(x,t)$: the **net** number of molecules crossing a unit area (a cross-section of the tube) per unit time, at position $x$. If $j > 0$, more molecules cross towards $+x$ than towards $-x$.
- **Diffusion coefficient** $D$, in m²/s: a property of the molecule and the medium it moves in.
- **Partial derivative** $\partial$: differentiate with respect to one variable while holding the others fixed. For example, $\partial c/\partial t$ is how fast $c$ changes in time **at a fixed place** $x$.

## 2. Conservation: what comes in minus what goes out

Take a thin slice of the tube between $x$ and $x + \Delta x$, with cross-section area $A$.

- Molecules inside the slice: $c\,A\,\Delta x$ (concentration × volume).
- They can only change by crossing the two faces: $j(x)\,A$ come in through the left face each unit of time, and $j(x+\Delta x)\,A$ leave through the right face.

$$\frac{\partial}{\partial t}\left(c\,A\,\Delta x\right) = j(x)\,A - j(x+\Delta x)\,A$$

**Example.** $j(x) = 5$ and $j(x + \Delta x) = 3$ molecules per Å² per ps. Then 5 come in, 3 leave: the slice **gains 2** per Å² per ps, and its concentration rises.

Divide both sides by $A\,\Delta x$ and let the slice become very thin ($\Delta x \to 0$; this limit is the definition of a derivative):

$$\boxed{\;\frac{\partial c}{\partial t} = -\frac{\partial j}{\partial x}\;}$$

**In words: the concentration rises wherever the flux decreases along $x$** (more comes in than goes out). There is no physics of diffusion here yet; this is just counting. It holds for anything that is conserved.

## 3. Fick's law: flow goes downhill

$$\boxed{\;j = -D\,\frac{\partial c}{\partial x}\;}$$

The flux is proportional to the slope of the concentration, and points **opposite** to it.

### 3.1 Why: a random-hop model

Suppose every time interval $t_h$, each molecule hops a distance $\ell$ to the left or to the right, with probability ½ each. Look at an imaginary plane at position $x$:

- **From the left:** the molecules in the strip of width $\ell$ just left of the plane number about $c(x - \ell/2)\,\ell$ per unit area (concentration at the middle of the strip × width). Half of them hop right and cross the plane.
- **From the right:** the strip just right of the plane holds $c(x + \ell/2)\,\ell$ per unit area. Half hop left and cross.

Net crossings to the right, per unit area, per unit time:

$$j = \frac{\ell}{2t_h}\Big[c(x - \tfrac{\ell}{2}) - c(x + \tfrac{\ell}{2})\Big]$$

For a small hop, $c(x+\ell/2) - c(x-\ell/2) \approx \ell\,\partial c/\partial x$ (the change in $c$ over a distance $\ell$ is slope × distance). So

$$j \approx -\frac{\ell^2}{2t_h}\,\frac{\partial c}{\partial x}, \qquad D = \frac{\ell^2}{2t_h}$$

**The key point: no molecule "knows" about the gradient.** Each one hops at random. The net flow appears only because there are more molecules on one side available to hop across.

**Numbers.** $\ell = 1$ Å $= 10^{-10}$ m, $t_h = 2$ ps $= 2\times10^{-12}$ s:

$$D = \frac{(10^{-10})^2}{2\times2\times10^{-12}} = 2.5\times10^{-9}\ \text{m}^2/\text{s}$$

close to water's self-diffusion coefficient at 25 °C, about $2.3\times10^{-9}$ m²/s. (A toy model, but the right size.)

### 3.2 The sign: two minus signs

**Rule: the flux always has the opposite sign to the slope.**

| Profile | Slope $\partial c/\partial x$ | Flux $j = -D\,\partial c/\partial x$ | Flow goes |
|---|---|---|---|
| $c$ falls with $x$ (8, 6, 4) | negative | **positive** | towards $+x$, from high to low |
| $c$ rises with $x$ (2, 5, 8) | positive | **negative** | towards $-x$, from high to low |

> **Common mistake (made while studying this).** For $c$ = 8, 6, 4 at $x$ = 0, 1, 2 nm with $D = 2$ nm²/ns, the answer "flux = −4" is wrong. Step by step: slope $= (4 - 8)/2 = -2$ nm⁻⁴; $j = -D\times\text{slope} = -(2)\times(-2) = +4$ molecules per nm² per ns. **Two minus signs cancel.** A flux of −4 would mean flow from 4 towards 8, uphill, which diffusion never does.

## 4. The diffusion equation and what the second derivative means

Put Fick's law into the conservation equation (with $D$ constant):

$$\frac{\partial c}{\partial t} = -\frac{\partial}{\partial x}\left(-D\,\frac{\partial c}{\partial x}\right) \quad\Rightarrow\quad \boxed{\;\frac{\partial c}{\partial t} = D\,\frac{\partial^2 c}{\partial x^2}\;}$$

### 4.1 The second derivative is curvature

Take three neighbouring grid points spaced $\Delta x$ apart, with concentrations $c_L$ (left), $c_M$ (middle), $c_R$ (right):

$$\frac{\partial^2 c}{\partial x^2} \approx \frac{c_L - 2c_M + c_R}{\Delta x^2}$$

It is **positive when the middle point is below the average of its neighbours** (a valley), and negative when it is above (a peak).

| Profile | $c_L - 2c_M + c_R$ | Middle point |
|---|---|---|
| 2, 5, 2 (peak) | −6 | **falls** |
| 3, 1, 3 (valley) | +4 | **rises** |
| 2, 2, 2 (flat) | 0 | no change |
| 1, 2, 3 (straight slope) | 0 | **no change, although there is a flux** |

**Diffusion flattens peaks and fills valleys.**

### 4.2 The hop picture of the second derivative

In a short time, let a small fraction $f$ of the molecules in each box hop left, and the same fraction hop right. The middle box **gains** $f\,c_L + f\,c_R$ from its neighbours and **loses** $2f\,c_M$:

$$\text{net change} = f\,(c_L + c_R - 2c_M)$$

That bracket is exactly the top of the second-derivative formula. **The second derivative literally is "what my neighbours send in minus what I send out".**

- Straight slope 8, 6, 4: gain $f(8 + 4) = 12f$, loss $2f\times6 = 12f$, net **0**. There is a flux (+4), but 4 flows in from the left and 4 flows out to the right, so nothing piles up. **Concentration changes only where the flux itself changes.**
- Valley 3, 1, 3: gain $f(3 + 3) = 6f$, loss $2f\times1 = 2f$, net $+4f$: the middle **rises**.

> **Common mistake (made while studying this).** For 3, 1, 3 with $D = 2$ nm²/ns at 1 nm spacing: $(3 - 2 + 3)/1^2 = +4$, so $\partial c/\partial t = 2\times4 = +8$ molecules/nm³ per ns. The **size** 8 is right, but the sign is **+**: the middle **rises**. It is a valley; molecules flow in from both sides.

## 5. The same thing in 3D: gradient, divergence, Laplacian

Now position is a vector $\vec r = (x, y, z)$ and the flux is a vector $\vec j = (j_x, j_y, j_z)$. Define the operator "del":

$$\nabla = \left(\frac{\partial}{\partial x},\ \frac{\partial}{\partial y},\ \frac{\partial}{\partial z}\right)$$

| Name | Formula | Meaning |
|---|---|---|
| gradient | $\nabla c = \left(\frac{\partial c}{\partial x}, \frac{\partial c}{\partial y}, \frac{\partial c}{\partial z}\right)$ | a vector pointing uphill (where $c$ increases fastest) |
| divergence | $\nabla\cdot\vec j = \frac{\partial j_x}{\partial x} + \frac{\partial j_y}{\partial y} + \frac{\partial j_z}{\partial z}$ | net **outflow** per unit volume (the "out minus in" of Section 2, done in each direction and added) |
| Laplacian | $\nabla^2 c = \nabla\cdot(\nabla c) = \frac{\partial^2 c}{\partial x^2} + \frac{\partial^2 c}{\partial y^2} + \frac{\partial^2 c}{\partial z^2}$ | curvature in all three directions, added |

With these, the lecture's three lines are exactly Sections 2, 3 and 4:

$$\frac{\partial c}{\partial t} = -\nabla\cdot\vec j \qquad\text{(conservation)}$$

$$\vec j = -D\,\nabla c \qquad\text{(Fick: flow opposite to the gradient, i.e. downhill)}$$

$$\frac{\partial c}{\partial t} = D\,\nabla^2 c \qquad\text{(diffusion equation)}$$

**Two things to remember for later:**

1. In an **equilibrium MD run of a pure liquid** the concentration is uniform, so nothing flows. Section 10 explains what diffuses then: one **tagged** molecule.
2. Fick's law assumed that $c$ hardly changes over one hop. **On scales shorter than a few hops (very short times) the diffusion equation is not valid.** This is the ballistic region of the MSD curve (Section 17).

---

# Part II. Spreading from a point

**The picture.** At $t = 0$, put all the ink at one single point. Later it is a cloud: densest in the middle, thinning outwards, and getting wider all the time. This part finds the exact shape of that cloud.

## 6. The starting condition: the Dirac delta

The lecture writes the starting condition as $c(\vec r, 0) = \delta(\vec r)$.

$\delta(\vec r)$, the **Dirac delta**, is a spike at the origin: zero everywhere except at $\vec r = 0$, infinitely tall there, with total volume integral exactly 1:

$$\int_{\text{all space}}\delta(\vec r)\,d\vec r = 1, \qquad d\vec r \equiv dx\,dy\,dz$$

($d\vec r$ is a tiny volume element, as in the lecture.)

- **Think of it as** a bell curve made narrower and narrower while its area is kept at 1.
- **Here it means:** one unit of stuff, all at one point.
- In 3D, $\delta(\vec r) = \delta(x)\,\delta(y)\,\delta(z)$.

## 7. The 1D solution is a Gaussian (and checking it)

### 7.1 The answer

$$c(x,t) = \frac{1}{\sqrt{4\pi Dt}}\exp\left(-\frac{x^2}{4Dt}\right)$$

where $\exp(u)$ means $e^u$.

**Compare it with a normal (Gaussian) distribution** with mean 0 and standard deviation $\sigma$:

$$f(x) = \frac{1}{\sqrt{2\pi\sigma^2}}\exp\left(-\frac{x^2}{2\sigma^2}\right)$$

They match if $2\sigma^2 = 4Dt$:

$$\boxed{\;\sigma^2 = 2Dt, \qquad \sigma = \sqrt{2Dt}\;}$$

**So the cloud is a normal distribution centred on the starting point, whose variance grows linearly with time.**

![Gaussian spreading and random walk](figures/fig_gaussian_randomwalk.png)

*Panel (a): the cloud for water ($D = 0.23$ Å²/ps) after 10, 100 and 1000 ps. Panel (b) is used in Section 9.*

### 7.2 Checking that it solves the diffusion equation

Only the chain rule and the product rule are needed. Write $c = (4\pi D)^{-1/2}\;t^{-1/2}\;\exp(-x^2/4Dt)$.

**Time derivative.** Two factors depend on $t$:

- $t^{-1/2}$ gives a factor $-\frac{1}{2t}$;
- the exponent $-x^2/(4Dt)$ has time derivative $+x^2/(4Dt^2)$.

$$\frac{\partial c}{\partial t} = c\left[-\frac{1}{2t} + \frac{x^2}{4Dt^2}\right]$$

**First x-derivative:**

$$\frac{\partial c}{\partial x} = c\cdot\left(-\frac{x}{2Dt}\right)$$

**Second x-derivative** (product rule: differentiate $c$, then differentiate $-x/2Dt$):

$$\frac{\partial^2 c}{\partial x^2} = c\cdot\frac{x^2}{4D^2t^2} + c\cdot\left(-\frac{1}{2Dt}\right) = c\left[\frac{x^2}{4D^2t^2} - \frac{1}{2Dt}\right]$$

**Multiply by $D$:**

$$D\,\frac{\partial^2 c}{\partial x^2} = c\left[\frac{x^2}{4Dt^2} - \frac{1}{2t}\right] = \frac{\partial c}{\partial t}\quad\checkmark$$

### 7.3 The amount stays 1

A standard integral: $\int_{-\infty}^{\infty}\exp(-x^2/2\sigma^2)\,dx = \sqrt{2\pi\sigma^2}$. So the area under $c$ is $\sqrt{2\pi\sigma^2}/\sqrt{2\pi\sigma^2} = 1$ at every time. **No molecules are lost**, as conservation (Section 2) demands. As $t\to0$, $\sigma\to0$ and the curve becomes the $\delta$ spike.

### 7.4 Where it falls and where it rises

The bracket in $\partial c/\partial t$ can be rewritten using $\sigma^2 = 2Dt$:

$$-\frac{1}{2t} + \frac{x^2}{4Dt^2} = \frac{x^2 - 2Dt}{4Dt^2} = \frac{x^2 - \sigma^2}{4Dt^2}$$

- Inside $\lvert x\rvert < \sigma$: negative, so $c$ **falls** (that part of the curve is a peak).
- Outside $\lvert x\rvert > \sigma$: positive, so $c$ **rises** (the tails curve upwards, like valleys).

**The peak slumps and the tails fill in: that is the spreading.** It is the "peaks fall, valleys rise" rule of Section 4.

## 8. The 3D solution

Because $r^2 = x^2 + y^2 + z^2$, the exponential splits into a product, and $(4\pi Dt)^{-3/2} = \left[(4\pi Dt)^{-1/2}\right]^3$:

$$c(\vec r,t) = \frac{1}{(4\pi Dt)^{3/2}}\,e^{-r^2/4Dt} = c_1(x,t)\;c_1(y,t)\;c_1(z,t)$$

where $c_1$ is the 1D solution of Section 7. **The three directions spread independently, each with variance $2Dt$.**

The power $3/2$ is $d/2$, with $d = 3$ the number of dimensions. On a surface ($d = 2$) the prefactor is $(4\pi Dt)^{-1}$; along a line ($d = 1$) it is $(4\pi Dt)^{-1/2}$.

## 9. Why the width grows as √t: the random walk

In the hop model of Section 3, $D = \ell^2/(2t_h)$. In a time $t$ a molecule makes $N = t/t_h$ hops. So

$$\sigma^2 = 2Dt = 2\cdot\frac{\ell^2}{2t_h}\cdot t = \ell^2\,\frac{t}{t_h} = N\ell^2, \qquad \sigma = \ell\sqrt N$$

This is the **block-averaging lecture (Lecture 9) again**: each hop is an independent step of $\pm\ell$ (mean 0, variance $\ell^2$), and **variances of independent steps add**, giving $N\ell^2$. The **Gaussian shape is the central limit theorem**: the sum of many independent hops is normally distributed.

Panel (b) of the figure shows the difference this makes: 400 steps all in one direction end 400 steps away; 400 steps in random directions end only about $\sqrt{400} = 20$ steps away (this particular walk: 28).

**Numbers for water**, $D = 2.3\times10^{-9}$ m²/s $= 0.23$ Å²/ps (because 1 m²/s $= 10^{20}$ Å² per $10^{12}$ ps $= 10^8$ Å²/ps):

| $t$ | $\sigma^2 = 2Dt$ | $\sigma$ (per direction) |
|---|---|---|
| 10 ps | 4.6 Å² | 2.1 Å |
| 100 ps | 46 Å² | 6.8 Å |
| 1 ns | 460 Å² | 21 Å |
| 10 ns | 4600 Å² | 68 Å |

**100 times longer only takes the molecule 10 times farther.**

## 10. What $c$ means in an MD simulation

In an equilibrium MD run the concentration is uniform, so there is nothing to flow (Section 5). Instead we **follow one tagged molecule**:

- at $t = 0$ we know exactly where it is, so we put the origin there; $\delta(\vec r)$ now means "certainly at zero displacement";
- at a later time $t$, $c(\vec r, t)\,d\vec r$ is the **probability that its displacement lies in the tiny volume $d\vec r$ around $\vec r$**.

The tagged molecule's spreading obeys the same diffusion equation (on scales longer than a few collisions). This is called **self-diffusion**, and its $D$ is the **self-diffusion coefficient**: the $D$ that MD simulations usually report.

---

# Part III. Mean squared displacement

## 11. What the MSD is, and why we square

The lecture defines

$$\langle r^2(t)\rangle = \int_{\text{all space}} r^2\,c(\vec r,t)\,d\vec r, \qquad r^2 = x^2 + y^2 + z^2$$

This is the ordinary "average = value × probability, integrated over all possibilities", like $\langle X\rangle = \int x\,f(x)\,dx$ in Lecture 9, with $r^2$ as the value and $c$ as the probability density (it integrates to 1, Section 7.3).

Reading the name word by word:

- **displacement:** $\vec r$ = where the molecule is now **minus where it started**;
- **squared:** we use its squared length $r^2$;
- **mean:** averaged over all the places it could be.

$t$ is the **time elapsed since the start** (the lag), not the clock time of the simulation.

### 11.1 Why not just the average displacement?

A 1D walker hops ±1 Å with equal probability. After 2 hops, the 4 equally likely paths end at:

| Path | Displacement | Displacement² |
|---|---|---|
| ++ | +2 | 4 |
| +− | 0 | 0 |
| −+ | 0 | 0 |
| −− | −2 | 4 |
| **Average** | **0** | **2** |

- The **mean displacement is always 0**: going left and going right are equally likely, so they cancel. It says nothing about how fast the walker spreads.
- The **mean squared displacement is 2 Å²** $= N\ell^2$ (2 hops of 1 Å), as Section 9 says. Squaring makes every excursion count positively, whichever way it went.

## 12. MSD = 6Dt straight from the Gaussian

Each direction of the cloud is a Gaussian with mean 0 and variance $2Dt$ (Section 8). For a quantity with mean 0, the average of its square **is** its variance ($\mathrm{Var}(X) = \langle X^2\rangle - \langle X\rangle^2$ and $\langle X\rangle = 0$). So $\langle x^2\rangle = \langle y^2\rangle = \langle z^2\rangle = 2Dt$, and

$$\langle r^2\rangle = \langle x^2\rangle + \langle y^2\rangle + \langle z^2\rangle = 3\times2Dt = 6Dt$$

| Dimensions $d$ | MSD |
|---|---|
| 1 (a line) | $2Dt$ |
| 2 (a surface) | $4Dt$ |
| 3 (a volume) | $6Dt$ |
| general | $2d\,Dt$ |

**This is why the slope of the MSD gives D:** in 3D, $D = \frac16\times$ slope.

## 13. The lecture's derivation of $d\langle r^2\rangle/dt = 6D$

Section 12 used the explicit Gaussian. The lecture proves the result **directly from the diffusion equation**, without solving it. The only tool is **integration by parts**, used twice. We do it in 1D first, where it is ordinary calculus, then repeat it in 3D.

### 13.1 Reminder: integration by parts

$$\int_a^b u\,\frac{dv}{dx}\,dx = \Big[u\,v\Big]_a^b - \int_a^b \frac{du}{dx}\,v\,dx$$

It moves a derivative from one factor ($v$) onto the other ($u$), at the price of a minus sign and a **boundary term** $[uv]$.

### 13.2 The 1D derivation

Start from $\langle x^2\rangle = \int x^2\,c\,dx$ (integrals from $-\infty$ to $+\infty$) and differentiate in time. Only $c$ depends on $t$:

$$\frac{d\langle x^2\rangle}{dt} = \int x^2\,\frac{\partial c}{\partial t}\,dx = D\int x^2\,\frac{\partial^2 c}{\partial x^2}\,dx$$

(the second step uses the diffusion equation).

**First integration by parts** with $u = x^2$ and $v = \partial c/\partial x$:

$$\int x^2\,\frac{\partial^2 c}{\partial x^2}\,dx = \Big[x^2\,\frac{\partial c}{\partial x}\Big]_{-\infty}^{\infty} - \int 2x\,\frac{\partial c}{\partial x}\,dx$$

The boundary term is **zero**: far away, the slope of $c$ dies off (like $e^{-x^2/4Dt}$) much faster than $x^2$ grows.

**Second integration by parts** with $u = 2x$ and $v = c$:

$$\int 2x\,\frac{\partial c}{\partial x}\,dx = \Big[2x\,c\Big]_{-\infty}^{\infty} - \int 2\,c\,dx = 0 - 2\times1 = -2$$

The boundary term is zero because $c \to 0$ far away, and $\int c\,dx = 1$ (Section 7.3).

**Put it together:**

$$\int x^2\,\frac{\partial^2 c}{\partial x^2}\,dx = 0 - (-2) = 2 \qquad\Rightarrow\qquad \boxed{\;\frac{d\langle x^2\rangle}{dt} = 2D\;}$$

which agrees with $\langle x^2\rangle = 2Dt$ from the Gaussian. Notice the signs: the first step gives a minus, the second gives another minus, and they cancel.

### 13.3 The 3D derivation (the lecture's version)

Exactly the same two steps, written with vectors. We need two tools.

**Tool 1: a product rule for the divergence.** For a function $f$ and a vector $\vec A$:

$$\nabla\cdot(f\vec A) = f\,\nabla\cdot\vec A + \nabla f\cdot\vec A$$

(Write out the x-part: $\partial(fA_x)/\partial x = f\,\partial A_x/\partial x + A_x\,\partial f/\partial x$; then add the y- and z-parts.)

**Tool 2: the divergence theorem.** The total outflow through a closed surface equals the sum of the divergence inside it:

$$\int_{\text{volume}}\nabla\cdot\vec F\,d\vec r = \oint_{\text{surface}}\vec F\cdot d\vec S$$

Here $d\vec S$ is a small piece of the surface, as a vector pointing outwards. When the volume is "all space", the surface is a huge sphere at infinity. **Tools 1 and 2 together are 3D integration by parts.**

**Step 1.** Differentiate in time and use the diffusion equation:

$$\frac{d\langle r^2\rangle}{dt} = \int r^2\,\frac{\partial c}{\partial t}\,d\vec r = D\int r^2\,\nabla^2 c\,d\vec r$$

**Step 2.** Tool 1 with $f = r^2$ and $\vec A = \nabla c$ (and $\nabla\cdot\nabla c = \nabla^2 c$):

$$\nabla\cdot(r^2\nabla c) = r^2\,\nabla^2 c + \nabla r^2\cdot\nabla c$$

so

$$\frac{d\langle r^2\rangle}{dt} = D\int\nabla\cdot(r^2\nabla c)\,d\vec r - D\int\nabla r^2\cdot\nabla c\,d\vec r$$

**Step 3.** The first integral becomes a surface integral (Tool 2), $D\oint r^2\,\nabla c\cdot d\vec S$, on a sphere at infinity. It is **zero** because $\nabla c \to 0$ far from the origin (the lecture's red note "$\nabla c = 0$ at surfaces far away from origin"). For the gradient of $r^2$:

$$\nabla r^2 = \left(\frac{\partial r^2}{\partial x}, \frac{\partial r^2}{\partial y}, \frac{\partial r^2}{\partial z}\right) = (2x, 2y, 2z) = 2\vec r$$

So

$$\frac{d\langle r^2\rangle}{dt} = -2D\int\vec r\cdot\nabla c\,d\vec r$$

**Step 4.** Tool 1 again, now with $f = c$ and $\vec A = \vec r$:

$$\nabla\cdot(c\,\vec r) = c\,\nabla\cdot\vec r + \nabla c\cdot\vec r \quad\Rightarrow\quad \vec r\cdot\nabla c = \nabla\cdot(c\,\vec r) - c\,\nabla\cdot\vec r$$

so

$$\frac{d\langle r^2\rangle}{dt} = -2D\int\nabla\cdot(c\,\vec r)\,d\vec r + 2D\int c\,\nabla\cdot\vec r\,d\vec r$$

**Step 5.** The first integral is again a surface integral, $-2D\oint c\,\vec r\cdot d\vec S$, which is **zero** because $c \to 0$ far away (the lecture's red note "$c = 0$ at far away distances"). For the second:

$$\nabla\cdot\vec r = \frac{\partial x}{\partial x} + \frac{\partial y}{\partial y} + \frac{\partial z}{\partial z} = 1 + 1 + 1 = 3 \;(= d,\ \text{the number of dimensions})$$

and $\int c\,d\vec r = 1$. Therefore

$$\boxed{\;\frac{d\langle r^2\rangle}{dt} = 2D\times3\times1 = 6D\;}$$

**Step 6.** Integrate in time. At $t = 0$ everything is at the origin, so $\langle r^2(0)\rangle = 0$, and

$$\langle r^2(t)\rangle = 6Dt$$

**Why the surface terms really vanish:** on a sphere of radius $R$, the surface area grows like $R^2$ and the factors $r^2$ or $\vec r$ grow like a power of $R$, but $c$ and $\nabla c$ fall off like $e^{-R^2/4Dt}$, which beats any power.

### 13.4 What this derivation tells you

- **The slope of the MSD is $6D$** for any $c$ that obeys the diffusion equation and vanishes far away; the Gaussian shape was never used.
- The lecture's "LHS ≡ slope of MSD" is this statement: $D = \frac16\,d\langle r^2\rangle/dt$.
- **It is only true at long times**, because the diffusion equation itself fails at very short times (Section 5, and Section 17). The safe statement is
$$D = \lim_{t\to\infty}\frac16\,\frac{d\langle r^2\rangle}{dt}$$

## 14. Computing the MSD from an MD trajectory

### 14.1 The formula

In MD there is no $c(\vec r,t)$ to integrate; there are trajectories. Write $\vec r_i(t)$ for the position of molecule $i$ ($i = 1, \dots, N$) at time $t$. **Every molecule is a tagged molecule**, so average over all of them:

$$\text{MSD}(t) = \frac1N\sum_{i=1}^{N}\left\lvert\vec r_i(t_0 + t) - \vec r_i(t_0)\right\rvert^2$$

$t_0$ is the **time origin**: the moment the stopwatch is started.

**Correction to the lecture's wording.** The notes write $r^2 = x^2 + y^2 + z^2$ as if $x$ were a position. In MD, $x$ must be the **displacement** $x_i(t_0 + t) - x_i(t_0)$. Squaring raw positions measures the distance from the corner of the box, which means nothing.

**Example** (positions in Å):

| Molecule | $\vec r_i(t_0)$ | $\vec r_i(t_0 + t)$ | Displacement | $r^2$ |
|---|---|---|---|---|
| 1 | (0, 0, 0) | (1, 2, 2) | (1, 2, 2) | 9 |
| 2 | (5, 5, 5) | (5, 2, 1) | (0, −3, −4) | 25 |
| 3 | (2, 8, 1) | (3, 8, 3) | (1, 0, 2) | 5 |

MSD $= (9 + 25 + 5)/3 =$ **13 Å²**.

> **Your question:** *"In this equation, at a time step, do we calculate it for all the molecules and then average?"*
>
> **Yes.** For one lag $t$: take each molecule's squared displacement over that lag, add them up, divide by $N$. That gives **one point** of the MSD curve. Then do the same for the next lag.
>
> The key idea: **the MSD is not a property of one time step** (like temperature is). Every MSD value needs **two frames**, a start and an end, and $t$ is the gap between them.

### 14.2 Way 1: one time origin

Start the stopwatch at frame 0. At every later frame, compute each molecule's squared displacement from frame 0, and average over molecules.

**Tiny example:** two molecules, 1D, frames 1 ps apart, positions in Å:

| Frame (ps) | 0 | 1 | 2 | 3 |
|---|---|---|---|---|
| molecule A | 0 | 1 | 1 | 3 |
| molecule B | 5 | 4 | 6 | 6 |

| Lag $t$ | A: $(\Delta x)^2$ | B: $(\Delta x)^2$ | MSD$(t)$ |
|---|---|---|---|
| 1 ps | $(1-0)^2 = 1$ | $(4-5)^2 = 1$ | 1.0 |
| 2 ps | $(1-0)^2 = 1$ | $(6-5)^2 = 1$ | 1.0 |
| 3 ps | $(3-0)^2 = 9$ | $(6-5)^2 = 1$ | 5.0 |

### 14.3 Way 2: many time origins

In equilibrium, frame 0 is not special: **any frame can be the start**.

> **Your question:** *"I don't understand the sum over k part."*
>
> $k$ just **numbers the different starting frames**.
>
> MSD at a lag of 1 ps asks: **"how far does a molecule move in 1 ps?"** Any 1-ps stretch of the trajectory can answer that, not only the first one. With frames 0, 1, 2, 3 there are three 1-ps windows:
>
> ```
> frame:     0     1     2     3
> k = 1:     |-----|                  start at frame 0, end at frame 1
> k = 2:           |-----|            start at frame 1, end at frame 2
> k = 3:                 |-----|      start at frame 2, end at frame 3
> ```
>
> $t_{0,k}$ is the **start time of window $k$**: $t_{0,1} = 0$, $t_{0,2} = 1$, $t_{0,3} = 2$ ps. $K = 3$ windows fit.
>
> **Every term written out** (lag 1 ps, $N = 2$, $K = 3$):
>
> | Window $k$ | Start → end | A: $(\Delta x)^2$ | B: $(\Delta x)^2$ | Inner sum (over molecules) | Window MSD (÷ N) |
> |---|---|---|---|---|---|
> | 1 | 0 → 1 | $(1-0)^2 = 1$ | $(4-5)^2 = 1$ | 2 | 1 |
> | 2 | 1 → 2 | $(1-1)^2 = 0$ | $(6-4)^2 = 4$ | 4 | 2 |
> | 3 | 2 → 3 | $(3-1)^2 = 4$ | $(6-6)^2 = 0$ | 4 | 2 |
>
> Outer sum over windows: $2 + 4 + 4 = 10$. Divide by $NK = 2\times3 = 6$: **MSD(1 ps) = 1.67 Å²**. (Same thing in two steps: the window MSDs are 1, 2, 2, and their average is 1.67.)
>
> **So the sum over $k$ means: do Way 1 several times, each time starting at a different frame, and average the results.**

The general formula, with $K$ origins $t_{0,1}, \dots, t_{0,K}$ that fit for this lag:

$$\text{MSD}(t) = \frac{1}{NK}\sum_{k=1}^{K}\sum_{i=1}^{N}\left\lvert\vec r_i(t_{0,k} + t) - \vec r_i(t_{0,k})\right\rvert^2$$

It is just an average over **all (molecule, origin) pairs** for that lag. Every origin contributes the same $N$ molecules, so averaging over molecules first or origins first gives the same answer.

**The full table for the example:**

| Lag $t$ | Windows (origins) | Squared displacements | MSD, many origins | MSD, one origin (Way 1) |
|---|---|---|---|---|
| 1 ps | 0→1, 1→2, 2→3 | A: 1, 0, 4; B: 1, 4, 0 | 10/6 = 1.67 | 1.0 |
| 2 ps | 0→2, 1→3 | A: 1, 4; B: 1, 4 | 10/4 = 2.5 | 1.0 |
| 3 ps | 0→3 only | A: 9; B: 1 | 10/2 = 5.0 | 5.0 |

**Why $K$ depends on the lag:** the window must fit inside the run. Using every frame as a start, $K$ = (number of frames) − (lag in frames).

**Why bother:** with one origin you get $N$ squared displacements per lag; with $K$ origins, $N\times K$, so the curve is much less noisy. (Like measuring how far someone walks in one minute: you don't time only the first minute of an hour's walk, you time every minute and average.)

**Two consequences:**

- **Long lags are noisy:** fewer windows fit (3 for 1 ps, only 1 for 3 ps). Trust the MSD only up to a fraction of the run length, and leave the last part out of fits.
- **Neighbouring origins are correlated** (Lecture 9): windows 0→1 and 1→2 share frame 1. So 1000 origins are worth fewer than 1000 independent measurements, and people often space the origins apart.

### 14.4 In Python

```python
import numpy as np

def msd(pos, max_lag):
    """pos: UNWRAPPED positions, shape (n_frames, N, 3).
    Returns the MSD for lags 0..max_lag (in frames), averaged over molecules and all origins."""
    out = np.zeros(max_lag + 1)
    for lag in range(1, max_lag + 1):
        d = pos[lag:] - pos[:-lag]             # displacement for every origin at once: (n_frames-lag, N, 3)
        out[lag] = (d**2).sum(axis=2).mean()   # |dr|^2, then average over origins AND molecules
    return out
```

**How the slicing makes the windows:** for `lag = 1`, `pos[1:]` is frames 1, 2, 3 (the window ends) and `pos[:-1]` is frames 0, 1, 2 (the window starts). Subtracting pairs them as 1−0, 2−1, 3−2: exactly windows $k = 1, 2, 3$. Then `.mean()` divides by $N\times K$. On the tiny example (as an array of shape (4, 2, 1)) it returns 1.67, 2.5, 5.0.

## 15. How LAMMPS computes the MSD

> **Your question:** *"So how does LAMMPS do it?"*
>
> LAMMPS's `compute msd` uses **only one time origin** (Way 1). There is no sum over $k$ ($K = 1$).

### 15.1 What `compute msd` does

1. **When the `compute msd` line is executed** in the input script, LAMMPS stores the position of every atom in the group as the reference. It stores them **unwrapped**, using the image flags (Section 16), so crossing the box boundary is handled.
2. **Every time the compute is evaluated later** (at a thermo output or a `fix ave/time` output), it calculates each atom's $\Delta x^2$, $\Delta y^2$, $\Delta z^2$ from that reference and **averages over the atoms** in the group.
3. It returns 4 numbers:

| Output | Meaning |
|---|---|
| `c_msd[1]` | $\langle\Delta x^2\rangle$ |
| `c_msd[2]` | $\langle\Delta y^2\rangle$ |
| `c_msd[3]` | $\langle\Delta z^2\rangle$ |
| `c_msd[4]` | total, $\langle\Delta x^2 + \Delta y^2 + \Delta z^2\rangle$: the MSD |

The lag is $t$ = (current step − step at which the compute was issued) × timestep.

> **Your question:** *"In which step is it divided by N?"*
>
> **Inside `compute msd`; you never divide by $N$ yourself.** The LAMMPS documentation says each value is "averaged over atoms in the group" and calls the output *intensive* (it does not grow with the number of atoms). The logic, simplified:
>
> ```
> # when the 'compute msd' line is executed
> ref = unwrapped positions of the group atoms      # r_i(t0)
> N   = number of atoms in the group
>
> # every time c_msd is evaluated
> sx = sy = sz = 0
> for i in group:
>     dx, dy, dz = unwrapped_pos[i] - ref[i]        # r_i(t0 + t) - r_i(t0)
>     sx += dx*dx;  sy += dy*dy;  sz += dz*dz       # the sum over i
> c_msd = [sx/N, sy/N, sz/N, (sx+sy+sz)/N]          # <-- the 1/N happens here
> ```
>
> (In a parallel run, each processor sums its own atoms, the partial sums are added across processors, then divided by $N$.)
>
> Do not confuse the different divisions:
>
> | Where | Divides by | Meaning |
> |---|---|---|
> | inside `compute msd` | $N$ | average over atoms |
> | `.mean()` in the `msd()` Python function | $N\times K$ | average over atoms **and** origins |
> | `slope / 6` when fitting | 6 | from MSD $= 6Dt$ (2 × 3 dimensions); nothing to do with $N$ |

### 15.2 A minimal input and the fit

```
# ... equilibration run above ...
reset_timestep  0                          # so that step number = lag
compute         msd all msd com yes        # reference stored NOW; com yes removes centre-of-mass drift
fix             msdout all ave/time 100 1 100 c_msd[1] c_msd[2] c_msd[3] c_msd[4] file msd.dat
run             500000
```

`ave/time 100 1 100` writes the instantaneous values every 100 steps. Then fit in Python:

```python
import numpy as np
d   = np.loadtxt("msd.dat")        # columns: step, msd_x, msd_y, msd_z, msd_total ('#' lines skipped)
dt  = 0.001                        # timestep in ps (metal units)
t   = d[:, 0] * dt                 # lag in ps (correct because of reset_timestep 0)
msd = d[:, 4]                      # Angstrom^2
mask = t > 10                      # fit only the straight (diffusive) part -- choose it from a log-log plot
slope = np.polyfit(t[mask], msd[mask], 1)[0]   # Angstrom^2/ps
D = slope / 6                      # Angstrom^2/ps; multiply by 1e-8 for m^2/s
```

You can also get the slope inside LAMMPS (`fix vector` plus a variable using `slope()`, as in `examples/DIFFUSE`), but the LAMMPS developers recommend writing the MSD out and fitting afterwards, because the early, non-linear part must be left out of the fit.

**Units:** in LAMMPS `metal` units, lengths are Å and times ps, so $D$ comes out in Å²/ps; $1$ Å²/ps $= 10^{-8}$ m²/s.

### 15.3 Getting many time origins with LAMMPS

**Option A: start several MSD computes during the run.** Each new `compute msd` takes the positions at that moment as its reference, so each is a new origin $k$:

```
compute  msd1 all msd com yes
fix      f1 all ave/time 100 1 100 c_msd1[4] file msd1.dat
run      100000
compute  msd2 all msd com yes          # origin k = 2
fix      f2 all ave/time 100 1 100 c_msd2[4] file msd2.dat
run      100000
# ... and so on
```

Afterwards shift each file so its own start is $t = 0$, and average the curves at equal lag. The origins are far apart, so they are nearly independent (good, by Lecture 9).

**Option B: dump unwrapped positions and post-process** with the `msd()` function of Section 14.4:

```
dump         traj all custom 1000 traj.lammpstrj id type xu yu zu
dump_modify  traj sort id      # essential: the same atom in the same row in every frame
```

Without `sort id`, atoms can appear in a different order in each frame, and `pos[frame, i]` would mix up different atoms. Dump often enough for the times you care about (every 1000 steps at 1 fs gives lags in multiples of 1 ps).

### 15.4 Things to watch with `compute msd`

- **Issue it after equilibration.** The reference is taken when the line is executed, not when the run starts.
- **The number of atoms in the group must not change** during the run. The documentation names `fix deposit` and `fix evaporate` specifically. In a deposition simulation, compute the MSD in a stage without deposition, or post-process a dump (Option B).
- **Do not use `average yes` for diffusion.** It measures displacements from a running-average position, meant for atoms that do not diffuse (crystals).
- **Restarts:** use the same compute ID so the stored reference positions are read back.

## 16. Practical traps: periodic boundaries, drift, species

**1. Periodic boundaries: use unwrapped coordinates.** When a molecule leaves through one face of the box, it re-enters through the opposite face, and its stored ("wrapped") coordinate jumps by the box length $L$.

Example: $L = 30$ Å. A molecule at $x = 29.5$ Å moves +1 Å and is wrapped to $x = 0.5$ Å.

- Naive displacement: $0.5 - 29.5 = -29$ Å, so $\Delta x^2 = 841$ Å².
- True displacement: $+1$ Å, so $\Delta x^2 = 1$ Å².

One crossing ruins the average. **Unwrapped** coordinates keep counting past the edge: $x_{\text{unwrapped}} = x + n_x L$, where the **image flag** $n_x$ is an integer counting net crossings in x. (LAMMPS: `xu yu zu` in a dump, or `x` plus `ix`; `compute msd` unwraps automatically.)

**2. Drift of the whole system.** If the centre of mass slowly drifts (momentum not kept at zero), every molecule picks up the same extra displacement, adding a spurious $\propto t^2$ term. Subtract it (`com yes`).

**3. One MSD per species.** Each component of a mixture has its own $D$. Use a group per species.

---

# Part IV. Reading the MSD curve

The lecture sketches the MSD of a fluid: a curved start near the origin ($\langle r^2\rangle \sim t^2$, "ballistic trajectory, occurs before collisions take place between molecules"), then a straight line ($\langle r^2\rangle \sim t$, "diffusive trajectory") with slope $6D$. Sections 17 and 18 explain both parts from scratch; Sections 19–24 look at real MSD curves of four phases.

## 17. Ballistic motion: what it means and why MSD ∝ t²

> **Your question:** *"What does ballistic basically mean? Teach it slowly, then explain why the MSD is proportional to $t^2$."*

### 17.1 What "ballistic" means

**The picture.** A billiard ball rolling across the table: until it hits another ball or the cushion, nothing pushes it sideways, so it keeps going **in a straight line at the same speed**. That is Newton's first law: **no force, no change of velocity.**

An atom does the same for a short while. Between collisions the forces on it are small, so it keeps the velocity it had and flies straight.

**"Ballistic" = free flight: the atom is still moving with the velocity it started with, because it has not collided with anything yet.** The word comes from *ballistics*, the study of projectiles (a bullet after it leaves the barrel, a thrown ball).

The opposite is **diffusive**: after many collisions the direction has been changed again and again at random, and the path is a zig-zag.

| | Ballistic | Diffusive |
|---|---|---|
| collisions so far | none | many |
| path | straight line | random zig-zag |
| velocity | the same as at the start | long forgotten |

### 17.2 One atom, one direction

Look along $x$ only. $v_x$ is the atom's x-velocity at the start, and $t$ the time since. With constant velocity:

$$\Delta x = v_x\,t \qquad\Rightarrow\qquad \Delta x^2 = v_x^2\,t^2$$

With $v_x = 2$ σ per time unit:

| $t$ | $\Delta x = 2t$ | $\Delta x^2$ |
|---|---|---|
| 1 | 2 | 4 |
| 2 | 4 | 16 |
| 3 | 6 | 36 |

The **distance** grows like $t$ (×2 when $t$ doubles); the **squared distance** grows like $t^2$ (×4 when $t$ doubles). That is all "∝ t²" means.

### 17.3 Many atoms: averaging keeps the t²

Different atoms move with different velocities, some right, some left. Four atoms with $v_x = +2, -1, +1, -2$:

| Atom | $v_x$ | $\Delta x = v_x t$ | $\Delta x^2$ |
|---|---|---|---|
| 1 | +2 | $+2t$ | $4t^2$ |
| 2 | −1 | $-t$ | $t^2$ |
| 3 | +1 | $+t$ | $t^2$ |
| 4 | −2 | $-2t$ | $4t^2$ |
| **Average** | | **0** | $\frac{4+1+1+4}{4}\,t^2 = 2.5\,t^2$ |

- The average $\Delta x$ is **0** (left and right cancel); that is why we square.
- **Every** atom's $\Delta x^2$ is (a number) × $t^2$, so the average is also (a number) × $t^2$. Averaging only fixes the number in front, the average of $v_x^2$:

$$\text{MSD}_x = \langle v_x^2\rangle\,t^2 \qquad(\text{here } 2.5\,t^2)$$

### 17.4 All three directions

$$\Delta r^2 = \Delta x^2 + \Delta y^2 + \Delta z^2 = (v_x^2 + v_y^2 + v_z^2)\,t^2 = v^2t^2 \qquad\Rightarrow\qquad \text{MSD} = \langle v^2\rangle\,t^2$$

### 17.5 What sets ⟨v²⟩: the temperature

**Equipartition** says that, in equilibrium, the average kinetic energy in each direction is $\frac12k_BT$ ($T$ the temperature, $m$ the mass, $k_B$ Boltzmann's constant):

$$\tfrac12 m\langle v_x^2\rangle = \tfrac12 k_BT \quad\Rightarrow\quad \langle v_x^2\rangle = \frac{k_BT}{m} \quad\Rightarrow\quad \langle v^2\rangle = \frac{3k_BT}{m}$$

$$\boxed{\;\text{MSD}_{\text{ballistic}} = \frac{3k_BT}{m}\,t^2\;} \qquad (= 3Tt^2 \text{ in Lennard-Jones reduced units, where } m = k_B = 1)$$

**Argon at 85 K** ($m = 6.63\times10^{-26}$ kg): $\langle v^2\rangle = 3\times1.381\times10^{-23}\times85/6.63\times10^{-26} = 5.3\times10^4$ m²/s², an rms speed of about 230 m/s.

| $t$ | MSD $= 5.3\times10^4\,t^2$ | rms distance |
|---|---|---|
| 0.1 ps | 0.053 Å² | 0.23 Å |
| 0.2 ps | 0.21 Å² (×4) | 0.46 Å |

On a log–log plot, $\log\text{MSD} = \log\langle v^2\rangle + 2\log t$: a straight line of **slope 2**.

### 17.6 Why it cannot last

The formula assumes the velocity **never changes**. The first collision breaks $\Delta x = v_x t$. In our liquid simulation ($T = 0.705$, so the ballistic prediction is $3Tt^2$):

| $t$ | MSD measured | Ballistic $3Tt^2$ | |
|---|---|---|---|
| 0.005 | 0.00005 | 0.00005 | match |
| 0.03 | 0.00187 | 0.00190 | match |
| 0.06 | 0.00705 | 0.00761 | starting to fall below |
| 0.125 | 0.0247 | 0.0330 | clearly below: collisions have begun |

A collision turns the atom sideways or back, so its later motion no longer adds up in line with its earlier motion, and it ends up closer to its start than straight flight would have taken it. **The MSD falls below the ballistic line.**

## 18. Diffusive motion: why MSD ∝ t

> **Your question:** *"I understand the ballistic derivation. Now tell me how it is proportional to $t$ in the diffusive region."*

### 18.1 The picture

After many collisions the path is a **chain of short flights**, each collision sending the atom off in a new direction unrelated to the old one.

- Call each flight a **step** $\vec s_k$ (a vector: length and direction; $k = 1, 2, \dots$).
- Each step takes a time $t_c$ (roughly the time between collisions).
- After a time $t$ the atom has made $N = t/t_c$ steps.

The total displacement is the sum of the steps:

$$\Delta\vec r = \vec s_1 + \vec s_2 + \dots + \vec s_N$$

### 18.2 Square two steps

$$\left\lvert\vec s_1 + \vec s_2\right\rvert^2 = s_1^2 + s_2^2 + 2\,\vec s_1\cdot\vec s_2$$

The last term is the **cross term**. ($\vec s_1\cdot\vec s_2 = s_1s_2\cos\theta$, with $\theta$ the angle between the steps: positive if they point roughly the same way, negative if opposite.)

**1D example, steps of ±1, all 4 combinations equally likely:**

| Step 1 | Step 2 | End point | End² | Cross term $2s_1s_2$ |
|---|---|---|---|---|
| +1 | +1 | +2 | 4 | +2 |
| +1 | −1 | 0 | 0 | −2 |
| −1 | +1 | 0 | 0 | −2 |
| −1 | −1 | −2 | 4 | +2 |
| **Average** | | | **2** | **0** |

The cross terms cancel on average, leaving MSD $= 1 + 1 = 2$: **just the sum of the squares.**

**Ballistic, for comparison:** the second step always points the same way as the first, so only $(+1, +1)$ happens: $2^2 = 4$. The cross term is always $+2$ and never cancels.

### 18.3 N steps: count the terms

Squaring a sum of $N$ steps gives $N^2$ terms:

- $N$ **square terms** $s_k^2$ (always positive);
- $N(N-1)$ **cross terms** $\vec s_j\cdot\vec s_k$ with $j \ne k$.

| | Square terms | Cross terms | MSD (steps of length $\ell$) | Grows like |
|---|---|---|---|---|
| **Ballistic** (all steps aligned) | $N\ell^2$ | each $= +\ell^2$, adding $N(N-1)\ell^2$ | $N^2\ell^2$ | $t^2$ |
| **Diffusive** (random directions) | $N\ell^2$ | each averages to 0 | $N\ell^2$ | $t$ |

**This is the whole story: the $t^2$ comes from the cross terms. Random directions wipe them out, and only the $N$ square terms are left.** With $N = t/t_c$:

$$\text{MSD} = N\ell^2 = \frac{\ell^2}{t_c}\,t \;\propto\; t$$

It is the rule of Lecture 9 again: **the variance of a sum of independent pieces is the sum of their variances.**

**Check by listing every path** (1D, steps ±1):

| $N$ | Random walk MSD | Ballistic MSD |
|---|---|---|
| 1 | 1 | 1 |
| 2 | 2 | 4 |
| 3 | 3 (8 paths: +3 once, +1 three times, −1 three times, −3 once; (9+3+3+9)/8) | 9 |
| 4 | 4 | 16 |

Panel (b) of the figure in Section 7 shows the same thing with 400 steps: 400 versus about 20.

### 18.4 The link to D

In 3D, MSD $= 6Dt$, so

$$6D = \frac{\ell^2}{t_c}, \qquad D = \frac{\ell^2}{6\,t_c}$$

**D is set by how long the steps are and how often they happen.** (In 1D the same argument gives $D = \ell^2/2t_h$, Section 3.)

**Numbers:** steps of 0.5 Å every 0.2 ps.

| $t$ | $N$ | Random walk MSD $= N\times0.25$ Å² | All steps aligned, $(N\times0.5)^2$ |
|---|---|---|---|
| 1 ps | 5 | 1.25 Å² | 6.25 Å² |
| 2 ps | 10 | 2.5 Å² (×2) | 25 Å² (×4) |
| 10 ps | 50 | 12.5 Å² | 625 Å² |

$D = 1.25/(6\times1) = 0.21$ Å²/ps $\approx 2\times10^{-9}$ m²/s, about water's value.

**In our liquid simulation:** from $t = 17.1$ to $t = 34.5$ (×2.02 in time) the MSD went from 3.17 to 6.40 σ² (×2.02). **Double the time, double the MSD.**

### 18.5 A real liquid is not a perfect random walk, and that does not matter

In a liquid an atom often **bounces back** off a neighbour, so a step is slightly anti-correlated with the next. The nearby cross terms are then not exactly zero. But each step is correlated only with its **few** neighbours, so the leftover cross terms add a **fixed amount per step**: still proportional to $N$, still linear. **Only the slope (D) changes.** This is exactly Lecture 9's statistical inefficiency, $g = 1 + 2\sum\rho_k$; the bounce-back gives negative $\rho$ and makes D smaller. Part V makes this exact.

**Why the curve switches from $t^2$ to $t$:** at short times the atom is still inside its first flight (all "steps" aligned: $t^2$); after many collisions the cross terms have died out ($t$).

## 19. The simulations used in these notes

To see real MSD curves, four small MD simulations were run (Python with Numba; Lennard-Jones atoms, truncated and shifted at $2.5\sigma$; velocity-Verlet with timestep 0.005; equilibrated with velocity rescaling, then **NVE** production; MSD from unwrapped positions, averaged over all atoms and many time origins).

| System | Model | Atoms | Density $\rho$ | $T$ |
|---|---|---|---|---|
| gas | Lennard-Jones (LJ) | 500 | 0.05 | 2.0 |
| liquid | LJ | 500 | 0.844 | 0.71 |
| supercooled liquid | Kob–Andersen 80:20 binary LJ mixture (a standard model that resists crystallizing; plain LJ would freeze) | 1000 | 1.2 | 0.60 |
| fcc crystal | LJ | 500 | 1.0 | 0.50 |

**Units (LJ reduced units):** lengths in $\sigma$, energies in $\epsilon$, temperature in $\epsilon/k_B$, mass in $m$, time in $\sigma\sqrt{m/\epsilon}$. For argon: $\sigma = 3.405$ Å, $\epsilon/k_B = 119.8$ K, time unit $= 2.16$ ps, unit of $D$ $= \sigma^2/\text{time} = 5.38\times10^{-8}$ m²/s.

![MSD of four phases](figures/fig_msd_phases.png)

*(a) MSD versus t on log–log axes. (b) Local slope of (a). (c) Liquid on linear axes with the fit. (d) Crystal on linear axes at short times.*

## 20. Log–log plots and the local slope α

If MSD $= A\,t^\alpha$ for some constant $A$, then $\ln\text{MSD} = \ln A + \alpha\ln t$: a straight line on log–log axes whose **slope is the exponent** $\alpha$. Panel (b) plots the **local slope**

$$\alpha(t) = \frac{d\ln\text{MSD}}{d\ln t}$$

| $\alpha$ | Meaning |
|---|---|
| 2 | ballistic: straight flight |
| 1 | diffusive: random walk |
| between 0 and 1 | subdiffusive: something holds the atom back |
| 0 | trapped: MSD constant |
| above 1 at long times | warning sign (Section 24) |

On linear axes the first 0.1 time units are squashed into the corner and a crystal's small plateau looks like zero. **Plot log–log to find the regimes; fit on linear axes** (panel c).

**The ballistic start is the same in every phase.** At $t = 0.02$:

| Phase | $T$ | MSD measured | $3Tt^2$ |
|---|---|---|---|
| gas | 2.016 | 0.002414 | 0.002419 |
| liquid | 0.705 | 0.000838 | 0.000846 |
| supercooled | 0.598 | 0.000709 | 0.000718 |
| crystal | 0.503 | 0.000599 | 0.000604 |

**The ballistic part depends only on $T$ and $m$, never on the phase.** It is a check of the temperature and says nothing about structure. (In panel (a) the gas starts higher only because it is hotter.)

**When the ballistic part ends** (where $\alpha$ first drops below 1.8):

| Phase | Ends at |
|---|---|
| gas | $t \approx 2.1$ |
| liquid, crystal, supercooled | $t \approx 0.05$–0.07 (about 0.15 ps for argon) |

In the dense phases the neighbours are already touching, so the first collision comes almost at once. For the gas, kinetic theory gives the time between collisions:

- **mean free path** (average distance between collisions): $\lambda = 1/(\sqrt2\,\pi\rho\sigma^2) = 1/(\sqrt2\,\pi\times0.05) = 4.5\sigma$;
- **mean speed:** $\bar v = \sqrt{8k_BT/\pi m} = \sqrt{8\times2.0/\pi} = 2.26$;
- **time between collisions:** $\lambda/\bar v \approx 2.0$, exactly where the gas leaves the ballistic regime. This is the lecture's red note: ballistic motion "occurs before collision takes place between molecules".

## 21. Gas versus liquid

> **Your question:** *"What is the difference between the gas and liquid MSD plots?"*

Both start ballistic and end linear. They differ in **how long the free flight lasts** and **how steep the line is**:

| | Gas ($\rho$ = 0.05) | Liquid ($\rho$ = 0.84) |
|---|---|---|
| free flight before a collision | long, about $4.5\sigma$ | almost none (neighbours already touching) |
| ballistic part ends at | $t \approx 2$ | $t \approx 0.07$ |
| change from $t^2$ to $t$ | slow, spread over $t \approx 1$–40 | quick, done by $t \approx 0.3$ |
| dip of $\alpha$ below 1 | none | yes, to 0.88 at $t \approx 0.7$ |
| $D$ (slope/6) | 5.72 | 0.0303 (about 190× smaller) |
| MSD at $t = 80$ | about 2,700 σ² (rms ≈ 52σ) | 14.5 σ² (rms ≈ 3.8σ) |

**In one line:** a gas atom takes **long, rare** steps (long free flights), so it stays ballistic longer and spreads fast; a liquid atom is boxed in by neighbours and takes **tiny, frequent** steps, so it becomes a random walk almost at once and spreads slowly. (The gas was also hotter, which adds to the difference, but density is the main cause.)

**More detail on each:**

- **Gas.** Fitting $t > 40$ gives $D = 5.72$. Hard-sphere kinetic theory (Chapman–Enskog), $D = \frac{3}{8\rho\sigma^2}\sqrt{k_BT/\pi m}$, gives 6.0: close, but LJ atoms are not hard spheres, so no exact match is expected. For argon, about $3\times10^{-7}$ m²/s. The bend from $\alpha = 2$ to 1 takes more than a decade of time because **one collision does not erase the velocity; it takes several.**
- **Fit the slope, do not divide MSD by 6t.** The gas fit line has intercept −85 σ², so MSD/(6t) at $t = 100$ gives 5.55, 3% low. The slope is not affected by the intercept.
- **Liquid.** The dip of $\alpha$ to 0.88 is weak **caging**: an atom bounces off its neighbours before it can slip past them. Fitting $t > 5$ gives slope 0.182, **$D = 0.0303$** (panel c), which for argon is $1.6\times10^{-9}$ m²/s at 84 K and 1419 kg/m³ (liquid-argon conditions). By $t = 80$ each atom has moved several neighbour distances: **that is what makes it a liquid.**
- **Crossover time.** The ballistic line $3Tt^2$ and the diffusive line $6Dt$ meet where $3Tt^2 = 6Dt$, at $t_c = 2D/T$ (in reduced units; in general $2mD/k_BT$): 5.7 for the gas, 0.086 for the liquid. It measures how long an atom remembers its velocity (Part V).

## 22. Supercooled liquid: the cage

Three stages are visible:

1. **Ballistic** up to $t \approx 0.05$.
2. **Caged:** $\alpha$ falls to **0.29** at $t \approx 1.2$, where MSD ≈ 0.05 σ² (rms 0.22σ). The atom rattles inside a **cage** formed by its neighbours.
3. **Escape:** the neighbours slowly rearrange, the atom gets out, and $\alpha$ climbs back, reaching 0.98 only at $t \approx 300$. $D \approx 1.2\times10^{-3}$.

**What you read off:** the plateau height gives the **cage size**; the time at which $\alpha$ climbs back gives the **cage lifetime** (the structural relaxation time). On cooling, the plateau gets longer and flatter. In a **glass** it lasts longer than any simulation, and the MSD looks like a crystal's. **The MSD alone cannot tell a glass from a crystal**; that needs structure (the radial distribution function, or Steinhardt order parameters).

**Practical rule:** in a slow system, check that $\alpha$ has really reached about 1 before fitting D. (Do not compare this D with the liquid's: different model, different density.)

## 23. Crystal: overshoot and plateau

Panel (d): ballistic until $t \approx 0.07$, a **peak of 0.0261 at $t \approx 0.26$**, then a constant **plateau of 0.0212 σ²**. No diffusion: the atoms only vibrate about their lattice sites.

### 23.1 Why the plateau is $2\langle u^2\rangle$

Write the position as $\vec r(t) = \vec R + \vec u(t)$, where $\vec R$ is the fixed lattice site and $\vec u$ the vibration displacement from it (average zero). Then $\Delta\vec r = \vec u(t) - \vec u(0)$ and

$$\left\langle\left\lvert\Delta\vec r\right\rvert^2\right\rangle = \langle u^2\rangle + \langle u^2\rangle - 2\langle\vec u(t)\cdot\vec u(0)\rangle$$

This is Lecture 9's $\mathrm{Var}(X - Y) = \mathrm{Var}X + \mathrm{Var}Y - 2\,\mathrm{Cov}(X, Y)$.

- **$t = 0$:** $\vec u(t)\cdot\vec u(0)$ averages to $\langle u^2\rangle$, so MSD = 0 ✓.
- **Long $t$:** the vibration has been forgotten, the cross term is 0, so **MSD $= 2\langle u^2\rangle$**: the plateau.
- **About half a vibration later:** the atom is on the opposite side of its site, the cross term is **negative**, and MSD $> 2\langle u^2\rangle$: **the overshoot.**

### 23.2 What you read off the plateau

1. **Vibration amplitude:** $\langle u^2\rangle = 0.0212/2 = 0.0106$ σ², rms amplitude $0.103\sigma$ (0.35 Å for argon).
2. **Lindemann ratio** = rms amplitude / nearest-neighbour distance. For fcc at $\rho = 1$ the nearest-neighbour distance is $2^{1/6}\sigma = 1.122\sigma$, so the ratio is $0.103/1.122 = 0.092$. Rule of thumb: **crystals melt when this ratio reaches roughly 0.1–0.15** (the Lindemann criterion). This crystal is below that range.
3. **The plateau is proportional to $T$.** At $T = 0.252$ the plateau was 0.0106, half the value at $T = 0.503$. Model: each atom tied to its site by a spring of stiffness $k_s$. Equipartition per direction, $\frac12k_s\langle u_x^2\rangle = \frac12k_BT$, gives $\langle u^2\rangle = 3k_BT/k_s$.
4. **The vibration time from the plateau alone:**
   - stiffness: $k_s = 3k_BT/\langle u^2\rangle = 3\times0.503/0.0106 = 142$ ε/σ²;
   - angular frequency: $\omega = \sqrt{k_s/m} = 11.9$;
   - period: $2\pi/\omega = 0.53$; half a period is **0.26**, exactly where the overshoot peaks.

   (Real crystals have a spread of vibration frequencies, so this is approximate in general; here it works well.)

Real crystals do diffuse, by atoms hopping into vacancies, but those hops are far too rare to appear in an MD run. A perfect crystal in MD shows a flat plateau.

## 24. Other things the MSD tells you, and warning signs

- **Per-direction MSD** (`c_msd[1]`, `c_msd[2]`, `c_msd[3]`) shows **anisotropy**. Atoms on a substrate: MSD$_x$ and MSD$_y$ grow linearly while MSD$_z$ stays flat (bound to the surface). The **surface diffusion coefficient** is $D_\parallel$ = slope of (MSD$_x$ + MSD$_y$)/4 (4 = 2 × the number of directions counted).
- **Confinement:** for atoms spread evenly across a slab of thickness $H$, the MSD in that direction levels off at $H^2/6$ (the variance of the difference of two independent uniform positions is $2\times H^2/12$).
- **Per-atom MSD:** each atom's squared displacement over a fixed window. In a solid–liquid coexistence box, liquid atoms move about $\sigma$ and solid atoms about $0.1\sigma$, so atoms can be labelled mobile or immobile: a **dynamic** complement to structural order parameters.
- **Warning signs:**
  - **$\alpha > 1$ at long times:** collective motion. Centre-of-mass drift, flow, or the system changing phase. A drift velocity $v_d$ adds $v_d^2t^2$, so the curve bends upwards.
  - **The MSD depends on the origin $t_0$:** not in equilibrium (aging, or melting/crystallizing in progress). Compare MSDs from early and late origins.
  - **Noisy tail:** few origins at long lags; leave it out of fits.
  - **$D$ depends on box size:** periodic images slow diffusion down; the correction needs the viscosity (Section 32).

**Summary of the phases:**

| Phase | $\alpha(t)$ | MSD at long times | What you read off |
|---|---|---|---|
| gas | stays 2 for a long time, then slowly bends to 1 | linear, steep | $D$, collision time |
| liquid | 2, small dip to ~0.9, then 1 | linear | $D$ |
| supercooled liquid | 2, plateau near 0.3, then 1 | linear after cage escape | cage size, cage lifetime, $D$ |
| glass | 2, then a plateau longer than the run | flat (on MD timescales) | cage size |
| crystal | 2, overshoot, then 0 | flat at $2\langle u^2\rangle$ | amplitude, Lindemann ratio, vibration time |

---

# Part V. Green–Kubo relation for diffusion

The lecture: "One can compute transport coefficients from integration of auto-correlation functions. Such relations are called Green–Kubo functions" (usually called Green–Kubo **relations**), with

$$D = \frac13\int_0^\infty\langle\vec v(t)\cdot\vec v(0)\rangle\,dt = \frac13\left[\int_0^\infty\langle v_x(t)v_x(0)\rangle dt + \int_0^\infty\langle v_y(t)v_y(0)\rangle dt + \int_0^\infty\langle v_z(t)v_z(0)\rangle dt\right]$$

This part explains what the integrand is, derives the formula from the MSD (using the same "squares and cross terms" idea as Section 18), and checks it with the simulations.

## 25. The velocity autocorrelation function (VACF)

$$C(t) = \langle\vec v(t)\cdot\vec v(0)\rangle$$

$\vec v(0)$ is an atom's velocity at a time origin, $\vec v(t)$ the **same** atom's velocity a time $t$ later; the average is over atoms and time origins, exactly like the MSD (Section 14). "Auto" means the atom is correlated with **itself**.

**What it asks: "If an atom is moving some way now, how much is it still moving that way a time $t$ later?"**

| $C(t)$ | Meaning |
|---|---|
| positive | still going the same way |
| zero | direction forgotten |
| negative | turned round (bounced back) |

At $t = 0$: $C(0) = \langle v^2\rangle = 3k_BT/m$ (equipartition, Section 17.5).

![Green-Kubo figure](figures/fig_green_kubo.png)

*(a) Normalized VACF, $C(t)/C(0)$, for the gas, liquid and crystal. (b) Running Green–Kubo integral for the liquid and crystal. (c) The same for the gas. (d) Running Green–Kubo integral for the viscosity of the liquid (Part VI).*

Panel (a):

- **Gas:** stays near 1 until $t \approx 0.5$, then fades smoothly; gone by $t \approx 15$. Long free flights, long memory.
- **Liquid:** drops fast, crosses zero at $t = 0.135$, and has a minimum of **−0.13** at $t = 0.19$: the atom has bounced off a neighbour and is moving **backwards**. Memory gone by $t \approx 1$.
- **Crystal:** a deep minimum of **−0.51** at $t = 0.175$ (the atom swings back in its cage), then oscillations.

## 26. Deriving $D = \frac13\int\langle\vec v(t)\cdot\vec v(0)\rangle dt$

**Step 1: chop the time into tiny pieces.** Split the time $t$ into $n$ small intervals of length $\Delta t$ ($n = t/\Delta t$). Over each one the velocity is nearly constant, so the step taken in interval $k$ is $v_{x,k}\,\Delta t$, where $v_{x,k}$ is the x-velocity then. The displacement is the sum of the steps (this is literally what an MD integrator does):

$$\Delta x = \sum_{k=1}^{n}v_{x,k}\,\Delta t$$

**Step 2: square and average** (all squares and all cross terms):

$$\langle\Delta x^2\rangle = \Delta t^2\sum_{j=1}^{n}\sum_{k=1}^{n}\langle v_{x,j}\,v_{x,k}\rangle$$

In equilibrium, $\langle v_{x,j}v_{x,k}\rangle$ depends only on how far apart the two intervals are, $q = \lvert j - k\rvert$. It equals $C_x(q\,\Delta t)$, where $C_x(s) = \langle v_x(s)\,v_x(0)\rangle = C(s)/3$ (the three directions are equivalent).

**Step 3: count the terms.** There are $n$ terms with $q = 0$ (the squares) and $2(n - q)$ terms with each separation $q \ge 1$ (cross terms; for 3 intervals: separation 1 occurs for the pairs (1,2), (2,3), (2,1), (3,2), which is $2\times2 = 4$):

$$\langle\Delta x^2\rangle = \Delta t^2\left[n\,C_x(0) + 2\sum_{q=1}^{n-1}(n - q)\,C_x(q\,\Delta t)\right]$$

**Step 4: the two limits.**

- **$t$ much shorter than the velocity memory:** every $C_x \approx C_x(0)$, so all $n^2$ terms survive, and $\langle\Delta x^2\rangle = n^2\Delta t^2\,C_x(0) = \langle v_x^2\rangle\,t^2$. **The ballistic result** (Section 17).
- **$t$ much longer than the memory:** $C_x$ is zero except for small $q$, where $(n - q) \approx n$:

$$\langle\Delta x^2\rangle \approx \underbrace{n\,\Delta t}_{=\,t}\;\Delta t\left[C_x(0) + 2\sum_{q\ge1}C_x(q\,\Delta t)\right] \;\longrightarrow\; t\cdot2\int_0^\infty C_x(s)\,ds \quad (\Delta t\to0)$$

because $\Delta t\sum_{q\ge1}C_x(q\Delta t)$ becomes the integral, and the single term $\Delta t\,C_x(0)$ vanishes. **Linear in $t$: the diffusive result** (Section 18).

**Step 5: compare with $\langle\Delta x^2\rangle = 2Dt$** (Section 12):

$$D = \int_0^\infty\langle v_x(s)\,v_x(0)\rangle\,ds \qquad\Rightarrow\qquad \boxed{\;D = \frac13\int_0^\infty\langle\vec v(t)\cdot\vec v(0)\rangle\,dt\;}$$

The second line of the lecture is just this with $\vec v\cdot\vec v$ split into its x, y, z parts. In an anisotropic system (atoms on a surface, layered materials) each direction gives its own coefficient, $D_x = \int_0^\infty\langle v_x(t)v_x(0)\rangle dt$, and so on.

**The link to Lecture 9.** This is the same mathematics as the variance of the mean of correlated samples: $\mathrm{Var} = g\,\sigma^2/M$ with $g = 1 + 2\sum\rho_k$. Here the "samples" are velocities, and $\rho_q = C_x(q\Delta t)/C_x(0)$. **Correlated velocities are why the MSD is not simply $n\Delta t^2\langle v_x^2\rangle$.**

**An exact relation.** Letting $\Delta t\to0$ in Step 3 without dropping anything gives, in 3D,

$$\text{MSD}(t) = 2\int_0^t (t - s)\,C(s)\,ds, \qquad \frac{d\,\text{MSD}}{dt} = 2\int_0^t C(s)\,ds$$

So the **running integral** $\frac13\int_0^t C(s)\,ds$ equals $\frac16\,d\text{MSD}/dt$ at the same $t$: it is the "D you would read off the MSD slope at time $t$". **The MSD and the VACF carry exactly the same information** (the MSD is the VACF integrated twice).

## 27. Checks with the simulations, and what D means

| Phase | D from MSD slope | D from Green–Kubo |
|---|---|---|
| gas | 5.72 | 5.53 |
| liquid | 0.0303 | 0.0304 |
| crystal | 0 | 0.0000 (wobbles within ±0.001) |

(The gas values come from separate short runs at slightly different temperatures, 2.016 and 1.987; a 3% gap is within the noise.)

### 27.1 D = (size of the velocity) × (how long it is remembered)

Define the **velocity memory time** as the area under the normalized VACF:

$$\tau_v = \frac{1}{C(0)}\int_0^\infty C(s)\,ds \qquad\Rightarrow\qquad D = \frac{C(0)}{3}\,\tau_v = \frac{k_BT}{m}\,\tau_v$$

- **Gas:** $\tau_v = 2.77$ (the collision time is about 2.0: the velocity survives roughly 1.4 collisions). $D = 1.99\times2.77 = 5.5$.
- **Liquid:** $\tau_v = 0.042$. $D = 0.72\times0.042 = 0.030$. In panel (b) the running integral first climbs to **0.049**, then the negative bounce-back lobe pulls it down to **0.030**. **Bounce-back reduces D by about 40%.** It is the same effect as the dip of $\alpha$ below 1 in the MSD (Section 21).
- **Crystal:** the positive and negative lobes cancel exactly, so $\tau_v = 0$ and $D = 0$. An atom that always returns to its site has zero net memory area. (Equivalently: the displacement of a crystal atom stays bounded, so the MSD is flat, its slope is zero, and by the exact relation above the integral of $C$ must go to zero.)

The crossover time of Section 21 is $t_c = 2D/T = 2\tau_v$ (with $m = k_B = 1$): for the gas, 5.7 against $2\times2.77 = 5.5$ ✓.

**Example: an exponential VACF.** If $C(t) = \langle v^2\rangle\,e^{-t/\tau_v}$ (a reasonable model for a gas), then $\int_0^\infty C\,dt = \langle v^2\rangle\tau_v$ and $D = \frac13\langle v^2\rangle\tau_v = \frac{k_BT}{m}\tau_v$.

## 28. MSD or VACF in practice

| | MSD (Einstein) | VACF (Green–Kubo) |
|---|---|---|
| needs | unwrapped positions, long times | velocities saved every few steps |
| robustness | robust, simple | the noisy tail can make the integral wander |
| physics shown | regimes (ballistic, cage, diffusive) | bounce-back, oscillations, memory time |

**Computing both is a good check.** In LAMMPS, `compute vacf` works like `compute msd`: the reference velocities are stored when the compute is defined (one time origin), and it outputs 4 values ($\langle v_xv_{x,0}\rangle$, $\langle v_yv_{y,0}\rangle$, $\langle v_zv_{z,0}\rangle$ and their total), averaged over the atoms in the group. The documentation's example stores `c_ID[4]` every step with `fix vector` and integrates it with `variable ... trap()`. Remember to divide by 3 to get $D$, and to sample often enough (the liquid VACF changes within 0.05 time units).

---

# Part VI. Shear viscosity

## 29. What viscosity is

Put a fluid between two flat plates. Hold the bottom plate still and slide the top plate sideways. The fluid layers slide over one another, fastest near the top, slowest near the bottom. The force per unit area needed to keep the top plate moving is the **shear stress**, and for most liquids

$$\text{shear stress} = \eta\times(\text{rate at which the fluid velocity changes across the gap})$$

**$\eta$ measures how strongly the layers resist sliding past each other** (honey: large $\eta$; water: small).

At the molecular level, sideways (x) momentum is passed from layer to layer in two ways:

1. molecules that **cross** from one layer to the next carry their x-momentum with them;
2. **forces** between molecules in neighbouring layers pull them along.

## 30. The pressure tensor

$P_{\alpha\beta}$, with $\alpha, \beta \in \lbrace x, y, z\rbrace$, is the flow of $\alpha$-momentum across a plane facing the $\beta$ direction, per unit area per unit time. The lecture writes (calling it the "stress tensor")

$$P_{\alpha\beta} = \frac1V\left(\sum_i\frac{p_{i\alpha}\,p_{i\beta}}{m} + \sum_i r_{i\alpha}\,f_{i\beta}\right)$$

- $V$: volume of the box. $p_{i\alpha} = m\,v_{i\alpha}$: the $\alpha$-component of the momentum of atom $i$. $f_{i\beta}$: the $\beta$-component of the total force on atom $i$. $r_{i\alpha}$: a position component.
- **Kinetic term** $\sum_i p_{i\alpha}p_{i\beta}/m$: momentum **carried** by atoms moving across a plane (way 1 above).
- **Virial term** $\sum_i r_{i\alpha}f_{i\beta}$: momentum **passed across** a plane by forces (way 2).

**Diagonal elements** ($xx$, $yy$, $zz$): their average is the ordinary pressure,
$$P = \tfrac13\left(P_{xx} + P_{yy} + P_{zz}\right)$$
(the virial equation of state). In our liquid all three averaged 0.89: isotropic, as a fluid at rest must be.

**Off-diagonal elements** ($xy$, $xz$, $yz$): **shear stresses**. In a fluid at rest they average to zero (we got 0.006, −0.002, −0.001), but they **fluctuate**, by about ±0.17. Those fluctuations are what the Green–Kubo formula uses.

**Three corrections to the lecture's formula:**

1. **Periodic boundaries.** $\sum_i r_{i\alpha}f_{i\beta}$ is ambiguous in a periodic box (which periodic image of atom $i$ is meant?). Use the **pair form**
$$\sum_{i<j} r_{ij,\alpha}\,f_{ij,\beta}$$
where $\vec r_{ij} = \vec r_i - \vec r_j$ (minimum image) and $\vec f_{ij}$ is the force on $i$ due to $j$. This is what LAMMPS computes.
2. **Masses:** use $m_i$ for each atom if the masses differ.
3. **Sign:** "stress" is conventionally $-P$. The sign cancels in the product $P_{xy}(t)P_{xy}(0)$, so it does not matter for the viscosity.

## 31. The Green–Kubo formula for viscosity

The lecture writes $\eta = \frac{V}{k_BT}\int_0^\infty\langle P_{\alpha\beta}(t)P_{\alpha\beta}(0)\rangle\,dt$ with "$\alpha, \beta \in \lbrace x, y, z\rbrace$". **Correction: $\alpha \ne \beta$.** The diagonal elements have a nonzero average (the pressure), and their fluctuations measure a mixture of shear and bulk viscosity. The correct formula uses an off-diagonal element:

$$\boxed{\;\eta = \frac{V}{k_BT}\int_0^\infty\langle P_{xy}(t)\,P_{xy}(0)\rangle\,dt\;}$$

and in practice the three off-diagonal estimates ($xy$, $xz$, $yz$) are averaged.

**Same structure as D: (size of the fluctuation) × (how long it is remembered):**

$$\eta = \underbrace{\frac{V\langle P_{xy}^2\rangle}{k_BT}}_{\text{size: } 24.4} \times \underbrace{\tau_s}_{\text{memory: } 0.136} = 3.32$$

where $\tau_s$ is the area under the normalized stress autocorrelation. The general principle is **fluctuation–dissipation**: the way spontaneous equilibrium fluctuations die away tells you how the system responds to being pushed. (Deriving the $V/k_BT$ prefactor needs linear-response theory, which is beyond this lecture. Roughly: the factor $V$ is needed because fluctuations of a whole-box average shrink as $1/V$, and $1/k_BT$ turns a fluctuation into a response.)

## 32. Our result, its error bar, and two cross-checks

**The run:** the LJ liquid ($\rho = 0.844$, $T \approx 0.72$, 500 atoms), 200,000 steps (1000 time units), $P_{xy}$, $P_{xz}$, $P_{yz}$ recorded every step, autocorrelation averaged over the three.

**Panel (d)** of the figure: the running integral rises quickly, levels off around $t \approx 1.5$–2.5, and then wanders because of noise. The value is read from the level part.

**Error bar by block averaging (Lecture 9).** Splitting the run into 10 blocks and computing $\eta$ in each:

block values: 3.58, 2.82, 3.86, 3.84, 2.85, 2.92, 3.50, 2.71, 3.50, 3.61

$$\eta = 3.32 \pm 0.26 \quad(90\%,\ t\text{-value } 1.833 \text{ for } 9 \text{ degrees of freedom})$$

For argon, the unit of viscosity is $\sqrt{m\epsilon}/\sigma^2 = 9.0\times10^{-5}$ Pa·s, so $\eta \approx 3.0\times10^{-4}$ Pa·s (0.3 mPa·s), the right size for liquid argon.

**Cross-check 1: Stokes–Einstein.** A sphere of radius $R$ diffusing in a liquid of viscosity $\eta$ has $D = k_BT/(c\,\pi\eta R)$, with $c = 6$ for a "sticky" surface and $c = 4$ for a "slippery" one. With $R \approx 0.5\sigma$ and our $D = 0.0304$, $T = 0.72$:

$$\eta = \frac{k_BT}{c\,\pi DR} = \frac{0.72}{c\,\pi\times0.0304\times0.5} = 2.5\ (c = 6)\ \text{to}\ 3.8\ (c = 4)$$

Our 3.32 sits inside that range. **A thick liquid diffuses slowly.**

**Why $\eta$ is so much noisier than $D$.** $D$ averages the velocities of 500 atoms; $P_{xy}$ is **one number for the whole box**. Our $D$ converged in 60 time units; $\eta$ needed 1000 and is still uncertain by ±8%.

**Cross-check 2: the box-size correction to D.** In a periodic box, a diffusing atom drags fluid that interacts with its own periodic images, so $D$ comes out too low. The Yeh–Hummer correction for a cubic box of side $L$ is

$$D_\infty = D_{\text{box}} + \frac{2.837\,k_BT}{6\pi\eta L}$$

For our box ($L = 8.4$): $2.837\times0.72/(6\pi\times3.32\times8.4) = 0.0039$, so $D_\infty \approx 0.0343$. **The small-box D was about 13% too low**, and correcting it requires $\eta$.

## 33. Viscosity in LAMMPS

The LAMMPS documentation ("Calculate viscosity") gives a Green–Kubo example. Its production part (the variables `s`, `p`, `d`, `convert`, `kB`, `T`, `V`, `dt` are set earlier in that script: `s` is the sample interval, `p` the correlation length, `d = p*s`):

```
reset_timestep 0
variable     pxy equal pxy
variable     pxz equal pxz
variable     pyz equal pyz
fix          SS all ave/correlate $s $p $d &
             v_pxy v_pxz v_pyz type auto file S0St.dat ave running
variable     scale equal ${convert}/(${kB}*$T)*$V*$s*${dt}
variable     v11 equal trap(f_SS[3])*${scale}
variable     v22 equal trap(f_SS[4])*${scale}
variable     v33 equal trap(f_SS[5])*${scale}
thermo_style custom step temp press v_pxy v_pxz v_pyz v_v11 v_v22 v_v33
run          100000
variable     v equal (v_v11+v_v22+v_v33)/3.0
```

- `fix ave/correlate ... type auto` computes $\langle P_{xy}(t)P_{xy}(0)\rangle$ (and $xz$, $yz$) using **many time origins** (`ave running` keeps accumulating).
- `trap()` integrates it (trapezoidal rule); `scale` multiplies by $V/k_BT$, the sampling interval, and a unit conversion.
- The three results are averaged, exactly as in Section 31.

The same documentation page lists non-equilibrium alternatives (shearing the box with `fix deform` and `fix nvt/sllod`, moving walls, the reverse method `fix viscosity`, a periodic-perturbation method) and an **Einstein form** for viscosity, the analogue of the MSD for diffusion.

---

# Part VII. Wrap-up

## 34. Corrections and clarifications to the lecture notes

| The notes say | Correct version | Section |
|---|---|---|
| $r^2 = x^2 + y^2 + z^2$ in an MD simulation | in MD, $\vec r$ is the **displacement** from the time origin, $\vec r_i(t_0 + t) - \vec r_i(t_0)$, from **unwrapped** coordinates | 14, 16 |
| $\frac{d}{dt}\langle r^2\rangle = 6D$ | true only at **long times**; the diffusion equation fails in the ballistic regime. Safe form: $D = \lim_{t\to\infty}\frac16\,d\langle r^2\rangle/dt$ | 13.4 |
| the factor 6 | it is $2d$ for $d$ dimensions: $2D$ along a line, $4D$ on a surface | 12 |
| $\eta = \frac{V}{k_BT}\int\langle P_{\alpha\beta}(t)P_{\alpha\beta}(0)\rangle dt$, $\alpha, \beta \in \lbrace x, y, z\rbrace$ | **$\alpha \ne \beta$**: use $xy$, $xz$, $yz$ and average them | 31 |
| $\sum_i r_{i\alpha}f_{i\beta}$ | under periodic boundaries use the pair form $\sum_{i<j}r_{ij,\alpha}f_{ij,\beta}$; use $m_i$ if masses differ | 30 |
| "stress tensor" $P$ | conventionally this is the **pressure** tensor (stress = −P); the sign cancels in the autocorrelation | 30 |
| "Green–Kubo functions" | usually called Green–Kubo **relations** | 25 |

## 35. Your questions, in one place

| Question | Short answer | Section |
|---|---|---|
| "At a time step, do we compute $\lvert\vec r_i(t_0+t) - \vec r_i(t_0)\rvert^2$ for all molecules and then average?" | Yes, for each lag: square each molecule's displacement, add, divide by $N$. One lag = one point of the curve. The MSD needs two frames, not one time step. | 14.1 |
| "I don't understand the sum over $k$." | $k$ numbers the starting frames (windows). Do the one-origin calculation from several starts and average. | 14.3 |
| "So how does LAMMPS do it?" | `compute msd`: one origin (when the command runs), unwrapped positions, averaged over the group's atoms, 4 outputs. Many origins: several computes, or post-process a dump. | 15 |
| "In which step is it divided by $N$?" | Inside `compute msd`; `c_msd[4]` is already the per-atom average. The `/6` in the fit is the dimension factor, not $N$. | 15.1 |
| "What does ballistic mean, and why MSD ∝ $t^2$?" | Free flight before any collision: $\Delta x = v_xt$, so $\Delta x^2 \propto t^2$; averaged, MSD $= (3k_BT/m)\,t^2$. | 17 |
| "Why is it proportional to $t$ in the diffusive region?" | Random step directions make all cross terms cancel; only the $N \propto t$ square terms remain. | 18 |
| "What is the difference between the gas and liquid MSD plots?" | Gas: long free flights, ballistic until $t \approx 2$, slow bend, large D. Liquid: ballistic only until $t \approx 0.07$, small cage dip, D about 190× smaller. | 21 |
| (Signs, from the first check questions) | Flux has the opposite sign to the slope; a straight profile does not change; a valley rises. | 3.2, 4 |

## 36. Practice questions with answers

### 36.1 The diffusion equation

1. **$c$ = 8, 6, 4 molecules/nm³ at $x$ = 0, 1, 2 nm; $D = 2$ nm²/ns (= $2\times10^{-9}$ m²/s). Flux? Does $c$ at 1 nm change?**
   Slope $= (4 - 8)/2 = -2$ nm⁻⁴; $j = -2\times(-2) = $ **+4** molecules/(nm² ns), towards $+x$. No change: $8 - 12 + 4 = 0$ (as much flows in as out).
2. **$c$ = 3, 1, 3 at 1 nm spacing, same $D$. Initial rate of change of the middle?**
   $(3 - 2 + 3)/1^2 = 4$; $\partial c/\partial t = 2\times4 = $ **+8** molecules/(nm³ ns): it **rises** (a valley fills).
3. **$c$ = 2, 5, 8 at $x$ = 0, 1, 2 nm, same $D$. Flux? Change?**
   Slope $+3$; $j = -6$ molecules/(nm² ns), towards $-x$ (from the 8 side to the 2 side). No change: $2 - 10 + 8 = 0$.

### 36.2 Spreading from a point

4. **Water, $D = 0.23$ Å²/ps. How long until $\sigma$ (per direction) reaches 10 Å?**
   $\sigma^2 = 100 = 2Dt$, so $t = 100/0.46 = $ **217 ps**.
5. **Each direction has mean 0 and variance $2Dt$. What is $\langle r^2\rangle$ in 3D? On a surface?**
   $3\times2Dt = $ **6Dt**; on a surface $2\times2Dt = $ **4Dt**.

### 36.3 Computing the MSD

6. **Displacements (3, 0, 4) Å and (−1, 2, 2) Å. MSD?**
   $(25 + 9)/2 = $ **17 Å²**.
7. **$L = 40$ Å; wrapped $x$ is 39.0 Å at $t_0$ and 1.0 Å at $t_0 + t$; the image flag went from 0 to 1. True and naive $\Delta x$?**
   Unwrapped end $= 1.0 + 40 = 41.0$: true $\Delta x = $ **+2.0 Å** ($\Delta x^2 = 4$). Naive: $1.0 - 39.0 = -38.0$ Å ($\Delta x^2 = 1444$, 361 times too big).
8. **1D walker, ±1 Å, 3 hops. MSD?**
   Paths: +3 (1), +1 (3), −1 (3), −3 (1). MSD $= (9 + 3 + 3 + 9)/8 = $ **3 Å²** $= N\ell^2$.
9. **In the example of Section 14.3, compute the lag of 2 ps window by window.**
   Window 0→2: A $(1-0)^2 = 1$, B $(6-5)^2 = 1$, window MSD 1. Window 1→3: A $(3-1)^2 = 4$, B $(6-4)^2 = 4$, window MSD 4. Average: **2.5**.

### 36.4 Ballistic motion

10. **$v_x = 3$, no collisions: $\Delta x$ and $\Delta x^2$ at $t = 0.01$ and 0.02?**
    0.03 and 0.0009; 0.06 and 0.0036. Doubling $t$ multiplies $\Delta x^2$ by **4**.
11. **$v_x = +1, -3, +2$: average $\Delta x$ and MSD$_x$ at $t = 0.1$?**
    $\Delta x$ = 0.1, −0.3, 0.2, average **0**. Squares 0.01, 0.09, 0.04: MSD$_x = 0.14/3 = $ **0.0467** $= \langle v_x^2\rangle t^2 = (14/3)\times0.01$.
12. **Why does a collision make the MSD fall below the ballistic line?**
    It turns the atom sideways or back, so the later motion no longer adds up in line with the earlier motion; the atom ends up closer to its start than straight flight would take it.

### 36.5 Diffusive motion

13. **Random steps of ±0.5 Å every 0.1 ps (1D). MSD after 1 ps and 4 ps? D?**
    $N = 10$: $10\times0.25 = $ **2.5 Å²**; $N = 40$: **10 Å²** (4× time, 4× MSD). $D = 2.5/(2\times1) = $ **1.25 Å²/ps** $= 1.25\times10^{-8}$ m²/s.
14. **The same steps all in one direction: MSD after 1 ps?**
    $(10\times0.5)^2 = $ **25 Å²** (10 times the random walk; after 4 ps, 400 Å² against 10 Å²).
15. **For $N = 100$ steps: how many square and cross terms? Why does killing the cross terms turn $t^2$ into $t$?**
    100 squares, $100\times99 = 9900$ cross terms. All $N^2$ terms grow like $t^2$; only the $N$ squares (∝ $t$) survive random directions.

### 36.6 Phases

16. **LJ liquid at $T = 0.705$: predicted MSD at $t = 0.01$? A crystal at the same $T$?**
    $3\times0.705\times0.01^2 = $ **$2.1\times10^{-4}$ σ²**, the **same** for the crystal (the ballistic part depends only on $T$ and $m$).
17. **A crystal's plateau is 0.030 σ², nearest-neighbour distance $1.122\sigma$. Lindemann ratio?**
    $\langle u^2\rangle = 0.015$, rms $= 0.1225\sigma$, ratio $= $ **0.109**: inside 0.1–0.15, **close to melting**. (Since the plateau ∝ $T$, this corresponds to $T \approx 0.5\times0.030/0.0212 \approx 0.71$, right at the LJ triple point.)
18. **Surface atoms: MSD$_x$ + MSD$_y$ goes from 2.0 to 6.0 Å² between 10 and 30 ps; MSD$_z$ stays at 0.4 Å². $D_\parallel$?**
    Slope $= 4.0/20 = 0.2$ Å²/ps; $D_\parallel = 0.2/4 = $ **0.05 Å²/ps** $= 5\times10^{-10}$ m²/s. The flat MSD$_z$ means the atoms are bound to the surface.
19. **Log–log slope 1.0 from 10 to 100 ps, then 1.6 beyond 200 ps. Suspect?**
    Collective motion. Check: centre-of-mass drift (`com yes`, total momentum), too few origins at long lags, and non-equilibrium (MSD from early versus late origins; is the system still changing?).

### 36.7 Green–Kubo

20. **$C(t) = \langle v^2\rangle e^{-t/\tau_v}$. Show $D = (k_BT/m)\,\tau_v$.**
    $D = \frac13\int_0^\infty\langle v^2\rangle e^{-t/\tau_v}dt = \frac13\langle v^2\rangle\tau_v = \frac13\cdot\frac{3k_BT}{m}\tau_v = \frac{k_BT}{m}\tau_v$.
21. **Why does a crystal's VACF integrate to zero, although $C(t)$ is far from zero at short times?**
    The positive and negative lobes cancel: the atom keeps coming back. Equivalently, its displacement stays bounded, the MSD is flat, $d\text{MSD}/dt = 2\int_0^t C\,ds \to 0$.
22. **A Green–Kubo viscosity from one long run keeps drifting as the integration limit grows. What do you do?**
    Stop the integral where the plateau starts (before noise dominates), and get an error bar from independent blocks or runs (block averaging), averaging $xy$, $xz$, $yz$ and using many time origins.
23. **$V\langle P_{xy}^2\rangle/k_BT = 24.4$ and the stress memory time is 0.136. $\eta$?**
    $24.4\times0.136 = $ **3.32**.
24. **$D_{\text{box}} = 0.0304$, $\eta = 3.32$, $T = 0.72$, $L = 8.4$ (reduced units). $D_\infty$?**
    $0.0304 + 2.837\times0.72/(6\pi\times3.32\times8.4) = 0.0304 + 0.0039 = $ **0.0343**.

## 37. Python code used in these notes

**MSD with many time origins:** see Section 14.4.

**Autocorrelation with all time origins at once (FFT trick).** Computing $\langle x(t_0 + t)\,x(t_0)\rangle$ for every origin by brute force is slow; a Fourier transform does all origins together. Zero-padding to length $2n$ prevents the end of the series wrapping round onto the start, and dividing by $(n - \text{lag})$ turns the sum into an average over the origins that fit.

```python
import numpy as np

def acf_fft(x, nmax):
    """Autocorrelation of x along axis 0 (x has shape (n_frames, ...)) for lags 0..nmax-1,
    averaged over all time origins."""
    n = x.shape[0]
    F = np.fft.rfft(x, n=2 * n, axis=0)                 # zero-pad to 2n
    c = np.fft.irfft(F * np.conj(F), axis=0)[:nmax]     # sum over origins, for every lag
    norm = (n - np.arange(nmax)).reshape((-1,) + (1,) * (x.ndim - 1))
    return c / norm                                      # number of origins for each lag
```

**Green–Kubo D from saved velocities** (`vel` has shape (n_frames, N, 3), frames `dt_frame` apart):

```python
C = acf_fft(vel, nlag).sum(axis=2).mean(axis=1)      # <v(t).v(0)>, averaged over atoms
t = np.arange(nlag) * dt_frame
D_running = np.concatenate([[0], np.cumsum(0.5 * (C[1:] + C[:-1]) * dt_frame)]) / 3   # trapezoid rule
# read D from the level part of D_running (panel b of the figure)
```

**Green–Kubo viscosity with a block-averaged error bar** (`P_off` has shape (n_steps, 3) holding $P_{xy}, P_{xz}, P_{yz}$ every step `dt`; `V` the volume, `T` the temperature, reduced units with $k_B = 1$):

```python
from scipy import stats

def eta_running(P_off, nlag, dt, V, T):
    C = acf_fft(P_off, nlag).mean(axis=1)                # average over xy, xz, yz
    return V / T * np.concatenate([[0], np.cumsum(0.5 * (C[1:] + C[:-1]) * dt)])

nlag, nb = 1000, 10                                      # integrate up to t = 5; 10 blocks
t = np.arange(nlag) * dt
window = (t >= 1.5) & (t <= 2.5)                         # the level part of the running integral
B = len(P_off) // nb
etas = np.array([eta_running(P_off[k*B:(k+1)*B], nlag, dt, V, T)[window].mean() for k in range(nb)])
mean, se = etas.mean(), etas.std(ddof=1) / np.sqrt(nb)
e = stats.t.ppf(0.95, nb - 1)                            # 90% confidence, Student t
print(f"eta = {mean:.2f} +/- {e*se:.2f} (90%)")         # our run: 3.32 +/- 0.26
```

## 38. Summary: the whole chain

1. **Diffusion equation.** Conservation, $\partial c/\partial t = -\nabla\cdot\vec j$, plus Fick's law, $\vec j = -D\nabla c$ (flow goes downhill because random hops carry more molecules from where there are more), give $\partial c/\partial t = D\nabla^2c$: peaks fall, valleys rise.
2. **Point source.** The solution is a Gaussian with variance $2Dt$ in each direction: the central limit theorem for a random walk, width $\propto\sqrt t$.
3. **MSD.** $\langle r^2\rangle = \int r^2c\,d\vec r = 6Dt$; the lecture's two integrations by parts give $d\langle r^2\rangle/dt = 6D$ ($2dD$ in $d$ dimensions), valid at long times.
4. **In MD.** MSD$(t) = \frac{1}{NK}\sum_k\sum_i\lvert\vec r_i(t_{0,k} + t) - \vec r_i(t_{0,k})\rvert^2$ with unwrapped coordinates; LAMMPS `compute msd` uses one origin and already divides by $N$.
5. **Ballistic.** Before collisions, $\Delta x = v_xt$, so MSD $= (3k_BT/m)\,t^2$ (slope 2 on log–log), the same in every phase.
6. **Diffusive.** After many collisions the cross terms cancel and MSD $= N\ell^2 \propto t$ (slope 1), with slope $6D$; fit the slope, not MSD/6t.
7. **Phases.** Gas: long ballistic part, large D. Liquid: quick switch, small cage dip. Supercooled: a cage plateau, then diffusion. Crystal: overshoot, then a plateau at $2\langle u^2\rangle$ (amplitude, Lindemann ratio, vibration time).
8. **Green–Kubo for D.** $D = \frac13\int_0^\infty\langle\vec v(t)\cdot\vec v(0)\rangle dt = (k_BT/m)\,\tau_v$; the MSD is the VACF integrated twice. Bounce-back lowers D; in a crystal the integral is zero.
9. **Viscosity.** $\eta = \frac{V}{k_BT}\int_0^\infty\langle P_{xy}(t)P_{xy}(0)\rangle dt$ with $\alpha \ne \beta$ and the pair form of the virial; much noisier than D, so block-average it. Our LJ liquid: $\eta = 3.32 \pm 0.26$ (90%), consistent with Stokes–Einstein.
10. **The pattern:** every Green–Kubo transport coefficient = **(size of a fluctuation) × (how long it is remembered)**.

