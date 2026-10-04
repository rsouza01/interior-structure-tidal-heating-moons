# Interior Structure & Tidal Heating of the Moons: A Hands-On Study Guide

Oct 4, 2026 · @Rodrigo de Souza

## How to use this guide

The goal is a Python code base that takes a moon's mass, radius and moment of inertia, builds layered interior models, computes its tidal response and heating, and reproduces published numbers for Ganymede, Europa, Io and Enceladus. You do it in ten phases, easiest first, and each phase ends with a number you must reproduce before moving on.

Every phase has the same six parts, so you always know where you are:

1. **Goal**: what you will be able to do afterwards.
2. **Concepts**: the physics, in the order you need it.
3. **Math**: the equations, with the derivation steps you should do on paper.
4. **Hands-on**: what to code.
5. **Check**: a known result to reproduce. If you can't hit it, don't move on.
6. **Gate**: a one-line self-test of understanding.

| Phase | Topic                                                      | Difficulty     | Rough time (5 to 6 h/week) |
| ----- | ---------------------------------------------------------- | -------------- | -------------------------- |
| 0     | Setup, repo, units                                         | Easy           | 1 week                     |
| 1     | Physical picture and vocabulary                            | Easy           | 1 to 2 weeks               |
| 2     | Observables: gravity, J2, C22, k2, libration               | Easy to medium | 2 weeks                    |
| 3     | Homogeneous and two-layer hydrostatic bodies, Radau-Darwin | Medium         | 2 weeks                    |
| 4     | Multi-layer models with equations of state                 | Medium         | 3 weeks                    |
| 5     | Degeneracy and Bayesian inversion                          | Medium to hard | 3 weeks                    |
| 6     | Tides and Love numbers, analytic case                      | Medium         | 2 weeks                    |
| 7     | Layered viscoelastic response, propagator matrix           | Hard           | 4 to 5 weeks               |
| 8     | Tidal heating, rheology, thermal feedback                  | Hard           | 4 weeks                    |
| 9     | Validation, data, write-up                                 | Medium         | 2 weeks                    |

That is roughly four to five months. It will run longer if you pause to derive things, which is the point.

### Ground rules

- **Derive before you code.** Every equation in Phases 3, 6 and 7 should be derived once by hand. You already have the habits for this from the TOV and Kaluza-Klein work.
- **SI units everywhere** in the code. Convert only when printing (km, g/cm³, GW).
- **One test per phase** in `tests/`, asserting the check value within a stated tolerance.
- **Start with Ganymede and Europa.** Their structure is the best constrained, and Ganymede's data are the cleanest. Io and Enceladus come in Phase 8.
- **Citations in this guide are from my memory of the literature**, not freshly opened sources. Treat every number and reference as "approximate, verify before relying on it", and always check the primary paper when a test depends on it.

### Suggested repo layout

```text
moon-interiors/
  README.md
  pyproject.toml
  data/            # constants, published gravity solutions (as YAML)
  moons/           # one module per body: mass, radius, C/MR2, k2, refs
  interior/        # structure: eos.py, hydrostatic.py, radau.py
  tides/           # love.py, propagator.py, rheology.py, heating.py
  inversion/       # priors.py, likelihood.py, sampler.py
  notebooks/       # one per phase
  tests/           # one per phase, reproducing the check values
  notes/           # your derivations (LaTeX or markdown)
```

## Phase 0: Setup, tools and constants

**Goal:** a clean Python environment and a constants module you will trust for the next several months.

### Python stack

- `numpy`, `scipy` (ODE solvers, root finding, interpolation, linear algebra)
- `matplotlib` for plots
- `emcee` for MCMC (Phase 5), and optionally `corner` for posterior plots
- `pyyaml` for the data files
- `pytest` for the checks
- Worth knowing about, **not** to be used as a black box until you have built your own: `burnman` (mineral equations of state), `SeaFreeze` (water and ice thermodynamics), `PlanetProfile` (a published ocean-world interior code), `TidalPy` (tidal and thermal-orbital code). Their value to you is as a cross-check in Phases 4, 7 and 9. Verify they are still maintained before you depend on them.

### Constants and body data

Put these in `data/constants.yaml` and `moons/*.py`, with a source comment on every line. The values below are approximate, from memory; replace them with the figures in the papers you cite.

| Body      | Mass (kg) | Mean radius (km) | Mean density (kg/m³) | C/MR²            | Orbital period (days) |
| --------- | --------- | ---------------- | -------------------- | ---------------- | --------------------- |
| Ganymede  | 1.482e23  | 2634.1           | 1942                 | 0.3115           | 7.155                 |
| Europa    | 4.80e22   | 1560.8           | 3013                 | 0.346            | 3.551                 |
| Io        | 8.93e22   | 1821.5           | 3528                 | 0.377            | 1.769                 |
| Enceladus | 1.08e20   | 252.1            | 1609                 | not firmly known | 1.370                 |
| Titan     | 1.345e23  | 2574.7           | 1881                 | about 0.34       | 15.945                |

Other constants: G = 6.674e-11 m³ kg⁻¹ s⁻²; Jupiter's GM is about 1.267e17 m³ s⁻², Saturn's about 3.793e16 m³ s⁻². Work with GM products where you can, since they are measured far more precisely than G or M separately. That is why published moon masses carry a GM.

### Hands-on

1. Create the repo with the layout from the overview and an editable install (`pip install -e .`).
2. Write `moons/ganymede.py` as a small dataclass: `GM`, `R`, `n` (mean motion, rad/s), `e`, `a`, `C_MR2`, `J2`, `C22`, and a `refs` string.
3. Write a helper `mean_motion(period_days)` returning `2*pi/(period*86400)`.
4. Write a plotting helper that makes a standard "profile" figure: density, gravity, pressure and mass versus radius, stacked, sharing the radius axis. You will reuse it every phase.

### Check

Compute the mean density of each body from M and R and match the table to three digits. For Ganymede you should get 1942 kg/m³.

### Gate

Can you explain why a moon's mean density alone does not tell you whether it has a core? (If not, you are ready for Phase 1.)

## Phase 1: The physical picture and the vocabulary

**Goal:** be able to sketch the interior of each moon from memory and say what keeps it warm, before touching any equation.

### Concepts, in the order you need them

1. **Differentiation.** A body that got warm enough early on lets dense metal sink and light material rise. A fully differentiated moon has a metal core, a rock mantle and an outer ice or water layer. An undifferentiated one is a uniform rock and ice mixture. Callisto is the classic partial case, which is why it matters as a contrast.
2. **The layers and what they are made of.**
   - **Iron or iron-sulfide core** (Fe, FeS), density roughly 5000 to 8000 kg/m³ depending on sulfur content.
   - **Silicate mantle**, density roughly 3200 to 3600 kg/m³ (olivine and pyroxene, possibly hydrated in the ice-rich moons).
   - **High-pressure ice layers** (ice III, V, VI) in the largest moons, where the pressure at the base of the water layer is several hundred MPa to over 1 GPa. These are denser than liquid water, so on Ganymede the ocean can be sandwiched between two ice layers.
   - **Liquid water ocean**, with salts and possibly ammonia. Density about 1000 to 1200 kg/m³.
   - **Ice Ih shell**, density about 920 kg/m³, rigid and cold at the top.
3. **Heat sources.** Three matter: leftover heat from accretion and differentiation, radioactive decay in rock (long-lived U, Th and K isotopes, of order a few pW per kg of rock today), and **tidal dissipation**. For Io the third dominates by a huge margin.
4. **Tidal locking.** All these moons keep one face toward their planet, so their rotation rate equals their mean motion n. That fixes the shape of the permanent tidal and rotational bulge, which you will exploit in Phase 3.
5. **The Laplace resonance.** Io, Europa and Ganymede have orbital periods in a 1:2:4 ratio. The resonance keeps their orbits slightly eccentric, which means the tide they feel changes through each orbit, which flexes them, which heats them.
6. **Hydrostatic equilibrium.** For a body that behaves like a fluid on long timescales, pressure balances gravity at every depth. This single idea builds the whole structure model in Phases 3 and 4.

### Pocket glossary

| Term         | Meaning                                                                                                           |
| ------------ | ----------------------------------------------------------------------------------------------------------------- |
| C/MR²        | Normalized polar moment of inertia. 0.4 for a uniform sphere, lower when mass is concentrated toward the center   |
| J2, C22      | Degree-2 gravity coefficients: how flattened the gravity field is, and how much it is stretched toward the planet |
| k2           | Tidal Love number: how much the gravity field changes when the body is tidally deformed                           |
| h2           | Love number for the radial displacement of the surface                                                            |
| Maxwell time | Viscosity divided by shear modulus. Below it the material acts elastic, above it viscous                          |
| Q            | Quality factor, inversely related to how much energy is lost per tidal cycle                                      |

### The bodies at a glance

| Body      | Likely structure (top to bottom)                                   | Evidence quality                               |
| --------- | ------------------------------------------------------------------ | ---------------------------------------------- |
| Ganymede  | Ice Ih, ocean, high-pressure ice, rock mantle, metallic core       | Strong (gravity plus intrinsic magnetic field) |
| Europa    | Ice Ih, ocean, rock mantle, metallic core                          | Good for ocean, weak for core                  |
| Io        | Rock mantle with partial melt, iron-rich core, no ice              | Good for core, open debate on mantle melt      |
| Enceladus | Ice shell, regional or global ocean, low-density porous rocky core | Good for ocean, open for core                  |
| Titan     | Ice, ocean or slush, high-pressure ice, hydrated rock interior     | Moderate                                       |

Treat this table as a hypothesis to test in later phases, not as settled fact. The Io mantle question in particular has been revised recently, so check the latest papers before you quote it.

### Hands-on

No heavy code yet. In a notebook, compute for each moon in the Phase 0 table:

- surface gravity g = GM/R². You should get about 1.43 (Ganymede), 1.32 (Europa), 1.80 (Io) and 0.113 (Enceladus) m/s².
- the central pressure of a **uniform-density** body, P_c = (2π/3) G ρ² R². For Ganymede this gives about 3.7 GPa. The real central pressure is roughly twice that, because mass is concentrated toward the center. Keep this in mind, it is your first clue that structure matters.
- the ratio of the tidal-heating timescale you would need (say 1e14 W for Io) to the body's total thermal energy content, to see why Io cannot be passively cooling.

### Reading

Skim a review chapter on satellite interiors (the Treatise on Geophysics volume on planets and moons has one, and Nimmo and Pappalardo 2016 is a good ocean-worlds review). You only need the pictures and the vocabulary at this stage.

### Check

The three surface gravities and the Ganymede central pressure above.

### Gate

Without notes, can you explain in three sentences why a tidally locked moon on an eccentric orbit heats up, and why a perfectly circular orbit would not?

## Phase 2: Observables, what a spacecraft actually measures

**Goal:** know every number that constrains an interior, where it comes from, and what it is sensitive to. Then turn a measured C22 into a moment of inertia, which is your first real result.

### The measurements

| Observable                     | How it is measured                                                    | What it constrains                                                                                                      |
| ------------------------------ | --------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------- |
| Mass (GM) and radius R         | Doppler tracking of flybys, imaging and limb fits                     | Mean density, the starting point                                                                                        |
| J2 and C22                     | Doppler tracking over several flybys at different geometries          | Degree-2 shape of the gravity field, which gives C/MR² if the body is in hydrostatic equilibrium                        |
| Tidal Love number k2           | Gravity change around the orbit (needs an orbiter or repeated flybys) | Whether a liquid layer decouples the shell, and the interior's softness                                                 |
| Tidal Love number h2           | Laser altimetry of the surface tide                                   | Same as k2, from displacement instead of gravity. k2 and h2 together break ambiguities                                  |
| Physical libration             | Tracking surface features in images across years                      | Whether the shell is decoupled from the interior (ocean)                                                                |
| Magnetic field                 | Magnetometer flybys                                                   | An intrinsic field implies a dynamo in a liquid core. An induced field implies a conducting layer such as a salty ocean |
| Heat flow and thermal emission | Infrared mapping                                                      | Total tidal and radiogenic output, particularly for Io and Enceladus                                                    |

### Gravity field essentials

The external potential of a body is expanded in spherical harmonics. At degree 2, the part you need is

$$
U(r,\theta,\phi) = -\frac{GM}{r}\left[1 - \left(\frac{R}{r}\right)^2 \left( J_2 P_{20}(\cos\theta) + \left(C_{22}\cos 2\phi + S_{22}\sin 2\phi\right) P_{22}(\cos\theta)\right) + \dots\right]
$$

J2 measures the polar flattening, and C22 measures the equatorial elongation toward the planet. Normalization conventions differ between papers (unnormalized versus fully normalized coefficients). **Always check which one a paper uses**. A factor of √5 or √(5/12) is the usual trap.

### From C22 to moment of inertia

For a body in hydrostatic equilibrium that is tidally locked, the quadrupole field comes only from its rotation and the planet's tide. Define the rotation parameter

$$
q = \frac{n^2 R^3}{GM}
$$

Then the gravity coefficients are tied to the **fluid Love number** k_f of the body:

$$
J_2 = \frac{5}{6}\,k_f\, q, \qquad C_{22} = \frac{1}{4}\,k_f\, q, \qquad \frac{J_2}{C_{22}} = \frac{10}{3}
$$

and k_f is related to the moment of inertia by the Radau-Darwin relation (you will derive it in Phase 3, use it on trust for now):

$$
\frac{C}{MR^2} = \frac{2}{3}\left[1 - \frac{2}{5}\sqrt{\frac{4 - k_f}{1 + k_f}}\right]
$$

Sanity check on the formula: a uniform sphere has k_f = 3/2, and the formula gives exactly 0.4.

Two caveats you must carry through every later phase:

- This only works if the body is **in hydrostatic equilibrium**. A rigid, cold shell can hold some non-hydrostatic anomalies, and for some bodies that biases C/MR². The J2/C22 = 10/3 ratio is your test: if the measured ratio is far off, the assumption is in trouble.
- Even when valid, C/MR² has an uncertainty from the measured C22, which you should propagate (see Phase 5).

### Hands-on

1. Write `radau.py` with two functions: `kf_from_C22(C22, q)` and `moi_from_kf(kf)`.
2. Feed in the Ganymede and Europa numbers. Use the published gravity solutions from the Galileo papers (approximately: Ganymede J2 = 127.5e-6, C22 = 38.3e-6; Europa J2 = 435.5e-6, C22 = 131.0e-6; verify these in the primary papers before you rely on them).
3. Also check the ratio J2/C22 for both.
4. Make a plot of C/MR² against k_f from 0 to 1.5 and mark the five moons on it.

### Check

You should get C/MR² close to 0.311 for Ganymede and close to 0.346 for Europa. The ratio J2/C22 should be near 3.33 for both. If your Ganymede value is far off, the usual culprit is the wrong mean motion or radius, or a normalization mismatch in C22.

### Gate

Why does a _smaller_ C/MR² mean a more centrally condensed body? Explain it using the definition of the moment of inertia, then say what C/MR² of 0.31 versus 0.35 implies about whether Ganymede or Europa has the larger iron core fraction. (Be careful: the answer depends on the ice and ocean thickness too, which is exactly the degeneracy of Phase 5.)

## Phase 3: Your first structure models, homogeneous and two-layer

**Goal:** write a hydrostatic structure solver, test it against analytic cases, and meet degeneracy for the first time. This is the Newtonian cousin of your TOV solver, and it should feel familiar.

### The equations

For a spherically symmetric, non-rotating, static body, with density ρ(r), enclosed mass m(r), gravity g(r) and pressure P(r):

$$
\frac{dm}{dr} = 4\pi r^2 \rho, \qquad g(r) = \frac{G\,m(r)}{r^2}, \qquad \frac{dP}{dr} = -\rho\, g
$$

The quantities you extract are the total mass, and the moment of inertia about the polar axis:

$$
C = \frac{8\pi}{3}\int_0^R \rho(r)\, r^4\, dr, \qquad \frac{C}{MR^2}\ \text{(dimensionless)}
$$

Compare with TOV: you drop the general-relativistic corrections, and the density comes from a prescribed or tabulated law instead of an EOS closure you integrate. In this phase the density is just a **step function** of radius, so P does not feed back on ρ.

### Step 1: the uniform body (analytic test)

For constant ρ, derive on paper:

$$
m(r) = \tfrac{4\pi}{3}\rho r^3, \qquad P(r) = \tfrac{2\pi}{3}\, G \rho^2 \left(R^2 - r^2\right), \qquad \frac{C}{MR^2} = \frac{2}{5}
$$

Code it with a numerical integrator (`scipy.integrate.solve_ivp`, integrating P inward from the surface where P = 0) and compare. Your numerical C/MR² must equal 0.4 to better than 1e-6, and P(0) must match the formula.

### Step 2: the two-layer body

Take a dense inner region of radius r_c and density ρ_c, and an outer region of density ρ_o. Mass conservation fixes one of the unknowns. Derive:

$$
\frac{C}{MR^2} = \frac{2}{5}\,\frac{\rho_c x^5 + \rho_o\,(1 - x^5)}{\bar\rho}, \qquad x = \frac{r_c}{R}, \qquad \bar\rho = \rho_o + (\rho_c - \rho_o)\,x^3
$$

Then explore, for Ganymede (mean density 1942 kg/m³):

1. Rock plus ice only: ρ_c = 3300, ρ_o = 1000. Solve for x from the mean density. You should find x of about 0.74 and C/MR² of about 0.31.
2. Compare that with the measured 0.3115. **A two-layer rock and ice body already matches the moment of inertia.** This is the key lesson of the phase: C/MR² alone cannot tell you whether a metal core exists. You need extra information such as the magnetic field, and a physically reasoned composition.
3. Now fix a metal core at ρ_c = 5150 kg/m³ (FeS) and 8000 kg/m³ (Fe) and make a plot of C/MR² against core radius fraction for several mantle densities. You will see families of models giving the same C/MR².

### Step 3: where the Radau-Darwin relation comes from

In Phase 2 you used the Radau-Darwin relation on trust. Now derive it, in four steps. Take your time, this is the most conceptually dense derivation in the first half.

1. **Clairaut's equation.** For a rotating fluid body in hydrostatic equilibrium, the surface flattening ε(r) of each equipotential obeys a second-order ODE involving the density profile and the mean density ρ̄(r) inside radius r. Derive it, to first order in the small rotation parameter, from the condition that the total potential is constant on each level surface.
2. **Radau's transformation.** Define η(r) = r ε'(r)/ε(r). Show that Clairaut's equation becomes a first-order nonlinear ODE for η, which is much better behaved numerically.
3. **Boundary condition and the Love number.** At the surface, η_s is tied to the fluid Love number k_f by an algebraic relation. Derive it. It must give k_f = 3/2 for the uniform body (where ε is constant, so η = 0).
4. **The Darwin-Radau approximation.** The moment of inertia is expressed through an integral weighted by ε. Using a smooth approximation for how η varies, it reduces to the closed form you used, C/MR² = (2/3)\[1 − (2/5)√((4 − k*f)/(1 + k_f))\]. Follow Hubbard's \_Planetary Interiors* or Zharkov and Trubitsyn for the detailed steps. The standard literature states the relation as accurate to a small fraction of a percent for bodies like these, but you will test that yourself below.

### Step 4: test the approximation numerically

Implement Clairaut's equation (or the Radau form) and integrate it through your two-layer Ganymede model with a shooting method: integrate from the center with ε(0) = 1 (arbitrary normalization), read off η_s, convert to k_f, then compute C/MR² from the Radau-Darwin formula. Compare with the exact C/MR² from your direct integral in Step 2. The difference tells you how much the Radau-Darwin shortcut costs. Keep that number, it is a systematic error you will add to your Phase 5 uncertainties.

### Check

- Uniform sphere: C/MR² = 0.4 and P(0) = (2π/3) G ρ² R² to numerical precision.
- Ganymede rock plus ice: x about 0.74, C/MR² about 0.31.
- Your Radau-Darwin versus direct integral difference, reported as a percentage for at least three different two-layer models.

### Gate

In one paragraph: why is it a feature, not a bug, that two very different interiors can have the same C/MR²? (Hint: the integral C only cares about the _weighted_ density distribution, and many distributions share the same weighted sum.)

## Phase 4: Multi-layer models with equations of state

**Goal:** a general layered solver where density depends on pressure (and later temperature), used to build models of Io, Europa and Ganymede. After this phase you can answer "what would the interior look like if it had X?" for any X.

### What changes from Phase 3

In Phase 3 each layer had a constant density. Now each layer carries an **equation of state** ρ(P, T) and the structure equations are coupled: pressure sets density, density sets gravity, gravity sets pressure. This is exactly the structure of the TOV problem, so keep your TOV code open as a reference.

### Building blocks

**1. Equations of state.** The workhorse for solids is the third-order Birch-Murnaghan form, written in terms of the compression ratio ρ/ρ₀:

$$
P = \frac{3}{2} K_0 \left[\left(\frac{\rho}{\rho_0}\right)^{7/3} - \left(\frac{\rho}{\rho_0}\right)^{5/3}\right]\left[1 + \frac{3}{4}\left(K_0' - 4\right)\left(\left(\frac{\rho}{\rho_0}\right)^{2/3} - 1\right)\right]
$$

It gives P(ρ) explicitly, so you invert it for ρ(P) with a root finder, or pre-tabulate ρ on a pressure grid and interpolate. Illustrative starting parameters (approximate, take proper values from the literature or from the `burnman` database):

| Material                       | ρ₀ (kg/m³)         | K₀ (GPa)        | K₀'     |
| ------------------------------ | ------------------ | --------------- | ------- |
| Silicate mantle (olivine-like) | about 3300         | about 130       | about 4 |
| Iron                           | about 7800         | about 160       | about 5 |
| Troilite (FeS)                 | about 4800 to 5000 | about 60 to 100 | about 4 |
| Ice Ih                         | about 920          | about 9         | about 5 |

For liquid water and the high-pressure ices, begin with a **simple linear compressibility** law and graduate to tabulated thermodynamics (the `SeaFreeze` package, or the IAPWS tables) once your solver is validated. Remember that the real water phase diagram switches ice phases at specific pressures and temperatures. That is a modeling choice you make explicitly, not something that falls out of the equations.

**2. Temperature.** Start isothermal (T doesn't enter ρ). Then add a simple profile: conductive in the cold ice shell, nearly constant in the ocean, adiabatic in the convecting mantle. Thermal expansion changes densities by roughly a percent over these ranges. That is small next to the compositional uncertainty for C/MR², but it matters a lot for the _viscosity_ you need in Phase 7, which is exponentially sensitive to temperature.

**3. Layer bookkeeping.** A model is a list of layers from the surface down, each with a material and a thickness.

### The algorithm: integrate inward and close with a root find

Integrating from the center needs an unknown central pressure. It is neater to integrate **inward from the surface**, where P = 0 and the enclosed mass is the known total M:

$$
\frac{dm}{dr} = 4\pi r^2\rho(P,T), \qquad \frac{dP}{dr} = -\rho\,\frac{G\,m}{r^2}, \qquad m(R) = M,\ \ P(R) = 0
$$

Integrate toward r = 0 (so `dr` is negative). The model is self-consistent only if m(0) = 0. That residual is one equation, so you use it to solve for **one** unknown, usually the core radius or the core density, with `scipy.optimize.brentq`. Everything else (ice thickness, ocean thickness, mantle composition) is an input you chose. C/MR² is then a _prediction_ to compare against the measurement.

```python
@dataclass
class Layer:
    name: str
    eos: Callable[[float, float], float]   # rho(P, T)
    thickness: float | None                # None = solved for (the core)

def integrate(layers, M, R, T_profile):
    # state y = [m, P]; walk inward layer by layer
    # at each boundary, switch EOS; pressure and mass are continuous
    # return r, m, P, rho, g arrays
    ...

def residual(core_radius, layers, M, R, T_profile):
    r, m, P, rho, g = integrate(...)
    return m[-1]            # must be 0

core_radius = brentq(residual, 1e3, 0.95*R, args=...)
C_over_MR2 = (8*pi/3) * trapz(rho*r**4, r) / (M*R**2)
```

Use the integration with `dense_output=True` or a fixed fine grid, and make sure each layer boundary is a _grid node_, otherwise you smear the discontinuity and lose accuracy.

### Hands-on ladder (do these in order)

1. **Reproduce Phase 3.** With every EOS set to incompressible (K₀ huge), your new solver must give the Phase 3 numbers exactly. This is the unit test that guards everything after it.
2. **Io first: two compressible layers, no ice.** An iron-rich core under a silicate mantle. With incompressible densities of 5150 and 3200 kg/m³ you should find a core radius fraction near 0.55 and C/MR² near 0.374, against a measured 0.377. Then turn compressibility on and see how those numbers move. Report how the core radius changes.
3. **Europa: add the water layer.** Core, mantle, ocean, ice Ih shell. Fix a total water-layer thickness anywhere in the roughly 80 to 170 km range found in the literature, and solve for the core. Plot C/MR² against water-layer thickness.
4. **Ganymede: add high-pressure ice.** Core, mantle, a high-pressure ice layer below the ocean, ocean, ice Ih shell. Reproduce the standard textbook sandwich. Check the central pressure against published values (roughly several GPa, order 7 GPa in the literature I know, verify).
5. **Convergence test.** Halve the step size and confirm C/MR² changes by less than your target precision, say 1e-5. Make a convergence plot.

### Check

- Step 1 (incompressible limit) matches Phase 3 to 1e-8.
- Io two-layer: core fraction and C/MR² as above, with the compressibility shift quantified.
- Europa and Ganymede produce a smooth, monotonic gravity profile with a clear jump in density at each layer boundary and g peaking somewhere inside the body. Plot it with the helper from Phase 0.

### Gate

Predict, before running it: if you make the mantle denser while keeping the core radius fixed, what has to happen to the solved core density for the mass to stay the same, and what does that do to C/MR²? Then run it and see whether you were right.

## Phase 5: Degeneracy and Bayesian inversion

**Goal:** turn "the interior is probably like this" into a statement with error bars. You will learn exactly which structural parameters the data constrain and which they don't. This is the phase that separates a toy from a research-grade result.

### Why a single best-fit model is the wrong answer

You have, at best, three numbers (M, R and C/MR²), and the model has many free parameters: core radius, core composition, mantle density, ocean thickness, ice thicknesses. A one-line count shows that infinitely many models fit exactly. Phase 3 already showed it: rock plus ice and rock plus metal core both gave 0.31. So the honest output is a **probability distribution** over structures, not one picture.

### The framework

Let θ be the parameter vector and d the data. Bayes' theorem gives the posterior:

$$
p(\theta\,|\,d) \;\propto\; \mathcal{L}(d\,|\,\theta)\; p(\theta), \qquad \ln \mathcal{L} = -\frac{1}{2}\sum_i \frac{\left(d_i - f_i(\theta)\right)^2}{\sigma_i^2} + \text{const}
$$

Here f_i(θ) is your forward model from Phase 4 (and later Phase 7), and σ_i is the measurement uncertainty plus any systematic error you add in quadrature.

### Define the problem concretely

1. **Parameters (illustrative for Ganymede):** ice Ih shell thickness, ocean thickness, high-pressure ice thickness, mantle density at the reference pressure (a proxy for Fe/Si ratio and hydration), core composition (sulfur fraction, which sets core density). The **core radius is not free**, because the mass closure in Phase 4 solves for it. This is an important point: closure removes one dimension.
2. **Priors:** begin with uniform priors over physically motivated ranges, which you must write down and justify in a table in your notes. The ice Ih shell might be 0 to 150 km, the ocean 0 to 300 km, and so on. Use ranges from the literature rather than convenience.
3. **Data:** C/MR² with its published uncertainty. For Ganymede the Galileo solution is roughly 0.3115 ± 0.0028 (verify in the source). Later add k2 or h2 when you have them (Phase 7).
4. **Systematic error:** add in quadrature the Radau-Darwin accuracy you measured in Phase 3, and think about the non-hydrostatic issue. Studies of icy satellites show a rigid, non-hydrostatic shell can bias the inferred moment of inertia, so a conservative error bar is more honest than a tight one.

### Hands-on, in increasing sophistication

1. **Brute-force 2D map.** Pick two parameters (say core sulfur fraction and mantle density), grid them, compute C/MR² at every point, and contour the misfit. You will see a long, curved valley of acceptable models. That valley is the degeneracy, drawn.
2. **Rejection sampling.** Draw ten thousand parameter sets from the priors, keep those within 1σ or 2σ of the data. Histogram the surviving core radii. This is the cheapest posterior and is surprisingly informative.
3. **MCMC with `emcee`.** Set up `log_prob(theta)` which returns −∞ outside the priors, otherwise the log-likelihood. Run an ensemble of about 32 walkers for several thousand steps. Check convergence with the integrated autocorrelation time (`sampler.get_autocorr_time()`): the chain should be at least a few tens of times longer than it. Discard burn-in, thin, and make a corner plot.
4. **Summarize** with medians and 68 or 95 percent credible intervals for derived quantities you care about: core radius, core mass fraction, total water-layer thickness, central pressure.
5. **Synthetic recovery test.** Pick a known model, compute its C/MR², add noise at the published level, and run the inversion. Does the truth fall inside your credible intervals about 68 percent of the time over many trials? If not, your sampler or likelihood has a bug. Do this before you believe any real-data result.
6. **Prior sensitivity.** Rerun with a different but equally defensible prior (say uniform in log of ocean thickness). If the posterior for a parameter moves a lot, the data are not constraining it, and you must say so.

### What you will probably find

The moment of inertia pins down a _combination_ of parameters, the radial mass distribution, but leaves the partition between core, mantle and ice only loosely determined. Expect a broad core-radius posterior for Europa in particular. That is a real scientific statement, and it is exactly why future missions (Europa Clipper and JUICE) want k2, h2 and magnetometer data. They break degeneracies that gravity cannot.

### Check

- The brute-force valley, rejection samples and MCMC posterior agree for the same priors.
- The synthetic recovery test passes, with its coverage statistics recorded.
- A prior-sensitivity table showing which parameters are data-driven and which are prior-driven.

### Gate

What would you need to measure to shrink the core-radius posterior for Europa by a factor of two? Name the observable and the physical reason it helps, then check whether your forward model can already compute that observable (it can't yet, which motivates Phases 6 and 7).

## Phase 6: Tides and Love numbers, the analytic warm-up

**Goal:** understand how a body deforms under a tide, what k2 and h2 mean, and how dissipation comes out of a complex rigidity. You will do it all analytically for a uniform body first, so that Phase 7 (the layered version) has an exact limit to test against.

### The tidal potential

A planet of mass M_p at distance a raises, at a point on the moon at radius r and angle ψ from the sub-planet point, the degree-2 tidal potential

$$
W_2(r,\psi) = -\frac{G M_p\, r^2}{a^3}\, P_2(\cos\psi), \qquad P_2(x) = \tfrac{1}{2}\left(3x^2 - 1\right)
$$

The moon responds in three ways, each described by a dimensionless **Love number**, evaluated at the surface r = R:

- **k2**: the extra gravitational potential produced by the deformation is k2 times the applied one.
- **h2**: the radial surface displacement is h2 W₂/g.
- **l2** (the Shida number): the tangential displacement, which you will need for the strain energy in Phase 7.

Two limits to know cold. A perfectly **rigid** body has k2 = h2 = 0. A uniform **fluid** body has k2 = 3/2 and h2 = 5/2, which is the same k_f that appeared in the J2 and C22 relations of Phase 2. That is not a coincidence: the permanent tide of a synchronous moon is just the fluid limit of the same response.

### The uniform elastic body

For a homogeneous incompressible elastic sphere of rigidity μ, density ρ and surface gravity g, derive (or follow Love's classic treatment, or Munk and MacDonald) the closed forms:

$$
k_2 = \frac{3/2}{1 + A\mu}, \qquad h_2 = \frac{5/2}{1 + A\mu}, \qquad A \equiv \frac{19}{2\rho g R}
$$

The combination Aμ compares the rigidity with the body's self-gravitational stress scale ρgR. A big Aμ means "stiff compared with gravity", so k2 is small.

### Adding viscoelasticity: the correspondence principle

For a linear viscoelastic body forced at angular frequency ω, you can reuse the elastic solution with the rigidity replaced by a complex, frequency-dependent one. The simplest rheology is **Maxwell** (a spring and dashpot in series), with viscosity η and Maxwell time τ = η/μ:

$$
\mu^*(\omega) = \frac{i\omega\eta\,\mu}{\mu + i\omega\eta} = \mu\,\frac{i x}{1 + i x}, \qquad x = \omega\tau = \frac{\omega\eta}{\mu}
$$

Put it into the elastic formula and the Love number becomes complex. Do the algebra, with B = 1 + Aμ:

$$
k_2^*(x) = \frac{3}{2}\,\frac{1 + i x}{1 + i x B}, \qquad \mathrm{Im}\,k_2 = -\frac{3}{2}\,\frac{x\,(B - 1)}{1 + x^2 B^2}
$$

The real part is the usual tidal response. The imaginary part is the **phase lag** of the deformation behind the forcing, and it is what dissipates energy. Two results to derive yourself:

- Dissipation is maximal when x = 1/B, and then |Im k2| = (3/4)(B − 1)/B = (3/4) Aμ/(1 + Aμ).
- For very high viscosity (x large) the body is elastic and dissipates little. For very low viscosity (x small) it behaves like a fluid and also dissipates little. **In between there is a peak.** This single curve is the physical heart of Io's thermal story, and you will use it in Phase 8.

### Tidal heating formula

For a synchronously rotating moon on an orbit of eccentricity e, the total dissipation (from the classical Peale and Cassen result, to lowest order in e) is

$$
\dot{E} = -\frac{21}{2}\,\mathrm{Im}(k_2)\,\frac{G M_p^2\, n\, R^5\, e^2}{a^6}
$$

with n the orbital mean motion. Because Im(k2) is negative in this convention the heating is positive. Try to derive it: the e-dependent part of the tide at frequency n has a time-varying amplitude proportional to e, and the dissipated power is the product of the strain-rate and the out-of-phase stress, integrated over the body. Peale and Cassen (1978) is the original, and Segatz et al. (1988) is the layered version for Io.

### Hands-on

1. Write `tides/love.py` with `k2_h2_uniform(mu, rho, g, R)` and the complex Maxwell version.
2. **Limits test.** μ → 0 must give 3/2 and 5/2, and μ → ∞ must give 0.
3. **Moon check.** A uniform Moon with a rock-like rigidity of about 6.5e10 Pa should give a k2 of order 0.02, close to the roughly 0.024 measured by lunar laser ranging. A uniform Moon is a crude model, so the agreement is a sanity check, not a result.
4. **Dissipation curve.** For a uniform Io-like body, plot |Im k2| against log viscosity at Io's orbital frequency (ω = n). Mark the peak and compare its position and height with your analytic prediction. The peak height for rock-like rigidity comes out near 0.7, which is far above what Io actually shows.
5. **Io power budget.** Compute the prefactor G M_p² n R⁵ e²/a⁶ for Io (M_p = 1.898e27 kg, a = 4.217e8 m, e = 0.0041). It should come out near 6e14 W. Multiplied by 21/2, that gives roughly 6e15 W × |Im k2|. The observed dissipation of about 1e14 W then needs |Im k2| near 0.016, which agrees with the astrometric estimate k2/Q of about 0.015 (Lainey et al. 2009, verify). Notice what this means: Io sits **far from the dissipation peak** that a uniform model would predict, which tells you the real interior is very different from uniform.
6. Repeat the exercise for Europa (e about 0.009, a = 6.71e8 m) with a uniform model, and compute the prefactor, which comes out near 4e14 W before the |Im k2| factor. Europa's total heating is usually estimated at roughly 1e11 to 1e12 W (verify), so the needed |Im k2| is of order a few 1e-3. A uniform model cannot say where the dissipation happens, for instance in a soft ice shell above an ocean, which is the motivation for the layered model.

### Check

- Fluid and rigid limits, to numerical precision.
- The moon rigidity check, within a factor of order unity of the measured k2.
- Your analytic peak location and height match the numerical curve to better than 1 percent.
- The Io prefactor of about 6e14 W and the implied |Im k2| near 0.016.

### Gate

Explain, without equations, why a body that is too cold (high viscosity) and a body that is too hot (low viscosity) both dissipate little tidal energy, and why that implies a thermostat. Hold that thought for Phase 8.

## Phase 7: Layered viscoelastic response, the propagator-matrix method

**Goal:** compute complex k2, h2 and l2 for a stratified body that mixes solid and liquid layers, and use them to ask whether an ocean exists. This is the hardest phase, so it is split into seven small steps, each with a test, and it relies on your Phase 4 structure as input.

### What problem you are solving

The moon is subjected to a small, slowly varying degree-2 tidal potential. You want its linear response: displacements, stresses, and the extra gravitational potential, as a function of radius. Three physical ingredients, all linearized about the hydrostatic state of Phase 4:

1. **Momentum balance** (quasi-static, because tidal periods of days are very long compared with seismic ones), including the perturbation of gravity and the advection of the background pressure by the displacement.
2. **Poisson's equation** for the perturbed potential, with the density perturbation caused by the displacement, δρ = −∇·(ρ₀u).
3. **A constitutive law** relating stress to strain, taken as complex via the correspondence principle: replace μ by μ\*(ω) from Phase 6.

Assuming a degree-2 spheroidal solution, the three-dimensional problem collapses to **six first-order ordinary differential equations in radius** for six radial functions, the standard "y-variables" of Alterman, Jarosch and Pekeris (1959) and Takeuchi and Saito (1972). Roughly: radial displacement, radial stress, tangential displacement, tangential stress, perturbed potential, and a potential-gradient term. **The ordering and normalization differ between papers.** Pick one reference (Sabadini and Vermeersen's textbook, Tobie et al. 2005, or Roberts and Nimmo 2008) and follow its definitions strictly, including its sign conventions. Mixing two sources is the single most common source of bugs.

$$
\frac{d\,\mathbf{y}}{dr} = \mathbf{A}(r;\ \rho,\ g,\ \mu^*,\ \lambda)\ \mathbf{y}, \qquad \mathbf{y} = (y_1,\dots,y_6)^T
$$

The 6×6 matrix **A** depends on local density, gravity, rigidity and the Lamé parameter, and you copy its entries from your chosen reference. Derive at least the incompressible case yourself to see where each term comes from.

### Boundary and interface conditions

- **Center.** The solution must be regular, which leaves **three** independent solutions. The total solution is a combination of them with three unknown constants.
- **Solid to solid interface.** All six y-variables are continuous.
- **Solid to liquid interface.** The liquid cannot support shear, so the tangential stress is zero on the solid side, the tangential displacement may jump, and the radial stress, potential and potential gradient stay continuous. Inside the liquid layer the system is reduced to the variables that remain meaningful, which you derive from the hydrostatic response of an inviscid fluid. This step is where most implementations stumble, and it deserves a careful derivation in your notes.
- **Surface.** Radial and tangential stress vanish (y2 = y4 = 0), and the potential-gradient variable is set by the applied tidal forcing. With unit potential normalization, the Love numbers are read off as k2 = y5(R) − 1, h2 = g y1(R), l2 = g y3(R). Check these conventions in your reference.

The three constants follow from three surface conditions, a 3×3 complex linear solve.

### Propagator matrix versus direct integration

In a layer of constant properties the ODE system can be solved analytically, and the solution at one radius is a **propagator matrix** times the solution at another, so a stack of layers is a product of matrices (Sabadini and Vermeersen). With thin layers or numerically stratified profiles, simple direct integration with a complex-valued ODE solver also works for degree 2 on these bodies. Do it in this order: start with direct integration (easier to debug), and move to a propagator matrix only if you need speed for the Phase 5 and Phase 8 sampling. If numbers grow unreasonably between layers, orthonormalize the three solution vectors as you go (Gram-Schmidt), which suppresses the dominant growing mode.

```python
def love_numbers(r, rho, g, mu_c, lam, omega, layers):
    # nondimensionalise: r/R, mu/(rho_s g_s R), etc.
    # for each of 3 start vectors regular at the center:
    #     integrate dy/dr = A(r) y through the layers (dtype=complex)
    #     apply interface conditions at every boundary
    # solve the 3x3 system for the surface conditions
    return k2, h2, l2    # complex
```

### Build it in seven small steps

1. **7a. Homogeneous incompressible elastic sphere.** Test: k2 = (3/2)/(1 + 19μ/(2ρgR)) and h2 = (5/2)/(1 + 19μ/(2ρgR)) to 1e-6.
2. **7b. Homogeneous Maxwell body.** Test: reproduce the complex k2 and the dissipation peak from Phase 6. This confirms your complex arithmetic and sign conventions.
3. **7c. Two solid layers** (core and mantle, different μ and ρ). Test: if both layers are given identical properties you recover 7a. A thin shell with ten times the mantle rigidity should lower k2 smoothly.
4. **7d. Fluid core under a solid mantle.** This is the first interface with a liquid. Test: shrink the solid mantle's rigidity to zero and recover the fluid limit k2 = 3/2.
5. **7e. Europa with an ocean.** Ice Ih shell, ocean, rock mantle and core from your Phase 4 model. Plot k2 and h2 against ice shell thickness. The values expected with an ocean are of order 0.2 to 0.3 for k2 and about 1 to 1.3 for h2, and without an ocean about 0.02 for k2 (Moore and Schubert 2000, and Wahr et al. 2006, verify). The contrast, a factor of ten, is why measuring k2 settles the ocean question.
6. **7f. Compressibility.** Add the bulk modulus and the full compressible matrix. Quantify how much k2 and h2 change, usually a few percent.
7. **7g. Cross-check against an independent code.** Run one case through `TidalPy` or a published tool (check its current status and conventions), and show agreement to about 1 percent. If you disagree, the disagreement is nearly always a convention or a layer-ordering mistake.

A useful physical anchor for step 7e: ice near its melting point has a viscosity of roughly 1e14 to 1e15 Pa s and a rigidity of about 3.5e9 Pa, giving a Maxwell time of order 3e4 to 3e5 s. Europa's orbital period is 3.55 days, about 3e5 s, so the dimensionless ωτ is of order 1. The warm base of the shell sits near the dissipation peak. Check that your own numbers reproduce this.

### Check

- 7a and 7b analytic tests pass to 1e-6.
- Fluid limit and rigid limit behave correctly in 7c and 7d.
- 7e: an ocean raises k2 by about an order of magnitude over the no-ocean case, and h2 depends strongly on shell thickness.
- Titan sanity check: a model with a global ocean gives a k2 of order 0.5 (the Cassini value is about 0.6, verify), while a fully solid model gives something close to zero.

### Gate

Why do k2 and h2 _together_ constrain the ice shell thickness better than k2 alone? Think about what each one measures (a potential perturbation versus a surface displacement), and which one saturates once the ocean is present.

## Phase 8: Tidal heating, rheology and thermal feedback (Io and Enceladus)

**Goal:** go from "how much does it deform" to "how hot does it get, where, and is that state stable". This is where the interior model becomes a story about volcanoes, oceans and thermostats, and where you can address Io's magma-ocean debate and Enceladus's missing heat directly.

### Step 1: where the heat is deposited

The surface formula from Phase 6 gives the _total_ power. The layered model also tells you the **volumetric heating rate** as a function of depth and position, from the strain amplitudes in each layer and the imaginary part of the local complex rigidity:

$$
h(r,\theta,\phi) \;=\; \frac{\omega}{2}\,\mathrm{Im}\!\left[\mu^*(r)\right]\ \left|\varepsilon(r,\theta,\phi)\right|^2 \ \ (\text{schematically, with the strain invariant built from the } y\text{-variables})
$$

Segatz et al. (1988) derived the expressions for Io, Tobie et al. (2005) and Roberts and Nimmo (2008) for ice shells, and Beuthe (2013) for the spatial patterns. Copy the exact strain expressions from one of these, with its conventions. Your **energy-conservation test** is powerful:

$$
\int_V h\ dV \;=\; -\frac{21}{2}\,\mathrm{Im}(k_2)\,\frac{G M_p^2\, n R^5 e^2}{a^6}
$$

If the volume integral of your dissipation does not match the surface formula to better than about 1 percent, there is a bug in either the strain expressions or the Love numbers. This single check catches most mistakes in Phase 7 too.

### Step 2: better rheologies

The Maxwell model is the simplest, and for tides it tends to **underestimate dissipation** away from the peak, because real rock and ice show a broad, frequency-dependent anelastic response. Two upgrades, in order:

- **Andrade rheology**, adding a transient creep term to the Maxwell compliance, with an exponent α typically around 0.2 to 0.4:

$$
J^{*(\omega)} = J_U + \beta\,\Gamma(1+\alpha)\,(i\omega)^{-\alpha} - \frac{i}{\omega\eta}, \qquad \mu^{*} = \frac{1}{J^{*}}
$$

Use the parameterization of Efroimsky (2012) or Renaud and Henning (2018) for β in terms of η and μ, and cite which you chose.

- **Burgers or Sundberg-Cooper**, which adds a second relaxation and is standard in recent Io work.

Also make viscosity depend on temperature and melt:

$$
\eta(T) = \eta_0\,\exp\!\left[\frac{E_a}{R_g}\left(\frac{1}{T} - \frac{1}{T_m}\right)\right]
$$

with activation energy E*a, and add a sharp drop in both rigidity and viscosity once the melt fraction passes a critical value (a \_rheological transition*, around 40 percent in much of the Io literature; verify). Viscosity changes by many orders of magnitude over a few hundred kelvin, so heating is exquisitely sensitive to temperature.

### Step 3: the thermal thermostat (Io)

Combine the dissipation curve from Phase 6 with a heat-loss law. A zero-dimensional mantle model is enough to show the idea:

$$
M c_p\,\frac{dT}{dt} = H_{\mathrm{tide}}(T) + H_{\mathrm{rad}} - Q_{\mathrm{loss}}(T)
$$

where H_tide(T) comes from your Love-number code with η(T), and Q_loss(T) is a convective or conductive loss, for example scaling with the temperature contrast to a power near 4/3 for convection. For Io, mantle melt transport ("heat pipes") is an important extra loss channel (O'Reilly and Davies 1981, Moore 2001).

Do this:

1. Plot H_tide(T) and Q_loss(T) on the same axes. H_tide **rises, peaks, then falls** as the body warms (the Phase 6 curve). Q_loss rises monotonically.
2. Every intersection is an equilibrium. It is stable where the loss curve is steeper than the heating curve. With a peaked H_tide you can get one, two or three equilibria depending on the parameters, which is the classic origin of a tidally heated body's "high-dissipation" and "low-dissipation" states.
3. Integrate dT/dt from several initial temperatures and watch which equilibrium each reaches. Report honestly whether any stable state delivers the observed output of about 1e14 W. Simple models struggle to match Io, which is partly why the interior structure and the melt distribution remain debated.

### Step 4: Io's magma ocean question, using your own code

A shallow global magma ocean would make Io _very soft_, with a large k2 (the models predict values several times larger than for a solid mantle). Recent analysis of Juno data (Park et al., published 2025, verify) argues that the measured k2 is much smaller than a shallow magma-ocean model gives. You can test the logic yourself:

1. Build Io with a liquid core, a mostly solid mantle, and a thin crust (Phase 7 machinery).
2. Insert a molten layer a few tens of km thick under a thin lithosphere. Compute k2.
3. Compare with the same Io without the molten layer, and with the measured value. Report which hypotheses the measurement rules out, and under which assumptions.

That last sentence is what real papers do, and it is a valid, publishable-style exercise.

### Step 5: Enceladus, a different puzzle

Enceladus radiates roughly 10 to 16 GW from its south pole (Howett et al. 2011, verify), with a tiny orbital eccentricity of about 0.005. Simple models that dissipate only in the ice shell give a few GW at most, which is well short of that. The proposed fix is dissipation in a **porous, fractured rocky core** (Choblet et al. 2017), which is warm, water-filled and soft. To explore it:

1. Build a layered Enceladus: ice shell of about 20 km (thinner near the south pole), a global or regional ocean, and a low-density rocky core, with density near 2400 kg/m³ (verify).
2. Compute the heating from the shell alone, then add a core with a soft (low-μ, viscous) rheology and see how total dissipation changes.
3. Discuss: which parameters of the core rheology are unconstrained, and what range of heating they allow?

### Step 6 (optional): coupled thermal and orbital evolution

Tidal dissipation also drains the orbit. Add an equation for e(t) and for the semimajor axis, include the Laplace resonance forcing, and integrate for billions of years. This is the territory of Fischer and Spohn (1990), Hussmann and Spohn (2004) and the thermal-orbital literature, and it can produce oscillations rather than a single steady state. It is a stretch goal, not a requirement.

### Check

- Volume-integrated dissipation matches the surface formula to better than 1 percent.
- Andrade versus Maxwell: you quantify the difference in |Im k2| at Io-like parameters.
- The H_tide and Q_loss intersection plot, with equilibria classified as stable or unstable.
- A short table of what k2 each Io hypothesis predicts, and which are disfavored by the observation.
- The Enceladus shell-only versus porous-core dissipation comparison.

### Gate

In plain words, explain why a tidally heated moon can have an unstable middle equilibrium, and what physical process would let Io switch from a low-dissipation state to a high-dissipation one.

## Phase 9: Data, validation and write-up

**Goal:** confront your code with real numbers, learn where the public data live, and package the work so that someone else can reproduce it and a specialist would take it seriously.

### Where the public data are

| Data                                                                                       | Where                                                                                                                                                                    | What you use it for                                                                                 |
| ------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | --------------------------------------------------------------------------------------------------- |
| Spacecraft radio science, magnetometer, imaging, spectra (Galileo, Cassini, Juno, Voyager) | NASA's Planetary Data System (PDS) nodes: Geosciences for radio science, Planetary Plasma Interactions for magnetometer data, Imaging and Atmospheres nodes for the rest | Raw material. Mostly you will not process it yourself at first                                      |
| Trajectories, orientations, geometry                                                       | NASA NAIF SPICE kernels, accessed from Python with `spiceypy`                                                                                                            | Computing orbital positions, tidal geometry, flyby timing                                           |
| Published gravity solutions (J2, C22, GM, Love numbers)                                    | Tables in the papers: Anderson et al. for Galileo (Io, Europa, Ganymede, Callisto), Iess et al. for Titan and Enceladus, Park et al. for Io from Juno                    | Your inputs and your validation targets                                                             |
| Open-source codes                                                                          | `PlanetProfile` (ocean-world structure), `TidalPy` (tidal and thermal-orbital), `burnman` (mineral physics), `SeaFreeze` (water and ice)                                 | Cross-checks, never replacements for understanding. Verify each is maintained and check conventions |

A note on expectations: turning raw Doppler tracking data into a gravity field is a full orbit-determination problem and is research-level work in itself. For this guide, start from the published gravity coefficients. If that later interests you, SPICE plus the PDS radio science archives are the entry point.

### Validation targets

Run these as automated tests and report every result, including failures. The reference values are approximate and from memory, so verify against the primary papers before you cite them.

| Quantity                       | Body      | Reference (approx.)     | Your phase |
| ------------------------------ | --------- | ----------------------- | ---------- |
| C/MR²                          | Ganymede  | 0.3115                  | 2, 3, 4    |
| C/MR²                          | Europa    | 0.346                   | 2, 4       |
| C/MR²                          | Io        | 0.377                   | 4          |
| k2 with global ocean / without | Europa    | about 0.25 / about 0.02 | 7          |
| k2                             | Titan     | about 0.6               | 7          |
| k2                             | Moon      | about 0.024             | 6, 7       |
| k2/Q                           | Io        | about 0.015             | 6, 8       |
| Total tidal heating            | Io        | about 1e14 W            | 8          |
| South polar heat output        | Enceladus | about 10 to 16 GW       | 8          |

### Writing it up

A good write-up is the best check that you understand your own results. Aim for a short technical note, structured like a paper:

1. **Question.** One sentence: what are you determining, for which body?
2. **Model.** Layers, EOS, rheology, with a table of every parameter, its prior and its source.
3. **Method.** Structure solver, Love-number solver, inversion, with the tests that validate each.
4. **Results.** Posterior on structure, predicted k2 and h2, predicted heating, each with uncertainties.
5. **Discussion.** What is robust, what is prior-driven, what a future measurement (Europa Clipper's gravity and magnetometer data, JUICE at Ganymede) would change.
6. **Reproducibility.** A README so that `pip install -e . && pytest && python run_all.py` rebuilds every figure.

Mission timing, as I know it: Europa Clipper arrives at Jupiter in 2030, JUICE reaches Jupiter in 2031, and Dragonfly launches in 2028 for Titan. Check the current schedules before you quote any of them, since launch and arrival dates move.

### Capstone options

Pick one to finish, and treat it as your first real project:

- **Europa:** a combined structure and Love-number inversion that predicts what k2 and h2 Europa Clipper should see, as a function of ice shell thickness and ocean salinity.
- **Ganymede:** how much the high-pressure ice layers change the predicted k2 and h2, and whether a sandwiched ocean is distinguishable from a simple one.
- **Io:** a systematic scan of mantle melt fraction against k2 and dissipation, in the Park et al. framework.
- **Enceladus:** the shell-only versus porous-core heating comparison, with a posterior over core rheology.
- **Titan:** how a slushy rather than liquid ocean changes the k2 prediction (an active recent question, verify the literature).

### Check

- Every row of the validation table has a recorded pass, near miss with explanation, or failure with diagnosis.
- `pytest` passes from a clean checkout.
- A written note of about ten pages that a colleague with a physics background could follow.

### Gate

If you showed your note to a planetary scientist, what is the first question they would ask about your priors or your rheology choice? Write that question, and your answer, into the discussion before they have to ask.

## Reading list, pitfalls and milestone checklist

### Reading list

All references below are from my memory of the literature, not freshly opened sources. Look up each one before you rely on it, and prefer the primary paper when a number matters.

| Reference                                                                                                       | Use it for                                                    | Phase   |
| --------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------- | ------- |
| Turcotte and Schubert, _Geodynamics_                                                                            | Elasticity, heat transfer, viscoelasticity, the basic toolkit | 3, 6, 8 |
| Murray and Dermott, _Solar System Dynamics_                                                                     | Tides, resonances, orbital elements                           | 1, 6    |
| Hubbard, _Planetary Interiors_, and Zharkov and Trubitsyn, _Physics of Planetary Interiors_                     | Clairaut and Radau theory, hydrostatic figures                | 3       |
| Nimmo and Pappalardo (2016), review of ocean worlds in the outer solar system                                   | The big picture and vocabulary                                | 1       |
| Schubert, Anderson, Spohn and McKinnon (2004), chapter on the Galilean satellites' interiors                    | Standard interior models and data                             | 4, 5    |
| Anderson et al. (1996, 1998, 2001), Galileo gravity papers for Ganymede, Europa, Io                             | Your published J2, C22 and C/MR² inputs                       | 2, 4    |
| Gao and Stevenson (2013), non-hydrostatic effects on moment of inertia                                          | Limits of the Radau-Darwin approach                           | 2, 5    |
| Vance et al. (2018), ocean-world interior modeling (PlanetProfile)                                              | A complete reference implementation to compare against        | 4, 5    |
| Foreman-Mackey et al. (2013), the `emcee` paper, and Sivia and Skilling, _Data Analysis: A Bayesian Tutorial_   | MCMC practice and Bayesian reasoning                          | 5       |
| Takeuchi and Saito (1972); Sabadini and Vermeersen, _Global Dynamics of the Earth_                              | The y-variable equations and propagator matrices              | 7       |
| Moore and Schubert (2000), the tidal response of Europa; Wahr et al. (2006), tides and Europa's shell thickness | Expected k2 and h2 with and without an ocean                  | 7       |
| Tobie, Mocquet and Sotin (2005); Roberts and Nimmo (2008)                                                       | Dissipation in icy bodies with layered models                 | 7, 8    |
| Peale, Cassen and Reynolds (1979); Segatz et al. (1988)                                                         | Io's tidal heating, from prediction to layered model          | 6, 8    |
| Efroimsky (2012); Renaud and Henning (2018)                                                                     | Andrade and other anelastic rheologies in tidal codes         | 8       |
| Lainey et al. (2009); Park et al. (2025)                                                                        | Io's k2/Q from astrometry and its k2 from Juno                | 6, 8    |
| Iess et al. (2010, 2012, 2014)                                                                                  | Titan and Enceladus gravity, Titan's k2                       | 2, 7, 8 |
| Beuthe (2013, 2016); Choblet et al. (2017)                                                                      | Spatial heating patterns and Enceladus's porous core          | 8       |

### Pitfalls that will cost you days

1. **Normalization of gravity coefficients.** Unnormalized versus fully normalized J2 and C22 differ by factors of √5 or √(5/12). Check which one each paper uses.
2. **Mean motion versus rotation rate.** For a synchronous moon they are equal, but make sure the code uses n in rad/s from the period in days, not in revolutions per day.
3. **Sign conventions for complex quantities.** The choice between e^{iωt} and e^{−iωt} flips the sign of every imaginary part. Fix one convention in `notes/` and never mix it with a source that uses the other.
4. **Mixed y-variable conventions** between textbooks and papers. Follow a single reference throughout.
5. **Layer boundaries not on grid nodes**, which smears density jumps and ruins convergence tests.
6. **Unit slips** (km versus m, GPa versus Pa). Work in SI internally, convert at the edges.
7. **Assuming hydrostatic equilibrium blindly.** Always report the J2/C22 ratio and add a systematic error to C/MR².
8. **Believing a posterior that your prior produced.** Run the prior-sensitivity test every time.
9. **Using a constant Q** as though it were a material property. Dissipation depends on frequency, temperature and structure. Use it only as a quick estimate.
10. **Forgetting the right forcing frequency.** The eccentricity tide in a synchronous moon forces at ω = n.

### Milestone checklist

- [ ] Phase 0: repo, constants file, mean-density check passes
- [ ] Phase 1: surface gravities and uniform-body central pressures reproduced
- [ ] Phase 2: C/MR² from C22 for Ganymede and Europa, J2/C22 ratio checked
- [ ] Phase 3: homogeneous and two-layer solver validated, Radau-Darwin error measured
- [ ] Phase 4: compressible layered solver with Io, Europa and Ganymede models and convergence plot
- [ ] Phase 5: MCMC with synthetic recovery test and prior-sensitivity table
- [ ] Phase 6: analytic Love numbers, dissipation peak and Io power budget reproduced
- [ ] Phase 7: layered viscoelastic Love numbers, 7a to 7g tests passing, cross-check with an independent code
- [ ] Phase 8: energy-conservation test, Andrade versus Maxwell, thermal equilibria, Io and Enceladus analyses
- [ ] Phase 9: validation table complete, write-up and reproducible repo

If you get stuck on any phase, bring me the specific equation or test that fails and what you got instead of the expected value. That is usually enough to locate the bug quickly.
