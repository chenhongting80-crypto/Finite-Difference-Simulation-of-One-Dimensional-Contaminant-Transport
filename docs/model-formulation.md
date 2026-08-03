# Model formulation

## Conceptual domain

The model represents a 65 m one-dimensional groundwater-flow path.

### Method 1: distributed treatment

- Zone 1: 35 m TCE source zone with mass transfer and transport.
- Zone 2: 30 m treatment zone with transport and sequential biodegradation.

### Method 2: permeable reactive biobarrier

- Zone 1: 35 m source zone.
- Zone 2: 20 m nonreactive transition zone.
- Zone 3: 10 m permeable reactive barrier containing the reactive microbial population.

## Governing equation

For each dissolved species,

$$
\frac{\partial C}{\partial t}
= -v_x\frac{\partial C}{\partial x}
+ D_x\frac{\partial^2C}{\partial x^2}
+ R.
$$

The pore-water velocity is

$$
v_x = \frac{Ki}{n},
$$

where $K$ is hydraulic conductivity, $i$ is hydraulic gradient, and $n$ is porosity.

The longitudinal hydrodynamic dispersion coefficient is

$$
D_x = v_x\alpha_L + n\tau D_0,
$$

where $\alpha_L$ is longitudinal dispersivity, $\tau$ is tortuosity, and $D_0$ is molecular diffusion in water.

## Source-zone mass transfer

TCE release in Zone 1 is represented as

$$
R_{source} = \left(\frac{K_TA}{V}\right)(C_{sat}^{w} - C_i),
$$

where $C_{sat}^{w}$ is the aqueous TCE solubility and $C_i$ is the local dissolved concentration.

The source is finite. Release stops when the remaining source mass falls below 0.1% of its initial value.

## Sequential biodegradation

The reaction sequence is

$$
\mathrm{TCE} \rightarrow \mathrm{cDCE} \rightarrow \mathrm{VC} \rightarrow \mathrm{ethene/ethane}.
$$

The substrate-consumption term is represented by

$$
R = -\frac{\mu_{max}CX}{Y(K_s+C)},
$$

where $\mu_{max}$ is the maximum specific rate, $C$ is substrate concentration, $X$ is cell density, $Y$ is yield, and $K_s$ is the half-saturation constant.

For cDCE and VC, the numerical update includes production from the preceding compound and consumption to the following compound.

## Numerical method

The model uses an explicit time update, a centered first derivative for advection, and a centered second derivative for dispersion:

$$
\left.\frac{\partial C}{\partial x}\right|_i
\approx \frac{C_{i+1}-C_{i-1}}{2\Delta x},
$$

$$
\left.\frac{\partial^2C}{\partial x^2}\right|_i
\approx \frac{C_{i+1}-2C_i+C_{i-1}}{\Delta x^2},
$$

$$
C_i^{n+1} = C_i^n + \Delta t
\left(
-v_x\frac{\partial C}{\partial x}
+D_x\frac{\partial^2C}{\partial x^2}
+R
\right).
$$

The report applies the following step-size checks:

$$
\Delta x < \frac{2D_x}{v_x},
\qquad
\Delta t < \frac{\Delta x^2}{2D_x}.
$$

## Boundary conditions

- At the upstream boundary, a Dirichlet condition sets concentration to zero.
- At the downstream boundary of each zone, a ghost node approximates a zero-gradient Neumann condition.
- At zone interfaces, the downstream zone inlet receives the TCE concentration from the preceding zone.
- cDCE and VC enter each reactive zone at zero concentration and are generated internally by sequential transformation.
