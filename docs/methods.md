# Methods

## Model concept

The model represents a two-dimensional North-European clayey glacial-till
slope subjected to seasonal unsaturated hydraulic forcing and prescribed
coastal toe recession.

The modelling chain is:

1. coupled seasonal Richards-mechanics simulation;
2. extraction of dry, reference, and wet antecedent hydraulic states;
3. transfer to a Mohr–Coulomb coupled hydro-mechanical slope model;
4. intact equilibrium verification;
5. prescribed toe recession using a moving element-deactivation front;
6. distributed slope-body deformation analysis;
7. numerical robustness and attribution diagnostics.

## Governing processes

The OpenGeoSys RichardsMechanics process couples unsaturated flow and
deformation.

The mechanical screening description uses a local Mohr–Coulomb constitutive
model.

The analysis deliberately distinguishes:

- pore-pressure / saturation state;
- inherited effective-stress state;
- deformation response induced after toe recession.

## Antecedent states

Three seasonal antecedent states are analysed:

- dry;
- reference;
- wet.

After transfer to the production medium mesh, all states pass an intact
Mohr–Coulomb equilibrium hold with zero EquivalentPlasticStrain.

## Toe recession

Coastal erosion is represented as prescribed basal toe-notch recession using
OpenGeoSys element deactivation.

The continuation coordinate is numerical rather than physical time.

Actual removed element area is retained alongside nominal horizontal recession
because the deactivation geometry is discrete.

## Distributed deformation metric

The primary response is the RMS norm of displacement increment over a fixed
slope-body monitoring zone:

```math
R_\mathrm{RMS}
=
\sqrt{
\frac{1}{N}
\sum_i
\|\Delta\mathbf u_i\|^2
}.
```

The monitoring zone excludes the immediate coastal notch-front region.

A 95th-percentile displacement magnitude is retained as a complementary
distributed response diagnostic.

## Pressure/stress decomposition

To distinguish hydraulic-state effects from inherited stress-history effects,
five numerical combinations are evaluated:

- reference pressure + reference stress;
- dry pressure + reference stress;
- wet pressure + reference stress;
- reference pressure + dry stress;
- reference pressure + wet stress.

These mixed states are numerical attribution experiments. They are not claimed
to correspond to realizable physical loading histories.

## Deformation-mode similarity

Displacement-field similarity is quantified by cosine similarity:

```math
C =
\frac{
\Delta\mathbf u_i
\cdot
\Delta\mathbf u_\mathrm{ref}
}{
\|\Delta\mathbf u_i\|
\|\Delta\mathbf u_\mathrm{ref}\|
}.
```

A best-fit scalar amplitude and normalized residual are also calculated.

## Signal-floor diagnostic

The final intact-hold displacement increments are used as a numerical drift
proxy.

For each erosion state,

```math
\mathrm{SNR}
=
\frac{
R_\mathrm{erosion}
}{
R_\mathrm{intact\,drift}
}.
```

The minimum wet-state ratio exceeded $3\times10^3$, supporting interpretation
of the observed wet/reference mode difference as numerically resolved.

## Mesh-sensitivity safeguard

Coarse, medium, and fine local toe meshes were tested.

Although intact medium/fine solutions were closely consistent, moving
element-deactivation events produced mesh-dependent local plasticity and
nonconvergence.

Plastic-strain hotspots remained adjacent to the active/eroded interface.

Consequently, neither local plastic strain nor solver nonconvergence is used as
a physical landslide-initiation threshold.
