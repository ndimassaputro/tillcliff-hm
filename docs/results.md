# Results

## Antecedent-state response

At identical removed toe geometry, dry, reference, and wet hydraulic states
produce strongly different incremental distributed deformation.

At nominal recession \(E=0.25\) m:

- dry/reference RMS response ratio: **1.2456**;
- wet/reference RMS response ratio: **0.1558**.

The monitored slope body remains elastic over the common comparison range.

## Hydraulic versus inherited-stress contribution

Holding inherited reference effective stress fixed while changing hydraulic
state reproduces almost the entire dry/wet response separation.

At \(E=0.25\) m:

- dry pressure + reference stress:
  **1.2456 × reference**;
- wet pressure + reference stress:
  **0.1558 × reference**.

By comparison, changing inherited stress while holding reference hydraulic
state gives:

- reference pressure + dry stress:
  **1.0027 × reference**;
- reference pressure + wet stress:
  **0.9876 × reference**.

The tested response is therefore dominated by the antecedent hydraulic-state
component.

## Spatial deformation mode

Dry and reference responses retain almost identical displacement-vector shape:

- dry/reference cosine similarity remains above **0.99968**.

Wet/reference similarity decreases with erosion:

- \(E=0.05\) m: **0.92488**;
- \(E=0.15\) m: **0.92561**;
- \(E=0.25\) m: **0.76888**.

Thus the wet-state response is not simply a scalar reduction of the reference
response.

## Signal-to-drift verification

The smallest wet-state signal-to-intact-drift ratio is approximately:

\[
3.07\times10^3.
\]

The resolved mode-shape difference is therefore far above the measured
late-intact numerical displacement drift.

## Negative result: no defensible erosion failure threshold

Attempts to define an erosion threshold from solver nonconvergence were
rejected after mesh refinement.

Moving-front nonconvergence and local EquivalentPlasticStrain were strongly
dependent on the discrete element-removal event sequence.

The local plastic hotspot followed the notch interface on all tested meshes.

Therefore no physical \(E_\mathrm{crit}\), factor of safety, or landslide
initiation threshold is reported.
