from pathlib import Path

import numpy as np
import ogstools as ot


ROOT = Path.cwd()

RUN = (
    ROOT
    / "results"
    / "phase04d_mc_hm"
)

PVD = sorted(
    RUN.glob("*.pvd")
)

if not PVD:
    raise SystemExit(
        "FAIL: no Phase 04D PVD"
    )

series = ot.MeshSeries(
    str(PVD[0])
)

times = np.asarray(
    series.timevalues,
    dtype=float,
)

if len(times) < 2:
    raise SystemExit(
        "FAIL: insufficient output timesteps"
    )

mesh = series.mesh(
    len(times) - 1
)


# ============================================================
# FIELD HELPERS
# ============================================================

def field(mesh, name):

    if name in mesh.point_data:

        return (
            np.asarray(
                mesh.point_data[name],
                dtype=float,
            ),
            "point",
        )

    if name in mesh.cell_data:

        return (
            np.asarray(
                mesh.cell_data[name],
                dtype=float,
            ),
            "cell",
        )

    raise KeyError(
        f"{name!r} missing. "
        f"Point fields: "
        f"{list(mesh.point_data.keys())}; "
        f"Cell fields: "
        f"{list(mesh.cell_data.keys())}"
    )


def scalar_field(mesh, name):

    arr, location = field(
        mesh,
        name,
    )

    arr = np.asarray(
        arr,
        dtype=float,
    ).squeeze()

    return arr, location


def surface_height(x):

    x = np.asarray(
        x,
        dtype=float,
    )

    h = np.full_like(
        x,
        10.0,
    )

    mid = (
        (x > 8.0)
        & (x < 22.0)
    )

    h[mid] = (
        10.0
        - 8.0
        * (x[mid] - 8.0)
        / 14.0
    )

    h[x >= 22.0] = 2.0

    return h


# ============================================================
# DISPLACEMENT
# ============================================================

u, u_location = field(
    mesh,
    "displacement",
)

if u_location != "point":
    raise SystemExit(
        "FAIL: displacement is not point data"
    )

u = np.asarray(
    u,
    dtype=float,
)

u_mag = np.linalg.norm(
    u,
    axis=1,
)

u_idx = int(
    np.nanargmax(
        u_mag
    )
)

u_max = float(
    u_mag[u_idx]
)

u_xyz = mesh.points[
    u_idx
]

u_vec = u[
    u_idx
]


# ============================================================
# SLOPE-FACE DISPLACEMENT
# ============================================================

pts = np.asarray(
    mesh.points,
    dtype=float,
)

local_surface = surface_height(
    pts[:, 0]
)

depth = (
    local_surface
    - pts[:, 1]
)

slope_zone = (
    (pts[:, 0] >= 8.0)
    & (pts[:, 0] <= 22.0)
    & (depth >= -1e-8)
    & (depth <= 2.0)
)

if not np.any(
    slope_zone
):
    raise SystemExit(
        "FAIL: slope zone empty"
    )

slope_u_max = float(
    np.nanmax(
        u_mag[slope_zone]
    )
)


# ============================================================
# EQUIVALENT PLASTIC STRAIN
# ============================================================

epsp, epsp_location = scalar_field(
    mesh,
    "EquivalentPlasticStrain",
)

if epsp_location == "point":

    epsp_coords = pts

else:

    epsp_coords = (
        mesh
        .cell_centers()
        .points
    )

epsp_idx = int(
    np.nanargmax(
        epsp
    )
)

epsp_max = float(
    epsp[epsp_idx]
)

epsp_min = float(
    np.nanmin(
        epsp
    )
)

epsp_xyz = epsp_coords[
    epsp_idx
]


# Diagnostic thresholds only.
fractions = {}

for threshold in [
    1e-8,
    1e-6,
    1e-4,
    1e-3,
]:

    fractions[threshold] = float(
        np.mean(
            epsp > threshold
        )
    )


# ============================================================
# HYDRAULIC FIELDS
# ============================================================

pressure, _ = scalar_field(
    mesh,
    "pressure",
)

saturation, _ = scalar_field(
    mesh,
    "saturation",
)

porosity, _ = scalar_field(
    mesh,
    "porosity",
)


# ============================================================
# LATE-TIME DISPLACEMENT INCREMENT
# ============================================================

prev_mesh = series.mesh(
    len(times) - 2
)

u_prev, prev_location = field(
    prev_mesh,
    "displacement",
)

if prev_location != "point":
    raise SystemExit(
        "FAIL: previous displacement not point data"
    )

u_prev = np.asarray(
    u_prev,
    dtype=float,
)

du = (
    u
    - u_prev
)

du_mag = np.linalg.norm(
    du,
    axis=1,
)

final_du_max = float(
    np.nanmax(
        du_mag
    )
)

increment_ratio = (
    final_du_max / u_max
    if u_max > 0.0
    else np.nan
)


# ============================================================
# FINITE CHECK
# ============================================================

finite = all(
    [
        np.all(
            np.isfinite(u)
        ),
        np.all(
            np.isfinite(epsp)
        ),
        np.all(
            np.isfinite(pressure)
        ),
        np.all(
            np.isfinite(saturation)
        ),
        np.all(
            np.isfinite(porosity)
        ),
    ]
)


# ============================================================
# SCREENING CLASSIFICATION
# ============================================================

if not finite:

    numerical_status = (
        "FAIL: NON-FINITE FIELD"
    )

elif epsp_min < -1e-8:

    numerical_status = (
        "REVIEW: NEGATIVE PLASTIC STRAIN"
    )

else:

    numerical_status = (
        "PASS"
    )


if epsp_max <= 1e-8:

    plastic_pattern = (
        "NO RESOLVED YIELDING UNDER BASELINE"
    )

elif fractions[1e-4] <= 0.20:

    plastic_pattern = (
        "LIMITED / LOCALIZED PLASTICITY"
    )

else:

    plastic_pattern = (
        "WIDESPREAD PLASTICITY — REVIEW"
    )


# ============================================================
# PRINT
# ============================================================

print(
    "========================================"
)

print(
    "PHASE 04D COUPLED HM + MOHR-COULOMB"
)

print(
    "========================================"
)

print()

print(
    f"Outputs: {len(times)}"
)

print(
    "Final time [days]: "
    f"{times[-1]/86400.0:.6f}"
)

print()

print(
    "=== DEFORMATION ==="
)

print(
    "Global max |u| [mm]: "
    f"{u_max*1000.0:.8f}"
)

print(
    "Global max location [m]: "
    f"({u_xyz[0]:.6f}, "
    f"{u_xyz[1]:.6f})"
)

print(
    "Global max components [mm]: "
    f"({u_vec[0]*1000.0:.8f}, "
    f"{u_vec[1]*1000.0:.8f})"
)

print(
    "Slope-zone max |u| [mm]: "
    f"{slope_u_max*1000.0:.8f}"
)

print(
    "Final-step max Δu [mm]: "
    f"{final_du_max*1000.0:.10f}"
)

print(
    "Final-step Δu / total max |u|: "
    f"{increment_ratio:.10e}"
)

print()

print(
    "=== PLASTICITY ==="
)

print(
    "EquivalentPlasticStrain data location: "
    f"{epsp_location}"
)

print(
    "EquivalentPlasticStrain min: "
    f"{epsp_min:.10e}"
)

print(
    "EquivalentPlasticStrain max: "
    f"{epsp_max:.10e}"
)

print(
    "Max plastic-strain location [m]: "
    f"({epsp_xyz[0]:.6f}, "
    f"{epsp_xyz[1]:.6f})"
)

for threshold, frac in fractions.items():

    print(
        f"Fraction eps_p > {threshold:.0e}: "
        f"{100.0*frac:.4f}%"
    )

print()

print(
    "Plasticity screening: "
    f"{plastic_pattern}"
)

print()

print(
    "=== HYDRAULICS ==="
)

print(
    "Pressure range [kPa]: "
    f"{np.nanmin(pressure)/1000.0:.6f} "
    f"to "
    f"{np.nanmax(pressure)/1000.0:.6f}"
)

print(
    "Saturation range [-]: "
    f"{np.nanmin(saturation):.8f} "
    f"to "
    f"{np.nanmax(saturation):.8f}"
)

print(
    "Porosity range [-]: "
    f"{np.nanmin(porosity):.8f} "
    f"to "
    f"{np.nanmax(porosity):.8f}"
)

print()

print(
    "NUMERICAL STATUS: "
    f"{numerical_status}"
)


# ============================================================
# SUMMARY FILE
# ============================================================

summary = (
    RUN
    / "phase04d_summary.txt"
)

lines = [
    "PHASE 04D COUPLED HM + MOHR-COULOMB",
    "",
    f"Final time [days]: {times[-1]/86400.0:.8f}",
    (
        "Global max displacement [mm]: "
        f"{u_max*1000.0:.10f}"
    ),
    (
        "Global max location [m]: "
        f"{u_xyz[0]:.8f}, "
        f"{u_xyz[1]:.8f}"
    ),
    (
        "Slope-zone max displacement [mm]: "
        f"{slope_u_max*1000.0:.10f}"
    ),
    (
        "Final-step max delta-u [mm]: "
        f"{final_du_max*1000.0:.12f}"
    ),
    (
        "Final-step delta-u / total max-u: "
        f"{increment_ratio:.12e}"
    ),
    "",
    (
        "EquivalentPlasticStrain min: "
        f"{epsp_min:.12e}"
    ),
    (
        "EquivalentPlasticStrain max: "
        f"{epsp_max:.12e}"
    ),
    (
        "Max plastic-strain location [m]: "
        f"{epsp_xyz[0]:.8f}, "
        f"{epsp_xyz[1]:.8f}"
    ),
]

for threshold, frac in fractions.items():

    lines.append(
        f"Fraction eps_p > {threshold:.0e}: "
        f"{100.0*frac:.8f}%"
    )

lines += [
    "",
    (
        "Plasticity screening: "
        f"{plastic_pattern}"
    ),
    "",
    (
        "Pressure range [kPa]: "
        f"{np.nanmin(pressure)/1000.0:.8f} "
        f"to "
        f"{np.nanmax(pressure)/1000.0:.8f}"
    ),
    (
        "Saturation range [-]: "
        f"{np.nanmin(saturation):.10f} "
        f"to "
        f"{np.nanmax(saturation):.10f}"
    ),
    (
        "Porosity range [-]: "
        f"{np.nanmin(porosity):.10f} "
        f"to "
        f"{np.nanmax(porosity):.10f}"
    ),
    "",
    (
        "NUMERICAL STATUS: "
        f"{numerical_status}"
    ),
    "",
    (
        "NOTE: c=25 kPa, phi=28 deg, psi=0 deg "
        "are provisional screening values only."
    ),
]

summary.write_text(
    "\n".join(lines) + "\n",
    encoding="utf-8",
)

print()
print(
    f"PASS: {summary}"
)
