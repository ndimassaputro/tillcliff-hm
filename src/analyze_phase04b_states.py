from pathlib import Path

import numpy as np
import pyvista as pv
import ogstools as ot


ROOT = Path.cwd()

RUN_ROOT = (
    ROOT
    / "results"
    / "phase04b_states"
)

OUT = (
    ROOT
    / "results"
    / "phase04b_analysis"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)

STATES = [
    "dry",
    "reference",
    "wet",
]


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


def load_final(state):

    folder = RUN_ROOT / state

    pvds = sorted(
        folder.glob("*.pvd")
    )

    if not pvds:
        raise SystemExit(
            f"FAIL: no PVD for {state}"
        )

    series = ot.MeshSeries(
        str(pvds[0])
    )

    times = np.asarray(
        series.timevalues,
        dtype=float,
    )

    mesh = series.mesh(
        len(times) - 1
    )

    return (
        mesh,
        times[-1],
    )


def prepare_cell_fields(mesh):

    m = mesh.copy()

    # Convert point fields to cell averages when needed.
    if (
        "pressure" in m.point_data
        or "saturation" in m.point_data
    ):
        m = m.point_data_to_cell_data(
            pass_point_data=True
        )

    m = m.compute_cell_sizes(
        length=False,
        area=True,
        volume=False,
    )

    return m


def get_cell_array(mesh, name):

    if name in mesh.cell_data:
        return np.asarray(
            mesh.cell_data[name],
            dtype=float,
        ).squeeze()

    raise KeyError(
        f"{name!r} not found in cell data. "
        f"Available: {list(mesh.cell_data.keys())}"
    )


def get_point_array(mesh, name):

    if name not in mesh.point_data:
        raise KeyError(
            f"{name!r} not found in point data"
        )

    return np.asarray(
        mesh.point_data[name],
        dtype=float,
    )


records = {}

for state in STATES:

    mesh, final_time = load_final(
        state
    )

    cell_mesh = prepare_cell_fields(
        mesh
    )

    centers = (
        cell_mesh
        .cell_centers()
        .points
    )

    x = centers[:, 0]
    y = centers[:, 1]

    local_surface = surface_height(
        x
    )

    depth = (
        local_surface
        - y
    )

    active = (
        (depth >= -1e-8)
        & (depth <= 2.0)
    )

    if np.sum(active) < 10:
        raise SystemExit(
            f"FAIL: too few active-zone cells for {state}"
        )

    area = get_cell_array(
        cell_mesh,
        "Area",
    )

    pressure = get_cell_array(
        cell_mesh,
        "pressure",
    )

    saturation = get_cell_array(
        cell_mesh,
        "saturation",
    )

    suction = np.maximum(
        -pressure,
        0.0,
    )

    # Bishop exponent = 1 in this model:
    # chi = Sr.
    effective_suction = (
        saturation
        * suction
    )

    w = area[active]

    active_mean_suction = np.average(
        suction[active],
        weights=w,
    )

    active_mean_sat = np.average(
        saturation[active],
        weights=w,
    )

    active_effective_suction = np.average(
        effective_suction[active],
        weights=w,
    )

    # --------------------------------------------------------
    # Mechanical metrics
    # --------------------------------------------------------

    disp = get_point_array(
        mesh,
        "displacement",
    )

    mag = np.linalg.norm(
        disp,
        axis=1,
    )

    pts = mesh.points

    max_idx = int(
        np.nanargmax(mag)
    )

    # Slope-face zone:
    # x=8..22 and points close to local surface.
    point_surface = surface_height(
        pts[:, 0]
    )

    point_depth = (
        point_surface
        - pts[:, 1]
    )

    slope_zone = (
        (pts[:, 0] >= 8.0)
        & (pts[:, 0] <= 22.0)
        & (point_depth >= -1e-8)
        & (point_depth <= 2.0)
    )

    slope_max_disp = float(
        np.nanmax(
            mag[slope_zone]
        )
    )

    records[state] = {
        "final_day":
            final_time / 86400.0,

        "active_mean_suction_kpa":
            active_mean_suction / 1000.0,

        "active_mean_saturation":
            active_mean_sat,

        "effective_suction_pa":
            active_effective_suction,

        "max_disp_mm":
            float(
                mag[max_idx] * 1000.0
            ),

        "max_disp_x":
            float(
                pts[max_idx, 0]
            ),

        "max_disp_y":
            float(
                pts[max_idx, 1]
            ),

        "slope_zone_max_disp_mm":
            slope_max_disp * 1000.0,

        "active_cells":
            int(
                np.sum(active)
            ),
    }


# ============================================================
# Candidate Antecedent Suction-Loss Index
# ============================================================

q_ref = records[
    "reference"
]["effective_suction_pa"]

for state in STATES:

    q = records[
        state
    ]["effective_suction_pa"]

    As = (
        1.0
        - q / q_ref
    )

    records[
        state
    ]["As"] = As


# Incremental displacement relative to reference.
ref_slope_u = records[
    "reference"
]["slope_zone_max_disp_mm"]

for state in STATES:

    records[
        state
    ]["delta_slope_disp_vs_ref_mm"] = (
        records[state][
            "slope_zone_max_disp_mm"
        ]
        - ref_slope_u
    )


# ============================================================
# Diagnostics
# ============================================================

print(
    "========================================"
)
print(
    "PHASE 04B ANTECEDENT HYDRAULIC STATES"
)
print(
    "========================================"
)

print()

for state in STATES:

    r = records[state]

    print(
        f"{state.upper():9s} | "
        f"As={r['As']:+.6f} | "
        f"mean suction="
        f"{r['active_mean_suction_kpa']:.3f} kPa | "
        f"mean Sr="
        f"{r['active_mean_saturation']:.6f} | "
        f"slope-zone umax="
        f"{r['slope_zone_max_disp_mm']:.4f} mm | "
        f"Δu(ref)="
        f"{r['delta_slope_disp_vs_ref_mm']:+.4f} mm"
    )


# ============================================================
# Ordering checks
# ============================================================

dry_As = records[
    "dry"
]["As"]

ref_As = records[
    "reference"
]["As"]

wet_As = records[
    "wet"
]["As"]

ordering_pass = (
    dry_As < ref_As
    and ref_As < wet_As
)

ref_zero_pass = (
    abs(ref_As) < 1e-10
)


print()
print("=== CANDIDATE As CHECK ===")

print(
    f"Dry As       = {dry_As:+.8f}"
)

print(
    f"Reference As = {ref_As:+.8f}"
)

print(
    f"Wet As       = {wet_As:+.8f}"
)

if ordering_pass:
    print(
        "STATE ORDERING: PASS"
    )
else:
    print(
        "STATE ORDERING: REVIEW"
    )

if ref_zero_pass:
    print(
        "REFERENCE NORMALIZATION: PASS"
    )
else:
    print(
        "REFERENCE NORMALIZATION: REVIEW"
    )


# ============================================================
# CSV
# ============================================================

csv = (
    OUT
    / "antecedent_state_metrics.csv"
)

with csv.open(
    "w",
    encoding="utf-8",
) as f:

    f.write(
        "state,"
        "As,"
        "active_mean_suction_kpa,"
        "active_mean_saturation,"
        "effective_suction_pa,"
        "slope_zone_max_disp_mm,"
        "delta_slope_disp_vs_ref_mm,"
        "global_max_disp_mm,"
        "global_max_x_m,"
        "global_max_y_m,"
        "active_cells\n"
    )

    for state in STATES:

        r = records[state]

        f.write(
            f"{state},"
            f"{r['As']:.10f},"
            f"{r['active_mean_suction_kpa']:.10f},"
            f"{r['active_mean_saturation']:.10f},"
            f"{r['effective_suction_pa']:.10f},"
            f"{r['slope_zone_max_disp_mm']:.10f},"
            f"{r['delta_slope_disp_vs_ref_mm']:.10f},"
            f"{r['max_disp_mm']:.10f},"
            f"{r['max_disp_x']:.10f},"
            f"{r['max_disp_y']:.10f},"
            f"{r['active_cells']}\n"
        )


# ============================================================
# Summary
# ============================================================

summary = (
    OUT
    / "phase04b_summary.txt"
)

lines = [
    "PHASE 04B ANTECEDENT HYDRAULIC STATES",
    "",
]

for state in STATES:

    r = records[state]

    lines.append(
        f"{state.upper():9s} | "
        f"As={r['As']:+.6f} | "
        f"mean_suction="
        f"{r['active_mean_suction_kpa']:.3f} kPa | "
        f"mean_Sr="
        f"{r['active_mean_saturation']:.6f} | "
        f"slope_u="
        f"{r['slope_zone_max_disp_mm']:.4f} mm | "
        f"delta_u_ref="
        f"{r['delta_slope_disp_vs_ref_mm']:+.4f} mm"
    )

lines += [
    "",
    (
        "STATE ORDERING: "
        + (
            "PASS"
            if ordering_pass
            else "REVIEW"
        )
    ),
    (
        "REFERENCE NORMALIZATION: "
        + (
            "PASS"
            if ref_zero_pass
            else "REVIEW"
        )
    ),
    "",
    (
        "STATUS: PASS"
        if ordering_pass
        else "STATUS: REVIEW"
    ),
]

summary.write_text(
    "\n".join(lines) + "\n",
    encoding="utf-8",
)

print()
print("PASS:", csv)
print("PASS:", summary)

print()
print(
    "========================================"
)
print(
    "PHASE 04B ANALYSIS COMPLETE"
)
print(
    "========================================"
)
