from pathlib import Path
import csv
import re

import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot
import pyvista as pv


ROOT = Path.cwd()

RUN = (
    ROOT
    / "results"
    / "phase05b_v2_moving_front"
)

MODEL = (
    ROOT
    / "model"
    / "phase05b_v2_moving_front"
)

OUT = (
    ROOT
    / "results"
    / "phase05b_v2_analysis"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# READ FAILURE TIME FROM LOG
#
# Critical correction relative to Phase 05B-v1:
#
# A failed timestep is NOT a valid equilibrium state.
# Any output at or after failure time is excluded.
# ============================================================

log_path = (
    RUN
    / "ogs.log"
)

failure_time = None

if log_path.exists():

    log_text = log_path.read_text(
        encoding="utf-8",
        errors="ignore",
    )

    matches = re.findall(
        r"nonlinear solver failed .*?"
        r"at t = ([0-9eE+.\-]+)",
        log_text,
        flags=re.IGNORECASE,
    )

    if matches:

        failure_time = float(
            matches[-1]
        )


# ============================================================
# LOAD OUTPUT
# ============================================================

pvds = sorted(
    RUN.glob("*.pvd")
)

if not pvds:

    raise SystemExit(
        "FAIL: no PVD output"
    )


series = ot.MeshSeries(
    str(pvds[0])
)

times = np.asarray(
    series.timevalues,
    dtype=float,
)


if failure_time is None:

    valid_time_mask = np.ones(
        len(times),
        dtype=bool,
    )

else:

    # STRICTLY before failed timestep.
    valid_time_mask = (
        times
        < failure_time - 1e-10
    )


valid_indices = np.where(
    valid_time_mask
)[0]


if len(valid_indices) == 0:

    raise SystemExit(
        "FAIL: no accepted output "
        "before solver failure"
    )


# ============================================================
# TOPOLOGY
# ============================================================

topology = pv.read(
    MODEL
    / "slope_toe_notch_ready.vtu"
)

material_ids = np.asarray(
    topology.cell_data[
        "MaterialIDs"
    ],
    dtype=int,
)

centers = (
    topology
    .cell_centers()
    .points
)

points = np.asarray(
    topology.points,
    dtype=float,
)


# ============================================================
# EROSION GEOMETRY
# ============================================================

TOE_X = 22.0


candidate = (
    material_ids > 0
)

recession_cell = (
    TOE_X
    - centers[:, 0]
)


def erosion_from_time(
    t,
):

    if t <= 10.0:

        return 0.0

    if t >= 110.0:

        return 2.0

    return (
        2.0
        * (
            t - 10.0
        )
        / 100.0
    )


def active_masks(
    E,
):

    removed = (
        candidate
        & (
            recession_cell
            <= E + 1e-10
        )
    )

    active_cells = ~removed

    active_nodes = np.zeros(
        topology.n_points,
        dtype=bool,
    )

    for cell_id in np.where(
        active_cells
    )[0]:

        cell = topology.get_cell(
            int(cell_id)
        )

        active_nodes[
            cell.point_ids
        ] = True

    return (
        active_cells,
        active_nodes,
    )


# ============================================================
# SLOPE MONITORING ZONE
# ============================================================

def surface_height(
    x,
):

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
        * (
            x[mid] - 8.0
        )
        / 14.0
    )

    h[
        x >= 22.0
    ] = 2.0

    return h


surface = surface_height(
    points[:, 0]
)

depth = (
    surface
    - points[:, 1]
)

slope_zone = (
    (points[:, 0] >= 8.0)
    & (points[:, 0] <= 22.0)
    & (depth >= -1e-8)
    & (depth <= 4.0)
)


# ============================================================
# FIELD HELPERS
# ============================================================

def point_field(
    mesh,
    name,
):

    if name not in mesh.point_data:

        raise KeyError(
            f"{name!r} missing "
            "from point data"
        )

    return np.asarray(
        mesh.point_data[name],
        dtype=float,
    )


def any_field(
    mesh,
    name,
):

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
        f"{name!r} missing"
    )


def latest_valid_index_at_or_before(
    target_time,
):

    candidates = [
        i
        for i in valid_indices
        if times[i]
        <= target_time + 1e-9
    ]

    if not candidates:

        return None

    return int(
        candidates[-1]
    )


# ============================================================
# BASELINE
#
# t=10 s is intact hold state.
# ============================================================

i0 = latest_valid_index_at_or_before(
    10.0
)

if i0 is None:

    raise SystemExit(
        "FAIL: no valid intact "
        "baseline near t=10 s"
    )


mesh0 = series.mesh(
    i0
)

u0 = point_field(
    mesh0,
    "displacement",
)


# ============================================================
# TARGET E VALUES
#
# We sample every 0.2 m.
# If solver failed before a target,
# that target is explicitly INVALID.
# ============================================================

target_E = np.arange(
    0.0,
    2.0001,
    0.2,
)


records = []


for E_target in target_E:

    target_time = (
        10.0
        + 50.0
        * E_target
    )

    idx = latest_valid_index_at_or_before(
        target_time
    )


    if idx is None:

        records.append(
            {
                "E_target_m":
                    E_target,

                "status":
                    "UNAVAILABLE",
            }
        )

        continue


    actual_time = float(
        times[idx]
    )

    actual_E = erosion_from_time(
        actual_time
    )


    # If the accepted state is too far
    # behind the requested erosion target,
    # do not pretend that it represents it.
    if (
        actual_E
        < E_target - 0.051
    ):

        records.append(
            {
                "E_target_m":
                    E_target,

                "time_s":
                    actual_time,

                "E_actual_m":
                    actual_E,

                "status":
                    "NOT_REACHED",
            }
        )

        continue


    mesh = series.mesh(
        idx
    )

    active_cells, active_nodes = (
        active_masks(
            actual_E
        )
    )


    # --------------------------------------------------------
    # DISPLACEMENT RESPONSE
    # --------------------------------------------------------

    u = point_field(
        mesh,
        "displacement",
    )

    du = (
        u
        - u0
    )

    du_mag = np.linalg.norm(
        du,
        axis=1,
    )

    monitor_nodes = (
        active_nodes
        & slope_zone
    )

    if not np.any(
        monitor_nodes
    ):

        raise SystemExit(
            "FAIL: empty slope monitor zone"
        )


    R = float(
        np.nanmax(
            du_mag[
                monitor_nodes
            ]
        )
        * 1000.0
    )


    # --------------------------------------------------------
    # PLASTICITY — ACTIVE DOMAIN ONLY
    # --------------------------------------------------------

    epsp, location = any_field(
        mesh,
        "EquivalentPlasticStrain",
    )

    epsp = np.asarray(
        epsp,
        dtype=float,
    ).squeeze()


    if location == "cell":

        epsp_active = epsp[
            active_cells
        ]

    else:

        epsp_active = epsp[
            active_nodes
        ]


    epsp_active = epsp_active[
        np.isfinite(
            epsp_active
        )
    ]


    epsp_max = float(
        np.nanmax(
            epsp_active
        )
    )

    f6 = float(
        np.mean(
            epsp_active > 1e-6
        )
    )

    f4 = float(
        np.mean(
            epsp_active > 1e-4
        )
    )


    # --------------------------------------------------------
    # HYDRAULICS — ACTIVE NODES ONLY
    # --------------------------------------------------------

    pressure = point_field(
        mesh,
        "pressure",
    ).squeeze()

    p = pressure[
        active_nodes
    ]

    p = p[
        np.isfinite(
            p
        )
    ]


    records.append(
        {
            "E_target_m":
                E_target,

            "time_s":
                actual_time,

            "E_actual_m":
                actual_E,

            "R_slope_mm":
                R,

            "epsp_max":
                epsp_max,

            "frac_epsp_gt_1e6":
                f6,

            "frac_epsp_gt_1e4":
                f4,

            "pressure_min_kpa":
                float(
                    np.min(p)
                    / 1000.0
                ),

            "pressure_max_kpa":
                float(
                    np.max(p)
                    / 1000.0
                ),

            "active_cells":
                int(
                    np.sum(
                        active_cells
                    )
                ),

            "status":
                "VALID",
        }
    )


# ============================================================
# COMPLIANCE
#
# Only between consecutive VALID states.
# ============================================================

previous_valid = None

for r in records:

    r[
        "C_E_mm_per_m"
    ] = np.nan

    r[
        "compliance_ratio"
    ] = np.nan


    if (
        r.get("status")
        != "VALID"
    ):

        continue


    if previous_valid is not None:

        dE = (
            r["E_actual_m"]
            - previous_valid[
                "E_actual_m"
            ]
        )

        dR = (
            r["R_slope_mm"]
            - previous_valid[
                "R_slope_mm"
            ]
        )


        if dE > 1e-10:

            r[
                "C_E_mm_per_m"
            ] = (
                dR / dE
            )


            prev_C = previous_valid.get(
                "C_E_mm_per_m",
                np.nan,
            )

            if (
                np.isfinite(
                    prev_C
                )
                and abs(
                    prev_C
                ) > 1e-12
            ):

                r[
                    "compliance_ratio"
                ] = (
                    r[
                        "C_E_mm_per_m"
                    ]
                    / prev_C
                )


    previous_valid = r


# ============================================================
# PRINT
# ============================================================

print(
    "========================================"
)

print(
    "PHASE 05B-V2 MOVING-FRONT EROSION"
)

print(
    "========================================"
)

print()


if failure_time is None:

    print(
        "Solver failure time: NONE"
    )

else:

    print(
        "Solver failure time [s]: "
        f"{failure_time:.6f}"
    )


print(
    "Last accepted output time [s]: "
    f"{times[valid_indices[-1]]:.6f}"
)

print(
    "Last accepted nominal E [m]: "
    f"{erosion_from_time(times[valid_indices[-1]]):.6f}"
)

print()


for r in records:

    if (
        r.get("status")
        != "VALID"
    ):

        print(
            f"E_target={r['E_target_m']:.1f} m | "
            f"{r.get('status')}"
        )

        continue


    C = r[
        "C_E_mm_per_m"
    ]

    Ctxt = (
        "-"
        if not np.isfinite(C)
        else f"{C:.6f}"
    )


    print(
        f"E={r['E_actual_m']:.3f} m | "
        f"R={r['R_slope_mm']:.6f} mm | "
        f"C_E={Ctxt} mm/m | "
        f"epsp_max={r['epsp_max']:.3e} | "
        f"frac>1e-4="
        f"{100*r['frac_epsp_gt_1e4']:.4f}% | "
        f"active_cells="
        f"{r['active_cells']}"
    )


# ============================================================
# RESPONSE TRANSITION DIAGNOSTICS
# ============================================================

valid_records = [
    r
    for r in records
    if r.get("status") == "VALID"
]


print()
print(
    "=== COMPLIANCE RATIOS ==="
)


for r in valid_records:

    ratio = r[
        "compliance_ratio"
    ]

    if np.isfinite(
        ratio
    ):

        print(
            f"E={r['E_actual_m']:.3f} m | "
            f"C_i/C_prev="
            f"{ratio:.6f}"
        )


# ============================================================
# CSV
# ============================================================

csv_path = (
    OUT
    / "reference_moving_front_response.csv"
)

fieldnames = [
    "E_target_m",
    "time_s",
    "E_actual_m",
    "R_slope_mm",
    "C_E_mm_per_m",
    "compliance_ratio",
    "epsp_max",
    "frac_epsp_gt_1e6",
    "frac_epsp_gt_1e4",
    "pressure_min_kpa",
    "pressure_max_kpa",
    "active_cells",
    "status",
]


with csv_path.open(
    "w",
    newline="",
    encoding="utf-8",
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=fieldnames,
    )

    writer.writeheader()

    for r in records:

        writer.writerow(
            {
                key:
                    r.get(
                        key,
                        ""
                    )
                for key in fieldnames
            }
        )


# ============================================================
# FIGURES
# ============================================================

if len(valid_records) >= 2:

    E = np.array(
        [
            r["E_actual_m"]
            for r in valid_records
        ]
    )

    R = np.array(
        [
            r["R_slope_mm"]
            for r in valid_records
        ]
    )


    fig, ax = plt.subplots(
        figsize=(7.2, 4.8)
    )

    ax.plot(
        E,
        R,
        marker="o",
    )

    ax.set_xlabel(
        "Toe recession E [m]"
    )

    ax.set_ylabel(
        "Slope response R(E) [mm]"
    )

    ax.set_title(
        "Moving-front toe-erosion response"
    )

    ax.grid(
        alpha=0.25
    )

    fig.tight_layout()

    fig.savefig(
        OUT
        / "figure_01_moving_front_response.png",
        dpi=220,
    )

    plt.close(
        fig
    )


    valid_C = [
        r
        for r in valid_records
        if np.isfinite(
            r["C_E_mm_per_m"]
        )
    ]


    if valid_C:

        fig, ax = plt.subplots(
            figsize=(7.2, 4.8)
        )

        ax.plot(
            [
                r["E_actual_m"]
                for r in valid_C
            ],
            [
                r["C_E_mm_per_m"]
                for r in valid_C
            ],
            marker="o",
        )

        ax.set_xlabel(
            "Toe recession E [m]"
        )

        ax.set_ylabel(
            "Erosion compliance C_E [mm/m]"
        )

        ax.set_title(
            "Incremental erosion compliance"
        )

        ax.grid(
            alpha=0.25
        )

        fig.tight_layout()

        fig.savefig(
            OUT
            / "figure_02_moving_front_compliance.png",
            dpi=220,
        )

        plt.close(
            fig
        )


# ============================================================
# SUMMARY
# ============================================================

summary = (
    OUT
    / "phase05b_v2_summary.txt"
)


lines = [
    "PHASE 05B-V2 MOVING-FRONT TOE EROSION",
    "",
]


if failure_time is None:

    lines.append(
        "Solver failure time: NONE"
    )

else:

    lines.append(
        "Solver failure time [s]: "
        f"{failure_time:.8f}"
    )


lines += [
    (
        "Last accepted output time [s]: "
        f"{times[valid_indices[-1]]:.8f}"
    ),
    (
        "Last accepted nominal E [m]: "
        f"{erosion_from_time(times[valid_indices[-1]]):.8f}"
    ),
    "",
]


for r in records:

    if (
        r.get("status")
        != "VALID"
    ):

        lines.append(
            f"E_target={r['E_target_m']:.1f} m | "
            f"{r.get('status')}"
        )

        continue


    C = (
        "-"
        if not np.isfinite(
            r["C_E_mm_per_m"]
        )
        else
        f"{r['C_E_mm_per_m']:.8f}"
    )


    lines.append(
        f"E={r['E_actual_m']:.3f} m | "
        f"R={r['R_slope_mm']:.8f} mm | "
        f"C_E={C} mm/m | "
        f"epsp_max="
        f"{r['epsp_max']:.8e} | "
        f"frac_epsp_gt_1e4="
        f"{100*r['frac_epsp_gt_1e4']:.6f}% | "
        f"active_cells="
        f"{r['active_cells']}"
    )


lines += [
    "",
    (
        "IMPORTANT: failed timesteps and "
        "unreached erosion targets are "
        "excluded from interpretation."
    ),
    (
        "No Ecrit is assigned automatically."
    ),
]


summary.write_text(
    "\n".join(lines) + "\n",
    encoding="utf-8",
)


print()
print(
    "PASS:",
    csv_path,
)

print(
    "PASS:",
    summary,
)

print()
print(
    "PHASE 05B-V2 ANALYSIS COMPLETE"
)
