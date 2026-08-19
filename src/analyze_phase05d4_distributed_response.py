from pathlib import Path
import csv
import re

import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot
import pyvista as pv


ROOT = Path.cwd()

MODEL_ROOT = (
    ROOT
    / "model"
    / "phase05d2_mesh_erosion"
)

RUN_ROOT = (
    ROOT
    / "results"
    / "phase05d2_mesh_erosion"
)

OUT = (
    ROOT
    / "results"
    / "phase05d4_distributed_response"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)

CASES = [
    "coarse",
    "medium",
    "fine",
]

TOE_X = 22.0


# ============================================================
# NUMERICAL EROSION LAW
# ============================================================

def erosion_from_time(t):

    if t <= 20.0:
        return 0.0

    if t >= 100.0:
        return 0.8

    return (
        0.01
        * (
            t - 20.0
        )
    )


# ============================================================
# FAILURE TIME
# ============================================================

def failure_time_from_log(path):

    if not path.exists():
        return None

    text = path.read_text(
        encoding="utf-8",
        errors="ignore",
    )

    matches = re.findall(
        r"nonlinear solver failed .*?"
        r"at t = ([0-9eE+.\-]+)",
        text,
        flags=re.IGNORECASE,
    )

    if not matches:
        return None

    return float(
        matches[-1]
    )


# ============================================================
# SLOPE SURFACE
# ============================================================

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
        * (
            x[mid] - 8.0
        )
        / 14.0
    )

    h[
        x >= 22.0
    ] = 2.0

    return h


# ============================================================
# FIELD HELPERS
# ============================================================

def point_field(mesh, name):

    if name not in mesh.point_data:

        raise KeyError(
            f"{name!r} missing "
            "from point data"
        )

    return np.asarray(
        mesh.point_data[name],
        dtype=float,
    )


def epsp_to_cells(
    output_mesh,
    bulk,
):

    name = (
        "EquivalentPlasticStrain"
    )

    if name in output_mesh.cell_data:

        return np.asarray(
            output_mesh.cell_data[name],
            dtype=float,
        ).squeeze()


    if name in output_mesh.point_data:

        p = np.asarray(
            output_mesh.point_data[name],
            dtype=float,
        ).squeeze()

        values = np.empty(
            bulk.n_cells,
            dtype=float,
        )

        for cid in range(
            bulk.n_cells
        ):

            ids = (
                bulk
                .get_cell(cid)
                .point_ids
            )

            values[cid] = np.nanmax(
                p[ids]
            )

        return values


    raise KeyError(
        "EquivalentPlasticStrain missing"
    )


# ============================================================
# ANALYZE CASE
# ============================================================

def analyze_case(case):

    model_dir = (
        MODEL_ROOT
        / case
    )

    run_dir = (
        RUN_ROOT
        / case
    )

    bulk = pv.read(
        model_dir
        / "bulk.vtu"
    )

    material_ids = np.asarray(
        bulk.cell_data[
            "MaterialIDs"
        ],
        dtype=int,
    )

    centers = (
        bulk
        .cell_centers()
        .points
    )

    points = np.asarray(
        bulk.points,
        dtype=float,
    )

    sizes = bulk.compute_cell_sizes(
        length=False,
        area=True,
        volume=False,
    )

    cell_area = np.asarray(
        sizes.cell_data[
            "Area"
        ],
        dtype=float,
    )

    candidate = (
        material_ids > 0
    )

    recession = (
        TOE_X
        - centers[:, 0]
    )


    # --------------------------------------------------------
    # Fixed BODY monitoring zone.
    #
    # Intentionally excludes the coastal notch tip.
    #
    # x <= 20.5 means:
    # even at E = 0.8 m,
    # the front is x = 21.2 m,
    # leaving >=0.7 m separation.
    # --------------------------------------------------------

    top_point = surface_height(
        points[:, 0]
    )

    point_depth = (
        top_point
        - points[:, 1]
    )

    body_nodes = (
        (points[:, 0] >= 8.0)
        & (points[:, 0] <= 20.5)
        & (point_depth >= -1e-8)
        & (point_depth <= 4.0)
    )


    top_cell = surface_height(
        centers[:, 0]
    )

    cell_depth = (
        top_cell
        - centers[:, 1]
    )

    body_cells = (
        (centers[:, 0] >= 8.0)
        & (centers[:, 0] <= 20.5)
        & (cell_depth >= -1e-8)
        & (cell_depth <= 4.0)
    )


    if not np.any(
        body_nodes
    ):
        raise RuntimeError(
            f"Empty body node zone: {case}"
        )

    if not np.any(
        body_cells
    ):
        raise RuntimeError(
            f"Empty body cell zone: {case}"
        )


    # --------------------------------------------------------
    # Load accepted output.
    # --------------------------------------------------------

    fail_t = failure_time_from_log(
        run_dir
        / "ogs.log"
    )

    pvds = sorted(
        run_dir.glob(
            "*.pvd"
        )
    )

    if not pvds:
        raise RuntimeError(
            f"No PVD for {case}"
        )

    series = ot.MeshSeries(
        str(
            pvds[0]
        )
    )

    times = np.asarray(
        series.timevalues,
        dtype=float,
    )

    if fail_t is None:

        valid_indices = np.arange(
            len(times)
        )

    else:

        valid_indices = np.where(
            times
            < fail_t
            - 1e-10
        )[0]


    # --------------------------------------------------------
    # Baseline at final intact hold state.
    # --------------------------------------------------------

    baseline_ids = [
        i
        for i in valid_indices
        if times[i]
        <= 20.0 + 1e-9
    ]

    if not baseline_ids:
        raise RuntimeError(
            f"No baseline for {case}"
        )

    i0 = int(
        baseline_ids[-1]
    )

    baseline = series.mesh(
        i0
    )

    u0 = point_field(
        baseline,
        "displacement",
    )


    records = []

    previous_removed_area = None


    for idx in valid_indices:

        t = float(
            times[idx]
        )

        if t < 20.0 - 1e-9:
            continue

        E = erosion_from_time(
            t
        )

        removed = (
            candidate
            & (
                recession
                <= E + 1e-10
            )
        )

        active_cells = ~removed


        # ----------------------------------------------------
        # Actual geometric perturbation.
        # ----------------------------------------------------

        removed_area = float(
            np.sum(
                cell_area[
                    removed
                ]
            )
        )

        if previous_removed_area is None:

            topology_event = 0

        else:

            topology_event = int(
                removed_area
                > previous_removed_area
                + 1e-12
            )


        output = series.mesh(
            int(idx)
        )


        # ----------------------------------------------------
        # DISTRIBUTED DISPLACEMENT.
        # ----------------------------------------------------

        u = point_field(
            output,
            "displacement",
        )

        du = (
            u - u0
        )

        du_mag = np.linalg.norm(
            du,
            axis=1,
        )

        body_values = du_mag[
            body_nodes
        ]

        body_values = body_values[
            np.isfinite(
                body_values
            )
        ]


        R_rms = float(
            np.sqrt(
                np.mean(
                    body_values ** 2
                )
            )
            * 1000.0
        )

        R_mean = float(
            np.mean(
                body_values
            )
            * 1000.0
        )

        R95 = float(
            np.percentile(
                body_values,
                95.0,
            )
            * 1000.0
        )

        Rmax_body = float(
            np.max(
                body_values
            )
            * 1000.0
        )


        # ----------------------------------------------------
        # BODY PLASTICITY ONLY.
        #
        # Notch-front cells are excluded geometrically.
        # ----------------------------------------------------

        epsp = epsp_to_cells(
            output,
            bulk,
        )

        valid_body = (
            body_cells
            & active_cells
            & np.isfinite(
                epsp
            )
        )

        body_epsp = epsp[
            valid_body
        ]

        body_area = cell_area[
            valid_body
        ]


        epsp_body_max = float(
            np.max(
                body_epsp
            )
        )

        plastic_area_1e6 = float(
            np.sum(
                body_area[
                    body_epsp > 1e-6
                ]
            )
        )

        plastic_area_1e4 = float(
            np.sum(
                body_area[
                    body_epsp > 1e-4
                ]
            )
        )

        plastic_integral = float(
            np.sum(
                body_epsp
                * body_area
            )
        )


        records.append(
            {
                "case":
                    case,

                "time_s":
                    t,

                "E_nominal_m":
                    E,

                "removed_area_m2":
                    removed_area,

                "topology_event":
                    topology_event,

                "R_body_rms_mm":
                    R_rms,

                "R_body_mean_mm":
                    R_mean,

                "R_body_p95_mm":
                    R95,

                "R_body_max_mm":
                    Rmax_body,

                "epsp_body_max":
                    epsp_body_max,

                "body_plastic_area_gt_1e6_m2":
                    plastic_area_1e6,

                "body_plastic_area_gt_1e4_m2":
                    plastic_area_1e4,

                "body_integral_epsp_dA":
                    plastic_integral,
            }
        )

        previous_removed_area = (
            removed_area
        )


    return {
        "case":
            case,

        "failure_E_m":
            (
                erosion_from_time(
                    fail_t
                )
                if fail_t
                is not None
                else np.nan
            ),

        "records":
            records,
    }


# ============================================================
# ANALYZE
# ============================================================

results = [
    analyze_case(
        case
    )
    for case in CASES
]


all_records = []

for result in results:

    all_records.extend(
        result[
            "records"
        ]
    )


# ============================================================
# FULL CSV
# ============================================================

history_csv = (
    OUT
    / "distributed_response_history.csv"
)

fieldnames = list(
    all_records[0].keys()
)

with history_csv.open(
    "w",
    newline="",
    encoding="utf-8",
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=fieldnames,
    )

    writer.writeheader()
    writer.writerows(
        all_records
    )


# ============================================================
# COMMON VALID E WINDOW
#
# Fine fails first.
# We compare meshes only where all three
# have accepted output.
# ============================================================

max_common_E = min(
    result[
        "records"
    ][-1][
        "E_nominal_m"
    ]
    for result in results
)


targets = np.arange(
    0.0,
    max_common_E
    + 1e-9,
    0.01,
)


def record_at_E(
    result,
    E,
):

    records = result[
        "records"
    ]

    values = np.asarray(
        [
            r[
                "E_nominal_m"
            ]
            for r in records
        ],
        dtype=float,
    )

    idx = int(
        np.argmin(
            np.abs(
                values - E
            )
        )
    )

    if abs(
        values[idx] - E
    ) > 0.0051:

        return None

    return records[idx]


common_rows = []


for E in targets:

    row = {
        "E_nominal_m":
            E,
    }

    ok = True

    for result in results:

        r = record_at_E(
            result,
            E,
        )

        if r is None:
            ok = False
            break

        case = result[
            "case"
        ]

        for key in [
            "removed_area_m2",
            "R_body_rms_mm",
            "R_body_p95_mm",
            "epsp_body_max",
            "body_plastic_area_gt_1e4_m2",
            "body_integral_epsp_dA",
        ]:

            row[
                f"{case}_{key}"
            ] = r[key]


    if not ok:
        continue


    # --------------------------------------------------------
    # Medium / fine relative differences.
    #
    # For near-zero responses, relative difference is
    # intentionally left NaN.
    # --------------------------------------------------------

    for metric in [
        "R_body_rms_mm",
        "R_body_p95_mm",
    ]:

        m = row[
            f"medium_{metric}"
        ]

        f = row[
            f"fine_{metric}"
        ]

        absolute = abs(
            m - f
        )

        scale = max(
            abs(m),
            abs(f),
        )

        relative = (
            absolute / scale
            if scale > 1e-8
            else np.nan
        )

        row[
            f"medium_fine_absdiff_{metric}"
        ] = absolute

        row[
            f"medium_fine_reldiff_{metric}"
        ] = relative


    common_rows.append(
        row
    )


common_csv = (
    OUT
    / "common_window_mesh_comparison.csv"
)

if common_rows:

    with common_csv.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=list(
                common_rows[0].keys()
            ),
        )

        writer.writeheader()
        writer.writerows(
            common_rows
        )


# ============================================================
# PRINT CASE HISTORY AT TOPOLOGY EVENTS
# ============================================================

print(
    "========================================"
)

print(
    "PHASE 05D-4 DISTRIBUTED RESPONSE AUDIT"
)

print(
    "========================================"
)

print(
    f"Common accepted E range: "
    f"0 to {max_common_E:.3f} m"
)


for result in results:

    print()
    print(
        f"=== {result['case'].upper()} ==="
    )

    print(
        "Numerical nonconvergence E [m]: "
        f"{result['failure_E_m']}"
    )

    print()
    print(
        "TOPOLOGY STATES"
    )


    topology_records = [
        r
        for r in result[
            "records"
        ]
        if (
            r[
                "topology_event"
            ] == 1
        )
    ]


    for r in topology_records:

        print(
            f"E={r['E_nominal_m']:.3f} | "
            f"Arem={r['removed_area_m2']:.6f} m2 | "
            f"Rrms={r['R_body_rms_mm']:.8f} mm | "
            f"R95={r['R_body_p95_mm']:.8f} mm | "
            f"epsp_body_max="
            f"{r['epsp_body_max']:.3e} | "
            f"Aplastic_body="
            f"{r['body_plastic_area_gt_1e4_m2']:.8f} m2"
        )


# ============================================================
# MEDIUM/FINE COMMON-WINDOW ROBUSTNESS
# ============================================================

valid_rel_rms = np.asarray(
    [
        r[
            "medium_fine_reldiff_R_body_rms_mm"
        ]
        for r in common_rows
        if np.isfinite(
            r[
                "medium_fine_reldiff_R_body_rms_mm"
            ]
        )
    ],
    dtype=float,
)


valid_rel_p95 = np.asarray(
    [
        r[
            "medium_fine_reldiff_R_body_p95_mm"
        ]
        for r in common_rows
        if np.isfinite(
            r[
                "medium_fine_reldiff_R_body_p95_mm"
            ]
        )
    ],
    dtype=float,
)


print()
print(
    "========================================"
)

print(
    "MEDIUM / FINE DISTRIBUTED RESPONSE"
)

print(
    "========================================"
)


if len(valid_rel_rms):

    print(
        "R_body RMS relative difference:"
    )

    print(
        "  median = "
        f"{100*np.median(valid_rel_rms):.4f}%"
    )

    print(
        "  max    = "
        f"{100*np.max(valid_rel_rms):.4f}%"
    )

else:

    print(
        "R_body RMS relative difference:"
        " no non-zero common states"
    )


if len(valid_rel_p95):

    print(
        "R_body P95 relative difference:"
    )

    print(
        "  median = "
        f"{100*np.median(valid_rel_p95):.4f}%"
    )

    print(
        "  max    = "
        f"{100*np.max(valid_rel_p95):.4f}%"
    )

else:

    print(
        "R_body P95 relative difference:"
        " no non-zero common states"
    )


# ============================================================
# BODY PLASTICITY CHECK
# ============================================================

print()
print(
    "=== FAR-FIELD/BODY PLASTICITY ==="
)


for result in results:

    max_epsp_body = max(
        r[
            "epsp_body_max"
        ]
        for r in result[
            "records"
        ]
    )

    max_area_body = max(
        r[
            "body_plastic_area_gt_1e4_m2"
        ]
        for r in result[
            "records"
        ]
    )

    print(
        f"{result['case']:6s} | "
        f"max epsp_body="
        f"{max_epsp_body:.3e} | "
        f"max A(epsp>1e-4)="
        f"{max_area_body:.8f} m2"
    )


# ============================================================
# DECISION
# ============================================================

median_rms = (
    float(
        np.median(
            valid_rel_rms
        )
    )
    if len(
        valid_rel_rms
    )
    else np.nan
)

median_p95 = (
    float(
        np.median(
            valid_rel_p95
        )
    )
    if len(
        valid_rel_p95
    )
    else np.nan
)


print()
print(
    "========================================"
)

print(
    "DISTRIBUTED-METRIC DECISION"
)

print(
    "========================================"
)


if (
    np.isfinite(
        median_rms
    )
    and np.isfinite(
        median_p95
    )
    and median_rms <= 0.15
    and median_p95 <= 0.15
):

    decision = (
        "DISTRIBUTED BODY RESPONSE "
        "IS PROMISING"
    )

    print(
        decision
    )

    print(
        "Medium/fine distributed displacement "
        "is substantially more consistent than "
        "notch-tip plasticity/nonconvergence."
    )

else:

    decision = (
        "DISTRIBUTED BODY RESPONSE "
        "STILL NEEDS REVIEW"
    )

    print(
        decision
    )

    print(
        "Do not define a transition metric yet."
    )


# ============================================================
# FIGURE 1 — RMS vs nominal E
# ============================================================

fig, ax = plt.subplots(
    figsize=(7.5, 4.8)
)


for result in results:

    records = result[
        "records"
    ]

    ax.plot(
        [
            r[
                "E_nominal_m"
            ]
            for r in records
        ],
        [
            r[
                "R_body_rms_mm"
            ]
            for r in records
        ],
        marker="o",
        markersize=3,
        label=result[
            "case"
        ],
    )


ax.set_xlabel(
    "Nominal toe recession E [m]"
)

ax.set_ylabel(
    "Body RMS erosion-induced displacement [mm]"
)

ax.set_title(
    "Distributed slope-body response"
)

ax.legend()

ax.grid(
    alpha=0.25
)

fig.tight_layout()

fig.savefig(
    OUT
    / "figure_01_body_rms_vs_E.png",
    dpi=220,
)

plt.close(
    fig
)


# ============================================================
# FIGURE 2 — RMS vs actual removed area
# ============================================================

fig, ax = plt.subplots(
    figsize=(7.5, 4.8)
)


for result in results:

    records = result[
        "records"
    ]

    ax.plot(
        [
            r[
                "removed_area_m2"
            ]
            for r in records
        ],
        [
            r[
                "R_body_rms_mm"
            ]
            for r in records
        ],
        marker="o",
        markersize=3,
        label=result[
            "case"
        ],
    )


ax.set_xlabel(
    "Actual removed notch area [m²]"
)

ax.set_ylabel(
    "Body RMS erosion-induced displacement [mm]"
)

ax.set_title(
    "Distributed response versus actual erosion geometry"
)

ax.legend()

ax.grid(
    alpha=0.25
)

fig.tight_layout()

fig.savefig(
    OUT
    / "figure_02_body_rms_vs_removed_area.png",
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
    / "phase05d4_summary.txt"
)


lines = [
    "PHASE 05D-4 DISTRIBUTED RESPONSE AUDIT",
    "",
    (
        "Body monitoring zone: "
        "x=8.0 to 20.5 m and "
        "0 to 4 m below local slope surface."
    ),
    (
        "The coastal notch-tip region is "
        "excluded from the body metric."
    ),
    "",
    (
        "Common accepted nominal E range [m]: "
        f"0 to {max_common_E:.3f}"
    ),
    "",
    (
        "Medium/fine median relative "
        "difference R_body_RMS: "
        f"{median_rms}"
    ),
    (
        "Medium/fine median relative "
        "difference R_body_P95: "
        f"{median_p95}"
    ),
    "",
    (
        "DECISION: "
        f"{decision}"
    ),
    "",
    (
        "Notch-tip epsp_max and solver "
        "nonconvergence remain excluded "
        "from any Ecrit definition."
    ),
    (
        "Actual removed area is retained as "
        "a geometry diagnostic because "
        "element deactivation produces "
        "mesh-dependent staircase recession."
    ),
]


summary.write_text(
    "\n".join(
        lines
    )
    + "\n",
    encoding="utf-8",
)


print()
print(
    "PASS:",
    history_csv,
)

print(
    "PASS:",
    common_csv,
)

print(
    "PASS:",
    OUT
    / "figure_01_body_rms_vs_E.png"
)

print(
    "PASS:",
    OUT
    / "figure_02_body_rms_vs_removed_area.png"
)

print(
    "PASS:",
    summary,
)

print()
print(
    "PHASE 05D-4 ANALYSIS COMPLETE"
)
