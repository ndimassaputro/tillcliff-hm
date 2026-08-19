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
    / "phase05e3_antecedent_erosion"
)

RUN_ROOT = (
    ROOT
    / "results"
    / "phase05e3_antecedent_erosion"
)

OUT = (
    ROOT
    / "results"
    / "phase05e3_analysis"
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

TOE_X = 22.0

HOLD_END = 20.0

EROSION_RATE = 0.01

E_MAX = 0.55


def erosion_from_time(t):

    if t <= HOLD_END:

        return 0.0

    E = (
        EROSION_RATE
        * (
            t - HOLD_END
        )
    )

    return float(
        np.clip(
            E,
            0.0,
            E_MAX,
        )
    )


def failure_time(path):

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


def point_field(
    mesh,
    name,
):

    if name not in mesh.point_data:

        raise KeyError(
            f"{name} missing "
            "from point data"
        )

    return np.asarray(
        mesh.point_data[name],
        dtype=float,
    )


def epsp_cells(
    mesh,
    bulk,
):

    name = (
        "EquivalentPlasticStrain"
    )

    if name in mesh.cell_data:

        return np.asarray(
            mesh.cell_data[name],
            dtype=float,
        ).squeeze()


    if name in mesh.point_data:

        values = np.asarray(
            mesh.point_data[name],
            dtype=float,
        ).squeeze()

        cell_values = np.empty(
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

            cell_values[cid] = np.nanmax(
                values[ids]
            )

        return cell_values


    raise KeyError(
        "EquivalentPlasticStrain missing"
    )


def analyze_state(state):

    model_dir = (
        MODEL_ROOT
        / state
    )

    run_dir = (
        RUN_ROOT
        / state
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

    points = np.asarray(
        bulk.points,
        dtype=float,
    )

    centers = (
        bulk
        .cell_centers()
        .points
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


    # ========================================================
    # FIXED DISTRIBUTED SLOPE-BODY ZONE
    #
    # Excludes notch-front localization.
    # ========================================================

    top_nodes = surface_height(
        points[:, 0]
    )

    depth_nodes = (
        top_nodes
        - points[:, 1]
    )

    body_nodes = (
        (points[:, 0] >= 8.0)
        & (points[:, 0] <= 20.5)
        & (depth_nodes >= -1e-8)
        & (depth_nodes <= 4.0)
    )


    top_cells = surface_height(
        centers[:, 0]
    )

    depth_cells = (
        top_cells
        - centers[:, 1]
    )

    body_cells = (
        (centers[:, 0] >= 8.0)
        & (centers[:, 0] <= 20.5)
        & (depth_cells >= -1e-8)
        & (depth_cells <= 4.0)
    )


    fail_t = failure_time(
        run_dir
        / "ogs.log"
    )

    fail_E = (
        erosion_from_time(
            fail_t
        )
        if fail_t is not None
        else np.nan
    )


    pvds = sorted(
        run_dir.glob(
            "*.pvd"
        )
    )

    if not pvds:

        raise RuntimeError(
            f"No PVD for {state}"
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


    # ========================================================
    # BASELINE = t=20 s AFTER INTACT HOLD
    # ========================================================

    baseline_ids = [
        i
        for i in valid_indices
        if times[i]
        <= HOLD_END + 1e-9
    ]

    if not baseline_ids:

        raise RuntimeError(
            f"No baseline for {state}"
        )


    i0 = int(
        baseline_ids[-1]
    )

    base = series.mesh(
        i0
    )

    u0 = point_field(
        base,
        "displacement",
    )


    records = []

    previous_area = None


    for idx in valid_indices:

        t = float(
            times[idx]
        )

        if t < HOLD_END - 1e-9:
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


        removed_area = float(
            np.sum(
                cell_area[
                    removed
                ]
            )
        )


        if previous_area is None:

            topology_event = 0

        else:

            topology_event = int(
                removed_area
                > previous_area
                + 1e-12
            )


        mesh = series.mesh(
            int(idx)
        )


        # ----------------------------------------------------
        # DISTRIBUTED BODY DISPLACEMENT
        # ----------------------------------------------------

        u = point_field(
            mesh,
            "displacement",
        )

        du = (
            u - u0
        )

        du_mag = np.linalg.norm(
            du,
            axis=1,
        )

        body_u = du_mag[
            body_nodes
        ]

        body_u = body_u[
            np.isfinite(
                body_u
            )
        ]


        R_rms = float(
            np.sqrt(
                np.mean(
                    body_u ** 2
                )
            )
            * 1000.0
        )

        R_mean = float(
            np.mean(
                body_u
            )
            * 1000.0
        )

        R95 = float(
            np.percentile(
                body_u,
                95.0,
            )
            * 1000.0
        )

        Rmax = float(
            np.max(
                body_u
            )
            * 1000.0
        )


        # ----------------------------------------------------
        # HYDRAULIC STATE
        # ----------------------------------------------------

        pressure = point_field(
            mesh,
            "pressure",
        ).squeeze()

        saturation = point_field(
            mesh,
            "saturation",
        ).squeeze()


        p_body = pressure[
            body_nodes
        ]

        sr_body = saturation[
            body_nodes
        ]


        # ----------------------------------------------------
        # FAR-FIELD PLASTICITY
        # ----------------------------------------------------

        epsp = epsp_cells(
            mesh,
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

        area_body_1e4 = float(
            np.sum(
                body_area[
                    body_epsp > 1e-4
                ]
            )
        )


        records.append(
            {
                "state":
                    state,

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
                    Rmax,

                "body_pressure_mean_kpa":
                    float(
                        np.mean(
                            p_body
                        )
                        / 1000.0
                    ),

                "body_saturation_mean":
                    float(
                        np.mean(
                            sr_body
                        )
                    ),

                "epsp_body_max":
                    epsp_body_max,

                "body_plastic_area_gt_1e4_m2":
                    area_body_1e4,
            }
        )


        previous_area = (
            removed_area
        )


    return {
        "state":
            state,

        "failure_E_m":
            fail_E,

        "records":
            records,
    }


# ============================================================
# ANALYZE ALL STATES
# ============================================================

results = [
    analyze_state(
        state
    )
    for state in STATES
]


all_records = []

for result in results:

    all_records.extend(
        result[
            "records"
        ]
    )


history_csv = (
    OUT
    / "antecedent_erosion_history.csv"
)


with history_csv.open(
    "w",
    newline="",
    encoding="utf-8",
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=list(
            all_records[0].keys()
        ),
    )

    writer.writeheader()
    writer.writerows(
        all_records
    )


# ============================================================
# COMMON ACCEPTED E RANGE
# ============================================================

common_max_E = min(
    result[
        "records"
    ][-1][
        "E_nominal_m"
    ]
    for result in results
)


def nearest_record(
    result,
    E,
):

    records = result[
        "records"
    ]

    e = np.asarray(
        [
            r["E_nominal_m"]
            for r in records
        ],
        dtype=float,
    )

    idx = int(
        np.argmin(
            np.abs(
                e - E
            )
        )
    )

    if abs(
        e[idx] - E
    ) > 0.0051:

        return None

    return records[idx]


# ============================================================
# TOPOLOGY / GEOMETRY STATES
#
# Medium mesh has actual events approximately:
# 0.05, 0.15, 0.25, 0.35, 0.45, 0.55
# ============================================================

candidate_E = [
    0.05,
    0.15,
    0.25,
    0.35,
    0.45,
    0.55,
]


comparison = []


for target_E in candidate_E:

    if target_E > common_max_E + 1e-8:
        continue


    state_records = {}

    valid = True


    for result in results:

        record = nearest_record(
            result,
            target_E,
        )

        if record is None:

            valid = False
            break

        state_records[
            result["state"]
        ] = record


    if not valid:
        continue


    dry = state_records[
        "dry"
    ]

    ref = state_records[
        "reference"
    ]

    wet = state_records[
        "wet"
    ]


    # Geometry must be identical between states.
    areas = [
        dry[
            "removed_area_m2"
        ],
        ref[
            "removed_area_m2"
        ],
        wet[
            "removed_area_m2"
        ],
    ]


    geometry_match = (
        max(areas)
        - min(areas)
        <= 1e-12
    )


    row = {
        "E_nominal_m":
            target_E,

        "removed_area_m2":
            ref[
                "removed_area_m2"
            ],

        "geometry_match":
            int(
                geometry_match
            ),
    }


    for state, r in (
        state_records.items()
    ):

        for metric in [
            "R_body_rms_mm",
            "R_body_p95_mm",
            "body_pressure_mean_kpa",
            "body_saturation_mean",
            "epsp_body_max",
            "body_plastic_area_gt_1e4_m2",
        ]:

            row[
                f"{state}_{metric}"
            ] = r[
                metric
            ]


    ref_rms = ref[
        "R_body_rms_mm"
    ]

    ref_p95 = ref[
        "R_body_p95_mm"
    ]


    row[
        "wet_minus_ref_Rrms_mm"
    ] = (
        wet[
            "R_body_rms_mm"
        ]
        - ref_rms
    )

    row[
        "dry_minus_ref_Rrms_mm"
    ] = (
        dry[
            "R_body_rms_mm"
        ]
        - ref_rms
    )


    row[
        "wet_over_ref_Rrms"
    ] = (
        wet[
            "R_body_rms_mm"
        ]
        / ref_rms
        if abs(
            ref_rms
        ) > 1e-10
        else np.nan
    )

    row[
        "dry_over_ref_Rrms"
    ] = (
        dry[
            "R_body_rms_mm"
        ]
        / ref_rms
        if abs(
            ref_rms
        ) > 1e-10
        else np.nan
    )


    row[
        "wet_minus_ref_Rp95_mm"
    ] = (
        wet[
            "R_body_p95_mm"
        ]
        - ref_p95
    )

    row[
        "dry_minus_ref_Rp95_mm"
    ] = (
        dry[
            "R_body_p95_mm"
        ]
        - ref_p95
    )


    comparison.append(
        row
    )


comparison_csv = (
    OUT
    / "antecedent_state_comparison.csv"
)


if comparison:

    with comparison_csv.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=list(
                comparison[0].keys()
            ),
        )

        writer.writeheader()
        writer.writerows(
            comparison
        )


# ============================================================
# PRINT
# ============================================================

print(
    "========================================"
)

print(
    "PHASE 05E-3 ANTECEDENT × EROSION"
)

print(
    "========================================"
)

print(
    "Common accepted E range [m]: "
    f"0 to {common_max_E:.3f}"
)


for result in results:

    print()
    print(
        f"=== {result['state'].upper()} ==="
    )

    print(
        "Numerical nonconvergence E [m]: "
        f"{result['failure_E_m']}"
    )


    event_records = [
        r
        for r in result[
            "records"
        ]
        if r[
            "topology_event"
        ] == 1
    ]


    print(
        "TOPOLOGY STATES"
    )


    for r in event_records:

        print(
            f"E={r['E_nominal_m']:.3f} | "
            f"Arem="
            f"{r['removed_area_m2']:.6f} m2 | "
            f"Rrms="
            f"{r['R_body_rms_mm']:.8f} mm | "
            f"R95="
            f"{r['R_body_p95_mm']:.8f} mm | "
            f"pbody="
            f"{r['body_pressure_mean_kpa']:.3f} kPa | "
            f"Srbody="
            f"{r['body_saturation_mean']:.6f} | "
            f"epsp_body="
            f"{r['epsp_body_max']:.3e}"
        )


print()
print(
    "========================================"
)

print(
    "DRY / REFERENCE / WET COMPARISON"
)

print(
    "========================================"
)


for r in comparison:

    print(
        f"E={r['E_nominal_m']:.2f} | "
        f"Arem={r['removed_area_m2']:.6f} | "
        f"Rrms dry/ref/wet="
        f"{r['dry_R_body_rms_mm']:.8f}/"
        f"{r['reference_R_body_rms_mm']:.8f}/"
        f"{r['wet_R_body_rms_mm']:.8f} mm | "
        f"wet/ref="
        f"{r['wet_over_ref_Rrms']:.6f} | "
        f"dry/ref="
        f"{r['dry_over_ref_Rrms']:.6f}"
    )


# ============================================================
# FAR-FIELD PLASTICITY
# ============================================================

print()
print(
    "========================================"
)

print(
    "BODY PLASTICITY CHECK"
)

print(
    "========================================"
)


body_plasticity_zero = True


for result in results:

    max_epsp = max(
        r[
            "epsp_body_max"
        ]
        for r in result[
            "records"
        ]
    )

    max_area = max(
        r[
            "body_plastic_area_gt_1e4_m2"
        ]
        for r in result[
            "records"
        ]
    )


    print(
        f"{result['state']:9s} | "
        f"max epsp_body="
        f"{max_epsp:.3e} | "
        f"max plastic area="
        f"{max_area:.8f} m2"
    )


    if (
        max_epsp > 1e-8
        or max_area > 0
    ):

        body_plasticity_zero = False


# ============================================================
# ORDERING DIAGNOSTIC
#
# No assumption about which response SHOULD be largest.
# ============================================================

wet_minus_ref = np.asarray(
    [
        r[
            "wet_minus_ref_Rrms_mm"
        ]
        for r in comparison
        if np.isfinite(
            r[
                "wet_minus_ref_Rrms_mm"
            ]
        )
    ],
    dtype=float,
)


dry_minus_ref = np.asarray(
    [
        r[
            "dry_minus_ref_Rrms_mm"
        ]
        for r in comparison
        if np.isfinite(
            r[
                "dry_minus_ref_Rrms_mm"
            ]
        )
    ],
    dtype=float,
)


if len(
    wet_minus_ref
):

    wet_sign_fraction = float(
        np.mean(
            wet_minus_ref > 0
        )
    )

else:

    wet_sign_fraction = np.nan


if len(
    dry_minus_ref
):

    dry_sign_fraction = float(
        np.mean(
            dry_minus_ref > 0
        )
    )

else:

    dry_sign_fraction = np.nan


print()
print(
    "========================================"
)

print(
    "ANTECEDENT-STATE RESPONSE DIAGNOSTIC"
)

print(
    "========================================"
)

print(
    "Fraction of comparison states with "
    "Rwet > Rreference: "
    f"{wet_sign_fraction}"
)

print(
    "Fraction of comparison states with "
    "Rdry > Rreference: "
    f"{dry_sign_fraction}"
)

print(
    "No physical ordering is imposed "
    "a priori."
)


# ============================================================
# FIGURE 1
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
            "state"
        ],
    )


ax.set_xlabel(
    "Actual removed notch area [m²]"
)

ax.set_ylabel(
    "Body RMS erosion-induced displacement [mm]"
)

ax.set_title(
    "Antecedent hydraulic state × toe erosion"
)

ax.legend()

ax.grid(
    alpha=0.25
)

fig.tight_layout()

fig.savefig(
    OUT
    / "figure_01_antecedent_Rrms_vs_removed_area.png",
    dpi=220,
)

plt.close(
    fig
)


# ============================================================
# FIGURE 2
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
                "R_body_p95_mm"
            ]
            for r in records
        ],
        marker="o",
        markersize=3,
        label=result[
            "state"
        ],
    )


ax.set_xlabel(
    "Actual removed notch area [m²]"
)

ax.set_ylabel(
    "Body P95 erosion-induced displacement [mm]"
)

ax.set_title(
    "Distributed upper-tail deformation response"
)

ax.legend()

ax.grid(
    alpha=0.25
)

fig.tight_layout()

fig.savefig(
    OUT
    / "figure_02_antecedent_R95_vs_removed_area.png",
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
    / "phase05e3_summary.txt"
)


lines = [
    "PHASE 05E-3 ANTECEDENT × TOE EROSION",
    "",
    (
        "Production mesh: medium "
        "toe refinement."
    ),
    (
        "Maximum prescribed nominal "
        f"toe recession: {E_MAX} m."
    ),
    (
        "Common accepted E range [m]: "
        f"0 to {common_max_E:.3f}"
    ),
    "",
]


for result in results:

    lines.append(
        f"{result['state']} | "
        f"nonconv_E="
        f"{result['failure_E_m']} | "
        f"last_E="
        f"{result['records'][-1]['E_nominal_m']}"
    )


lines += [
    "",
    (
        "Body plasticity zero across "
        "accepted histories: "
        f"{body_plasticity_zero}"
    ),
    (
        "Fraction Rwet > Rref: "
        f"{wet_sign_fraction}"
    ),
    (
        "Fraction Rdry > Rref: "
        f"{dry_sign_fraction}"
    ),
    "",
    (
        "Interpretation rule: compare "
        "distributed deformation at the "
        "same actual removed geometry."
    ),
    (
        "Notch-tip plastic strain and solver "
        "nonconvergence are not treated as "
        "physical failure metrics."
    ),
    (
        "No wet/dry response ordering is "
        "assumed before observing results."
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
    history_csv
)

print(
    "PASS:",
    comparison_csv
)

print(
    "PASS:",
    OUT
    / "figure_01_antecedent_Rrms_vs_removed_area.png"
)

print(
    "PASS:",
    OUT
    / "figure_02_antecedent_R95_vs_removed_area.png"
)

print(
    "PASS:",
    summary
)

print()
print(
    "PHASE 05E-3 ANALYSIS COMPLETE"
)
