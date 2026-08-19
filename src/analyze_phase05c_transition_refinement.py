from pathlib import Path
import csv
import re

import numpy as np
import ogstools as ot
import pyvista as pv


ROOT = Path.cwd()

MODEL = (
    ROOT
    / "model"
    / "phase05c_transition_refinement"
)

RUN_ROOT = (
    ROOT
    / "results"
    / "phase05c_transition_refinement"
)

OUT = (
    ROOT
    / "results"
    / "phase05c_analysis"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# TOPOLOGY
# ============================================================

mesh = pv.read(
    MODEL
    / "slope_toe_notch_ready.vtu"
)

material_ids = np.asarray(
    mesh.cell_data[
        "MaterialIDs"
    ],
    dtype=int,
)

centers = (
    mesh
    .cell_centers()
    .points
)

points = np.asarray(
    mesh.points,
    dtype=float,
)

TOE_X = 22.0

recession = (
    TOE_X
    - centers[:, 0]
)

candidate = (
    material_ids > 0
)


cell_points = [
    np.asarray(
        mesh.get_cell(i).point_ids,
        dtype=int,
    )
    for i in range(
        mesh.n_cells
    )
]


def erosion_from_time(t):

    if t <= 10.0:
        return 0.0

    if t >= 110.0:
        return 2.0

    return (
        0.02
        * (
            t - 10.0
        )
    )


def active_masks(E):

    removed = (
        candidate
        & (
            recession
            <= E + 1e-10
        )
    )

    active_cells = ~removed

    active_nodes = np.zeros(
        mesh.n_points,
        dtype=bool,
    )

    for cell_id in np.where(
        active_cells
    )[0]:

        active_nodes[
            cell_points[
                int(cell_id)
            ]
        ] = True

    return (
        active_cells,
        active_nodes,
    )


# ============================================================
# SLOPE MONITOR ZONE
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
# HELPERS
# ============================================================

def point_field(
    m,
    name,
):

    if name not in m.point_data:

        raise KeyError(
            f"{name} missing from point data"
        )

    return np.asarray(
        m.point_data[name],
        dtype=float,
    )


def any_field(
    m,
    name,
):

    if name in m.cell_data:

        return (
            np.asarray(
                m.cell_data[name],
                dtype=float,
            ),
            "cell",
        )

    if name in m.point_data:

        return (
            np.asarray(
                m.point_data[name],
                dtype=float,
            ),
            "point",
        )

    raise KeyError(
        f"{name} missing"
    )


def parse_failure_time(
    log_path,
):

    if not log_path.exists():
        return None

    text = log_path.read_text(
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
# ANALYZE CASE
# ============================================================

def analyze_case(
    case,
):

    folder = (
        RUN_ROOT
        / case
    )

    failure_time = parse_failure_time(
        folder / "ogs.log"
    )

    pvds = sorted(
        folder.glob("*.pvd")
    )

    if not pvds:

        return {
            "case":
                case,

            "failure_time_s":
                failure_time,

            "failure_E_m":
                (
                    erosion_from_time(
                        failure_time
                    )
                    if failure_time
                    is not None
                    else np.nan
                ),

            "records":
                [],
        }


    series = ot.MeshSeries(
        str(pvds[0])
    )

    times = np.asarray(
        series.timevalues,
        dtype=float,
    )


    if failure_time is None:

        valid_indices = np.arange(
            len(times)
        )

    else:

        valid_indices = np.where(
            times
            < failure_time
            - 1e-10
        )[0]


    baseline_candidates = [
        i
        for i in valid_indices
        if times[i] <= 10.0 + 1e-9
    ]

    if not baseline_candidates:

        raise SystemExit(
            f"FAIL: no baseline for {case}"
        )

    i0 = int(
        baseline_candidates[-1]
    )

    base = series.mesh(
        i0
    )

    u0 = point_field(
        base,
        "displacement",
    )


    records = []

    previous_active_count = None
    previous_R95 = None
    previous_event_E = None


    for idx in valid_indices:

        t = float(
            times[idx]
        )

        if t < 10.0 - 1e-9:
            continue

        E = erosion_from_time(
            t
        )

        m = series.mesh(
            int(idx)
        )

        active_cells, active_nodes = (
            active_masks(
                E
            )
        )

        monitor = (
            active_nodes
            & slope_zone
        )

        u = point_field(
            m,
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

        vals = du_mag[
            monitor
        ]

        vals = vals[
            np.isfinite(
                vals
            )
        ]

        Rmax = float(
            np.max(vals)
            * 1000.0
        )

        R95 = float(
            np.percentile(
                vals,
                95.0,
            )
            * 1000.0
        )

        Rrms = float(
            np.sqrt(
                np.mean(
                    vals ** 2
                )
            )
            * 1000.0
        )


        epsp, loc = any_field(
            m,
            "EquivalentPlasticStrain",
        )

        epsp = np.asarray(
            epsp,
            dtype=float,
        ).squeeze()


        if loc == "cell":

            ids = np.where(
                active_cells
                & np.isfinite(
                    epsp
                )
            )[0]

            local = epsp[
                ids
            ]

            j = int(
                ids[
                    np.argmax(
                        local
                    )
                ]
            )

            epsp_xyz = centers[
                j
            ]

        else:

            ids = np.where(
                active_nodes
                & np.isfinite(
                    epsp
                )
            )[0]

            local = epsp[
                ids
            ]

            j = int(
                ids[
                    np.argmax(
                        local
                    )
                ]
            )

            epsp_xyz = points[
                j
            ]


        epsp_max = float(
            np.max(
                local
            )
        )

        frac_1e4 = float(
            np.mean(
                local > 1e-4
            )
        )

        active_count = int(
            np.sum(
                active_cells
            )
        )


        if previous_active_count is None:

            newly_removed = 0

        else:

            newly_removed = (
                previous_active_count
                - active_count
            )


        event_compliance = np.nan

        if (
            newly_removed > 0
            and previous_event_E
            is not None
            and previous_R95
            is not None
            and E > previous_event_E
        ):

            event_compliance = (
                (
                    R95
                    - previous_R95
                )
                / (
                    E
                    - previous_event_E
                )
            )


        record = {
            "time_s":
                t,

            "E_m":
                E,

            "active_cells":
                active_count,

            "newly_removed_cells":
                newly_removed,

            "Rmax_mm":
                Rmax,

            "R95_mm":
                R95,

            "Rrms_mm":
                Rrms,

            "epsp_max":
                epsp_max,

            "frac_epsp_gt_1e4":
                frac_1e4,

            "epsp_x_m":
                float(
                    epsp_xyz[0]
                ),

            "epsp_y_m":
                float(
                    epsp_xyz[1]
                ),

            "event_C95_mm_per_m":
                event_compliance,
        }

        records.append(
            record
        )


        if newly_removed > 0:

            previous_event_E = E
            previous_R95 = R95


        previous_active_count = (
            active_count
        )


    return {
        "case":
            case,

        "failure_time_s":
            failure_time,

        "failure_E_m":
            (
                erosion_from_time(
                    failure_time
                )
                if failure_time
                is not None
                else np.nan
            ),

        "records":
            records,
    }


# ============================================================
# RUN BOTH
# ============================================================

cases = [
    analyze_case(
        "dt0p10"
    ),
    analyze_case(
        "dt0p05"
    ),
]


print(
    "========================================"
)

print(
    "PHASE 05C TRANSITION REFINEMENT"
)

print(
    "========================================"
)


summary_lines = [
    "PHASE 05C TRANSITION REFINEMENT",
    "",
]


for result in cases:

    case = result[
        "case"
    ]

    records = result[
        "records"
    ]

    print()
    print(
        f"=== {case} ==="
    )


    if result[
        "failure_time_s"
    ] is None:

        print(
            "Solver failure: NONE"
        )

    else:

        print(
            "Failure time [s]: "
            f"{result['failure_time_s']:.6f}"
        )

        print(
            "Failure nominal E [m]: "
            f"{result['failure_E_m']:.6f}"
        )


    if not records:

        print(
            "NO VALID RECORDS"
        )

        continue


    first_yield = next(
        (
            r
            for r in records
            if r["epsp_max"] > 1e-8
        ),
        None,
    )

    first_sig_yield = next(
        (
            r
            for r in records
            if r["epsp_max"] > 1e-4
        ),
        None,
    )


    last = records[-1]

    print(
        "Last saved valid E [m]: "
        f"{last['E_m']:.6f}"
    )

    print(
        "Last saved Rmax [mm]: "
        f"{last['Rmax_mm']:.8f}"
    )

    print(
        "Last saved R95 [mm]: "
        f"{last['R95_mm']:.8f}"
    )

    print(
        "Last saved epsp_max: "
        f"{last['epsp_max']:.8e}"
    )

    print(
        "Last epsp location [m]: "
        f"({last['epsp_x_m']:.5f}, "
        f"{last['epsp_y_m']:.5f})"
    )


    if first_yield:

        print(
            "First resolved yielding E [m]: "
            f"{first_yield['E_m']:.6f}"
        )

    else:

        print(
            "First resolved yielding: NONE"
        )


    if first_sig_yield:

        print(
            "First epsp_max > 1e-4 E [m]: "
            f"{first_sig_yield['E_m']:.6f}"
        )

    else:

        print(
            "First epsp_max > 1e-4: NONE"
        )


    print()
    print(
        "LAST 12 VALID SAVED STATES"
    )

    for r in records[-12:]:

        print(
            f"E={r['E_m']:.3f} | "
            f"removed+={r['newly_removed_cells']:2d} | "
            f"Rmax={r['Rmax_mm']:.6f} mm | "
            f"R95={r['R95_mm']:.6f} mm | "
            f"epsp={r['epsp_max']:.3e}"
        )


    summary_lines += [
        f"CASE: {case}",
        (
            "Failure E [m]: "
            f"{result['failure_E_m']}"
        ),
        (
            "Last valid saved E [m]: "
            f"{last['E_m']:.8f}"
        ),
        (
            "Last Rmax [mm]: "
            f"{last['Rmax_mm']:.10f}"
        ),
        (
            "Last R95 [mm]: "
            f"{last['R95_mm']:.10f}"
        ),
        (
            "Last epsp_max: "
            f"{last['epsp_max']:.12e}"
        ),
        (
            "First resolved yield E [m]: "
            + (
                f"{first_yield['E_m']:.8f}"
                if first_yield
                else "NONE"
            )
        ),
        (
            "First epsp>1e-4 E [m]: "
            + (
                f"{first_sig_yield['E_m']:.8f}"
                if first_sig_yield
                else "NONE"
            )
        ),
        "",
    ]


# ============================================================
# GEOMETRY EVENTS NEAR FAILURE
# ============================================================

audit_csv = (
    OUT
    / "deactivation_events_035_060.csv"
)

print()
print(
    "=== GEOMETRIC EVENTS 0.35–0.60 m ==="
)

if audit_csv.exists():

    with audit_csv.open(
        "r",
        encoding="utf-8",
    ) as f:

        rows = list(
            csv.DictReader(
                f
            )
        )

    for row in rows:

        print(
            f"E={float(row['E_event_m']):.6f} | "
            f"cells={row['cell_count']} | "
            f"area={float(row['area_m2']):.8f} m2 | "
            f"IDs={row['material_ids']}"
        )


# ============================================================
# TEMPORAL-RESOLUTION COMPARISON
# ============================================================

fail_E = [
    r["failure_E_m"]
    for r in cases
    if np.isfinite(
        r["failure_E_m"]
    )
]


print()
print(
    "=== FAILURE-E COMPARISON ==="
)


if len(fail_E) == 2:

    delta = abs(
        fail_E[0]
        - fail_E[1]
    )

    print(
        "dt=0.10 failure E [m]: "
        f"{fail_E[0]:.8f}"
    )

    print(
        "dt=0.05 failure E [m]: "
        f"{fail_E[1]:.8f}"
    )

    print(
        "|Delta E_fail| [m]: "
        f"{delta:.8f}"
    )


    if delta <= 0.01:

        print(
            "TEMPORAL-RESOLUTION CHECK: "
            "FAILURE LOCATION STABLE"
        )

    else:

        print(
            "TEMPORAL-RESOLUTION CHECK: "
            "TIME-STEP SENSITIVE"
        )

else:

    print(
        "Failure comparison unavailable "
        "because not both cases failed."
    )


# ============================================================
# WRITE SUMMARY
# ============================================================

summary = (
    OUT
    / "phase05c_summary.txt"
)

summary_lines += [
    (
        "Interpretation rule:"
    ),
    (
        "Solver nonconvergence alone is NOT "
        "defined as physical failure."
    ),
    (
        "Candidate transition requires "
        "response acceleration and/or plastic "
        "localization before nonconvergence."
    ),
]

summary.write_text(
    "\n".join(
        summary_lines
    )
    + "\n",
    encoding="utf-8",
)

print()
print(
    "PASS:",
    summary,
)

print()
print(
    "PHASE 05C ANALYSIS COMPLETE"
)
