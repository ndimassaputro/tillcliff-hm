from pathlib import Path
import csv
import re

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
    / "phase05d2_analysis"
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
# Numerical erosion law.
# ============================================================

def erosion_from_time(
    t,
):

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
# Parse OGS failure time.
# ============================================================

def failure_time_from_log(
    path,
):

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
# Surface geometry.
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
            x[mid]
            - 8.0
        )
        / 14.0
    )


    h[
        x >= 22.0
    ] = 2.0


    return h


# ============================================================
# Field helpers.
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
        mesh.point_data[
            name
        ],
        dtype=float,
    )


def any_field(
    mesh,
    name,
):

    if name in mesh.cell_data:

        return (
            np.asarray(
                mesh.cell_data[
                    name
                ],
                dtype=float,
            ),
            "cell",
        )


    if name in mesh.point_data:

        return (
            np.asarray(
                mesh.point_data[
                    name
                ],
                dtype=float,
            ),
            "point",
        )


    raise KeyError(
        f"{name!r} missing"
    )


# ============================================================
# Analyze one mesh.
# ============================================================

def analyze_case(
    case,
):

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


    candidate = (
        material_ids > 0
    )


    recession = (
        TOE_X
        - centers[:, 0]
    )


    # Cache cell->node connectivity.
    cell_nodes = [
        np.asarray(
            bulk
            .get_cell(i)
            .point_ids,
            dtype=int,
        )
        for i in range(
            bulk.n_cells
        )
    ]


    def active_masks(
        E,
    ):

        removed = (
            candidate
            & (
                recession
                <= E + 1e-10
            )
        )


        active_cells = ~removed


        active_nodes = np.zeros(
            bulk.n_points,
            dtype=bool,
        )


        for cid in np.where(
            active_cells
        )[0]:

            active_nodes[
                cell_nodes[
                    int(cid)
                ]
            ] = True


        return (
            active_cells,
            active_nodes,
        )


    # --------------------------------------------------------
    # Slope monitoring zone.
    # --------------------------------------------------------

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


    # A smaller toe/localization window.
    local_zone = (
        (points[:, 0] >= 18.0)
        & (points[:, 0] <= 22.0)
        & (points[:, 1] >= 1.0)
        & (points[:, 1] <= 5.0)
    )


    fail_t = failure_time_from_log(
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

        return {
            "case":
                case,

            "failure_time_s":
                fail_t,

            "failure_E_m":
                fail_E,

            "records":
                [],
        }


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


    if len(valid_indices) == 0:

        raise SystemExit(
            f"FAIL: no valid output "
            f"for {case}"
        )


    # --------------------------------------------------------
    # Intact baseline at t=20 s.
    # --------------------------------------------------------

    baseline_ids = [
        i
        for i in valid_indices
        if (
            times[i]
            <= 20.0 + 1e-9
        )
    ]


    if not baseline_ids:

        raise SystemExit(
            f"FAIL: no intact baseline "
            f"for {case}"
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

    previous_active_cells = None


    for idx in valid_indices:

        t = float(
            times[idx]
        )


        if t < 20.0 - 1e-9:
            continue


        E = erosion_from_time(
            t
        )


        mesh = series.mesh(
            int(idx)
        )


        active_cells, active_nodes = (
            active_masks(
                E
            )
        )


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


        slope_nodes = (
            active_nodes
            & slope_zone
        )


        local_nodes = (
            active_nodes
            & local_zone
        )


        Rmax = float(
            np.max(
                du_mag[
                    slope_nodes
                ]
            )
            * 1000.0
        )


        R95 = float(
            np.percentile(
                du_mag[
                    slope_nodes
                ],
                95.0,
            )
            * 1000.0
        )


        Rlocal = float(
            np.max(
                du_mag[
                    local_nodes
                ]
            )
            * 1000.0
        )


        # ----------------------------------------------------
        # Plastic strain.
        # ----------------------------------------------------

        epsp, location = any_field(
            mesh,
            "EquivalentPlasticStrain",
        )


        epsp = np.asarray(
            epsp,
            dtype=float,
        ).squeeze()


        if location == "cell":

            ids = np.where(
                active_cells
                & np.isfinite(
                    epsp
                )
            )[0]


            values = epsp[
                ids
            ]


            peak_id = int(
                ids[
                    np.argmax(
                        values
                    )
                ]
            )


            peak_xyz = centers[
                peak_id
            ]


        else:

            ids = np.where(
                active_nodes
                & np.isfinite(
                    epsp
                )
            )[0]


            values = epsp[
                ids
            ]


            peak_id = int(
                ids[
                    np.argmax(
                        values
                    )
                ]
            )


            peak_xyz = points[
                peak_id
            ]


        epsp_max = float(
            np.max(
                values
            )
        )


        frac_1e4 = float(
            np.mean(
                values
                > 1e-4
            )
        )


        active_count = int(
            np.sum(
                active_cells
            )
        )


        if previous_active_cells is None:

            newly_removed = 0

        else:

            newly_removed = (
                previous_active_cells
                - active_count
            )


        records.append(
            {
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

                "Rlocal_mm":
                    Rlocal,

                "epsp_max":
                    epsp_max,

                "frac_epsp_gt_1e4":
                    frac_1e4,

                "epsp_x_m":
                    float(
                        peak_xyz[0]
                    ),

                "epsp_y_m":
                    float(
                        peak_xyz[1]
                    ),
            }
        )


        previous_active_cells = (
            active_count
        )


    # --------------------------------------------------------
    # Event table from exact mesh geometry.
    # --------------------------------------------------------

    events = []

    with (
        model_dir
        / "erosion_events.csv"
    ).open(
        "r",
        encoding="utf-8",
    ) as f:

        reader = csv.DictReader(
            f
        )

        for row in reader:

            events.append(
                {
                    "E_m":
                        float(
                            row["E_m"]
                        ),

                    "cell_count":
                        int(
                            row[
                                "cell_count"
                            ]
                        ),

                    "area_m2":
                        float(
                            row[
                                "area_m2"
                            ]
                        ),
                }
            )


    if fail_t is not None:

        nearest_failure_event = min(
            events,
            key=lambda r:
                abs(
                    r["E_m"]
                    - fail_E
                ),
        )

    else:

        nearest_failure_event = None


    return {
        "case":
            case,

        "failure_time_s":
            fail_t,

        "failure_E_m":
            fail_E,

        "failure_event":
            nearest_failure_event,

        "events":
            events,

        "records":
            records,
    }


# ============================================================
# Analyze all.
# ============================================================

results = [
    analyze_case(
        case
    )
    for case in CASES
]


# ============================================================
# Summary diagnostics.
# ============================================================

comparison = []


print(
    "========================================"
)

print(
    "PHASE 05D-2 MESH EROSION ROBUSTNESS"
)

print(
    "========================================"
)


for result in results:

    case = result[
        "case"
    ]

    records = result[
        "records"
    ]


    print()
    print(
        f"=== {case.upper()} ==="
    )


    if not records:

        print(
            "NO VALID RECORDS"
        )

        comparison.append(
            {
                "case":
                    case,

                "failure_E_m":
                    result[
                        "failure_E_m"
                    ],

                "last_valid_E_m":
                    np.nan,

                "first_yield_E_m":
                    np.nan,

                "first_epsp1e4_E_m":
                    np.nan,

                "last_Rmax_mm":
                    np.nan,

                "last_epsp_max":
                    np.nan,
            }
        )

        continue


    first_yield = next(
        (
            r
            for r in records
            if r["epsp_max"]
            > 1e-8
        ),
        None,
    )


    first_1e4 = next(
        (
            r
            for r in records
            if r["epsp_max"]
            > 1e-4
        ),
        None,
    )


    last = records[-1]


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


        event = result[
            "failure_event"
        ]


        if event is not None:

            print(
                "Nearest topology event:"
            )

            print(
                f"  E     = "
                f"{event['E_m']:.6f} m"
            )

            print(
                f"  cells = "
                f"{event['cell_count']}"
            )

            print(
                f"  area  = "
                f"{event['area_m2']:.10f} m2"
            )


    print(
        "Last valid saved E [m]: "
        f"{last['E_m']:.6f}"
    )

    print(
        "Last Rmax [mm]: "
        f"{last['Rmax_mm']:.8f}"
    )

    print(
        "Last Rlocal [mm]: "
        f"{last['Rlocal_mm']:.8f}"
    )

    print(
        "Last epsp_max: "
        f"{last['epsp_max']:.8e}"
    )


    if first_yield is None:

        first_yield_E = np.nan

        print(
            "First yielding: NONE"
        )

    else:

        first_yield_E = (
            first_yield[
                "E_m"
            ]
        )

        print(
            "First resolved yielding E [m]: "
            f"{first_yield_E:.6f}"
        )


    if first_1e4 is None:

        first_1e4_E = np.nan

        print(
            "First epsp_max > 1e-4: NONE"
        )

    else:

        first_1e4_E = (
            first_1e4[
                "E_m"
            ]
        )

        print(
            "First epsp_max > 1e-4 E [m]: "
            f"{first_1e4_E:.6f}"
        )


    print()
    print(
        "TOPOLOGY-CHANGE STATES"
    )


    event_records = [
        r
        for r in records
        if r[
            "newly_removed_cells"
        ] > 0
    ]


    for r in event_records:

        print(
            f"E={r['E_m']:.3f} | "
            f"removed+="
            f"{r['newly_removed_cells']:2d} | "
            f"Rmax="
            f"{r['Rmax_mm']:.6f} mm | "
            f"Rlocal="
            f"{r['Rlocal_mm']:.6f} mm | "
            f"epsp="
            f"{r['epsp_max']:.3e}"
        )


    comparison.append(
        {
            "case":
                case,

            "failure_E_m":
                result[
                    "failure_E_m"
                ],

            "last_valid_E_m":
                last[
                    "E_m"
                ],

            "first_yield_E_m":
                first_yield_E,

            "first_epsp1e4_E_m":
                first_1e4_E,

            "last_Rmax_mm":
                last[
                    "Rmax_mm"
                ],

            "last_epsp_max":
                last[
                    "epsp_max"
                ],
        }
    )


# ============================================================
# Mesh-robustness decision.
#
# Medium vs fine are the key convergence pair.
# ============================================================

medium = next(
    r
    for r in comparison
    if r["case"] == "medium"
)

fine = next(
    r
    for r in comparison
    if r["case"] == "fine"
)


print()
print(
    "========================================"
)

print(
    "MESH-ROBUSTNESS DECISION"
)

print(
    "========================================"
)


# ------------------------------------------------------------
# Yield onset.
# ------------------------------------------------------------

if (
    np.isfinite(
        medium[
            "first_yield_E_m"
        ]
    )
    and np.isfinite(
        fine[
            "first_yield_E_m"
        ]
    )
):

    delta_yield = abs(
        medium[
            "first_yield_E_m"
        ]
        - fine[
            "first_yield_E_m"
        ]
    )


    print(
        "Medium first yield E [m]: "
        f"{medium['first_yield_E_m']:.6f}"
    )

    print(
        "Fine first yield E [m]: "
        f"{fine['first_yield_E_m']:.6f}"
    )

    print(
        "|Delta E_yield| [m]: "
        f"{delta_yield:.6f}"
    )


    yield_robust = (
        delta_yield
        <= 0.05
    )

else:

    delta_yield = np.nan
    yield_robust = False

    print(
        "Yield-onset comparison unavailable."
    )


# ------------------------------------------------------------
# Nonconvergence.
# ------------------------------------------------------------

medium_fail = (
    medium[
        "failure_E_m"
    ]
)

fine_fail = (
    fine[
        "failure_E_m"
    ]
)


if (
    np.isfinite(
        medium_fail
    )
    and np.isfinite(
        fine_fail
    )
):

    delta_fail = abs(
        medium_fail
        - fine_fail
    )


    print()
    print(
        "Medium nonconvergence E [m]: "
        f"{medium_fail:.6f}"
    )

    print(
        "Fine nonconvergence E [m]: "
        f"{fine_fail:.6f}"
    )

    print(
        "|Delta E_nonconv| [m]: "
        f"{delta_fail:.6f}"
    )


    nonconv_robust = (
        delta_fail
        <= 0.05
    )

else:

    delta_fail = np.nan
    nonconv_robust = False


    print()
    print(
        "Medium/fine do not both "
        "show nonconvergence."
    )


# ============================================================
# Classification.
# ============================================================

print()
print(
    "=== INTERPRETATION ==="
)


if (
    yield_robust
    and nonconv_robust
):

    interpretation = (
        "MESH-CONVERGING TRANSITION CANDIDATE"
    )

    print(
        interpretation
    )

    print(
        "Yield onset and nonconvergence "
        "locations are similar on medium "
        "and fine toe meshes."
    )


elif (
    yield_robust
    and (
        not np.isfinite(
            fine_fail
        )
    )
):

    interpretation = (
        "YIELD ONSET ROBUST; "
        "COARSE NONCONVERGENCE NOT ROBUST"
    )

    print(
        interpretation
    )


elif (
    not np.isfinite(
        medium_fail
    )
    and not np.isfinite(
        fine_fail
    )
):

    interpretation = (
        "COARSE FAILURE IS MESH-DEPENDENT"
    )

    print(
        interpretation
    )

    print(
        "Medium and fine meshes remain "
        "convergent through the tested "
        "erosion range."
    )


else:

    interpretation = (
        "MESH SENSITIVITY REMAINS"
    )

    print(
        interpretation
    )

    print(
        "Do not define Ecrit yet."
    )


# ============================================================
# CSV.
# ============================================================

csv_path = (
    OUT
    / "mesh_erosion_comparison.csv"
)


with csv_path.open(
    "w",
    newline="",
    encoding="utf-8",
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "case",
            "failure_E_m",
            "last_valid_E_m",
            "first_yield_E_m",
            "first_epsp1e4_E_m",
            "last_Rmax_mm",
            "last_epsp_max",
        ],
    )

    writer.writeheader()
    writer.writerows(
        comparison
    )


# ============================================================
# Summary.
# ============================================================

summary = (
    OUT
    / "phase05d2_summary.txt"
)


lines = [
    "PHASE 05D-2 MESH EROSION ROBUSTNESS",
    "",
]


for r in comparison:

    lines.append(
        f"{r['case']} | "
        f"first_yield_E="
        f"{r['first_yield_E_m']} | "
        f"first_epsp1e4_E="
        f"{r['first_epsp1e4_E_m']} | "
        f"last_valid_E="
        f"{r['last_valid_E_m']} | "
        f"nonconv_E="
        f"{r['failure_E_m']}"
    )


lines += [
    "",
    (
        "Medium-fine Delta E_yield [m]: "
        f"{delta_yield}"
    ),
    (
        "Medium-fine Delta E_nonconv [m]: "
        f"{delta_fail}"
    ),
    "",
    (
        "INTERPRETATION: "
        f"{interpretation}"
    ),
    "",
    (
        "Solver nonconvergence is not by "
        "itself defined as physical failure."
    ),
    (
        "The medium/fine pair is used as "
        "the primary mesh-robustness test."
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
    csv_path,
)

print(
    "PASS:",
    summary,
)

print()
print(
    "PHASE 05D-2 ANALYSIS COMPLETE"
)
