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
    / "phase05d3_localization"
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
# NUMERICAL EROSION LAW FROM PHASE 05D-2
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
# FIELD -> CELL VALUES
#
# EquivalentPlasticStrain is expected to be cell data
# for the MFront material.
#
# Point-data fallback is included for robustness.
# ============================================================

def epsp_cell_values(
    output_mesh,
    bulk_mesh,
):

    name = "EquivalentPlasticStrain"

    if name in output_mesh.cell_data:

        arr = np.asarray(
            output_mesh.cell_data[name],
            dtype=float,
        ).squeeze()

        if len(arr) != bulk_mesh.n_cells:

            raise RuntimeError(
                "Cell EPSP size mismatch"
            )

        return arr


    if name in output_mesh.point_data:

        p = np.asarray(
            output_mesh.point_data[name],
            dtype=float,
        ).squeeze()

        values = np.empty(
            bulk_mesh.n_cells,
            dtype=float,
        )

        for cid in range(
            bulk_mesh.n_cells
        ):

            ids = (
                bulk_mesh
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
# CELL EDGE ADJACENCY
#
# QUAD mesh:
# two cells sharing one complete edge are neighbours.
# ============================================================

def build_neighbours(mesh):

    edge_map = {}

    for cid in range(
        mesh.n_cells
    ):

        ids = list(
            mesh
            .get_cell(cid)
            .point_ids
        )

        n = len(ids)

        for j in range(n):

            a = int(
                ids[j]
            )

            b = int(
                ids[
                    (j + 1) % n
                ]
            )

            edge = tuple(
                sorted(
                    (
                        a,
                        b,
                    )
                )
            )

            edge_map.setdefault(
                edge,
                [],
            ).append(
                cid
            )


    neighbours = [
        set()
        for _ in range(
            mesh.n_cells
        )
    ]


    for cells in (
        edge_map.values()
    ):

        if len(cells) != 2:
            continue

        a, b = cells

        neighbours[a].add(
            b
        )

        neighbours[b].add(
            a
        )


    return neighbours


# ============================================================
# MIN DISTANCE TO SET OF CELL CENTRES
# ============================================================

def min_distances(
    points,
    targets,
):

    if len(targets) == 0:

        return np.full(
            len(points),
            np.nan,
        )


    result = np.full(
        len(points),
        np.inf,
        dtype=float,
    )


    # Small target sets, so a simple loop is
    # memory-safe and fast enough.
    for target in targets:

        d = np.linalg.norm(
            points[:, :2]
            - target[:2],
            axis=1,
        )

        result = np.minimum(
            result,
            d,
        )


    return result


# ============================================================
# ANALYZE ONE MESH
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


    bulk_path = (
        model_dir
        / "bulk.vtu"
    )

    if not bulk_path.exists():

        raise RuntimeError(
            f"Missing {bulk_path}"
        )


    bulk = pv.read(
        bulk_path
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
    # Toe horizontal resolution.
    # --------------------------------------------------------

    x_unique = np.unique(
        np.round(
            bulk.points[:, 0],
            10,
        )
    )

    dx = np.diff(
        x_unique
    )

    xmid = (
        0.5
        * (
            x_unique[:-1]
            + x_unique[1:]
        )
    )

    toe_dx_values = dx[
        (xmid >= 19.8)
        & (xmid <= 22.2)
    ]

    toe_dx = float(
        np.median(
            toe_dx_values
        )
    )


    # --------------------------------------------------------
    # Cell neighbours.
    # --------------------------------------------------------

    neighbours = build_neighbours(
        bulk
    )


    # --------------------------------------------------------
    # Failure / valid outputs.
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


    records = []

    previous_active_count = None


    for idx in valid_indices:

        t = float(
            times[idx]
        )

        if t < 20.0 - 1e-10:
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

        active = ~removed


        active_count = int(
            np.sum(
                active
            )
        )


        if previous_active_count is None:

            newly_removed = 0

        else:

            newly_removed = (
                previous_active_count
                - active_count
            )


        # ----------------------------------------------------
        # Active cells sharing an EDGE with removed cells.
        #
        # These are the cells directly adjacent to the
        # numerical excavation/notch interface.
        # ----------------------------------------------------

        interface_active = set()


        removed_ids = np.where(
            removed
        )[0]


        for rid in removed_ids:

            for nid in neighbours[
                int(rid)
            ]:

                if active[nid]:

                    interface_active.add(
                        int(nid)
                    )


        interface_active = np.asarray(
            sorted(
                interface_active
            ),
            dtype=int,
        )


        # ----------------------------------------------------
        # EPSP.
        # ----------------------------------------------------

        output_mesh = series.mesh(
            int(idx)
        )


        epsp = epsp_cell_values(
            output_mesh,
            bulk,
        )


        finite_active = (
            active
            & np.isfinite(
                epsp
            )
        )


        active_ids = np.where(
            finite_active
        )[0]


        active_epsp = epsp[
            active_ids
        ]


        peak_local = int(
            np.argmax(
                active_epsp
            )
        )


        peak_id = int(
            active_ids[
                peak_local
            ]
        )


        peak_epsp = float(
            epsp[
                peak_id
            ]
        )


        peak_xyz = centers[
            peak_id
        ]


        # ----------------------------------------------------
        # Distance to active/eroded interface.
        # ----------------------------------------------------

        if len(
            interface_active
        ):

            interface_centers = centers[
                interface_active
            ]


            d_peak_interface = float(
                np.min(
                    np.linalg.norm(
                        interface_centers[:, :2]
                        - peak_xyz[:2],
                        axis=1,
                    )
                )
            )


            peak_is_interface = bool(
                peak_id
                in set(
                    interface_active.tolist()
                )
            )


            all_d_interface = (
                min_distances(
                    centers,
                    interface_centers,
                )
            )

        else:

            d_peak_interface = np.nan

            peak_is_interface = False

            all_d_interface = np.full(
                bulk.n_cells,
                np.nan,
            )


        d_peak_norm = (
            d_peak_interface
            / toe_dx
            if np.isfinite(
                d_peak_interface
            )
            else np.nan
        )


        # ----------------------------------------------------
        # Nominal moving plane position.
        # ----------------------------------------------------

        nominal_front_x = (
            TOE_X
            - E
        )


        dx_peak_front = abs(
            float(
                peak_xyz[0]
            )
            - nominal_front_x
        )


        # ----------------------------------------------------
        # Plastic area.
        #
        # These are integrated cell-area metrics,
        # much less vulnerable than epsp_max alone.
        # ----------------------------------------------------

        plastic_1e6 = (
            finite_active
            & (
                epsp > 1e-6
            )
        )


        plastic_1e4 = (
            finite_active
            & (
                epsp > 1e-4
            )
        )


        area_1e6 = float(
            np.sum(
                cell_area[
                    plastic_1e6
                ]
            )
        )


        area_1e4 = float(
            np.sum(
                cell_area[
                    plastic_1e4
                ]
            )
        )


        # ----------------------------------------------------
        # How much plastic area lies close to interface?
        #
        # "near" = within 1.5 local toe dx.
        #
        # This is a diagnostic geometry scale,
        # NOT a physical failure criterion.
        # ----------------------------------------------------

        if len(
            interface_active
        ):

            near_interface = (
                all_d_interface
                <= 1.5
                * toe_dx
            )


            plastic_near_1e6 = (
                plastic_1e6
                & near_interface
            )

            plastic_near_1e4 = (
                plastic_1e4
                & near_interface
            )


            area_near_1e6 = float(
                np.sum(
                    cell_area[
                        plastic_near_1e6
                    ]
                )
            )


            area_near_1e4 = float(
                np.sum(
                    cell_area[
                        plastic_near_1e4
                    ]
                )
            )

        else:

            area_near_1e6 = 0.0
            area_near_1e4 = 0.0


        frac_near_1e6 = (
            area_near_1e6
            / area_1e6
            if area_1e6 > 0
            else np.nan
        )


        frac_near_1e4 = (
            area_near_1e4
            / area_1e4
            if area_1e4 > 0
            else np.nan
        )


        records.append(
            {
                "case":
                    case,

                "time_s":
                    t,

                "E_m":
                    E,

                "toe_dx_m":
                    toe_dx,

                "active_cells":
                    active_count,

                "newly_removed_cells":
                    newly_removed,

                "interface_active_cells":
                    int(
                        len(
                            interface_active
                        )
                    ),

                "epsp_max":
                    peak_epsp,

                "epsp_peak_cell":
                    peak_id,

                "epsp_peak_x_m":
                    float(
                        peak_xyz[0]
                    ),

                "epsp_peak_y_m":
                    float(
                        peak_xyz[1]
                    ),

                "nominal_front_x_m":
                    nominal_front_x,

                "peak_dx_to_nominal_front_m":
                    dx_peak_front,

                "peak_is_interface_adjacent":
                    int(
                        peak_is_interface
                    ),

                "peak_distance_to_interface_m":
                    d_peak_interface,

                "peak_distance_to_interface_over_dx":
                    d_peak_norm,

                "plastic_area_gt_1e6_m2":
                    area_1e6,

                "plastic_area_gt_1e4_m2":
                    area_1e4,

                "plastic_area_near_interface_gt_1e6_m2":
                    area_near_1e6,

                "plastic_area_near_interface_gt_1e4_m2":
                    area_near_1e4,

                "fraction_plastic_area_near_interface_gt_1e6":
                    frac_near_1e6,

                "fraction_plastic_area_near_interface_gt_1e4":
                    frac_near_1e4,
            }
        )


        previous_active_count = (
            active_count
        )


    return {
        "case":
            case,

        "failure_time_s":
            fail_t,

        "failure_E_m":
            (
                erosion_from_time(
                    fail_t
                )
                if fail_t
                is not None
                else np.nan
            ),

        "toe_dx_m":
            toe_dx,

        "records":
            records,
    }


# ============================================================
# RUN ALL CASES
# ============================================================

results = [
    analyze_case(
        case
    )
    for case in CASES
]


# ============================================================
# WRITE FULL HISTORY
# ============================================================

all_records = []

for result in results:

    all_records.extend(
        result[
            "records"
        ]
    )


history_csv = (
    OUT
    / "localization_history.csv"
)


fieldnames = [
    "case",
    "time_s",
    "E_m",
    "toe_dx_m",
    "active_cells",
    "newly_removed_cells",
    "interface_active_cells",
    "epsp_max",
    "epsp_peak_cell",
    "epsp_peak_x_m",
    "epsp_peak_y_m",
    "nominal_front_x_m",
    "peak_dx_to_nominal_front_m",
    "peak_is_interface_adjacent",
    "peak_distance_to_interface_m",
    "peak_distance_to_interface_over_dx",
    "plastic_area_gt_1e6_m2",
    "plastic_area_gt_1e4_m2",
    "plastic_area_near_interface_gt_1e6_m2",
    "plastic_area_near_interface_gt_1e4_m2",
    "fraction_plastic_area_near_interface_gt_1e6",
    "fraction_plastic_area_near_interface_gt_1e4",
]


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
# CASE SUMMARIES
# ============================================================

summary_rows = []


print(
    "========================================"
)

print(
    "PHASE 05D-3 LOCALIZATION AUDIT"
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


    first_yield = next(
        (
            r
            for r in records
            if r[
                "epsp_max"
            ] > 1e-8
        ),
        None,
    )


    first_1e4 = next(
        (
            r
            for r in records
            if r[
                "epsp_max"
            ] > 1e-4
        ),
        None,
    )


    last = records[-1]


    print(
        "Toe dx [m]: "
        f"{result['toe_dx_m']:.6f}"
    )

    print(
        "Failure E [m]: "
        f"{result['failure_E_m']:.6f}"
    )

    print(
        "Last valid E [m]: "
        f"{last['E_m']:.6f}"
    )


    if first_yield is not None:

        print()
        print(
            "FIRST YIELD STATE"
        )

        print(
            f"E = "
            f"{first_yield['E_m']:.6f} m"
        )

        print(
            "epsp_max = "
            f"{first_yield['epsp_max']:.8e}"
        )

        print(
            "peak location = "
            f"("
            f"{first_yield['epsp_peak_x_m']:.6f}, "
            f"{first_yield['epsp_peak_y_m']:.6f}"
            f")"
        )

        print(
            "peak interface-adjacent = "
            f"{first_yield['peak_is_interface_adjacent']}"
        )

        print(
            "d_interface / dx = "
            f"{first_yield['peak_distance_to_interface_over_dx']:.6f}"
        )

        print(
            "plastic area >1e-4 [m2] = "
            f"{first_yield['plastic_area_gt_1e4_m2']:.10f}"
        )

        print(
            "fraction plastic area >1e-4 "
            "near interface = "
            f"{first_yield['fraction_plastic_area_near_interface_gt_1e4']}"
        )


    print()
    print(
        "LAST VALID STATE"
    )

    print(
        "epsp_max = "
        f"{last['epsp_max']:.8e}"
    )

    print(
        "peak location = "
        f"("
        f"{last['epsp_peak_x_m']:.6f}, "
        f"{last['epsp_peak_y_m']:.6f}"
        f")"
    )

    print(
        "nominal front x = "
        f"{last['nominal_front_x_m']:.6f}"
    )

    print(
        "peak interface-adjacent = "
        f"{last['peak_is_interface_adjacent']}"
    )

    print(
        "peak distance to interface [m] = "
        f"{last['peak_distance_to_interface_m']:.8f}"
    )

    print(
        "d_interface / dx = "
        f"{last['peak_distance_to_interface_over_dx']:.6f}"
    )

    print(
        "plastic area >1e-6 [m2] = "
        f"{last['plastic_area_gt_1e6_m2']:.10f}"
    )

    print(
        "plastic area >1e-4 [m2] = "
        f"{last['plastic_area_gt_1e4_m2']:.10f}"
    )

    print(
        "fraction plastic area >1e-4 "
        "near interface = "
        f"{last['fraction_plastic_area_near_interface_gt_1e4']}"
    )


    # --------------------------------------------------------
    # Topology events with plasticity.
    # --------------------------------------------------------

    print()
    print(
        "PLASTIC TOPOLOGY EVENTS"
    )


    plastic_events = [
        r
        for r in records
        if (
            r[
                "newly_removed_cells"
            ] > 0
            and r[
                "epsp_max"
            ] > 1e-8
        )
    ]


    for r in plastic_events:

        print(
            f"E={r['E_m']:.3f} | "
            f"removed+="
            f"{r['newly_removed_cells']:2d} | "
            f"epsp="
            f"{r['epsp_max']:.3e} | "
            f"adjacent="
            f"{r['peak_is_interface_adjacent']} | "
            f"d/dx="
            f"{r['peak_distance_to_interface_over_dx']:.3f} | "
            f"A1e-4="
            f"{r['plastic_area_gt_1e4_m2']:.6f}"
        )


    summary_rows.append(
        {
            "case":
                case,

            "toe_dx_m":
                result[
                    "toe_dx_m"
                ],

            "failure_E_m":
                result[
                    "failure_E_m"
                ],

            "first_yield_E_m":
                (
                    first_yield[
                        "E_m"
                    ]
                    if first_yield
                    is not None
                    else np.nan
                ),

            "first_yield_peak_adjacent":
                (
                    first_yield[
                        "peak_is_interface_adjacent"
                    ]
                    if first_yield
                    is not None
                    else np.nan
                ),

            "first_yield_d_over_dx":
                (
                    first_yield[
                        "peak_distance_to_interface_over_dx"
                    ]
                    if first_yield
                    is not None
                    else np.nan
                ),

            "last_valid_E_m":
                last[
                    "E_m"
                ],

            "last_peak_adjacent":
                last[
                    "peak_is_interface_adjacent"
                ],

            "last_d_over_dx":
                last[
                    "peak_distance_to_interface_over_dx"
                ],

            "last_plastic_area_1e4_m2":
                last[
                    "plastic_area_gt_1e4_m2"
                ],

            "last_fraction_plastic_near_interface_1e4":
                last[
                    "fraction_plastic_area_near_interface_gt_1e4"
                ],
        }
    )


# ============================================================
# MEDIUM/FINE DIAGNOSTIC
#
# This is a numerical localization diagnostic,
# NOT a physical failure criterion.
# ============================================================

medium = next(
    r
    for r in summary_rows
    if r["case"] == "medium"
)

fine = next(
    r
    for r in summary_rows
    if r["case"] == "fine"
)


def close_to_interface(row):

    first_ok = (
        row[
            "first_yield_peak_adjacent"
        ] == 1
        or (
            np.isfinite(
                row[
                    "first_yield_d_over_dx"
                ]
            )
            and row[
                "first_yield_d_over_dx"
            ] <= 1.5
        )
    )


    last_ok = (
        row[
            "last_peak_adjacent"
        ] == 1
        or (
            np.isfinite(
                row[
                    "last_d_over_dx"
                ]
            )
            and row[
                "last_d_over_dx"
            ] <= 1.5
        )
    )


    return (
        first_ok
        and last_ok
    )


medium_local = close_to_interface(
    medium
)

fine_local = close_to_interface(
    fine
)


print()
print(
    "========================================"
)

print(
    "LOCALIZATION DECISION"
)

print(
    "========================================"
)

print(
    "Medium first-yield d/dx: "
    f"{medium['first_yield_d_over_dx']}"
)

print(
    "Fine first-yield d/dx: "
    f"{fine['first_yield_d_over_dx']}"
)

print(
    "Medium last-valid d/dx: "
    f"{medium['last_d_over_dx']}"
)

print(
    "Fine last-valid d/dx: "
    f"{fine['last_d_over_dx']}"
)


if (
    medium_local
    and fine_local
):

    decision = (
        "NOTCH-FRONT LOCALIZATION SUPPORTED"
    )

    print()
    print(
        decision
    )

    print(
        "Plastic-strain hotspots remain "
        "within approximately 1.5 local toe "
        "cell widths of the active/eroded "
        "interface on medium and fine meshes."
    )

else:

    decision = (
        "NOTCH-FRONT LOCALIZATION "
        "NOT YET ESTABLISHED"
    )

    print()
    print(
        decision
    )


# ============================================================
# SUMMARY CSV
# ============================================================

summary_csv = (
    OUT
    / "mesh_localization_summary.csv"
)


with summary_csv.open(
    "w",
    newline="",
    encoding="utf-8",
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=list(
            summary_rows[0].keys()
        ),
    )

    writer.writeheader()
    writer.writerows(
        summary_rows
    )


# ============================================================
# FIGURE 01:
# hotspot distance normalized by local dx
# ============================================================

fig, ax = plt.subplots(
    figsize=(7.5, 4.8)
)


for result in results:

    rec = [
        r
        for r in result[
            "records"
        ]
        if (
            r["epsp_max"] > 1e-8
            and np.isfinite(
                r[
                    "peak_distance_to_interface_over_dx"
                ]
            )
        )
    ]

    if not rec:
        continue


    ax.plot(
        [
            r["E_m"]
            for r in rec
        ],
        [
            r[
                "peak_distance_to_interface_over_dx"
            ]
            for r in rec
        ],
        marker="o",
        label=result[
            "case"
        ],
    )


ax.axhline(
    1.5,
    linestyle="--",
)

ax.set_xlabel(
    "Toe recession E [m]"
)

ax.set_ylabel(
    "Hotspot distance / local toe dx"
)

ax.set_title(
    "Plastic hotspot proximity to erosion interface"
)

ax.legend()

ax.grid(
    alpha=0.25
)

fig.tight_layout()

fig.savefig(
    OUT
    / "figure_01_hotspot_interface_distance.png",
    dpi=220,
)

plt.close(
    fig
)


# ============================================================
# FIGURE 02:
# integrated plastic area
# ============================================================

fig, ax = plt.subplots(
    figsize=(7.5, 4.8)
)


for result in results:

    rec = result[
        "records"
    ]

    ax.plot(
        [
            r["E_m"]
            for r in rec
        ],
        [
            r[
                "plastic_area_gt_1e4_m2"
            ]
            for r in rec
        ],
        marker="o",
        label=result[
            "case"
        ],
    )


ax.set_xlabel(
    "Toe recession E [m]"
)

ax.set_ylabel(
    "Active-domain area with epsp > 1e-4 [m²]"
)

ax.set_title(
    "Integrated plastic-zone response"
)

ax.legend()

ax.grid(
    alpha=0.25
)

fig.tight_layout()

fig.savefig(
    OUT
    / "figure_02_integrated_plastic_area.png",
    dpi=220,
)

plt.close(
    fig
)


# ============================================================
# TEXT SUMMARY
# ============================================================

summary = (
    OUT
    / "phase05d3_summary.txt"
)


lines = [
    "PHASE 05D-3 LOCALIZATION AUDIT",
    "",
]


for row in summary_rows:

    lines += [
        (
            f"{row['case']} | "
            f"dx={row['toe_dx_m']:.6f} m | "
            f"yield_E={row['first_yield_E_m']} | "
            f"yield_d/dx="
            f"{row['first_yield_d_over_dx']} | "
            f"last_E={row['last_valid_E_m']} | "
            f"last_d/dx="
            f"{row['last_d_over_dx']} | "
            f"A_epsp>1e-4="
            f"{row['last_plastic_area_1e4_m2']}"
        ),
        "",
    ]


lines += [
    (
        "DECISION: "
        f"{decision}"
    ),
    "",
    (
        "Diagnostic interpretation only: "
        "distance to the erosion interface "
        "is not a physical failure criterion."
    ),
    (
        "If localization follows the notch "
        "front across refined meshes, "
        "solver nonconvergence should not "
        "be used as Ecrit."
    ),
]


summary.write_text(
    "\n".join(lines) + "\n",
    encoding="utf-8",
)


print()
print(
    "PASS:",
    history_csv,
)

print(
    "PASS:",
    summary_csv,
)

print(
    "PASS:",
    OUT
    / "figure_01_hotspot_interface_distance.png"
)

print(
    "PASS:",
    OUT
    / "figure_02_integrated_plastic_area.png"
)

print(
    "PASS:",
    summary,
)

print()
print(
    "PHASE 05D-3 ANALYSIS COMPLETE"
)
