from pathlib import Path
import csv
import shutil
import xml.etree.ElementTree as ET

import numpy as np
import pyvista as pv


ROOT = Path.cwd()

SRC_MESH = (
    ROOT
    / "model"
    / "phase05a_v2_toe_notch"
    / "slope_toe_notch_ready.vtu"
)

SRC_PRJ = (
    ROOT
    / "model"
    / "phase04d_v2"
    / "screening.prj"
)

MODEL_ROOT = (
    ROOT
    / "model"
    / "phase05d1_restart_mesh"
)

OUT = (
    ROOT
    / "results"
    / "phase05d1_analysis"
)

MODEL_ROOT.mkdir(
    parents=True,
    exist_ok=True,
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


for path in [
    SRC_MESH,
    SRC_PRJ,
]:
    if not path.exists():
        raise SystemExit(
            f"FAIL: missing {path}"
        )


# ============================================================
# SETTINGS
#
# Horizontal refinement only.
#
# The problematic erosion events occur in the coastal toe.
# We therefore refine x only over 19.8 <= x <= 22.2 m.
#
# Vertical eta levels are held fixed.
# ============================================================

REFINE_X_MIN = 19.8
REFINE_X_MAX = 22.2

CASES = {
    "coarse": 1,
    "medium": 2,
    "fine": 4,
}

TOE_X = 22.0
TOE_Y = 2.0
NOTCH_TOP = 3.2

BAND_WIDTH = 0.4
N_BANDS = 5


# ============================================================
# LOAD COARSE VERIFIED RESTART MESH
# ============================================================

src = pv.read(
    SRC_MESH
)

if "p0" not in src.point_data:
    raise SystemExit(
        "FAIL: p0 missing"
    )

if "sigma0" not in src.cell_data:
    raise SystemExit(
        "FAIL: sigma0 missing"
    )

if "MaterialIDs" not in src.cell_data:
    raise SystemExit(
        "FAIL: MaterialIDs missing"
    )


src_points = np.asarray(
    src.points,
    dtype=float,
)

src_centers = (
    src
    .cell_centers()
    .points
)


# ============================================================
# RECONSTRUCT STRUCTURED LOGICAL GRID
# ============================================================

x0 = np.unique(
    np.round(
        src_points[:, 0],
        10,
    )
)

x0.sort()

NX0 = len(x0) - 1


top0 = np.empty(
    len(x0),
    dtype=float,
)

point_indices_by_column = []


for i, xv in enumerate(x0):

    ids = np.where(
        np.isclose(
            src_points[:, 0],
            xv,
            atol=1e-9,
        )
    )[0]

    ids = ids[
        np.argsort(
            src_points[
                ids,
                1,
            ]
        )
    ]

    point_indices_by_column.append(
        ids
    )

    top0[i] = float(
        np.max(
            src_points[
                ids,
                1,
            ]
        )
    )


column_counts = [
    len(ids)
    for ids in point_indices_by_column
]

if len(set(column_counts)) != 1:
    raise SystemExit(
        "FAIL: nonuniform source column counts"
    )


NETA = column_counts[0]
NY = NETA - 1


first_ids = (
    point_indices_by_column[0]
)

eta = (
    src_points[
        first_ids,
        1,
    ]
    / top0[0]
)


if not np.all(
    np.diff(eta) > 0
):
    raise SystemExit(
        "FAIL: eta levels not monotonic"
    )


print(
    "SOURCE LOGICAL GRID"
)

print(
    f"NX={NX0}, NY={NY}"
)

print(
    f"x columns={len(x0)}"
)

print(
    f"eta levels={len(eta)}"
)


# ============================================================
# BUILD p0 LOGICAL GRID
# ============================================================

p0_src = np.asarray(
    src.point_data["p0"],
    dtype=float,
).squeeze()

p0_grid = np.empty(
    (
        len(x0),
        NETA,
    ),
    dtype=float,
)


for i, ids in enumerate(
    point_indices_by_column
):

    p0_grid[
        i,
        :,
    ] = p0_src[
        ids
    ]


# ============================================================
# BUILD sigma0 LOGICAL CELL GRID
#
# We do NOT assume original cell ordering.
# Each original cell is located from its centre.
# ============================================================

sigma_src = np.asarray(
    src.cell_data["sigma0"],
    dtype=float,
)

if sigma_src.ndim == 1:
    sigma_src = sigma_src[:, None]


NCOMP_SIGMA = sigma_src.shape[1]

sigma_grid = np.full(
    (
        NX0,
        NY,
        NCOMP_SIGMA,
    ),
    np.nan,
    dtype=float,
)


for cell_id, centre in enumerate(
    src_centers
):

    cx = float(
        centre[0]
    )

    cy = float(
        centre[1]
    )

    i = (
        np.searchsorted(
            x0,
            cx,
            side="right",
        )
        - 1
    )

    i = int(
        np.clip(
            i,
            0,
            NX0 - 1,
        )
    )

    h = float(
        np.interp(
            cx,
            x0,
            top0,
        )
    )

    eta_c = (
        cy / h
    )

    j = (
        np.searchsorted(
            eta,
            eta_c,
            side="right",
        )
        - 1
    )

    j = int(
        np.clip(
            j,
            0,
            NY - 1,
        )
    )

    if np.all(
        np.isfinite(
            sigma_grid[
                i,
                j,
                :,
            ]
        )
    ):
        raise SystemExit(
            "FAIL: duplicate logical sigma cell "
            f"at i={i}, j={j}"
        )

    sigma_grid[
        i,
        j,
        :,
    ] = sigma_src[
        cell_id,
        :,
    ]


if not np.all(
    np.isfinite(
        sigma_grid
    )
):
    raise SystemExit(
        "FAIL: incomplete logical sigma grid"
    )


# ============================================================
# REFINEMENT HELPER
# ============================================================

def make_refined_x(
    factor,
):

    values = []

    for a, b in zip(
        x0[:-1],
        x0[1:],
    ):

        midpoint = (
            0.5
            * (
                a + b
            )
        )

        if (
            midpoint
            >= REFINE_X_MIN
            - 1e-10
            and midpoint
            <= REFINE_X_MAX
            + 1e-10
        ):
            f = factor
        else:
            f = 1

        for k in range(f):

            values.append(
                a
                + (
                    b - a
                )
                * k
                / f
            )

    values.append(
        x0[-1]
    )

    result = np.asarray(
        values,
        dtype=float,
    )

    if not np.all(
        np.diff(result) > 0
    ):
        raise SystemExit(
            "FAIL: refined x not monotonic"
        )

    return result


# ============================================================
# PROJECT XML HELPER
# ============================================================

def build_project(
    case_dir,
    case,
):

    tree = ET.parse(
        SRC_PRJ
    )

    root = tree.getroot()


    # --------------------------------------------------------
    # Bulk mesh filename.
    # --------------------------------------------------------

    meshes = root.find(
        "meshes"
    )

    if meshes is None:
        raise SystemExit(
            "FAIL: meshes block missing"
        )

    replaced = False

    for m in meshes.findall(
        "mesh"
    ):

        text = (
            m.text.strip()
            if m.text
            else ""
        )

        if text in {
            "slope_initialized.vtu",
            "slope.vtu",
        }:

            m.text = "bulk.vtu"

            replaced = True

    if not replaced:
        raise SystemExit(
            f"FAIL: bulk mesh replacement "
            f"failed for {case}"
        )


    # --------------------------------------------------------
    # Constitutive relation valid for all erosion IDs.
    # --------------------------------------------------------

    process = None

    for p in root.findall(
        "./processes/process"
    ):

        if (
            p.findtext("type")
            == "RICHARDS_MECHANICS"
        ):
            process = p
            break

    if process is None:
        raise SystemExit(
            "FAIL: RM process missing"
        )

    cr = process.find(
        "constitutive_relation"
    )

    if cr is None:
        raise SystemExit(
            "FAIL: constitutive relation missing"
        )

    cr.set(
        "id",
        "0,1,2,3,4,5",
    )


    medium = root.find(
        "./media/medium"
    )

    if medium is None:
        raise SystemExit(
            "FAIL: medium missing"
        )

    medium.set(
        "id",
        "0,1,2,3,4,5",
    )


    # --------------------------------------------------------
    # Absolutely no erosion during restart transfer test.
    # --------------------------------------------------------

    root_pvs = root.find(
        "process_variables"
    )

    if root_pvs is None:
        raise SystemExit(
            "FAIL: root process variables missing"
        )

    for pv in root_pvs.findall(
        "process_variable"
    ):

        old = pv.find(
            "deactivated_subdomains"
        )

        if old is not None:
            pv.remove(
                old
            )


    # --------------------------------------------------------
    # Short numerical hold.
    #
    # This tests whether transferred p0/sigma0
    # are accepted by the refined discretization.
    # --------------------------------------------------------

    ts = root.find(
        "./time_loop/processes/"
        "process/time_stepping"
    )

    if ts is None:
        raise SystemExit(
            "FAIL: time stepping missing"
        )

    ts.clear()

    ET.SubElement(
        ts,
        "type",
    ).text = "FixedTimeStepping"

    ET.SubElement(
        ts,
        "t_initial",
    ).text = "0"

    ET.SubElement(
        ts,
        "t_end",
    ).text = "20"

    steps = ET.SubElement(
        ts,
        "timesteps",
    )

    pair = ET.SubElement(
        steps,
        "pair",
    )

    ET.SubElement(
        pair,
        "repeat",
    ).text = "20"

    ET.SubElement(
        pair,
        "delta_t",
    ).text = "1"


    # --------------------------------------------------------
    # Newton.
    # --------------------------------------------------------

    solver = root.find(
        "./nonlinear_solvers/"
        "nonlinear_solver"
    )

    if solver is None:
        raise SystemExit(
            "FAIL: nonlinear solver missing"
        )

    max_iter = solver.find(
        "max_iter"
    )

    if max_iter is None:

        max_iter = ET.SubElement(
            solver,
            "max_iter",
        )

    max_iter.text = "80"


    # --------------------------------------------------------
    # Output every step.
    # --------------------------------------------------------

    prefix = root.find(
        "./time_loop/output/prefix"
    )

    if prefix is None:
        raise SystemExit(
            "FAIL: output prefix missing"
        )

    prefix.text = (
        f"phase05d1_{case}"
    )


    out_steps = root.find(
        "./time_loop/output/timesteps"
    )

    if out_steps is None:
        raise SystemExit(
            "FAIL: output timesteps missing"
        )

    out_steps.clear()

    pair = ET.SubElement(
        out_steps,
        "pair",
    )

    ET.SubElement(
        pair,
        "repeat",
    ).text = "20"

    ET.SubElement(
        pair,
        "each_steps",
    ).text = "1"


    variables = root.find(
        "./time_loop/output/variables"
    )

    existing = {
        v.text
        for v in variables.findall(
            "variable"
        )
    }

    for name in [
        "displacement",
        "pressure",
        "saturation",
        "sigma",
        "EquivalentPlasticStrain",
    ]:

        if name not in existing:

            ET.SubElement(
                variables,
                "variable",
            ).text = name


    ET.indent(
        tree,
        space="    ",
    )

    dst = (
        case_dir
        / "restart_check.prj"
    )

    tree.write(
        dst,
        encoding="UTF-8",
        xml_declaration=True,
    )

    ET.parse(
        dst
    )


# ============================================================
# BUILD EACH MESH
# ============================================================

summary_rows = []


for case, factor in CASES.items():

    print()
    print(
        "========================================"
    )

    print(
        f"BUILD {case.upper()}"
    )

    print(
        "========================================"
    )


    case_dir = (
        MODEL_ROOT
        / case
    )

    case_dir.mkdir(
        parents=True,
        exist_ok=True,
    )


    x = make_refined_x(
        factor
    )

    nx = (
        len(x) - 1
    )


    top = np.interp(
        x,
        x0,
        top0,
    )


    # --------------------------------------------------------
    # POINTS
    # --------------------------------------------------------

    points = np.empty(
        (
            len(x)
            * NETA,
            3,
        ),
        dtype=float,
    )

    for i, xv in enumerate(x):

        h = top[i]

        for j, eta_j in enumerate(
            eta
        ):

            pid = (
                i * NETA
                + j
            )

            points[
                pid,
                :,
            ] = [
                xv,
                eta_j * h,
                0.0,
            ]


    # --------------------------------------------------------
    # QUAD CONNECTIVITY
    # --------------------------------------------------------

    ncells = (
        nx
        * NY
    )

    cells = np.empty(
        (
            ncells,
            5,
        ),
        dtype=np.int64,
    )

    cells[:, 0] = 4

    c = 0

    for i in range(nx):

        for j in range(NY):

            p00 = (
                i * NETA
                + j
            )

            p10 = (
                (i + 1)
                * NETA
                + j
            )

            p11 = (
                (i + 1)
                * NETA
                + j
                + 1
            )

            p01 = (
                i * NETA
                + j
                + 1
            )

            cells[
                c,
                1:,
            ] = [
                p00,
                p10,
                p11,
                p01,
            ]

            c += 1


    celltypes = np.full(
        ncells,
        9,
        dtype=np.uint8,
    )


    grid = pv.UnstructuredGrid(
        cells.ravel(),
        celltypes,
        points,
    )


    # --------------------------------------------------------
    # p0:
    #
    # eta levels unchanged.
    # Interpolate only along x.
    # --------------------------------------------------------

    p0_new = np.empty(
        (
            len(x),
            NETA,
        ),
        dtype=float,
    )

    for j in range(NETA):

        p0_new[
            :,
            j,
        ] = np.interp(
            x,
            x0,
            p0_grid[
                :,
                j,
            ],
        )


    grid.point_data[
        "p0"
    ] = p0_new.ravel(
        order="C"
    )


    # --------------------------------------------------------
    # sigma0:
    #
    # Piecewise-constant parent-cell transfer.
    #
    # This avoids artificial smoothing of the stress
    # tensor during restart interpolation.
    # --------------------------------------------------------

    sigma_new = np.empty(
        (
            ncells,
            NCOMP_SIGMA,
        ),
        dtype=float,
    )

    c = 0

    for i in range(nx):

        xmid = (
            0.5
            * (
                x[i]
                + x[i + 1]
            )
        )

        parent_i = (
            np.searchsorted(
                x0,
                xmid,
                side="right",
            )
            - 1
        )

        parent_i = int(
            np.clip(
                parent_i,
                0,
                NX0 - 1,
            )
        )

        for j in range(NY):

            sigma_new[
                c,
                :,
            ] = sigma_grid[
                parent_i,
                j,
                :,
            ]

            c += 1


    grid.cell_data[
        "sigma0"
    ] = sigma_new


    # --------------------------------------------------------
    # MATERIAL IDS:
    #
    # Same physical notch window.
    # --------------------------------------------------------

    centres = (
        grid
        .cell_centers()
        .points
    )

    cx = centres[:, 0]
    cy = centres[:, 1]

    E_cell = (
        TOE_X
        - cx
    )

    candidate = (
        (E_cell >= 0.0)
        & (
            E_cell
            < N_BANDS
            * BAND_WIDTH
        )
        & (cy >= TOE_Y)
        & (cy <= NOTCH_TOP)
    )


    material_ids = np.zeros(
        grid.n_cells,
        dtype=np.int32,
    )

    band = (
        np.floor(
            E_cell
            / BAND_WIDTH
        )
        .astype(int)
        + 1
    )

    valid = (
        candidate
        & (band >= 1)
        & (band <= N_BANDS)
    )

    material_ids[
        valid
    ] = band[
        valid
    ]


    grid.cell_data[
        "MaterialIDs"
    ] = material_ids


    # --------------------------------------------------------
    # SAVE
    # --------------------------------------------------------

    bulk_path = (
        case_dir
        / "bulk.vtu"
    )

    grid.save(
        bulk_path
    )


    # --------------------------------------------------------
    # EVENT AUDIT
    # --------------------------------------------------------

    sizes = grid.compute_cell_sizes(
        length=False,
        area=True,
        volume=False,
    )

    area = np.asarray(
        sizes.cell_data[
            "Area"
        ],
        dtype=float,
    )


    candidate_mask = (
        material_ids > 0
    )

    event_E = np.unique(
        np.round(
            E_cell[
                candidate_mask
                & (E_cell >= 0.0)
                & (E_cell <= 0.8)
            ],
            8,
        )
    )


    event_file = (
        case_dir
        / "erosion_events.csv"
    )

    event_rows = []


    for E in event_E:

        mask = (
            candidate_mask
            & np.isclose(
                E_cell,
                E,
                atol=1e-7,
            )
        )

        event_rows.append(
            {
                "E_m":
                    float(E),

                "cell_count":
                    int(
                        np.sum(mask)
                    ),

                "area_m2":
                    float(
                        np.sum(
                            area[mask]
                        )
                    ),
            }
        )


    with event_file.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=[
                "E_m",
                "cell_count",
                "area_m2",
            ],
        )

        writer.writeheader()
        writer.writerows(
            event_rows
        )


    # --------------------------------------------------------
    # Event nearest E=0.50.
    # --------------------------------------------------------

    if event_rows:

        near = min(
            event_rows,
            key=lambda r:
                abs(
                    r["E_m"]
                    - 0.5
                ),
        )

    else:

        near = {
            "E_m":
                np.nan,

            "cell_count":
                0,

            "area_m2":
                np.nan,
        }


    toe_dx = np.diff(
        x
    )

    toe_intervals = (
        (
            0.5
            * (
                x[:-1]
                + x[1:]
            )
            >= REFINE_X_MIN
        )
        & (
            0.5
            * (
                x[:-1]
                + x[1:]
            )
            <= REFINE_X_MAX
        )
    )


    summary_rows.append(
        {
            "case":
                case,

            "factor":
                factor,

            "points":
                grid.n_points,

            "cells":
                grid.n_cells,

            "toe_dx_min_m":
                float(
                    np.min(
                        toe_dx[
                            toe_intervals
                        ]
                    )
                ),

            "toe_dx_max_m":
                float(
                    np.max(
                        toe_dx[
                            toe_intervals
                        ]
                    )
                ),

            "nearest_event_E_m":
                near[
                    "E_m"
                ],

            "nearest_event_cells":
                near[
                    "cell_count"
                ],

            "nearest_event_area_m2":
                near[
                    "area_m2"
                ],
        }
    )


    print(
        f"points = {grid.n_points}"
    )

    print(
        f"cells  = {grid.n_cells}"
    )

    print(
        "toe dx [m] = "
        f"{np.min(toe_dx[toe_intervals]):.6f}"
    )

    print(
        "event nearest 0.50 m:"
    )

    print(
        f"  E     = "
        f"{near['E_m']:.6f} m"
    )

    print(
        f"  cells = "
        f"{near['cell_count']}"
    )

    print(
        f"  area  = "
        f"{near['area_m2']:.10f} m2"
    )


    # --------------------------------------------------------
    # Build intact restart test project.
    # --------------------------------------------------------

    build_project(
        case_dir,
        case,
    )


# ============================================================
# SUMMARY CSV
# ============================================================

summary_csv = (
    OUT
    / "mesh_refinement_summary.csv"
)

with summary_csv.open(
    "w",
    newline="",
    encoding="utf-8",
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "case",
            "factor",
            "points",
            "cells",
            "toe_dx_min_m",
            "toe_dx_max_m",
            "nearest_event_E_m",
            "nearest_event_cells",
            "nearest_event_area_m2",
        ],
    )

    writer.writeheader()
    writer.writerows(
        summary_rows
    )


print()
print(
    "========================================"
)

print(
    "PHASE 05D-1 MESH BUILD COMPLETE"
)

print(
    "========================================"
)

print(
    "PASS:",
    summary_csv,
)
