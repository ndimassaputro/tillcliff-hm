from pathlib import Path
import csv
import shutil

import numpy as np
import pyvista as pv


ROOT = Path.cwd()

MESH_PATH = (
    ROOT
    / "model"
    / "phase05a_v2_toe_notch"
    / "slope_toe_notch_ready.vtu"
)

MODEL_DIR = (
    ROOT
    / "model"
    / "phase05a_v2_toe_notch"
)

OUT = (
    ROOT
    / "results"
    / "phase05d_preflight"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


if not MESH_PATH.exists():
    raise SystemExit(
        f"FAIL: missing {MESH_PATH}"
    )


mesh = pv.read(
    MESH_PATH
)


# ============================================================
# BASIC TOPOLOGY
# ============================================================

print(
    "========================================"
)

print(
    "PHASE 05D-0 MESH REFINEMENT PREFLIGHT"
)

print(
    "========================================"
)

print()

print(
    "=== BULK MESH ==="
)

print(
    f"Points : {mesh.n_points}"
)

print(
    f"Cells  : {mesh.n_cells}"
)

print(
    "Bounds :"
)

print(
    f"  x = {mesh.bounds[0]:.8f} "
    f"to {mesh.bounds[1]:.8f}"
)

print(
    f"  y = {mesh.bounds[2]:.8f} "
    f"to {mesh.bounds[3]:.8f}"
)

print(
    f"  z = {mesh.bounds[4]:.8f} "
    f"to {mesh.bounds[5]:.8f}"
)


# ============================================================
# CELL TYPES
# ============================================================

vtk_names = {
    3: "LINE",
    5: "TRIANGLE",
    8: "PIXEL",
    9: "QUAD",
    10: "TETRA",
    12: "HEXAHEDRON",
}

celltypes = np.asarray(
    mesh.celltypes,
    dtype=int,
)

unique_types, counts = np.unique(
    celltypes,
    return_counts=True,
)

print()
print(
    "=== CELL TYPES ==="
)

for t, n in zip(
    unique_types,
    counts,
):

    print(
        f"VTK type {t:2d} "
        f"({vtk_names.get(int(t), 'OTHER')}): "
        f"{n}"
    )


# ============================================================
# ARRAYS / RESTART STATE
# ============================================================

print()
print(
    "=== BULK DATA ARRAYS ==="
)

print(
    "Point arrays:"
)

for key in mesh.point_data.keys():

    arr = np.asarray(
        mesh.point_data[key]
    )

    print(
        f"  {key}: "
        f"shape={arr.shape}"
    )


print(
    "Cell arrays:"
)

for key in mesh.cell_data.keys():

    arr = np.asarray(
        mesh.cell_data[key]
    )

    print(
        f"  {key}: "
        f"shape={arr.shape}"
    )


required_checks = {
    "MaterialIDs":
        "cell",

    "p0":
        "point",

    "sigma0":
        "cell",
}


print()
print(
    "=== REQUIRED RESTART FIELDS ==="
)

restart_ok = True

for name, location in (
    required_checks.items()
):

    if location == "point":

        present = (
            name
            in mesh.point_data
        )

    else:

        present = (
            name
            in mesh.cell_data
        )

    print(
        f"{name:12s} | "
        f"{location:5s} | "
        f"{'PRESENT' if present else 'MISSING'}"
    )

    if not present:

        restart_ok = False


# ============================================================
# STRUCTURED-COLUMN AUDIT
#
# If every x-column contains about the same number
# of points, the current slope mesh is likely based
# on a structured logical grid warped to the slope.
# ============================================================

xyz = np.asarray(
    mesh.points,
    dtype=float,
)

x_round = np.round(
    xyz[:, 0],
    8,
)

unique_x = np.unique(
    x_round
)

column_counts = []

for xv in unique_x:

    column_counts.append(
        int(
            np.sum(
                x_round == xv
            )
        )
    )

column_counts = np.asarray(
    column_counts,
    dtype=int,
)


print()
print(
    "=== STRUCTURED-COLUMN AUDIT ==="
)

print(
    f"Unique x-columns: "
    f"{len(unique_x)}"
)

print(
    "Points per x-column:"
)

print(
    f"  min    = "
    f"{column_counts.min()}"
)

print(
    f"  median = "
    f"{np.median(column_counts):.1f}"
)

print(
    f"  max    = "
    f"{column_counts.max()}"
)

uniform_columns = (
    column_counts.min()
    == column_counts.max()
)

print(
    "Uniform point count per column: "
    f"{'YES' if uniform_columns else 'NO'}"
)


dx = np.diff(
    unique_x
)

if len(dx):

    print(
        "x-column spacing [m]:"
    )

    print(
        f"  min    = "
        f"{np.min(dx):.8f}"
    )

    print(
        f"  median = "
        f"{np.median(dx):.8f}"
    )

    print(
        f"  max    = "
        f"{np.max(dx):.8f}"
    )


# ============================================================
# CELL AREA AUDIT
# ============================================================

sizes = mesh.compute_cell_sizes(
    length=False,
    area=True,
    volume=False,
)

areas = np.asarray(
    sizes.cell_data["Area"],
    dtype=float,
)

centers = (
    mesh
    .cell_centers()
    .points
)


print()
print(
    "=== CELL AREA ==="
)

print(
    "Whole domain [m2]:"
)

print(
    f"  min    = "
    f"{np.min(areas):.10f}"
)

print(
    f"  median = "
    f"{np.median(areas):.10f}"
)

print(
    f"  max    = "
    f"{np.max(areas):.10f}"
)


toe_region = (
    (centers[:, 0] >= 19.0)
    & (centers[:, 0] <= 22.2)
    & (centers[:, 1] >= 1.5)
    & (centers[:, 1] <= 3.5)
)

toe_areas = areas[
    toe_region
]

print(
    "Toe audit region [m2]:"
)

print(
    f"  cells  = "
    f"{len(toe_areas)}"
)

if len(toe_areas):

    print(
        f"  min    = "
        f"{np.min(toe_areas):.10f}"
    )

    print(
        f"  median = "
        f"{np.median(toe_areas):.10f}"
    )

    print(
        f"  max    = "
        f"{np.max(toe_areas):.10f}"
    )


# ============================================================
# DEACTIVATION EVENT AUDIT
#
# These are the exact geometric events seen by
# the moving-front erosion implementation.
# ============================================================

if "MaterialIDs" not in mesh.cell_data:

    raise SystemExit(
        "FAIL: MaterialIDs missing"
    )


material_ids = np.asarray(
    mesh.cell_data[
        "MaterialIDs"
    ],
    dtype=int,
)

candidate = (
    material_ids > 0
)

TOE_X = 22.0

E_cell = (
    TOE_X
    - centers[:, 0]
)


event_rows = []

event_values = np.unique(
    np.round(
        E_cell[
            candidate
            & (E_cell >= 0.0)
            & (E_cell <= 0.8)
        ],
        8,
    )
)


print()
print(
    "=== DISCRETE EROSION EVENTS "
    "E=0–0.8 m ==="
)

for E in event_values:

    mask = (
        candidate
        & np.isclose(
            E_cell,
            E,
            atol=1e-7,
        )
    )

    n_cells = int(
        np.sum(mask)
    )

    area = float(
        np.sum(
            areas[mask]
        )
    )

    mids = sorted(
        set(
            material_ids[
                mask
            ].tolist()
        )
    )

    print(
        f"E={E:.8f} m | "
        f"cells={n_cells:3d} | "
        f"area={area:.10f} m2 | "
        f"IDs={mids}"
    )

    event_rows.append(
        {
            "E_m":
                float(E),

            "cell_count":
                n_cells,

            "area_m2":
                area,

            "material_ids":
                " ".join(
                    str(v)
                    for v in mids
                ),
        }
    )


with (
    OUT
    / "erosion_event_audit.csv"
).open(
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
            "material_ids",
        ],
    )

    writer.writeheader()
    writer.writerows(
        event_rows
    )


# ============================================================
# BOUNDARY MESH AUDIT
# ============================================================

print()
print(
    "=== BOUNDARY MESHES ==="
)


boundary_names = [
    "slope_left.vtu",
    "slope_right.vtu",
    "slope_top.vtu",
    "slope_bottom.vtu",
]


boundary_ok = True

for name in boundary_names:

    path = (
        MODEL_DIR
        / name
    )

    if not path.exists():

        print(
            f"{name}: MISSING"
        )

        boundary_ok = False
        continue


    b = pv.read(
        path
    )

    print(
        f"{name}: "
        f"points={b.n_points}, "
        f"cells={b.n_cells}"
    )

    print(
        "  point arrays: "
        + (
            ", ".join(
                b.point_data.keys()
            )
            if len(
                b.point_data.keys()
            )
            else "(none)"
        )
    )

    print(
        "  cell arrays : "
        + (
            ", ".join(
                b.cell_data.keys()
            )
            if len(
                b.cell_data.keys()
            )
            else "(none)"
        )
    )


# ============================================================
# TOOL AVAILABILITY
# ============================================================

tools = [
    "generateStructuredMesh",
    "ExtractBoundary",
    "identifySubdomains",
    "GMSH2OGS",
    "msh2vtu",
    "gmsh",
]


print()
print(
    "=== AVAILABLE MESH TOOLS ==="
)

tool_rows = []

for tool in tools:

    path = shutil.which(
        tool
    )

    available = (
        path is not None
    )

    print(
        f"{tool:24s}: "
        f"{'YES' if available else 'NO'}"
        + (
            f" | {path}"
            if available
            else ""
        )
    )

    tool_rows.append(
        {
            "tool":
                tool,

            "available":
                int(
                    available
                ),

            "path":
                path or "",
        }
    )


with (
    OUT
    / "tool_availability.csv"
).open(
    "w",
    newline="",
    encoding="utf-8",
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "tool",
            "available",
            "path",
        ],
    )

    writer.writeheader()
    writer.writerows(
        tool_rows
    )


# ============================================================
# CURRENT FAILURE EVENT
# ============================================================

failure_event = [
    r
    for r in event_rows
    if abs(
        r["E_m"] - 0.5
    ) < 1e-6
]


print()
print(
    "=== CURRENT E=0.50 EVENT ==="
)

if failure_event:

    r = failure_event[0]

    print(
        f"Cells released : "
        f"{r['cell_count']}"
    )

    print(
        f"Area released  : "
        f"{r['area_m2']:.10f} m2"
    )

else:

    print(
        "No event exactly at E=0.50 m"
    )


# ============================================================
# PREFLIGHT CLASSIFICATION
# ============================================================

print()
print(
    "=== PREFLIGHT STATUS ==="
)


if (
    restart_ok
    and boundary_ok
):

    print(
        "RESTART / BOUNDARY DATA: PASS"
    )

else:

    print(
        "RESTART / BOUNDARY DATA: REVIEW"
    )


if uniform_columns:

    print(
        "MESH STRUCTURE: "
        "STRUCTURED-COLUMN CANDIDATE"
    )

else:

    print(
        "MESH STRUCTURE: "
        "GENERAL UNSTRUCTURED / WARPED"
    )


# ============================================================
# SUMMARY
# ============================================================

summary = (
    OUT
    / "phase05d_preflight_summary.txt"
)


available_tools = [
    r["tool"]
    for r in tool_rows
    if r["available"]
]


lines = [
    "PHASE 05D-0 MESH REFINEMENT PREFLIGHT",
    "",
    f"Points: {mesh.n_points}",
    f"Cells: {mesh.n_cells}",
    (
        "Unique x-columns: "
        f"{len(unique_x)}"
    ),
    (
        "Column count range: "
        f"{column_counts.min()} "
        f"to {column_counts.max()}"
    ),
    (
        "Uniform columns: "
        f"{uniform_columns}"
    ),
    "",
    (
        "Restart fields present: "
        f"{restart_ok}"
    ),
    (
        "Boundary meshes present: "
        f"{boundary_ok}"
    ),
    "",
    (
        "Available tools: "
        + (
            ", ".join(
                available_tools
            )
            if available_tools
            else "NONE"
        )
    ),
    "",
]


if failure_event:

    r = failure_event[0]

    lines += [
        (
            "E=0.50 event cells: "
            f"{r['cell_count']}"
        ),
        (
            "E=0.50 event area [m2]: "
            f"{r['area_m2']:.10f}"
        ),
        "",
    ]


lines += [
    (
        "NEXT DECISION:"
    ),
    (
        "Use this audit to choose the safest "
        "coarse/medium/fine mesh-generation "
        "and restart-transfer route."
    ),
]


summary.write_text(
    "\n".join(lines) + "\n",
    encoding="utf-8",
)


print()
print(
    "PASS:",
    OUT
    / "erosion_event_audit.csv"
)

print(
    "PASS:",
    OUT
    / "tool_availability.csv"
)

print(
    "PASS:",
    summary
)

print()
print(
    "PHASE 05D-0 PREFLIGHT COMPLETE"
)
