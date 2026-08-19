from pathlib import Path
import csv
import re

import numpy as np
import pyvista as pv


ROOT = Path.cwd()

RESULT_ROOT = (
    ROOT
    / "results"
    / "phase04b_states"
)

MODEL_MESH = (
    ROOT
    / "model"
    / "phase04b_states"
    / "slope.vtu"
)

OUT = (
    ROOT
    / "results"
    / "phase05e1_antecedent_audit"
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


# ============================================================
# HELPERS
# ============================================================

def time_from_name(path):

    match = re.search(
        r"_t_([0-9eE+.\-]+)\.vtu$",
        path.name,
    )

    if not match:

        return -np.inf

    return float(
        match.group(1)
    )


def final_vtu(state):

    folder = (
        RESULT_ROOT
        / state
    )

    files = list(
        folder.glob("*.vtu")
    )

    if not files:

        raise RuntimeError(
            f"No VTU files for {state}"
        )

    return max(
        files,
        key=time_from_name,
    )


def describe_array(
    mesh,
    name,
):

    if name in mesh.point_data:

        arr = np.asarray(
            mesh.point_data[name],
            dtype=float,
        )

        return (
            "point",
            arr,
        )

    if name in mesh.cell_data:

        arr = np.asarray(
            mesh.cell_data[name],
            dtype=float,
        )

        return (
            "cell",
            arr,
        )

    return (
        None,
        None,
    )


# ============================================================
# LOAD REFERENCE GEOMETRY
# ============================================================

if not MODEL_MESH.exists():

    raise RuntimeError(
        f"Missing {MODEL_MESH}"
    )

source_mesh = pv.read(
    MODEL_MESH
)

source_points = np.asarray(
    source_mesh.points,
    dtype=float,
)


# ============================================================
# AUDIT STATES
# ============================================================

records = []

state_meshes = {}


print(
    "========================================"
)

print(
    "PHASE 05E-1 ANTECEDENT STATE AUDIT"
)

print(
    "========================================"
)


for state in STATES:

    path = final_vtu(
        state
    )

    mesh = pv.read(
        path
    )

    state_meshes[
        state
    ] = mesh


    print()
    print(
        f"=== {state.upper()} ==="
    )

    print(
        "Final VTU:"
    )

    print(
        path
    )

    print(
        f"Time [s]: "
        f"{time_from_name(path):.6f}"
    )

    print(
        f"Points: {mesh.n_points}"
    )

    print(
        f"Cells : {mesh.n_cells}"
    )


    # --------------------------------------------------------
    # Geometry consistency with Phase04B slope mesh.
    # --------------------------------------------------------

    geometry_match = (
        mesh.n_points
        == source_mesh.n_points
        and mesh.n_cells
        == source_mesh.n_cells
    )


    if geometry_match:

        max_coord_diff = float(
            np.max(
                np.abs(
                    np.asarray(
                        mesh.points,
                        dtype=float,
                    )
                    - source_points
                )
            )
        )

    else:

        max_coord_diff = np.nan


    print(
        "Geometry dimensions match source: "
        f"{geometry_match}"
    )

    print(
        "Max coordinate difference [m]: "
        f"{max_coord_diff}"
    )


    # --------------------------------------------------------
    # Arrays.
    # --------------------------------------------------------

    print()
    print(
        "Point arrays:"
    )

    for name in mesh.point_data.keys():

        arr = np.asarray(
            mesh.point_data[name]
        )

        print(
            f"  {name}: "
            f"shape={arr.shape}"
        )


    print(
        "Cell arrays:"
    )

    for name in mesh.cell_data.keys():

        arr = np.asarray(
            mesh.cell_data[name]
        )

        print(
            f"  {name}: "
            f"shape={arr.shape}"
        )


    # --------------------------------------------------------
    # Expected coupled-HM fields.
    # --------------------------------------------------------

    candidates = {
        "pressure": [
            "pressure",
            "p",
            "liquid_pressure",
        ],
        "saturation": [
            "saturation",
            "Saturation",
            "liquid_saturation",
        ],
        "displacement": [
            "displacement",
        ],
        "sigma": [
            "sigma",
        ],
    }


    resolved = {}


    print()
    print(
        "Resolved fields:"
    )


    for logical_name, names in (
        candidates.items()
    ):

        found = None

        for name in names:

            location, arr = describe_array(
                mesh,
                name,
            )

            if location is not None:

                found = (
                    name,
                    location,
                    arr,
                )

                break


        resolved[
            logical_name
        ] = found


        if found is None:

            print(
                f"  {logical_name}: MISSING"
            )

            continue


        name, location, arr = found

        finite = arr[
            np.isfinite(
                arr
            )
        ]


        print(
            f"  {logical_name}: "
            f"{name} | "
            f"{location} | "
            f"shape={arr.shape}"
        )


        if finite.size:

            print(
                f"    min={np.min(finite):.10e}"
            )

            print(
                f"    mean={np.mean(finite):.10e}"
            )

            print(
                f"    max={np.max(finite):.10e}"
            )


    # --------------------------------------------------------
    # Pressure metrics.
    # --------------------------------------------------------

    pressure_entry = resolved[
        "pressure"
    ]


    if pressure_entry is not None:

        _, pressure_location, pressure = (
            pressure_entry
        )

        pressure = np.asarray(
            pressure,
            dtype=float,
        ).squeeze()

        p_finite = pressure[
            np.isfinite(
                pressure
            )
        ]

        p_min = float(
            np.min(
                p_finite
            )
        )

        p_mean = float(
            np.mean(
                p_finite
            )
        )

        p_max = float(
            np.max(
                p_finite
            )
        )

    else:

        pressure_location = ""
        p_min = np.nan
        p_mean = np.nan
        p_max = np.nan


    # --------------------------------------------------------
    # Saturation metrics.
    # --------------------------------------------------------

    sat_entry = resolved[
        "saturation"
    ]


    if sat_entry is not None:

        _, saturation_location, saturation = (
            sat_entry
        )

        saturation = np.asarray(
            saturation,
            dtype=float,
        ).squeeze()

        sat_finite = saturation[
            np.isfinite(
                saturation
            )
        ]

        sr_min = float(
            np.min(
                sat_finite
            )
        )

        sr_mean = float(
            np.mean(
                sat_finite
            )
        )

        sr_max = float(
            np.max(
                sat_finite
            )
        )

    else:

        saturation_location = ""
        sr_min = np.nan
        sr_mean = np.nan
        sr_max = np.nan


    records.append(
        {
            "state":
                state,

            "file":
                str(path),

            "time_s":
                time_from_name(
                    path
                ),

            "points":
                mesh.n_points,

            "cells":
                mesh.n_cells,

            "geometry_match":
                int(
                    geometry_match
                ),

            "max_coordinate_diff_m":
                max_coord_diff,

            "pressure_present":
                int(
                    pressure_entry
                    is not None
                ),

            "pressure_location":
                pressure_location,

            "pressure_min_pa":
                p_min,

            "pressure_mean_pa":
                p_mean,

            "pressure_max_pa":
                p_max,

            "saturation_present":
                int(
                    sat_entry
                    is not None
                ),

            "saturation_location":
                saturation_location,

            "saturation_min":
                sr_min,

            "saturation_mean":
                sr_mean,

            "saturation_max":
                sr_max,

            "displacement_present":
                int(
                    resolved[
                        "displacement"
                    ]
                    is not None
                ),

            "sigma_present":
                int(
                    resolved[
                        "sigma"
                    ]
                    is not None
                ),
        }
    )


# ============================================================
# CROSS-STATE GEOMETRY CHECK
# ============================================================

print()
print(
    "========================================"
)

print(
    "CROSS-STATE GEOMETRY"
)

print(
    "========================================"
)


reference = state_meshes[
    "reference"
]

reference_points = np.asarray(
    reference.points,
    dtype=float,
)


cross_geometry_ok = True


for state in STATES:

    mesh = state_meshes[
        state
    ]

    same_dims = (
        mesh.n_points
        == reference.n_points
        and mesh.n_cells
        == reference.n_cells
    )


    if same_dims:

        diff = float(
            np.max(
                np.abs(
                    np.asarray(
                        mesh.points,
                        dtype=float,
                    )
                    - reference_points
                )
            )
        )

    else:

        diff = np.nan


    ok = (
        same_dims
        and np.isfinite(
            diff
        )
        and diff <= 1e-10
    )


    cross_geometry_ok = (
        cross_geometry_ok
        and ok
    )


    print(
        f"{state:9s} | "
        f"same topology={same_dims} | "
        f"max coord diff="
        f"{diff}"
    )


# ============================================================
# HYDRAULIC ORDERING
# ============================================================

by_state = {
    r[
        "state"
    ]:
        r
    for r in records
}


dry_p = by_state[
    "dry"
][
    "pressure_mean_pa"
]

ref_p = by_state[
    "reference"
][
    "pressure_mean_pa"
]

wet_p = by_state[
    "wet"
][
    "pressure_mean_pa"
]


pressure_order_ok = (
    np.isfinite(
        dry_p
    )
    and np.isfinite(
        ref_p
    )
    and np.isfinite(
        wet_p
    )
    and dry_p
    < ref_p
    < wet_p
)


print()
print(
    "========================================"
)

print(
    "HYDRAULIC ORDERING"
)

print(
    "========================================"
)

print(
    "Mean pressure [kPa]:"
)

print(
    f"  DRY       = "
    f"{dry_p/1000:.6f}"
)

print(
    f"  REFERENCE = "
    f"{ref_p/1000:.6f}"
)

print(
    f"  WET       = "
    f"{wet_p/1000:.6f}"
)

print(
    "Expected dry < reference < wet: "
    f"{'PASS' if pressure_order_ok else 'REVIEW'}"
)


# ============================================================
# TRANSFER READINESS
# ============================================================

all_pressure_point = all(
    r[
        "pressure_present"
    ] == 1
    and r[
        "pressure_location"
    ] == "point"
    for r in records
)


transfer_ready = (
    cross_geometry_ok
    and all_pressure_point
    and pressure_order_ok
)


print()
print(
    "========================================"
)

print(
    "TRANSFER READINESS"
)

print(
    "========================================"
)

print(
    "Cross-state geometry identical: "
    f"{cross_geometry_ok}"
)

print(
    "Pressure available as point field "
    "for all states: "
    f"{all_pressure_point}"
)

print(
    "Hydraulic state ordering: "
    f"{pressure_order_ok}"
)

print()


if transfer_ready:

    status = (
        "PASS — READY FOR MEDIUM-MESH "
        "HYDRAULIC STATE TRANSFER"
    )

else:

    status = (
        "REVIEW — DO NOT TRANSFER YET"
    )


print(
    status
)


# ============================================================
# OUTPUT CSV
# ============================================================

csv_path = (
    OUT
    / "antecedent_state_audit.csv"
)


with csv_path.open(
    "w",
    newline="",
    encoding="utf-8",
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=list(
            records[0].keys()
        ),
    )

    writer.writeheader()
    writer.writerows(
        records
    )


# ============================================================
# SUMMARY
# ============================================================

summary = (
    OUT
    / "phase05e1_summary.txt"
)


lines = [
    "PHASE 05E-1 ANTECEDENT STATE AUDIT",
    "",
]


for r in records:

    lines.append(
        f"{r['state']} | "
        f"t={r['time_s']:.0f} s | "
        f"points={r['points']} | "
        f"cells={r['cells']} | "
        f"pressure={r['pressure_present']} "
        f"({r['pressure_location']}) | "
        f"pmean={r['pressure_mean_pa']/1000:.6f} kPa | "
        f"Srmean={r['saturation_mean']}"
    )


lines += [
    "",
    (
        "Cross-state geometry identical: "
        f"{cross_geometry_ok}"
    ),
    (
        "Pressure point-field available: "
        f"{all_pressure_point}"
    ),
    (
        "Hydraulic ordering pass: "
        f"{pressure_order_ok}"
    ),
    "",
    (
        "STATUS: "
        f"{status}"
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
    csv_path
)

print(
    "PASS:",
    summary
)

print()
print(
    "PHASE 05E-1 AUDIT COMPLETE"
)
