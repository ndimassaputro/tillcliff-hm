from pathlib import Path
import csv

import numpy as np
import ogstools as ot
import pyvista as pv


ROOT = Path.cwd()

RUN_ROOT = (
    ROOT
    / "results"
    / "phase05d1_restart_mesh"
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

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


CASES = [
    "coarse",
    "medium",
    "fine",
]


# ============================================================
# RUN STATUS
# ============================================================

run_status = {}

status_path = (
    RUN_ROOT
    / "run_status.csv"
)

if status_path.exists():

    with status_path.open(
        "r",
        encoding="utf-8",
    ) as f:

        reader = csv.DictReader(
            f
        )

        for row in reader:

            run_status[
                row["case"]
            ] = int(
                row["return_code"]
            )


def get_field(
    mesh,
    name,
):

    if name in mesh.point_data:

        return np.asarray(
            mesh.point_data[name],
            dtype=float,
        )

    if name in mesh.cell_data:

        return np.asarray(
            mesh.cell_data[name],
            dtype=float,
        )

    raise KeyError(
        f"{name!r} missing"
    )


records = []


print(
    "========================================"
)

print(
    "PHASE 05D-1 RESTART EQUILIBRIUM CHECK"
)

print(
    "========================================"
)


for case in CASES:

    print()
    print(
        f"=== {case.upper()} ==="
    )


    rc = run_status.get(
        case,
        999,
    )

    folder = (
        RUN_ROOT
        / case
    )

    pvds = sorted(
        folder.glob("*.pvd")
    )


    bulk = pv.read(
        MODEL_ROOT
        / case
        / "bulk.vtu"
    )


    if rc != 0 or not pvds:

        print(
            f"OGS return code: {rc}"
        )

        print(
            "STATUS: REVIEW"
        )

        records.append(
            {
                "case":
                    case,

                "return_code":
                    rc,

                "cells":
                    bulk.n_cells,

                "max_displacement_mm":
                    np.nan,

                "epsp_max":
                    np.nan,

                "pressure_min_kpa":
                    np.nan,

                "pressure_max_kpa":
                    np.nan,

                "status":
                    "REVIEW",
            }
        )

        continue


    series = ot.MeshSeries(
        str(
            pvds[0]
        )
    )

    times = np.asarray(
        series.timevalues,
        dtype=float,
    )

    final = series.mesh(
        len(times) - 1
    )


    u = get_field(
        final,
        "displacement",
    )

    umag = np.linalg.norm(
        u,
        axis=1,
    )

    umax = float(
        np.nanmax(
            umag
        )
    )


    epsp = get_field(
        final,
        "EquivalentPlasticStrain",
    ).squeeze()

    epsp = epsp[
        np.isfinite(
            epsp
        )
    ]

    epsp_max = float(
        np.nanmax(
            epsp
        )
    )


    pressure = get_field(
        final,
        "pressure",
    ).squeeze()

    pressure = pressure[
        np.isfinite(
            pressure
        )
    ]


    finite_ok = (
        np.isfinite(
            umax
        )
        and np.isfinite(
            epsp_max
        )
    )


    if (
        finite_ok
        and epsp_max <= 1e-8
    ):

        status = "PASS"

    else:

        status = "REVIEW"


    print(
        f"OGS return code: {rc}"
    )

    print(
        f"Cells: {bulk.n_cells}"
    )

    print(
        "Final time [s]: "
        f"{times[-1]:.6f}"
    )

    print(
        "Max displacement [mm]: "
        f"{umax*1000:.10f}"
    )

    print(
        "Max EquivalentPlasticStrain: "
        f"{epsp_max:.12e}"
    )

    print(
        "Pressure range [kPa]: "
        f"{np.min(pressure)/1000:.6f} "
        f"to "
        f"{np.max(pressure)/1000:.6f}"
    )

    print(
        f"STATUS: {status}"
    )


    records.append(
        {
            "case":
                case,

            "return_code":
                rc,

            "cells":
                bulk.n_cells,

            "max_displacement_mm":
                umax
                * 1000.0,

            "epsp_max":
                epsp_max,

            "pressure_min_kpa":
                float(
                    np.min(
                        pressure
                    )
                    / 1000.0
                ),

            "pressure_max_kpa":
                float(
                    np.max(
                        pressure
                    )
                    / 1000.0
                ),

            "status":
                status,
        }
    )


# ============================================================
# CROSS-MESH COMPARISON
# ============================================================

print()
print(
    "=== CROSS-MESH COMPARISON ==="
)


valid = [
    r
    for r in records
    if r[
        "status"
    ] == "PASS"
]


if valid:

    coarse = next(
        (
            r
            for r in valid
            if r["case"] == "coarse"
        ),
        None,
    )


    if coarse is not None:

        for r in valid:

            du = (
                r[
                    "max_displacement_mm"
                ]
                - coarse[
                    "max_displacement_mm"
                ]
            )

            print(
                f"{r['case']:6s} | "
                f"cells={r['cells']:5d} | "
                f"u={r['max_displacement_mm']:.10f} mm | "
                f"Delta u vs coarse="
                f"{du:+.10f} mm"
            )


all_pass = (
    len(records) == 3
    and all(
        r["status"] == "PASS"
        for r in records
    )
)


print()
print(
    "=== PHASE STATUS ==="
)

if all_pass:

    print(
        "REFINED RESTART EQUILIBRIUM: PASS"
    )

else:

    print(
        "REFINED RESTART EQUILIBRIUM: REVIEW"
    )


# ============================================================
# WRITE CSV + SUMMARY
# ============================================================

csv_path = (
    OUT
    / "restart_equilibrium_metrics.csv"
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
            "return_code",
            "cells",
            "max_displacement_mm",
            "epsp_max",
            "pressure_min_kpa",
            "pressure_max_kpa",
            "status",
        ],
    )

    writer.writeheader()
    writer.writerows(
        records
    )


summary = (
    OUT
    / "phase05d1_summary.txt"
)


lines = [
    "PHASE 05D-1 REFINED RESTART EQUILIBRIUM",
    "",
]


for r in records:

    lines.append(
        f"{r['case']} | "
        f"RC={r['return_code']} | "
        f"cells={r['cells']} | "
        f"u={r['max_displacement_mm']} mm | "
        f"epsp_max={r['epsp_max']} | "
        f"STATUS={r['status']}"
    )


lines += [
    "",
    (
        "REFINED RESTART EQUILIBRIUM: "
        + (
            "PASS"
            if all_pass
            else "REVIEW"
        )
    ),
    "",
    (
        "Interpretation: this phase tests "
        "restart-field transfer only. "
        "No erosion is applied."
    ),
]


summary.write_text(
    "\n".join(lines) + "\n",
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
