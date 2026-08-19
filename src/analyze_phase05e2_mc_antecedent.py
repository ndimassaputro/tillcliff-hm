from pathlib import Path
import csv

import numpy as np
import ogstools as ot
import pyvista as pv


ROOT = Path.cwd()

MODEL_ROOT = (
    ROOT
    / "model"
    / "phase05e2_mc_antecedent"
)

RUN_ROOT = (
    ROOT
    / "results"
    / "phase05e2_mc_antecedent"
)

OUT = (
    ROOT
    / "results"
    / "phase05e2_analysis"
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
# RUN STATUS
# ============================================================

status_map = {}

status_file = (
    RUN_ROOT
    / "run_status.csv"
)

if status_file.exists():

    with status_file.open(
        "r",
        encoding="utf-8",
    ) as f:

        for row in csv.DictReader(
            f
        ):

            status_map[
                row["state"]
            ] = int(
                row["return_code"]
            )


# ============================================================
# HELPERS
# ============================================================

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
    "PHASE 05E-2 MC ANTECEDENT EQUILIBRIUM"
)

print(
    "========================================"
)


for state in STATES:

    print()
    print(
        f"=== {state.upper()} ==="
    )


    rc = status_map.get(
        state,
        999,
    )

    model_bulk = pv.read(
        MODEL_ROOT
        / state
        / "bulk.vtu"
    )

    p0 = np.asarray(
        model_bulk.point_data[
            "p0"
        ],
        dtype=float,
    ).squeeze()

    sat0 = np.asarray(
        model_bulk.point_data[
            "antecedent_saturation"
        ],
        dtype=float,
    ).squeeze()


    run_dir = (
        RUN_ROOT
        / state
    )

    pvds = sorted(
        run_dir.glob(
            "*.pvd"
        )
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
                "state":
                    state,

                "return_code":
                    rc,

                "p0_mean_kpa":
                    np.mean(p0)
                    / 1000.0,

                "sat0_mean":
                    np.mean(sat0),

                "final_p_mean_kpa":
                    np.nan,

                "final_sat_mean":
                    np.nan,

                "p_change_rms_pa":
                    np.nan,

                "p_change_max_pa":
                    np.nan,

                "max_displacement_mm":
                    np.nan,

                "epsp_max":
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


    pressure = get_field(
        final,
        "pressure",
    ).squeeze()

    saturation = get_field(
        final,
        "saturation",
    ).squeeze()

    displacement = get_field(
        final,
        "displacement",
    )

    epsp = get_field(
        final,
        "EquivalentPlasticStrain",
    ).squeeze()


    du = np.linalg.norm(
        displacement,
        axis=1,
    )

    dp = (
        pressure
        - p0
    )


    p_change_rms = float(
        np.sqrt(
            np.mean(
                dp ** 2
            )
        )
    )

    p_change_max = float(
        np.max(
            np.abs(
                dp
            )
        )
    )

    umax = float(
        np.max(
            du
        )
        * 1000.0
    )

    epsp_finite = epsp[
        np.isfinite(
            epsp
        )
    ]

    epsp_max = float(
        np.max(
            epsp_finite
        )
    )


    finite_ok = (
        np.all(
            np.isfinite(
                pressure
            )
        )
        and np.all(
            np.isfinite(
                displacement
            )
        )
        and np.isfinite(
            epsp_max
        )
    )


    no_plasticity = (
        epsp_max <= 1e-8
    )


    status = (
        "PASS"
        if (
            rc == 0
            and finite_ok
            and no_plasticity
        )
        else "REVIEW"
    )


    print(
        f"OGS return code: {rc}"
    )

    print(
        f"Final time [s]: "
        f"{times[-1]:.6f}"
    )

    print(
        "Initial mean pressure [kPa]: "
        f"{np.mean(p0)/1000:.6f}"
    )

    print(
        "Final mean pressure [kPa]: "
        f"{np.mean(pressure)/1000:.6f}"
    )

    print(
        "Initial mean saturation: "
        f"{np.mean(sat0):.8f}"
    )

    print(
        "Final mean saturation: "
        f"{np.mean(saturation):.8f}"
    )

    print(
        "Pressure change RMS [Pa]: "
        f"{p_change_rms:.6f}"
    )

    print(
        "Pressure change max [Pa]: "
        f"{p_change_max:.6f}"
    )

    print(
        "Max displacement during hold [mm]: "
        f"{umax:.10f}"
    )

    print(
        "Max EquivalentPlasticStrain: "
        f"{epsp_max:.12e}"
    )

    print(
        f"STATUS: {status}"
    )


    records.append(
        {
            "state":
                state,

            "return_code":
                rc,

            "p0_mean_kpa":
                float(
                    np.mean(p0)
                    / 1000.0
                ),

            "sat0_mean":
                float(
                    np.mean(sat0)
                ),

            "final_p_mean_kpa":
                float(
                    np.mean(
                        pressure
                    )
                    / 1000.0
                ),

            "final_sat_mean":
                float(
                    np.mean(
                        saturation
                    )
                ),

            "p_change_rms_pa":
                p_change_rms,

            "p_change_max_pa":
                p_change_max,

            "max_displacement_mm":
                umax,

            "epsp_max":
                epsp_max,

            "status":
                status,
        }
    )


# ============================================================
# CROSS-STATE ORDERING
# ============================================================

by_state = {
    r["state"]:
        r
    for r in records
}


print()
print(
    "========================================"
)

print(
    "FINAL HYDRAULIC ORDERING"
)

print(
    "========================================"
)


ordering_ready = all(
    np.isfinite(
        by_state[s][
            "final_p_mean_kpa"
        ]
    )
    for s in STATES
)


if ordering_ready:

    dry_p = by_state[
        "dry"
    ][
        "final_p_mean_kpa"
    ]

    ref_p = by_state[
        "reference"
    ][
        "final_p_mean_kpa"
    ]

    wet_p = by_state[
        "wet"
    ][
        "final_p_mean_kpa"
    ]


    dry_s = by_state[
        "dry"
    ][
        "final_sat_mean"
    ]

    ref_s = by_state[
        "reference"
    ][
        "final_sat_mean"
    ]

    wet_s = by_state[
        "wet"
    ][
        "final_sat_mean"
    ]


    pressure_order = (
        dry_p
        < ref_p
        < wet_p
    )

    saturation_order = (
        dry_s
        < ref_s
        < wet_s
    )


    print(
        f"Pressure mean [kPa]: "
        f"dry={dry_p:.6f}, "
        f"reference={ref_p:.6f}, "
        f"wet={wet_p:.6f}"
    )

    print(
        f"Saturation mean: "
        f"dry={dry_s:.8f}, "
        f"reference={ref_s:.8f}, "
        f"wet={wet_s:.8f}"
    )

    print(
        "Pressure ordering: "
        f"{'PASS' if pressure_order else 'REVIEW'}"
    )

    print(
        "Saturation ordering: "
        f"{'PASS' if saturation_order else 'REVIEW'}"
    )

else:

    pressure_order = False
    saturation_order = False

    print(
        "Ordering unavailable"
    )


all_pass = (
    len(records) == 3
    and all(
        r["status"] == "PASS"
        for r in records
    )
    and pressure_order
    and saturation_order
)


print()
print(
    "========================================"
)

print(
    "PHASE STATUS"
)

print(
    "========================================"
)


if all_pass:

    phase_status = (
        "PASS — THREE MC ANTECEDENT "
        "BRANCHES READY"
    )

else:

    phase_status = (
        "REVIEW — DO NOT RUN EROSION YET"
    )


print(
    phase_status
)


# ============================================================
# OUTPUTS
# ============================================================

csv_path = (
    OUT
    / "mc_antecedent_equilibrium.csv"
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


summary = (
    OUT
    / "phase05e2_summary.txt"
)


lines = [
    "PHASE 05E-2 MC ANTECEDENT EQUILIBRIUM",
    "",
]


for r in records:

    lines.append(
        f"{r['state']} | "
        f"RC={r['return_code']} | "
        f"p0={r['p0_mean_kpa']} kPa | "
        f"pfinal={r['final_p_mean_kpa']} kPa | "
        f"Srfinal={r['final_sat_mean']} | "
        f"dp_rms={r['p_change_rms_pa']} Pa | "
        f"u={r['max_displacement_mm']} mm | "
        f"epsp={r['epsp_max']} | "
        f"{r['status']}"
    )


lines += [
    "",
    (
        "Pressure ordering: "
        f"{pressure_order}"
    ),
    (
        "Saturation ordering: "
        f"{saturation_order}"
    ),
    "",
    (
        "STATUS: "
        f"{phase_status}"
    ),
    "",
    (
        "NOTE: sigma0 is a bilinear projection "
        "of the Phase04B nodal sigma output "
        "to medium-mesh element centres."
    ),
    (
        "Displacement is reset; the hold "
        "therefore tests the consistency of "
        "the transferred hydraulic/stress state "
        "with the Mohr-Coulomb model."
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
    "PHASE 05E-2 ANALYSIS COMPLETE"
)
