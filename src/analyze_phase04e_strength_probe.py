from pathlib import Path
import csv

import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot


ROOT = Path.cwd()

MODEL = (
    ROOT
    / "model"
    / "phase04e_strength_probe"
)

RUN = (
    ROOT
    / "results"
    / "phase04e_strength_probe"
)

OUT = (
    ROOT
    / "results"
    / "phase04e_analysis"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# READ CASE MANIFEST
# ============================================================

cases = []

with (
    MODEL / "cases.csv"
).open(
    "r",
    encoding="utf-8",
) as f:

    reader = csv.DictReader(
        f
    )

    for row in reader:

        cases.append(
            {
                "case":
                    row["case"],

                "srf":
                    float(
                        row["srf"]
                    ),

                "cohesion_kpa":
                    float(
                        row[
                            "cohesion_kpa"
                        ]
                    ),

                "friction_angle_deg":
                    float(
                        row[
                            "friction_angle_deg"
                        ]
                    ),
            }
        )


# ============================================================
# READ SOLVER STATUS
# ============================================================

solver_rc = {}

status_csv = (
    RUN
    / "run_status.csv"
)

if status_csv.exists():

    with status_csv.open(
        "r",
        encoding="utf-8",
    ) as f:

        reader = csv.DictReader(
            f
        )

        for row in reader:

            solver_rc[
                row["case"]
            ] = int(
                row["return_code"]
            )


# ============================================================
# FIELD HELPER
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


# ============================================================
# ANALYZE EACH TARGET SRF
# ============================================================

records = []

for item in cases:

    case = item["case"]

    folder = (
        RUN
        / case
    )

    rc = solver_rc.get(
        case,
        999,
    )

    record = dict(
        item
    )

    record[
        "solver_return_code"
    ] = rc

    pvds = sorted(
        folder.glob("*.pvd")
    )

    if rc != 0 or not pvds:

        record.update(
            {
                "final_time_s":
                    np.nan,

                "max_displacement_mm":
                    np.nan,

                "epsp_max":
                    np.nan,

                "frac_epsp_gt_1e8":
                    np.nan,

                "frac_epsp_gt_1e6":
                    np.nan,

                "frac_epsp_gt_1e4":
                    np.nan,

                "pressure_min_kpa":
                    np.nan,

                "pressure_max_kpa":
                    np.nan,

                "classification":
                    "NONCONVERGED",
            }
        )

        records.append(
            record
        )

        continue

    series = ot.MeshSeries(
        str(pvds[0])
    )

    times = np.asarray(
        series.timevalues,
        dtype=float,
    )

    mesh = series.mesh(
        len(times) - 1
    )

    u = get_field(
        mesh,
        "displacement",
    )

    epsp = get_field(
        mesh,
        "EquivalentPlasticStrain",
    ).squeeze()

    pressure = get_field(
        mesh,
        "pressure",
    ).squeeze()

    umag = np.linalg.norm(
        u,
        axis=1,
    )

    max_u = float(
        np.nanmax(
            umag
        )
    )

    epsp_max = float(
        np.nanmax(
            epsp
        )
    )

    f8 = float(
        np.mean(
            epsp > 1e-8
        )
    )

    f6 = float(
        np.mean(
            epsp > 1e-6
        )
    )

    f4 = float(
        np.mean(
            epsp > 1e-4
        )
    )

    if epsp_max <= 1e-10:

        classification = (
            "ELASTIC"
        )

    elif f4 <= 0.01:

        classification = (
            "YIELD_ONSET"
        )

    elif f4 <= 0.20:

        classification = (
            "LOCALIZED_PLASTICITY"
        )

    else:

        classification = (
            "WIDESPREAD_PLASTICITY"
        )

    record.update(
        {
            "final_time_s":
                float(
                    times[-1]
                ),

            "max_displacement_mm":
                max_u * 1000.0,

            "epsp_max":
                epsp_max,

            "frac_epsp_gt_1e8":
                f8,

            "frac_epsp_gt_1e6":
                f6,

            "frac_epsp_gt_1e4":
                f4,

            "pressure_min_kpa":
                float(
                    np.nanmin(
                        pressure
                    )
                    / 1000.0
                ),

            "pressure_max_kpa":
                float(
                    np.nanmax(
                        pressure
                    )
                    / 1000.0
                ),

            "classification":
                classification,
        }
    )

    records.append(
        record
    )


# ============================================================
# PRINT TABLE
# ============================================================

print(
    "========================================"
)

print(
    "PHASE 04E STRENGTH-REDUCTION PROBE"
)

print(
    "========================================"
)

print()

for r in records:

    if (
        r["classification"]
        == "NONCONVERGED"
    ):

        print(
            f"SRF={r['srf']:.2f} | "
            f"c={r['cohesion_kpa']:.2f} kPa | "
            f"phi={r['friction_angle_deg']:.2f} deg | "
            f"NONCONVERGED"
        )

        continue

    print(
        f"SRF={r['srf']:.2f} | "
        f"c={r['cohesion_kpa']:.2f} kPa | "
        f"phi={r['friction_angle_deg']:.2f} deg | "
        f"epsp_max={r['epsp_max']:.3e} | "
        f"frac>1e-4="
        f"{100*r['frac_epsp_gt_1e4']:.3f}% | "
        f"u={r['max_displacement_mm']:.5f} mm | "
        f"{r['classification']}"
    )


# ============================================================
# FIND FIRST TRANSITION BRACKET
# ============================================================

last_elastic = None
first_transition = None

for r in records:

    if (
        r["classification"]
        == "ELASTIC"
    ):

        if first_transition is None:
            last_elastic = r

        continue

    first_transition = r
    break


print()
print(
    "=== TRANSITION BRACKET ==="
)

if (
    last_elastic is not None
    and first_transition is not None
):

    print(
        "Last elastic target SRF: "
        f"{last_elastic['srf']:.3f}"
    )

    print(
        "First yielding/nonconverged target SRF: "
        f"{first_transition['srf']:.3f}"
    )

    print(
        "MECHANICAL TRANSITION BRACKET: PASS"
    )

elif (
    first_transition is not None
    and last_elastic is None
):

    print(
        "Transition already occurs at "
        "the first tested SRF."
    )

    print(
        "MECHANICAL TRANSITION BRACKET: REVIEW"
    )

else:

    print(
        "No transition found through "
        f"SRF={records[-1]['srf']:.2f}."
    )

    print(
        "MECHANICAL TRANSITION BRACKET: NOT YET BRACKETED"
    )


# ============================================================
# WRITE CSV
# ============================================================

csv_path = (
    OUT
    / "strength_probe_metrics.csv"
)

fieldnames = [
    "case",
    "srf",
    "cohesion_kpa",
    "friction_angle_deg",
    "solver_return_code",
    "final_time_s",
    "max_displacement_mm",
    "epsp_max",
    "frac_epsp_gt_1e8",
    "frac_epsp_gt_1e6",
    "frac_epsp_gt_1e4",
    "pressure_min_kpa",
    "pressure_max_kpa",
    "classification",
]

with csv_path.open(
    "w",
    newline="",
    encoding="utf-8",
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=fieldnames,
    )

    writer.writeheader()

    for r in records:

        writer.writerow(
            {
                key:
                    r.get(
                        key,
                        ""
                    )
                for key in fieldnames
            }
        )


# ============================================================
# FIGURE 01 — PLASTICITY RESPONSE
# ============================================================

valid = [
    r
    for r in records
    if (
        r["classification"]
        != "NONCONVERGED"
        and np.isfinite(
            r["epsp_max"]
        )
    )
]

if valid:

    srf = np.array(
        [
            r["srf"]
            for r in valid
        ]
    )

    epsp = np.array(
        [
            max(
                r["epsp_max"],
                1e-14,
            )
            for r in valid
        ]
    )

    fig, ax = plt.subplots(
        figsize=(7.2, 4.8)
    )

    ax.semilogy(
        srf,
        epsp,
        marker="o",
    )

    ax.axhline(
        1e-8,
        linestyle="--",
        linewidth=1,
    )

    ax.set_xlabel(
        "Strength reduction factor, SRF [-]"
    )

    ax.set_ylabel(
        "Maximum equivalent plastic strain [-]"
    )

    ax.set_title(
        "Strength-reduction probe"
    )

    ax.grid(
        alpha=0.25
    )

    fig.tight_layout()

    fig.savefig(
        OUT
        / "figure_01_strength_probe_plasticity.png",
        dpi=220,
    )

    plt.close(
        fig
    )


# ============================================================
# FIGURE 02 — DISPLACEMENT RESPONSE
# ============================================================

if valid:

    u = np.array(
        [
            r[
                "max_displacement_mm"
            ]
            for r in valid
        ]
    )

    fig, ax = plt.subplots(
        figsize=(7.2, 4.8)
    )

    ax.plot(
        srf,
        u,
        marker="o",
    )

    ax.set_xlabel(
        "Strength reduction factor, SRF [-]"
    )

    ax.set_ylabel(
        "Incremental maximum displacement [mm]"
    )

    ax.set_title(
        "Mechanical response to strength reduction"
    )

    ax.grid(
        alpha=0.25
    )

    fig.tight_layout()

    fig.savefig(
        OUT
        / "figure_02_strength_probe_displacement.png",
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
    / "phase04e_summary.txt"
)

lines = [
    "PHASE 04E STRENGTH-REDUCTION PROBE",
    "",
]

for r in records:

    if (
        r["classification"]
        == "NONCONVERGED"
    ):

        lines.append(
            f"SRF={r['srf']:.2f} | "
            f"NONCONVERGED"
        )

    else:

        lines.append(
            f"SRF={r['srf']:.2f} | "
            f"c={r['cohesion_kpa']:.3f} kPa | "
            f"phi={r['friction_angle_deg']:.3f} deg | "
            f"epsp_max={r['epsp_max']:.6e} | "
            f"frac_epsp_gt_1e4="
            f"{100*r['frac_epsp_gt_1e4']:.5f}% | "
            f"u={r['max_displacement_mm']:.6f} mm | "
            f"{r['classification']}"
        )

lines.append(
    ""
)

if (
    last_elastic is not None
    and first_transition is not None
):

    lines += [
        (
            "Last elastic SRF: "
            f"{last_elastic['srf']:.6f}"
        ),
        (
            "First transition SRF: "
            f"{first_transition['srf']:.6f}"
        ),
        (
            "MECHANICAL TRANSITION BRACKET: PASS"
        ),
    ]

elif first_transition is not None:

    lines += [
        (
            "Transition occurs at first tested SRF."
        ),
        (
            "MECHANICAL TRANSITION BRACKET: REVIEW"
        ),
    ]

else:

    lines += [
        (
            "No transition bracketed within tested SRF range."
        ),
        (
            "MECHANICAL TRANSITION BRACKET: NOT YET BRACKETED"
        ),
    ]

lines += [
    "",
    (
        "NOTE: SRF probe is a numerical "
        "strength-proximity diagnostic."
    ),
    (
        "NOTE: It is not reported here as "
        "a formal slope factor of safety."
    ),
]

summary.write_text(
    "\n".join(lines) + "\n",
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
    "PHASE 04E ANALYSIS COMPLETE"
)
