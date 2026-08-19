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
    / "phase05e4_decomposition"
)

RUN = (
    ROOT
    / "results"
    / "phase05e4_decomposition"
)

OUT = (
    ROOT
    / "results"
    / "phase05e4_analysis"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


CASES = [
    "refP_refS",
    "dryP_refS",
    "wetP_refS",
    "refP_dryS",
    "refP_wetS",
]


HOLD_END = 20.0
EROSION_RATE = 0.01
E_MAX = 0.25
TOE_X = 22.0


def erosion_from_time(t):

    if t <= HOLD_END:
        return 0.0

    return float(
        np.clip(
            EROSION_RATE
            * (
                t - HOLD_END
            ),
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


def point_field(mesh, name):

    return np.asarray(
        mesh.point_data[
            name
        ],
        dtype=float,
    )


def epsp_max(mesh):

    name = (
        "EquivalentPlasticStrain"
    )

    if name in mesh.cell_data:

        values = np.asarray(
            mesh.cell_data[name],
            dtype=float,
        )

    elif name in mesh.point_data:

        values = np.asarray(
            mesh.point_data[name],
            dtype=float,
        )

    else:

        raise RuntimeError(
            "EquivalentPlasticStrain missing"
        )

    values = values[
        np.isfinite(
            values
        )
    ]

    return float(
        np.max(
            values
        )
    )


def analyze_case(case):

    model_dir = (
        MODEL
        / case
    )

    run_dir = (
        RUN
        / case
    )

    bulk = pv.read(
        model_dir
        / "bulk.vtu"
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

    material_ids = np.asarray(
        bulk.cell_data[
            "MaterialIDs"
        ],
        dtype=int,
    )

    sizes = bulk.compute_cell_sizes(
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

    candidate = (
        material_ids > 0
    )

    recession = (
        TOE_X
        - centers[:, 0]
    )


    top = surface_height(
        points[:, 0]
    )

    depth = (
        top
        - points[:, 1]
    )

    body_nodes = (
        (points[:, 0] >= 8.0)
        & (points[:, 0] <= 20.5)
        & (depth >= -1e-8)
        & (depth <= 4.0)
    )


    fail_t = failure_time(
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


    valid = (
        np.arange(
            len(times)
        )
        if fail_t is None
        else np.where(
            times < fail_t - 1e-10
        )[0]
    )


    baseline_ids = [
        i
        for i in valid
        if times[i]
        <= HOLD_END + 1e-9
    ]

    if not baseline_ids:

        raise RuntimeError(
            f"No baseline: {case}"
        )


    base_idx = int(
        baseline_ids[-1]
    )

    base = series.mesh(
        base_idx
    )

    u0 = point_field(
        base,
        "displacement",
    )


    hold_epsp = epsp_max(
        base
    )


    p0 = np.asarray(
        bulk.point_data[
            "p0"
        ],
        dtype=float,
    ).squeeze()

    p_hold = point_field(
        base,
        "pressure",
    ).squeeze()

    dp_hold_rms = float(
        np.sqrt(
            np.mean(
                (
                    p_hold - p0
                )
                ** 2
            )
        )
    )


    records = []


    for idx in valid:

        t = float(
            times[idx]
        )

        if t < HOLD_END - 1e-9:
            continue


        E = erosion_from_time(
            t
        )

        mesh = series.mesh(
            int(idx)
        )

        u = point_field(
            mesh,
            "displacement",
        )

        du = np.linalg.norm(
            u - u0,
            axis=1,
        )

        values = du[
            body_nodes
        ]

        values = values[
            np.isfinite(
                values
            )
        ]


        Rrms = float(
            np.sqrt(
                np.mean(
                    values ** 2
                )
            )
            * 1000.0
        )

        R95 = float(
            np.percentile(
                values,
                95
            )
            * 1000.0
        )


        removed = (
            candidate
            & (
                recession
                <= E + 1e-10
            )
        )

        removed_area = float(
            np.sum(
                area[
                    removed
                ]
            )
        )


        records.append(
            {
                "E_m":
                    E,

                "removed_area_m2":
                    removed_area,

                "Rrms_mm":
                    Rrms,

                "R95_mm":
                    R95,
            }
        )


    return {
        "case":
            case,

        "failure_E":
            (
                erosion_from_time(
                    fail_t
                )
                if fail_t is not None
                else np.nan
            ),

        "hold_epsp":
            hold_epsp,

        "hold_dp_rms_pa":
            dp_hold_rms,

        "records":
            records,
    }


results = [
    analyze_case(case)
    for case in CASES
]


TARGETS = [
    0.05,
    0.15,
    0.25,
]


def nearest(result, target):

    rec = result[
        "records"
    ]

    E = np.asarray(
        [
            r["E_m"]
            for r in rec
        ]
    )

    idx = int(
        np.argmin(
            np.abs(
                E - target
            )
        )
    )

    if abs(
        E[idx] - target
    ) > 0.0051:
        return None

    return rec[idx]


table = []


for target in TARGETS:

    row = {
        "E_m":
            target,
    }

    valid = True

    for result in results:

        r = nearest(
            result,
            target
        )

        if r is None:

            valid = False
            break

        case = result[
            "case"
        ]

        row[
            f"{case}_Rrms_mm"
        ] = r[
            "Rrms_mm"
        ]

        row[
            f"{case}_R95_mm"
        ] = r[
            "R95_mm"
        ]

        row[
            f"{case}_Arem_m2"
        ] = r[
            "removed_area_m2"
        ]


    if valid:
        table.append(
            row
        )


print(
    "========================================"
)

print(
    "PHASE 05E-4 DECOMPOSITION"
)

print(
    "========================================"
)


print()
print(
    "=== INTACT HOLD CHECK ==="
)


holds_ok = True


for result in results:

    ok = (
        result[
            "hold_epsp"
        ] <= 1e-8
    )

    holds_ok = (
        holds_ok
        and ok
    )

    print(
        f"{result['case']:12s} | "
        f"epsp={result['hold_epsp']:.3e} | "
        f"dp_rms="
        f"{result['hold_dp_rms_pa']:.3f} Pa | "
        f"{'PASS' if ok else 'REVIEW'}"
    )


print()
print(
    "=== RESPONSE DECOMPOSITION ==="
)


for r in table:

    print()
    print(
        f"E={r['E_m']:.2f} m"
    )

    print(
        "  baseline REF_P + REF_S : "
        f"{r['refP_refS_Rrms_mm']:.8f} mm"
    )

    print(
        "  DRY_P + REF_S          : "
        f"{r['dryP_refS_Rrms_mm']:.8f} mm"
    )

    print(
        "  WET_P + REF_S          : "
        f"{r['wetP_refS_Rrms_mm']:.8f} mm"
    )

    print(
        "  REF_P + DRY_S          : "
        f"{r['refP_dryS_Rrms_mm']:.8f} mm"
    )

    print(
        "  REF_P + WET_S          : "
        f"{r['refP_wetS_Rrms_mm']:.8f} mm"
    )


    ref = r[
        "refP_refS_Rrms_mm"
    ]


    print(
        "  pressure-only ratios:"
    )

    print(
        "    dryP/ref = "
        f"{r['dryP_refS_Rrms_mm']/ref:.6f}"
    )

    print(
        "    wetP/ref = "
        f"{r['wetP_refS_Rrms_mm']/ref:.6f}"
    )


    print(
        "  stress-only ratios:"
    )

    print(
        "    dryS/ref = "
        f"{r['refP_dryS_Rrms_mm']/ref:.6f}"
    )

    print(
        "    wetS/ref = "
        f"{r['refP_wetS_Rrms_mm']/ref:.6f}"
    )


# ============================================================
# HEADLINE DIAGNOSTIC AT E=0.25
# ============================================================

r25 = next(
    (
        r
        for r in table
        if abs(
            r["E_m"] - 0.25
        ) < 1e-8
    ),
    None,
)


if r25 is None:

    decision = (
        "REVIEW — E=0.25 NOT COMMON"
    )

else:

    ref = r25[
        "refP_refS_Rrms_mm"
    ]

    pressure_wet_ratio = (
        r25[
            "wetP_refS_Rrms_mm"
        ]
        / ref
    )

    stress_wet_ratio = (
        r25[
            "refP_wetS_Rrms_mm"
        ]
        / ref
    )


    pressure_effect = abs(
        pressure_wet_ratio
        - 1.0
    )

    stress_effect = abs(
        stress_wet_ratio
        - 1.0
    )


    if (
        pressure_effect
        > 1.5
        * stress_effect
    ):

        decision = (
            "HYDRAULIC-STATE EFFECT DOMINANT"
        )

    elif (
        stress_effect
        > 1.5
        * pressure_effect
    ):

        decision = (
            "INHERITED-STRESS EFFECT DOMINANT"
        )

    else:

        decision = (
            "BOTH PRESSURE AND STRESS "
            "CONTRIBUTE MATERIALLY"
        )


print()
print(
    "========================================"
)

print(
    "DECOMPOSITION DECISION"
)

print(
    "========================================"
)

print(
    decision
)

print(
    "This is a numerical attribution "
    "diagnostic, not a universal "
    "physical causal statement."
)


csv_path = (
    OUT
    / "decomposition_response.csv"
)


if table:

    with csv_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=list(
                table[0].keys()
            ),
        )

        writer.writeheader()
        writer.writerows(
            table
        )


summary = (
    OUT
    / "phase05e4_summary.txt"
)


lines = [
    "PHASE 05E-4 PRESSURE/STRESS DECOMPOSITION",
    "",
    (
        "All intact mixed-state holds pass: "
        f"{holds_ok}"
    ),
    "",
    (
        "DECISION: "
        f"{decision}"
    ),
    "",
    (
        "Purpose: separate direct hydraulic "
        "antecedent-state effects from "
        "inherited effective-stress-history "
        "effects."
    ),
    (
        "Production interpretation remains "
        "restricted to the present model "
        "and screening parameter set."
    ),
]


summary.write_text(
    "\n".join(lines)
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
    "PHASE 05E-4 ANALYSIS COMPLETE"
)
