from pathlib import Path
import csv

import numpy as np
import ogstools as ot


ROOT = Path.cwd()

RUN_ROOT = (
    ROOT
    / "results"
    / "phase05e4_decomposition"
)

OUT = (
    ROOT
    / "results"
    / "phase05e5_mode_audit"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


CASES = [
    "refP_refS",
    "dryP_refS",
    "wetP_refS",
]

TARGET_E = [
    0.05,
    0.15,
    0.25,
]

HOLD_TIME = 20.0
EROSION_RATE = 0.01


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


def nearest_index(
    times,
    target,
):

    idx = int(
        np.argmin(
            np.abs(
                times - target
            )
        )
    )

    error = abs(
        float(
            times[idx]
        )
        - target
    )

    if error > 0.051:

        raise RuntimeError(
            f"No output close to "
            f"t={target}; error={error}"
        )

    return idx


def load_case(case):

    folder = (
        RUN_ROOT
        / case
    )

    pvds = sorted(
        folder.glob(
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

    i0 = nearest_index(
        times,
        HOLD_TIME,
    )

    base = series.mesh(
        i0
    )

    points = np.asarray(
        base.points,
        dtype=float,
    )

    top = surface_height(
        points[:, 0]
    )

    depth = (
        top
        - points[:, 1]
    )

    body = (
        (points[:, 0] >= 8.0)
        & (points[:, 0] <= 20.5)
        & (depth >= -1e-8)
        & (depth <= 4.0)
    )

    if not np.any(
        body
    ):

        raise RuntimeError(
            f"Empty body mask: {case}"
        )

    u0 = np.asarray(
        base.point_data[
            "displacement"
        ],
        dtype=float,
    )

    return {
        "case":
            case,

        "series":
            series,

        "times":
            times,

        "body":
            body,

        "u0":
            u0,
    }


def state_at_E(
    data,
    E,
):

    target_time = (
        HOLD_TIME
        + E
        / EROSION_RATE
    )

    idx = nearest_index(
        data[
            "times"
        ],
        target_time,
    )

    mesh = data[
        "series"
    ].mesh(
        idx
    )

    u = np.asarray(
        mesh.point_data[
            "displacement"
        ],
        dtype=float,
    )

    du = (
        u
        - data[
            "u0"
        ]
    )

    body = data[
        "body"
    ]

    du_body = du[
        body
    ]

    flat = du_body.reshape(
        -1
    )

    mag = np.linalg.norm(
        du_body,
        axis=1,
    )

    Rrms = float(
        np.sqrt(
            np.mean(
                mag ** 2
            )
        )
        * 1000.0
    )

    R95 = float(
        np.percentile(
            mag,
            95.0,
        )
        * 1000.0
    )

    pressure = np.asarray(
        mesh.point_data[
            "pressure"
        ],
        dtype=float,
    ).squeeze()

    saturation = np.asarray(
        mesh.point_data[
            "saturation"
        ],
        dtype=float,
    ).squeeze()

    return {
        "time":
            float(
                data[
                    "times"
                ][idx]
            ),

        "vector":
            flat,

        "Rrms_mm":
            Rrms,

        "R95_mm":
            R95,

        "pbody_kpa":
            float(
                np.mean(
                    pressure[
                        body
                    ]
                )
                / 1000.0
            ),

        "Srbody":
            float(
                np.mean(
                    saturation[
                        body
                    ]
                )
            ),
    }


def comparison(
    case_vector,
    ref_vector,
):

    nr = float(
        np.linalg.norm(
            ref_vector
        )
    )

    nc = float(
        np.linalg.norm(
            case_vector
        )
    )

    if (
        nr <= 0.0
        or nc <= 0.0
    ):

        return {
            "cosine":
                np.nan,

            "amplitude":
                np.nan,

            "shape_residual":
                np.nan,
        }

    dot = float(
        np.dot(
            case_vector,
            ref_vector,
        )
    )

    cosine = (
        dot
        / (
            nr
            * nc
        )
    )

    amplitude = (
        dot
        / (
            nr ** 2
        )
    )

    residual = float(
        np.linalg.norm(
            case_vector
            - amplitude
            * ref_vector
        )
        / nc
    )

    return {
        "cosine":
            cosine,

        "amplitude":
            amplitude,

        "shape_residual":
            residual,
    }


print(
    "========================================"
)

print(
    "PHASE 05E-5 DEFORMATION-MODE AUDIT"
)

print(
    "========================================"
)


data = {
    case:
        load_case(
            case
        )
    for case in CASES
}


records = []


for E in TARGET_E:

    states = {
        case:
            state_at_E(
                data[
                    case
                ],
                E,
            )
        for case in CASES
    }

    ref = states[
        "refP_refS"
    ]

    dry = states[
        "dryP_refS"
    ]

    wet = states[
        "wetP_refS"
    ]


    dry_cmp = comparison(
        dry[
            "vector"
        ],
        ref[
            "vector"
        ],
    )

    wet_cmp = comparison(
        wet[
            "vector"
        ],
        ref[
            "vector"
        ],
    )


    record = {
        "E_m":
            E,

        "ref_Rrms_mm":
            ref[
                "Rrms_mm"
            ],

        "dry_Rrms_mm":
            dry[
                "Rrms_mm"
            ],

        "wet_Rrms_mm":
            wet[
                "Rrms_mm"
            ],

        "dry_Rrms_ratio":
            dry[
                "Rrms_mm"
            ]
            / ref[
                "Rrms_mm"
            ],

        "wet_Rrms_ratio":
            wet[
                "Rrms_mm"
            ]
            / ref[
                "Rrms_mm"
            ],

        "dry_cosine":
            dry_cmp[
                "cosine"
            ],

        "wet_cosine":
            wet_cmp[
                "cosine"
            ],

        "dry_bestfit_amplitude":
            dry_cmp[
                "amplitude"
            ],

        "wet_bestfit_amplitude":
            wet_cmp[
                "amplitude"
            ],

        "dry_shape_residual":
            dry_cmp[
                "shape_residual"
            ],

        "wet_shape_residual":
            wet_cmp[
                "shape_residual"
            ],

        "dry_pbody_kpa":
            dry[
                "pbody_kpa"
            ],

        "ref_pbody_kpa":
            ref[
                "pbody_kpa"
            ],

        "wet_pbody_kpa":
            wet[
                "pbody_kpa"
            ],

        "dry_Srbody":
            dry[
                "Srbody"
            ],

        "ref_Srbody":
            ref[
                "Srbody"
            ],

        "wet_Srbody":
            wet[
                "Srbody"
            ],
    }

    records.append(
        record
    )


    print()
    print(
        f"=== E = {E:.2f} m ==="
    )

    print(
        "Rrms dry/ref/wet [mm]: "
        f"{dry['Rrms_mm']:.8f} / "
        f"{ref['Rrms_mm']:.8f} / "
        f"{wet['Rrms_mm']:.8f}"
    )

    print(
        "DRY vs REF:"
    )

    print(
        "  cosine            = "
        f"{dry_cmp['cosine']:.8f}"
    )

    print(
        "  best-fit amplitude = "
        f"{dry_cmp['amplitude']:.8f}"
    )

    print(
        "  shape residual     = "
        f"{dry_cmp['shape_residual']:.8f}"
    )

    print(
        "WET vs REF:"
    )

    print(
        "  cosine            = "
        f"{wet_cmp['cosine']:.8f}"
    )

    print(
        "  best-fit amplitude = "
        f"{wet_cmp['amplitude']:.8f}"
    )

    print(
        "  shape residual     = "
        f"{wet_cmp['shape_residual']:.8f}"
    )


# ============================================================
# DECISION
# ============================================================

dry_cos = np.asarray(
    [
        r[
            "dry_cosine"
        ]
        for r in records
    ],
    dtype=float,
)

wet_cos = np.asarray(
    [
        r[
            "wet_cosine"
        ]
        for r in records
    ],
    dtype=float,
)

dry_res = np.asarray(
    [
        r[
            "dry_shape_residual"
        ]
        for r in records
    ],
    dtype=float,
)

wet_res = np.asarray(
    [
        r[
            "wet_shape_residual"
        ]
        for r in records
    ],
    dtype=float,
)


common_mode = (
    np.all(
        dry_cos >= 0.95
    )
    and np.all(
        wet_cos >= 0.95
    )
    and np.all(
        dry_res <= 0.30
    )
    and np.all(
        wet_res <= 0.30
    )
)


print()
print(
    "========================================"
)

print(
    "MODE-SHAPE DECISION"
)

print(
    "========================================"
)


if common_mode:

    decision = (
        "COMMON DEFORMATION MODE — "
        "HYDRAULIC STATE PRIMARILY "
        "SCALES RESPONSE AMPLITUDE"
    )

else:

    decision = (
        "HYDRAULIC STATE ALSO CHANGES "
        "DEFORMATION MODE SHAPE"
    )


print(
    decision
)


print()
print(
    "Minimum dry/ref cosine: "
    f"{np.min(dry_cos):.8f}"
)

print(
    "Minimum wet/ref cosine: "
    f"{np.min(wet_cos):.8f}"
)

print(
    "Maximum dry shape residual: "
    f"{np.max(dry_res):.8f}"
)

print(
    "Maximum wet shape residual: "
    f"{np.max(wet_res):.8f}"
)


# ============================================================
# CSV
# ============================================================

csv_path = (
    OUT
    / "mode_shape_audit.csv"
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
    / "phase05e5_summary.txt"
)


lines = [
    "PHASE 05E-5 DEFORMATION-MODE AUDIT",
    "",
    (
        "Pressure-only cases use a common "
        "reference inherited sigma0."
    ),
    (
        "Displacement increments are measured "
        "relative to each case's own intact "
        "t=20 s equilibrium state."
    ),
    "",
]


for r in records:

    lines.append(
        f"E={r['E_m']:.2f} | "
        f"dry/ref R={r['dry_Rrms_ratio']:.6f} | "
        f"wet/ref R={r['wet_Rrms_ratio']:.6f} | "
        f"dry cos={r['dry_cosine']:.6f} | "
        f"wet cos={r['wet_cosine']:.6f} | "
        f"dry residual={r['dry_shape_residual']:.6f} | "
        f"wet residual={r['wet_shape_residual']:.6f}"
    )


lines += [
    "",
    (
        "DECISION: "
        f"{decision}"
    ),
    "",
    (
        "Cosine similarity assesses vector-field "
        "directional similarity."
    ),
    (
        "Best-fit amplitude measures scaling "
        "relative to the reference displacement "
        "mode."
    ),
    (
        "This remains a model-specific numerical "
        "diagnostic, not a stability criterion."
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
    "PHASE 05E-5 ANALYSIS COMPLETE"
)
