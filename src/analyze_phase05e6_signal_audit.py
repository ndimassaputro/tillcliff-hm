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
    / "phase05e6_signal_audit"
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

HOLD_END = 20.0
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
            f"No output near t={target}; "
            f"nearest={times[idx]}, "
            f"error={error}"
        )

    return idx


def body_mask(mesh):

    points = np.asarray(
        mesh.points,
        dtype=float,
    )

    top = surface_height(
        points[:, 0]
    )

    depth = (
        top
        - points[:, 1]
    )

    mask = (
        (points[:, 0] >= 8.0)
        & (points[:, 0] <= 20.5)
        & (depth >= -1e-8)
        & (depth <= 4.0)
    )

    if not np.any(mask):

        raise RuntimeError(
            "Empty body monitoring zone"
        )

    return mask


def displacement(mesh):

    return np.asarray(
        mesh.point_data[
            "displacement"
        ],
        dtype=float,
    )


def rms_mm(
    u_a,
    u_b,
    mask,
):

    du = (
        u_a
        - u_b
    )[
        mask
    ]

    mag = np.linalg.norm(
        du,
        axis=1,
    )

    return float(
        np.sqrt(
            np.mean(
                mag ** 2
            )
        )
        * 1000.0
    )


def vector_body(
    u_a,
    u_b,
    mask,
):

    return (
        (
            u_a
            - u_b
        )[
            mask
        ]
        .reshape(-1)
    )


def cosine(
    a,
    b,
):

    na = float(
        np.linalg.norm(a)
    )

    nb = float(
        np.linalg.norm(b)
    )

    if (
        na <= 0.0
        or nb <= 0.0
    ):

        return np.nan

    return float(
        np.dot(
            a,
            b,
        )
        / (
            na
            * nb
        )
    )


def classify_snr(snr):

    if snr >= 20.0:

        return "VERY CLEAN"

    if snr >= 10.0:

        return "CLEAN"

    if snr >= 3.0:

        return "CAUTION"

    return "NOISE-SENSITIVE"


case_data = {}


print(
    "========================================"
)

print(
    "PHASE 05E-6 SIGNAL-TO-DRIFT AUDIT"
)

print(
    "========================================"
)


for case in CASES:

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


    indices = {
        t:
            nearest_index(
                times,
                t,
            )
        for t in [
            18.0,
            19.0,
            20.0,
        ]
    }


    meshes = {
        t:
            series.mesh(
                indices[t]
            )
        for t in indices
    }


    mask = body_mask(
        meshes[
            20.0
        ]
    )


    u18 = displacement(
        meshes[
            18.0
        ]
    )

    u19 = displacement(
        meshes[
            19.0
        ]
    )

    u20 = displacement(
        meshes[
            20.0
        ]
    )


    drift_18_19 = rms_mm(
        u19,
        u18,
        mask,
    )

    drift_19_20 = rms_mm(
        u20,
        u19,
        mask,
    )


    noise_floor = max(
        drift_18_19,
        drift_19_20,
    )


    case_data[
        case
    ] = {
        "series":
            series,

        "times":
            times,

        "mask":
            mask,

        "u20":
            u20,

        "drift_18_19_mm":
            drift_18_19,

        "drift_19_20_mm":
            drift_19_20,

        "noise_floor_mm":
            noise_floor,
    }


    print()
    print(
        f"=== {case} INTACT DRIFT ==="
    )

    print(
        "RMS drift 18->19 s [mm]: "
        f"{drift_18_19:.12e}"
    )

    print(
        "RMS drift 19->20 s [mm]: "
        f"{drift_19_20:.12e}"
    )

    print(
        "Noise-floor proxy [mm]: "
        f"{noise_floor:.12e}"
    )


records = []


for E in TARGET_E:

    print()
    print(
        "========================================"
    )

    print(
        f"E = {E:.2f} m"
    )

    print(
        "========================================"
    )


    vectors = {}


    for case in CASES:

        data = case_data[
            case
        ]

        target_time = (
            HOLD_END
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


        u = displacement(
            mesh
        )


        signal = rms_mm(
            u,
            data[
                "u20"
            ],
            data[
                "mask"
            ],
        )


        floor = data[
            "noise_floor_mm"
        ]


        if floor > 0.0:

            snr = (
                signal
                / floor
            )

        else:

            snr = np.inf


        category = classify_snr(
            snr
        )


        vec = vector_body(
            u,
            data[
                "u20"
            ],
            data[
                "mask"
            ],
        )


        vectors[
            case
        ] = vec


        records.append(
            {
                "E_m":
                    E,

                "case":
                    case,

                "signal_Rrms_mm":
                    signal,

                "drift_18_19_mm":
                    data[
                        "drift_18_19_mm"
                    ],

                "drift_19_20_mm":
                    data[
                        "drift_19_20_mm"
                    ],

                "noise_floor_mm":
                    floor,

                "signal_to_drift":
                    snr,

                "classification":
                    category,
            }
        )


        print(
            f"{case:12s} | "
            f"signal="
            f"{signal:.10f} mm | "
            f"floor="
            f"{floor:.3e} mm | "
            f"SNR="
            f"{snr:.2f} | "
            f"{category}"
        )


    dry_cos = cosine(
        vectors[
            "dryP_refS"
        ],
        vectors[
            "refP_refS"
        ],
    )

    wet_cos = cosine(
        vectors[
            "wetP_refS"
        ],
        vectors[
            "refP_refS"
        ],
    )


    print()

    print(
        "Mode cosine:"
    )

    print(
        "  dry/ref = "
        f"{dry_cos:.8f}"
    )

    print(
        "  wet/ref = "
        f"{wet_cos:.8f}"
    )


# ============================================================
# DECISION
# ============================================================

wet_records = [
    r
    for r in records
    if r[
        "case"
    ] == "wetP_refS"
]


wet_min_snr = float(
    min(
        r[
            "signal_to_drift"
        ]
        for r in wet_records
    )
)


print()
print(
    "========================================"
)

print(
    "FINAL SIGNAL-SIGNIFICANCE DECISION"
)

print(
    "========================================"
)

print(
    "Minimum WET signal/drift ratio: "
    f"{wet_min_snr:.4f}"
)


if wet_min_snr >= 10.0:

    decision = (
        "WET MODE-SHAPE DIFFERENCE "
        "IS NUMERICALLY RESOLVED"
    )

elif wet_min_snr >= 3.0:

    decision = (
        "WET MODE-SHAPE DIFFERENCE "
        "IS RESOLVED BUT SHOULD BE "
        "INTERPRETED CAUTIOUSLY"
    )

else:

    decision = (
        "WET MODE-SHAPE DIFFERENCE "
        "IS TOO CLOSE TO NUMERICAL "
        "DRIFT FOR A STRONG CLAIM"
    )


print(
    decision
)


# ============================================================
# CSV
# ============================================================

csv_path = (
    OUT
    / "signal_to_drift_audit.csv"
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
    / "phase05e6_summary.txt"
)


lines = [
    "PHASE 05E-6 SIGNAL-TO-DRIFT AUDIT",
    "",
]


for case in CASES:

    d = case_data[
        case
    ]

    lines.append(
        f"{case} | "
        f"drift18-19="
        f"{d['drift_18_19_mm']:.12e} mm | "
        f"drift19-20="
        f"{d['drift_19_20_mm']:.12e} mm | "
        f"floor="
        f"{d['noise_floor_mm']:.12e} mm"
    )


lines += [
    "",
]


for r in records:

    lines.append(
        f"E={r['E_m']:.2f} | "
        f"{r['case']} | "
        f"signal="
        f"{r['signal_Rrms_mm']:.10f} mm | "
        f"SNR="
        f"{r['signal_to_drift']:.4f} | "
        f"{r['classification']}"
    )


lines += [
    "",
    (
        "Minimum wet SNR: "
        f"{wet_min_snr:.4f}"
    ),
    (
        "DECISION: "
        f"{decision}"
    ),
    "",
    (
        "The intact-hold displacement drift "
        "is used only as a numerical "
        "signal-floor diagnostic."
    ),
    (
        "Solver nonconvergence remains "
        "excluded as a physical failure "
        "criterion."
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
    "PHASE 05E-6 ANALYSIS COMPLETE"
)
