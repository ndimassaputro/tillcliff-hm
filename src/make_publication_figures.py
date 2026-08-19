from pathlib import Path
import shutil

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path.cwd()

FIG = ROOT / "figures"
DATA = ROOT / "data" / "processed"

FIG.mkdir(
    parents=True,
    exist_ok=True,
)

DATA.mkdir(
    parents=True,
    exist_ok=True,
)


SOURCES = {
    "antecedent_state_comparison.csv":
        ROOT
        / "results"
        / "phase05e3_analysis"
        / "antecedent_state_comparison.csv",

    "decomposition_response.csv":
        ROOT
        / "results"
        / "phase05e4_analysis"
        / "decomposition_response.csv",

    "mode_shape_audit.csv":
        ROOT
        / "results"
        / "phase05e5_mode_audit"
        / "mode_shape_audit.csv",

    "signal_to_drift_audit.csv":
        ROOT
        / "results"
        / "phase05e6_signal_audit"
        / "signal_to_drift_audit.csv",
}


for name, src in SOURCES.items():

    if not src.exists():

        raise RuntimeError(
            f"Missing required result: {src}"
        )

    shutil.copy2(
        src,
        DATA / name,
    )


plt.rcParams.update(
    {
        "font.size": 10.5,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "legend.fontsize": 9.5,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "figure.dpi": 120,
        "savefig.dpi": 320,
        "axes.linewidth": 0.8,
        "lines.linewidth": 2.0,
        "lines.markersize": 6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


COLORS = {
    "dry": "#B36B21",
    "reference": "#404040",
    "wet": "#277DA1",
    "stress_dry": "#D9A066",
    "stress_wet": "#75B5CF",
}


def finish(ax):

    ax.spines["top"].set_visible(
        False
    )

    ax.spines["right"].set_visible(
        False
    )

    ax.grid(
        axis="y",
        alpha=0.18,
        linewidth=0.7,
    )


def save(fig, stem):

    fig.tight_layout()

    fig.savefig(
        FIG / f"{stem}.png",
        bbox_inches="tight",
    )

    fig.savefig(
        FIG / f"{stem}.pdf",
        bbox_inches="tight",
    )

    plt.close(
        fig
    )


# ============================================================
# FIGURE 1
# ANTECEDENT HYDRAULIC STATE × TOE RECESSION
# ============================================================

df = pd.read_csv(
    DATA
    / "antecedent_state_comparison.csv"
)

x = df[
    "removed_area_m2"
]


fig, ax = plt.subplots(
    figsize=(7.0, 4.5)
)


for state, label in [
    (
        "dry",
        "Dry antecedent",
    ),
    (
        "reference",
        "Reference",
    ),
    (
        "wet",
        "Wet antecedent",
    ),
]:

    ax.plot(
        x,
        df[
            f"{state}_R_body_rms_mm"
        ],
        marker="o",
        label=label,
        color=COLORS[
            state
        ],
    )


ax.set_xlabel(
    "Actual removed toe-notch area, "
    r"$A_{\mathrm{rem}}$ [m$^2$]"
)

ax.set_ylabel(
    "Distributed slope-body RMS "
    "displacement [mm]"
)

ax.set_title(
    "Antecedent hydraulic state controls "
    "erosion-induced deformation"
)

ax.legend(
    frameon=False
)

finish(
    ax
)

save(
    fig,
    "fig01_antecedent_hydraulic_response",
)


# ============================================================
# FIGURE 2
# PRESSURE / STRESS DECOMPOSITION
# ============================================================

df = pd.read_csv(
    DATA
    / "decomposition_response.csv"
)

E = df[
    "E_m"
]

ref = df[
    "refP_refS_Rrms_mm"
]


ratios = {
    "Dry pressure + reference stress":
        (
            df[
                "dryP_refS_Rrms_mm"
            ]
            / ref
        ),

    "Wet pressure + reference stress":
        (
            df[
                "wetP_refS_Rrms_mm"
            ]
            / ref
        ),

    "Reference pressure + dry stress":
        (
            df[
                "refP_dryS_Rrms_mm"
            ]
            / ref
        ),

    "Reference pressure + wet stress":
        (
            df[
                "refP_wetS_Rrms_mm"
            ]
            / ref
        ),
}


fig, ax = plt.subplots(
    figsize=(7.2, 4.6)
)


ax.axhline(
    1.0,
    color="0.55",
    linestyle="--",
    linewidth=1.1,
    label="Reference",
)


ax.plot(
    E,
    ratios[
        "Dry pressure + reference stress"
    ],
    marker="o",
    color=COLORS[
        "dry"
    ],
    label="Dry pressure + ref. stress",
)


ax.plot(
    E,
    ratios[
        "Wet pressure + reference stress"
    ],
    marker="o",
    color=COLORS[
        "wet"
    ],
    label="Wet pressure + ref. stress",
)


ax.plot(
    E,
    ratios[
        "Reference pressure + dry stress"
    ],
    marker="s",
    linestyle=":",
    color=COLORS[
        "stress_dry"
    ],
    label="Ref. pressure + dry stress",
)


ax.plot(
    E,
    ratios[
        "Reference pressure + wet stress"
    ],
    marker="s",
    linestyle=":",
    color=COLORS[
        "stress_wet"
    ],
    label="Ref. pressure + wet stress",
)


ax.set_xlabel(
    "Nominal toe recession, "
    r"$E$ [m]"
)

ax.set_ylabel(
    r"$R_{\mathrm{RMS}}/"
    r"R_{\mathrm{RMS,ref}}$"
)

ax.set_title(
    "Hydraulic-state effect dominates "
    "inherited effective-stress effect"
)

ax.legend(
    frameon=False,
    ncol=2,
)

finish(
    ax
)

save(
    fig,
    "fig02_pressure_stress_decomposition",
)


# ============================================================
# FIGURE 3
# DEFORMATION MODE SIMILARITY
# ============================================================

df = pd.read_csv(
    DATA
    / "mode_shape_audit.csv"
)


fig, ax = plt.subplots(
    figsize=(7.0, 4.5)
)


ax.plot(
    df[
        "E_m"
    ],
    df[
        "dry_cosine"
    ],
    marker="o",
    color=COLORS[
        "dry"
    ],
    label="Dry vs reference",
)


ax.plot(
    df[
        "E_m"
    ],
    df[
        "wet_cosine"
    ],
    marker="o",
    color=COLORS[
        "wet"
    ],
    label="Wet vs reference",
)


ax.axhline(
    0.95,
    color="0.55",
    linestyle="--",
    linewidth=1.0,
    label="0.95 similarity guide",
)


ax.set_ylim(
    0.70,
    1.01,
)

ax.set_xlabel(
    "Nominal toe recession, "
    r"$E$ [m]"
)

ax.set_ylabel(
    "Displacement-vector cosine similarity"
)

ax.set_title(
    "Wet antecedent state modifies "
    "the spatial deformation mode"
)

ax.legend(
    frameon=False
)

finish(
    ax
)

save(
    fig,
    "fig03_deformation_mode_similarity",
)


# ============================================================
# FIGURE 4
# SIGNAL TO NUMERICAL DRIFT
# ============================================================

df = pd.read_csv(
    DATA
    / "signal_to_drift_audit.csv"
)


fig, ax = plt.subplots(
    figsize=(7.0, 4.5)
)


case_style = {
    "dryP_refS":
        (
            "Dry pressure",
            COLORS[
                "dry"
            ],
        ),

    "refP_refS":
        (
            "Reference",
            COLORS[
                "reference"
            ],
        ),

    "wetP_refS":
        (
            "Wet pressure",
            COLORS[
                "wet"
            ],
        ),
}


for case, (
    label,
    color,
) in case_style.items():

    part = df[
        df[
            "case"
        ] == case
    ]

    ax.plot(
        part[
            "E_m"
        ],
        part[
            "signal_to_drift"
        ],
        marker="o",
        color=color,
        label=label,
    )


ax.set_yscale(
    "log"
)

ax.axhline(
    10.0,
    color="0.55",
    linestyle="--",
    linewidth=1.0,
    label="SNR = 10",
)

ax.set_xlabel(
    "Nominal toe recession, "
    r"$E$ [m]"
)

ax.set_ylabel(
    "Signal / intact-hold numerical drift"
)

ax.set_title(
    "Deformation signals remain far above "
    "the numerical drift floor"
)

ax.legend(
    frameon=False
)

finish(
    ax
)

save(
    fig,
    "fig04_signal_to_drift",
)


# ============================================================
# FIGURE 0
# WORKFLOW
# ============================================================

fig, ax = plt.subplots(
    figsize=(9.2, 2.7)
)

ax.axis(
    "off"
)


boxes = [
    (
        0.02,
        "Seasonal\nwetting–drying",
    ),
    (
        0.22,
        "Antecedent\nhydraulic state",
    ),
    (
        0.42,
        "Coupled Richards\nmechanics + MC",
    ),
    (
        0.64,
        "Prescribed coastal\ntoe recession",
    ),
    (
        0.84,
        "Distributed\nslope response",
    ),
]


for x0, text in boxes:

    ax.text(
        x0,
        0.50,
        text,
        ha="center",
        va="center",
        transform=ax.transAxes,
        fontsize=10.5,
        bbox={
            "boxstyle":
                "round,pad=0.55",

            "facecolor":
                "white",

            "edgecolor":
                "0.25",

            "linewidth":
                1.1,
        },
    )


for x0, x1 in [
    (
        0.09,
        0.15,
    ),
    (
        0.29,
        0.35,
    ),
    (
        0.50,
        0.56,
    ),
    (
        0.72,
        0.78,
    ),
]:

    ax.annotate(
        "",
        xy=(
            x1,
            0.50,
        ),
        xytext=(
            x0,
            0.50,
        ),
        xycoords=ax.transAxes,
        arrowprops={
            "arrowstyle":
                "->",

            "linewidth":
                1.4,

            "color":
                "0.3",
        },
    )


save(
    fig,
    "fig00_workflow",
)


print(
    "========================================"
)

print(
    "PUBLICATION FIGURES COMPLETE"
)

print(
    "========================================"
)


for path in sorted(
    FIG.glob(
        "*"
    )
):

    print(
        path
    )
