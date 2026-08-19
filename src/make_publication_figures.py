from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path.cwd()

DATA = (
    ROOT
    / "data"
    / "processed"
)

FIG = (
    ROOT
    / "figures"
)

FIG.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# FORMAL JOURNAL-LIKE STYLE
#
# Explicit requirement:
# top and right axis lines ARE visible.
# ============================================================

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": [
            "Times New Roman",
            "Times",
            "DejaVu Serif",
        ],
        "mathtext.fontset": "dejavuserif",
        "font.size": 10.5,
        "axes.titlesize": 11.5,
        "axes.labelsize": 10.5,
        "legend.fontsize": 8.8,
        "xtick.labelsize": 9.2,
        "ytick.labelsize": 9.2,
        "axes.linewidth": 0.9,
        "lines.linewidth": 1.8,
        "lines.markersize": 5.5,
        "savefig.dpi": 350,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


COLORS = {
    "dry": "#9A5A24",
    "reference": "#333333",
    "wet": "#1F617D",
    "dry_light": "#C99567",
    "wet_light": "#70A3B8",
    "accent": "#465C75",
}


def formal_axes(ax):

    for spine in [
        "left",
        "right",
        "bottom",
        "top",
    ]:

        ax.spines[
            spine
        ].set_visible(
            True
        )

        ax.spines[
            spine
        ].set_linewidth(
            0.9
        )

        ax.spines[
            spine
        ].set_color(
            "0.15"
        )


    ax.tick_params(
        axis="both",
        which="both",
        direction="in",
        top=True,
        right=True,
        length=4.0,
        width=0.8,
    )

    ax.grid(
        False
    )


def panel_label(
    ax,
    label,
):

    ax.text(
        -0.12,
        1.04,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontweight="bold",
        fontsize=11,
    )


def save(
    fig,
    stem,
):

    fig.savefig(
        FIG
        / f"{stem}.png",
        bbox_inches="tight",
        facecolor="white",
    )

    fig.savefig(
        FIG
        / f"{stem}.pdf",
        bbox_inches="tight",
        facecolor="white",
    )

    plt.close(
        fig
    )


# ============================================================
# LOAD VERIFIED PROCESSED RESULTS
# ============================================================

antecedent = pd.read_csv(
    DATA
    / "antecedent_state_comparison.csv"
)

decomp = pd.read_csv(
    DATA
    / "decomposition_response.csv"
)

mode = pd.read_csv(
    DATA
    / "mode_shape_audit.csv"
)

signal = pd.read_csv(
    DATA
    / "signal_to_drift_audit.csv"
)


# ============================================================
# FIGURE 0 — GRAPHICAL ABSTRACT
# ============================================================

fig, ax = plt.subplots(
    figsize=(10.0, 2.9)
)

ax.set_xlim(
    0,
    1,
)

ax.set_ylim(
    0,
    1,
)

ax.axis(
    "off"
)


boxes = [
    (
        0.09,
        "Seasonal\nwetting–drying",
        "Hydraulic forcing",
    ),
    (
        0.29,
        "Antecedent\nsoil-water state",
        "Dry / reference / wet",
    ),
    (
        0.50,
        "Coupled HM\nslope model",
        "Richards + deformation",
    ),
    (
        0.71,
        "Prescribed basal\ntoe recession",
        "Element deactivation",
    ),
    (
        0.91,
        "Distributed\nslope response",
        "Magnitude + mode",
    ),
]


for x, title, subtitle in boxes:

    ax.text(
        x,
        0.58,
        title,
        ha="center",
        va="center",
        fontsize=10.5,
        fontweight="bold",
        bbox={
            "boxstyle":
                "round,pad=0.6",

            "facecolor":
                "#F8F8F8",

            "edgecolor":
                "#333333",

            "linewidth":
                0.9,
        },
    )

    ax.text(
        x,
        0.22,
        subtitle,
        ha="center",
        va="center",
        fontsize=8.2,
        color="0.30",
    )


for x0, x1 in [
    (
        0.16,
        0.22,
    ),
    (
        0.36,
        0.43,
    ),
    (
        0.57,
        0.64,
    ),
    (
        0.78,
        0.84,
    ),
]:

    ax.annotate(
        "",
        xy=(
            x1,
            0.58,
        ),
        xytext=(
            x0,
            0.58,
        ),
        arrowprops={
            "arrowstyle":
                "->",

            "linewidth":
                1.2,

            "color":
                "0.25",
        },
    )


ax.text(
    0.50,
    0.96,
    (
        "TillCliff-HM: seasonal hydro-mechanical "
        "preconditioning under coastal toe recession"
    ),
    ha="center",
    va="top",
    fontsize=12.2,
    fontweight="bold",
)


save(
    fig,
    "fig00_graphical_abstract",
)


# ============================================================
# FIGURE 1 — MAIN RESPONSE + ATTRIBUTION
# ============================================================

fig, axes = plt.subplots(
    1,
    2,
    figsize=(10.6, 4.25),
)


# ------------------------------------------------------------
# Panel a:
# Production seasonal branches
# ------------------------------------------------------------

ax = axes[
    0
]

x = antecedent[
    "removed_area_m2"
]


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
        antecedent[
            f"{state}_R_body_rms_mm"
        ],
        marker="o",
        color=COLORS[
            state
        ],
        label=label,
    )


ax.set_xlabel(
    r"Actual removed area, "
    r"$A_{\mathrm{rem}}$ [m$^2$]"
)

ax.set_ylabel(
    r"Slope-body "
    r"$R_{\mathrm{RMS}}$ [mm]"
)

ax.set_title(
    "Seasonal antecedent branches"
)

ax.legend(
    frameon=False,
    loc="upper left",
)

formal_axes(
    ax
)

panel_label(
    ax,
    "(a)",
)


# ------------------------------------------------------------
# Panel b:
# Controlled pressure / inherited stress decomposition
# ------------------------------------------------------------

ax = axes[
    1
]

E = decomp[
    "E_m"
]

ref = decomp[
    "refP_refS_Rrms_mm"
]


ax.axhline(
    1.0,
    color="0.30",
    linestyle="--",
    linewidth=1.0,
    label="Reference",
)


ax.plot(
    E,
    decomp[
        "dryP_refS_Rrms_mm"
    ]
    / ref,
    marker="o",
    color=COLORS[
        "dry"
    ],
    label="Dry pressure + ref. stress",
)


ax.plot(
    E,
    decomp[
        "wetP_refS_Rrms_mm"
    ]
    / ref,
    marker="o",
    color=COLORS[
        "wet"
    ],
    label="Wet pressure + ref. stress",
)


ax.plot(
    E,
    decomp[
        "refP_dryS_Rrms_mm"
    ]
    / ref,
    marker="s",
    linestyle=":",
    color=COLORS[
        "dry_light"
    ],
    label="Ref. pressure + dry stress",
)


ax.plot(
    E,
    decomp[
        "refP_wetS_Rrms_mm"
    ]
    / ref,
    marker="s",
    linestyle=":",
    color=COLORS[
        "wet_light"
    ],
    label="Ref. pressure + wet stress",
)


ax.set_xlabel(
    r"Nominal toe recession, $E$ [m]"
)

ax.set_ylabel(
    r"$R_{\mathrm{RMS}}/"
    r"R_{\mathrm{RMS,ref}}$"
)

ax.set_title(
    "Hydraulic / stress-state attribution"
)

ax.legend(
    frameon=False,
    fontsize=7.8,
    loc="best",
)

formal_axes(
    ax
)

panel_label(
    ax,
    "(b)",
)


fig.subplots_adjust(
    left=0.085,
    right=0.985,
    bottom=0.16,
    top=0.87,
    wspace=0.28,
)


save(
    fig,
    "fig01_hydraulic_state_response",
)


# ============================================================
# FIGURE 2 — DEFORMATION MODE
# ============================================================

fig, axes = plt.subplots(
    1,
    2,
    figsize=(10.6, 4.25),
)


# ------------------------------------------------------------
# Panel a — cosine similarity
# ------------------------------------------------------------

ax = axes[
    0
]


ax.plot(
    mode[
        "E_m"
    ],
    mode[
        "dry_cosine"
    ],
    marker="o",
    color=COLORS[
        "dry"
    ],
    label="Dry vs reference",
)


ax.plot(
    mode[
        "E_m"
    ],
    mode[
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
    linestyle="--",
    linewidth=1.0,
    color="0.45",
    label="0.95 guide",
)


ax.set_ylim(
    0.70,
    1.01,
)

ax.set_xlabel(
    r"Nominal toe recession, $E$ [m]"
)

ax.set_ylabel(
    "Vector-field cosine similarity"
)

ax.set_title(
    "Displacement-mode similarity"
)

ax.legend(
    frameon=False,
)

formal_axes(
    ax
)

panel_label(
    ax,
    "(a)",
)


# ------------------------------------------------------------
# Panel b — residual after scalar best fit
# ------------------------------------------------------------

ax = axes[
    1
]


ax.plot(
    mode[
        "E_m"
    ],
    mode[
        "dry_shape_residual"
    ],
    marker="o",
    color=COLORS[
        "dry"
    ],
    label="Dry vs reference",
)


ax.plot(
    mode[
        "E_m"
    ],
    mode[
        "wet_shape_residual"
    ],
    marker="o",
    color=COLORS[
        "wet"
    ],
    label="Wet vs reference",
)


ax.set_xlabel(
    r"Nominal toe recession, $E$ [m]"
)

ax.set_ylabel(
    "Normalized mode-shape residual"
)

ax.set_title(
    "Residual after best-fit amplitude scaling"
)

ax.legend(
    frameon=False,
)

formal_axes(
    ax
)

panel_label(
    ax,
    "(b)",
)


fig.subplots_adjust(
    left=0.085,
    right=0.985,
    bottom=0.16,
    top=0.87,
    wspace=0.28,
)


save(
    fig,
    "fig02_deformation_mode",
)


# ============================================================
# FIGURE 3 — NUMERICAL SIGNAL QUALITY
# ============================================================

fig, ax = plt.subplots(
    figsize=(6.7, 4.5)
)


styles = {
    "dryP_refS": (
        "Dry hydraulic state",
        COLORS[
            "dry"
        ],
    ),
    "refP_refS": (
        "Reference",
        COLORS[
            "reference"
        ],
    ),
    "wetP_refS": (
        "Wet hydraulic state",
        COLORS[
            "wet"
        ],
    ),
}


for case, (
    label,
    color,
) in styles.items():

    part = signal[
        signal[
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


ax.axhline(
    10.0,
    linestyle="--",
    linewidth=1.0,
    color="0.45",
    label="SNR = 10 guide",
)


ax.set_yscale(
    "log"
)

ax.set_xlabel(
    r"Nominal toe recession, $E$ [m]"
)

ax.set_ylabel(
    "Erosion signal / intact-hold drift"
)

ax.set_title(
    "Resolved deformation signal relative "
    "to numerical drift"
)

ax.legend(
    frameon=False,
    loc="best",
)

formal_axes(
    ax
)


fig.subplots_adjust(
    left=0.14,
    right=0.97,
    bottom=0.15,
    top=0.88,
)


save(
    fig,
    "fig03_signal_quality",
)


print(
    "============================================================"
)

print(
    "FORMAL RESEARCH FIGURES GENERATED"
)

print(
    "============================================================"
)


for path in sorted(
    FIG.glob(
        "fig0*"
    )
):

    print(
        path
    )
