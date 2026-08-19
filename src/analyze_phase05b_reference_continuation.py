from pathlib import Path
import csv

import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot
import pyvista as pv


ROOT = Path.cwd()

RUN = (
    ROOT
    / "results"
    / "phase05b_reference_continuation"
)

MODEL = (
    ROOT
    / "model"
    / "phase05b_reference_continuation"
)

OUT = (
    ROOT
    / "results"
    / "phase05b_reference_analysis"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# LOAD SERIES
# ============================================================

pvds = sorted(
    RUN.glob("*.pvd")
)

if not pvds:
    raise SystemExit(
        "FAIL: no PVD"
    )

series = ot.MeshSeries(
    str(pvds[0])
)

times = np.asarray(
    series.timevalues,
    dtype=float,
)


# ============================================================
# TOPOLOGY / MATERIAL IDS
# ============================================================

topology = pv.read(
    MODEL
    / "slope_toe_notch_ready.vtu"
)

material_ids = np.asarray(
    topology.cell_data[
        "MaterialIDs"
    ],
    dtype=int,
)

points = np.asarray(
    topology.points,
    dtype=float,
)


def surface_height(
    x,
):

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


# ============================================================
# ACTIVE-NODE MASK
# ============================================================

def active_masks(
    level,
):

    if level == 0:

        active_cells = np.ones(
            topology.n_cells,
            dtype=bool,
        )

    else:

        active_cells = ~(
            (material_ids >= 1)
            & (
                material_ids
                <= level
            )
        )

    active_nodes = np.zeros(
        topology.n_points,
        dtype=bool,
    )

    active_cell_ids = np.where(
        active_cells
    )[0]

    for cell_id in active_cell_ids:

        cell = topology.get_cell(
            int(cell_id)
        )

        active_nodes[
            cell.point_ids
        ] = True

    return (
        active_cells,
        active_nodes,
    )


# ============================================================
# MONITORING ZONE
#
# Slope body:
# x = 8..22 m
# within 4 m below local surface.
# ============================================================

surface = surface_height(
    points[:, 0]
)

depth = (
    surface
    - points[:, 1]
)

slope_zone = (
    (points[:, 0] >= 8.0)
    & (points[:, 0] <= 22.0)
    & (depth >= -1e-8)
    & (depth <= 4.0)
)


# ============================================================
# FIELD HELPERS
# ============================================================

def index_nearest(
    target,
):

    return int(
        np.argmin(
            np.abs(
                times
                - target
            )
        )
    )


def point_field(
    mesh,
    name,
):

    if name in mesh.point_data:

        return np.asarray(
            mesh.point_data[name],
            dtype=float,
        )

    raise KeyError(
        f"{name!r} not in point data"
    )


def any_field(
    mesh,
    name,
):

    if name in mesh.point_data:

        return (
            np.asarray(
                mesh.point_data[name],
                dtype=float,
            ),
            "point",
        )

    if name in mesh.cell_data:

        return (
            np.asarray(
                mesh.cell_data[name],
                dtype=float,
            ),
            "cell",
        )

    raise KeyError(
        f"{name!r} absent"
    )


# ============================================================
# PRE-EROSION REFERENCE
#
# t=9 s:
# one second before first deactivation.
# ============================================================

i0 = index_nearest(
    9.0
)

mesh0 = series.mesh(
    i0
)

u0 = point_field(
    mesh0,
    "displacement",
)


# ============================================================
# RESPONSE LEVELS
#
# Sample immediately before the NEXT event:
#
# E=0.0 -> t=9
# E=0.4 -> t=19
# E=0.8 -> t=29
# E=1.2 -> t=39
# E=1.6 -> t=49
# E=2.0 -> t=60
# ============================================================

levels = [
    (0, 0.0, 9.0),
    (1, 0.4, 19.0),
    (2, 0.8, 29.0),
    (3, 1.2, 39.0),
    (4, 1.6, 49.0),
    (5, 2.0, 60.0),
]


records = []


for level, E, target_time in levels:

    idx = index_nearest(
        target_time
    )

    mesh = series.mesh(
        idx
    )

    active_cells, active_nodes = (
        active_masks(
            level
        )
    )

    u = point_field(
        mesh,
        "displacement",
    )

    du = (
        u
        - u0
    )

    du_mag = np.linalg.norm(
        du,
        axis=1,
    )

    active_global = (
        active_nodes
    )

    active_slope = (
        active_nodes
        & slope_zone
    )

    R_global = float(
        np.nanmax(
            du_mag[
                active_global
            ]
        )
    )

    R_slope = float(
        np.nanmax(
            du_mag[
                active_slope
            ]
        )
    )


    # --------------------------------------------------------
    # Plastic strain — only ACTIVE domain.
    # --------------------------------------------------------

    epsp, epsp_loc = any_field(
        mesh,
        "EquivalentPlasticStrain",
    )

    epsp = np.asarray(
        epsp,
        dtype=float,
    ).squeeze()

    if epsp_loc == "cell":

        epsp_active = epsp[
            active_cells
        ]

    else:

        epsp_active = epsp[
            active_nodes
        ]

    epsp_active = epsp_active[
        np.isfinite(
            epsp_active
        )
    ]

    epsp_max = float(
        np.nanmax(
            epsp_active
        )
    )

    frac_1e6 = float(
        np.mean(
            epsp_active > 1e-6
        )
    )

    frac_1e4 = float(
        np.mean(
            epsp_active > 1e-4
        )
    )


    # --------------------------------------------------------
    # Pressure — ACTIVE NODES ONLY.
    # --------------------------------------------------------

    pressure = point_field(
        mesh,
        "pressure",
    ).squeeze()

    p_active = pressure[
        active_nodes
    ]

    p_active = p_active[
        np.isfinite(
            p_active
        )
    ]


    records.append(
        {
            "level":
                level,

            "E_m":
                E,

            "time_s":
                float(
                    times[idx]
                ),

            "R_global_mm":
                R_global
                * 1000.0,

            "R_slope_mm":
                R_slope
                * 1000.0,

            "epsp_max":
                epsp_max,

            "frac_epsp_gt_1e6":
                frac_1e6,

            "frac_epsp_gt_1e4":
                frac_1e4,

            "pressure_min_kpa":
                float(
                    np.nanmin(
                        p_active
                    )
                    / 1000.0
                ),

            "pressure_max_kpa":
                float(
                    np.nanmax(
                        p_active
                    )
                    / 1000.0
                ),

            "active_cells":
                int(
                    np.sum(
                        active_cells
                    )
                ),

            "active_nodes":
                int(
                    np.sum(
                        active_nodes
                    )
                ),
        }
    )


# ============================================================
# EROSION COMPLIANCE
# ============================================================

for i, r in enumerate(
    records
):

    if i == 0:

        r[
            "C_E_mm_per_m"
        ] = np.nan

        r[
            "compliance_ratio"
        ] = np.nan

        continue

    dR = (
        r["R_slope_mm"]
        - records[
            i - 1
        ]["R_slope_mm"]
    )

    dE = (
        r["E_m"]
        - records[
            i - 1
        ]["E_m"]
    )

    C = (
        dR / dE
    )

    r[
        "C_E_mm_per_m"
    ] = C

    if i <= 1:

        r[
            "compliance_ratio"
        ] = np.nan

    else:

        prev = records[
            i - 1
        ][
            "C_E_mm_per_m"
        ]

        if (
            np.isfinite(
                prev
            )
            and abs(prev) > 1e-14
        ):

            r[
                "compliance_ratio"
            ] = (
                C / prev
            )

        else:

            r[
                "compliance_ratio"
            ] = np.nan


# ============================================================
# PRINT
# ============================================================

print(
    "========================================"
)

print(
    "PHASE 05B REFERENCE EROSION CONTINUATION"
)

print(
    "========================================"
)

print()

for r in records:

    if np.isfinite(
        r["C_E_mm_per_m"]
    ):

        Ctxt = (
            f"{r['C_E_mm_per_m']:.6f}"
        )

    else:

        Ctxt = "-"

    print(
        f"E={r['E_m']:.1f} m | "
        f"R_slope={r['R_slope_mm']:.6f} mm | "
        f"C_E={Ctxt} mm/m | "
        f"epsp_max={r['epsp_max']:.3e} | "
        f"frac>1e-4="
        f"{100*r['frac_epsp_gt_1e4']:.4f}% | "
        f"active_cells={r['active_cells']}"
    )


print()
print(
    "=== COMPLIANCE RATIOS ==="
)

for r in records[2:]:

    ratio = r[
        "compliance_ratio"
    ]

    if np.isfinite(
        ratio
    ):

        print(
            f"E={r['E_m']:.1f} m | "
            f"C_i/C_prev={ratio:.6f}"
        )


# ============================================================
# CSV
# ============================================================

csv_path = (
    OUT
    / "reference_erosion_response.csv"
)

fieldnames = [
    "level",
    "E_m",
    "time_s",
    "R_global_mm",
    "R_slope_mm",
    "C_E_mm_per_m",
    "compliance_ratio",
    "epsp_max",
    "frac_epsp_gt_1e6",
    "frac_epsp_gt_1e4",
    "pressure_min_kpa",
    "pressure_max_kpa",
    "active_cells",
    "active_nodes",
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
# FIGURE 01 — RESPONSE CURVE
# ============================================================

E = np.array(
    [
        r["E_m"]
        for r in records
    ]
)

R = np.array(
    [
        r["R_slope_mm"]
        for r in records
    ]
)


fig, ax = plt.subplots(
    figsize=(7.2, 4.8)
)

ax.plot(
    E,
    R,
    marker="o",
)

ax.set_xlabel(
    "Toe recession E [m]"
)

ax.set_ylabel(
    "Erosion-induced slope response R(E) [mm]"
)

ax.set_title(
    "Reference-state toe-erosion continuation"
)

ax.grid(
    alpha=0.25
)

fig.tight_layout()

fig.savefig(
    OUT
    / "figure_01_reference_response.png",
    dpi=220,
)

plt.close(
    fig
)


# ============================================================
# FIGURE 02 — EROSION COMPLIANCE
# ============================================================

E_c = np.array(
    [
        r["E_m"]
        for r in records[1:]
    ]
)

C = np.array(
    [
        r["C_E_mm_per_m"]
        for r in records[1:]
    ]
)


fig, ax = plt.subplots(
    figsize=(7.2, 4.8)
)

ax.plot(
    E_c,
    C,
    marker="o",
)

ax.set_xlabel(
    "Toe recession E [m]"
)

ax.set_ylabel(
    "Erosion compliance C_E [mm/m]"
)

ax.set_title(
    "Incremental sensitivity to toe erosion"
)

ax.grid(
    alpha=0.25
)

fig.tight_layout()

fig.savefig(
    OUT
    / "figure_02_reference_compliance.png",
    dpi=220,
)

plt.close(
    fig
)


# ============================================================
# FIGURE 03 — PLASTIC RESPONSE
# ============================================================

epsp = np.array(
    [
        max(
            r["epsp_max"],
            1e-14,
        )
        for r in records
    ]
)


fig, ax = plt.subplots(
    figsize=(7.2, 4.8)
)

ax.semilogy(
    E,
    epsp,
    marker="o",
)

ax.set_xlabel(
    "Toe recession E [m]"
)

ax.set_ylabel(
    "Maximum active-domain equivalent plastic strain"
)

ax.set_title(
    "Plastic response during toe erosion"
)

ax.grid(
    alpha=0.25
)

fig.tight_layout()

fig.savefig(
    OUT
    / "figure_03_reference_plasticity.png",
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
    / "phase05b_reference_summary.txt"
)

lines = [
    "PHASE 05B REFERENCE EROSION CONTINUATION",
    "",
    (
        "Response definition: "
        "R(E)=max slope-zone "
        "|u(E)-u(E=0)| over active nodes."
    ),
    (
        "Compliance: "
        "C_E=Delta R / Delta E."
    ),
    "",
]

for r in records:

    c = (
        "-"
        if not np.isfinite(
            r["C_E_mm_per_m"]
        )
        else
        f"{r['C_E_mm_per_m']:.8f}"
    )

    lines.append(
        f"E={r['E_m']:.1f} m | "
        f"R={r['R_slope_mm']:.8f} mm | "
        f"C_E={c} mm/m | "
        f"epsp_max={r['epsp_max']:.8e} | "
        f"frac>1e-4="
        f"{100*r['frac_epsp_gt_1e4']:.6f}%"
    )

lines += [
    "",
    (
        "NOTE: hydraulic extrema are "
        "computed over active nodes only."
    ),
    (
        "NOTE: no Ecrit is assigned at "
        "this stage; the response curve "
        "must first be inspected for a "
        "transition/change point."
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
    "PHASE 05B REFERENCE ANALYSIS: PASS"
)
