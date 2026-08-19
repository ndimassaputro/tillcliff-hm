from pathlib import Path
import csv
import shutil

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv


ROOT = Path.cwd()

SRC = (
    ROOT
    / "model"
    / "phase04d_v2"
    / "slope_initialized.vtu"
)

MODEL = (
    ROOT
    / "model"
    / "phase05a_v2_toe_notch"
)

OUT = (
    ROOT
    / "results"
    / "phase05a_v2_toe_notch"
)

MODEL.mkdir(
    parents=True,
    exist_ok=True,
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)

if not SRC.exists():
    raise SystemExit(
        f"FAIL: missing {SRC}"
    )


# ============================================================
# PHYSICAL TOE-NOTCH DEFINITION
# ============================================================

TOE_X = 22.0
TOE_ELEVATION = 2.0

BAND_WIDTH = 0.4
N_BANDS = 5

# Coastal toe attack window:
# 2.0 m to 3.2 m elevation.
NOTCH_HEIGHT = 1.2
NOTCH_TOP = (
    TOE_ELEVATION
    + NOTCH_HEIGHT
)


mesh = pv.read(
    SRC
)

centers = (
    mesh
    .cell_centers()
    .points
)

x = centers[:, 0]
y = centers[:, 1]


# ============================================================
# MATERIAL IDS
#
# 0 = retained material
# 1 = E 0.0 -> 0.4 m
# 2 = E 0.4 -> 0.8 m
# ...
# 5 = E 1.6 -> 2.0 m
# ============================================================

material_ids = np.zeros(
    mesh.n_cells,
    dtype=np.int32,
)

recession = (
    TOE_X
    - x
)

candidate = (
    (recession >= 0.0)
    & (
        recession
        < N_BANDS * BAND_WIDTH
    )
    & (y >= TOE_ELEVATION)
    & (y <= NOTCH_TOP)
)

band = (
    np.floor(
        recession
        / BAND_WIDTH
    )
    .astype(int)
    + 1
)

valid = (
    candidate
    & (band >= 1)
    & (band <= N_BANDS)
)

material_ids[
    valid
] = band[
    valid
]

mesh.cell_data[
    "MaterialIDs"
] = material_ids


# ============================================================
# AREAS
# ============================================================

mesh_area = mesh.compute_cell_sizes(
    length=False,
    area=True,
    volume=False,
)

areas = np.asarray(
    mesh_area.cell_data["Area"],
    dtype=float,
)


# ============================================================
# BAND METRICS
# ============================================================

rows = []

print(
    "========================================"
)

print(
    "PHASE 05A-V2 BASAL TOE NOTCH"
)

print(
    "========================================"
)

print()

print(
    f"Toe coordinate x = {TOE_X:.3f} m"
)

print(
    f"Toe elevation = {TOE_ELEVATION:.3f} m"
)

print(
    f"Notch elevation window = "
    f"{TOE_ELEVATION:.3f} to "
    f"{NOTCH_TOP:.3f} m"
)

print(
    f"Maximum recession E = "
    f"{N_BANDS*BAND_WIDTH:.3f} m"
)

print()

all_present = True

for mid in range(
    0,
    N_BANDS + 1,
):

    mask = (
        material_ids
        == mid
    )

    count = int(
        np.sum(mask)
    )

    area = float(
        np.sum(
            areas[mask]
        )
    )

    if (
        mid > 0
        and count == 0
    ):
        all_present = False

    if mid == 0:

        label = "retained"
        e = 0.0

    else:

        label = (
            f"toe_notch_band_{mid}"
        )

        e = (
            mid
            * BAND_WIDTH
        )

    rows.append(
        {
            "material_id":
                mid,

            "label":
                label,

            "cell_count":
                count,

            "area_m2":
                area,

            "cumulative_E_m":
                e,
        }
    )

    print(
        f"MaterialID={mid} | "
        f"cells={count} | "
        f"area={area:.6f} m2 | "
        f"E={e:.3f} m"
    )


if not all_present:

    raise SystemExit(
        "FAIL: one or more notch bands empty"
    )


# ============================================================
# CUMULATIVE REMOVAL
# ============================================================

print()
print(
    "=== CUMULATIVE TOE REMOVAL ==="
)

levels = []

for level in range(
    0,
    N_BANDS + 1,
):

    if level == 0:

        mask = np.zeros(
            mesh.n_cells,
            dtype=bool,
        )

    else:

        mask = (
            (material_ids >= 1)
            & (
                material_ids
                <= level
            )
        )

    count = int(
        np.sum(mask)
    )

    area = float(
        np.sum(
            areas[mask]
        )
    )

    E = (
        level
        * BAND_WIDTH
    )

    levels.append(
        {
            "level":
                level,

            "E_m":
                E,

            "removed_cells":
                count,

            "removed_area_m2":
                area,
        }
    )

    print(
        f"E={E:.3f} m | "
        f"removed_cells={count} | "
        f"removed_area={area:.6f} m2"
    )


# ============================================================
# SANITY CHECK:
# cumulative area must increase monotonically.
# ============================================================

cum_areas = np.array(
    [
        r["removed_area_m2"]
        for r in levels
    ]
)

if not np.all(
    np.diff(cum_areas)
    > 0.0
):

    raise SystemExit(
        "FAIL: cumulative removed area "
        "is not strictly increasing"
    )


# ============================================================
# SAVE MESH
# ============================================================

DST = (
    MODEL
    / "slope_toe_notch_ready.vtu"
)

mesh.save(
    DST
)


# ============================================================
# COPY BOUNDARIES
# ============================================================

for name in [
    "slope_left.vtu",
    "slope_right.vtu",
    "slope_top.vtu",
    "slope_bottom.vtu",
]:

    src = (
        ROOT
        / "model"
        / "phase04d_v2"
        / name
    )

    if not src.exists():

        raise SystemExit(
            f"FAIL: missing {src}"
        )

    shutil.copy2(
        src,
        MODEL / name,
    )


# ============================================================
# WRITE CSVs
# ============================================================

with (
    OUT
    / "toe_notch_bands.csv"
).open(
    "w",
    newline="",
    encoding="utf-8",
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "material_id",
            "label",
            "cell_count",
            "area_m2",
            "cumulative_E_m",
        ],
    )

    writer.writeheader()
    writer.writerows(
        rows
    )


with (
    OUT
    / "toe_notch_levels.csv"
).open(
    "w",
    newline="",
    encoding="utf-8",
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "level",
            "E_m",
            "removed_cells",
            "removed_area_m2",
        ],
    )

    writer.writeheader()
    writer.writerows(
        levels
    )


# ============================================================
# FIGURE
# ============================================================

fig, ax = plt.subplots(
    figsize=(9.5, 4.8)
)

scatter = ax.scatter(
    centers[:, 0],
    centers[:, 1],
    c=material_ids,
    s=7,
)

ax.axhline(
    TOE_ELEVATION,
    linestyle="--",
    linewidth=1,
)

ax.axhline(
    NOTCH_TOP,
    linestyle="--",
    linewidth=1,
)

ax.set_xlim(
    18.5,
    23.5,
)

ax.set_ylim(
    0.5,
    4.5,
)

ax.set_aspect(
    "equal"
)

ax.set_xlabel(
    "x [m]"
)

ax.set_ylabel(
    "Elevation [m]"
)

ax.set_title(
    "Basal coastal toe-notch discretization"
)

cbar = fig.colorbar(
    scatter,
    ax=ax,
)

cbar.set_label(
    "Material ID"
)

fig.tight_layout()

FIG = (
    OUT
    / "figure_01_basal_toe_notch.png"
)

fig.savefig(
    FIG,
    dpi=220,
)

plt.close(
    fig
)


# ============================================================
# SUMMARY
# ============================================================

candidate_mask = (
    material_ids > 0
)

candidate_area = float(
    np.sum(
        areas[
            candidate_mask
        ]
    )
)

domain_area = float(
    np.sum(
        areas
    )
)

summary = (
    OUT
    / "phase05a_v2_summary.txt"
)

lines = [
    "PHASE 05A-V2 BASAL COASTAL TOE NOTCH",
    "",
    f"Toe x [m]: {TOE_X:.3f}",
    (
        "Toe elevation [m]: "
        f"{TOE_ELEVATION:.3f}"
    ),
    (
        "Notch top elevation [m]: "
        f"{NOTCH_TOP:.3f}"
    ),
    (
        "Notch height [m]: "
        f"{NOTCH_HEIGHT:.3f}"
    ),
    (
        "Band width [m]: "
        f"{BAND_WIDTH:.3f}"
    ),
    (
        "Maximum recession E [m]: "
        f"{N_BANDS*BAND_WIDTH:.3f}"
    ),
    "",
]

for r in rows:

    lines.append(
        f"MaterialID={r['material_id']} | "
        f"cells={r['cell_count']} | "
        f"area={r['area_m2']:.6f} m2 | "
        f"E={r['cumulative_E_m']:.3f} m"
    )

lines += [
    "",
    (
        "Total candidate notch area [m2]: "
        f"{candidate_area:.8f}"
    ),
    (
        "Notch/domain area [%]: "
        f"{100*candidate_area/domain_area:.6f}"
    ),
    "",
    "BASAL TOE-NOTCH TOPOLOGY: PASS",
    "",
    (
        "Interpretation: erosion is confined "
        "to the basal coastal toe elevation window."
    ),
    (
        "E denotes horizontal inland recession."
    ),
    (
        "This remains process-informed geometric "
        "erosion, not wave-resolved sediment transport."
    ),
]

summary.write_text(
    "\n".join(lines) + "\n",
    encoding="utf-8",
)


print()
print(
    "PASS:",
    DST,
)

print(
    "PASS:",
    FIG,
)

print()
print(
    "BASAL TOE-NOTCH TOPOLOGY: PASS"
)
