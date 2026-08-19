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

OUT_MODEL = (
    ROOT
    / "model"
    / "phase05a_erosion_mesh"
)

OUT_RESULTS = (
    ROOT
    / "results"
    / "phase05a_erosion_mesh"
)

OUT_MODEL.mkdir(
    parents=True,
    exist_ok=True,
)

OUT_RESULTS.mkdir(
    parents=True,
    exist_ok=True,
)


if not SRC.exists():
    raise SystemExit(
        f"FAIL: missing {SRC}"
    )


# ============================================================
# SETTINGS
# ============================================================

TOE_X = 22.0

# Coarse erosion continuation:
#
# E = 0.0 m
# E = 0.4 m
# E = 0.8 m
# E = 1.2 m
# E = 1.6 m
# E = 2.0 m
#
# Each material band represents another
# 0.4 m horizontal recession increment.

BAND_WIDTH = 0.4
N_BANDS = 5

# Only the shallow near-surface toe zone is removed.
# This represents a process-informed toe notch/recession,
# NOT wave-resolved sediment transport.
NOTCH_DEPTH = 1.2


def surface_height(x):

    x = np.asarray(
        x,
        dtype=float,
    )

    h = np.full_like(
        x,
        10.0,
    )

    slope = (
        (x > 8.0)
        & (x < 22.0)
    )

    h[slope] = (
        10.0
        - 8.0
        * (x[slope] - 8.0)
        / 14.0
    )

    h[x >= 22.0] = 2.0

    return h


# ============================================================
# LOAD INITIALIZED HM MESH
# ============================================================

mesh = pv.read(
    SRC
)

print(
    "========================================"
)

print(
    "PHASE 05A EROSION-READY TOE MESH"
)

print(
    "========================================"
)

print()

print(
    f"Input cells : {mesh.n_cells}"
)

print(
    f"Input points: {mesh.n_points}"
)


# ============================================================
# CELL GEOMETRY
# ============================================================

centers = (
    mesh
    .cell_centers()
    .points
)

x = centers[:, 0]
y = centers[:, 1]

surface = surface_height(
    x
)

depth = (
    surface
    - y
)

recession_from_toe = (
    TOE_X
    - x
)


# ============================================================
# MATERIAL IDS
#
# 0 = retained soil
# 1 = first erosion increment
# 2 = second
# ...
# 5 = fifth
# ============================================================

material_ids = np.zeros(
    mesh.n_cells,
    dtype=np.int32,
)


toe_candidate = (
    (recession_from_toe >= 0.0)
    & (
        recession_from_toe
        < N_BANDS * BAND_WIDTH
    )
    & (depth >= 0.0)
    & (depth <= NOTCH_DEPTH)
)


band = (
    np.floor(
        recession_from_toe
        / BAND_WIDTH
    )
    .astype(int)
    + 1
)


valid_band = (
    toe_candidate
    & (band >= 1)
    & (band <= N_BANDS)
)

material_ids[
    valid_band
] = band[
    valid_band
]


mesh.cell_data[
    "MaterialIDs"
] = material_ids


# ============================================================
# CELL AREAS
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
# VALIDATION
# ============================================================

print()
print(
    "=== MATERIAL-ID CHECK ==="
)

rows = []

all_bands_present = True

for material_id in range(
    0,
    N_BANDS + 1,
):

    mask = (
        material_ids
        == material_id
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
        material_id > 0
        and count == 0
    ):

        all_bands_present = False

    if material_id == 0:

        label = "retained"

        recession = 0.0

    else:

        label = (
            f"erosion_band_{material_id}"
        )

        recession = (
            material_id
            * BAND_WIDTH
        )

    rows.append(
        {
            "material_id":
                material_id,

            "label":
                label,

            "cell_count":
                count,

            "area_m2":
                area,

            "cumulative_E_m":
                recession,
        }
    )

    print(
        f"Material ID {material_id}: "
        f"cells={count:5d} | "
        f"area={area:.6f} m2"
    )


erosion_mask = (
    material_ids > 0
)

erosion_cells = int(
    np.sum(
        erosion_mask
    )
)

erosion_area = float(
    np.sum(
        areas[
            erosion_mask
        ]
    )
)

domain_area = float(
    np.sum(
        areas
    )
)


print()
print(
    "=== EROSION REGION ==="
)

print(
    f"Toe recession range: "
    f"0 to {N_BANDS*BAND_WIDTH:.2f} m"
)

print(
    f"Nominal notch depth: "
    f"{NOTCH_DEPTH:.2f} m"
)

print(
    f"Total candidate erosion cells: "
    f"{erosion_cells}"
)

print(
    f"Total candidate erosion area: "
    f"{erosion_area:.6f} m2"
)

print(
    "Candidate erosion area / domain: "
    f"{100*erosion_area/domain_area:.4f}%"
)


if not all_bands_present:

    raise SystemExit(
        "FAIL: one or more erosion bands "
        "contain zero cells"
    )

if erosion_cells <= 0:

    raise SystemExit(
        "FAIL: no erosion cells selected"
    )


# ============================================================
# SAVE EROSION-READY MESH
# ============================================================

DST = (
    OUT_MODEL
    / "slope_erosion_ready.vtu"
)

mesh.save(
    DST
)

print()
print(
    "PASS:",
    DST,
)


# ============================================================
# COPY BOUNDARY MESHES
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

    dst = (
        OUT_MODEL
        / name
    )

    if not src.exists():

        raise SystemExit(
            f"FAIL: missing {src}"
        )

    shutil.copy2(
        src,
        dst,
    )


# ============================================================
# MANIFEST
# ============================================================

manifest = (
    OUT_RESULTS
    / "toe_erosion_bands.csv"
)

with manifest.open(
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


# ============================================================
# CUMULATIVE REMOVAL TABLE
# ============================================================

cumulative_path = (
    OUT_RESULTS
    / "toe_erosion_levels.csv"
)

with cumulative_path.open(
    "w",
    newline="",
    encoding="utf-8",
) as f:

    writer = csv.writer(
        f
    )

    writer.writerow(
        [
            "level",
            "E_m",
            "material_ids_removed",
            "removed_cells",
            "removed_area_m2",
        ]
    )

    # E = 0
    writer.writerow(
        [
            0,
            0.0,
            "",
            0,
            0.0,
        ]
    )

    for level in range(
        1,
        N_BANDS + 1,
    ):

        mask = (
            (material_ids >= 1)
            & (
                material_ids
                <= level
            )
        )

        removed_cells = int(
            np.sum(mask)
        )

        removed_area = float(
            np.sum(
                areas[mask]
            )
        )

        ids = ",".join(
            str(i)
            for i in range(
                1,
                level + 1,
            )
        )

        writer.writerow(
            [
                level,
                level * BAND_WIDTH,
                ids,
                removed_cells,
                removed_area,
            ]
        )


# ============================================================
# FIGURE
# ============================================================

fig, ax = plt.subplots(
    figsize=(10, 4.5)
)

scatter = ax.scatter(
    centers[:, 0],
    centers[:, 1],
    c=material_ids,
    s=6,
)

ax.set_aspect(
    "equal"
)

ax.set_xlim(
    7.0,
    24.0,
)

ax.set_ylim(
    0.0,
    11.0,
)

ax.set_xlabel(
    "x [m]"
)

ax.set_ylabel(
    "Elevation [m]"
)

ax.set_title(
    "Phase 05A — Toe erosion material bands"
)

cbar = fig.colorbar(
    scatter,
    ax=ax,
)

cbar.set_label(
    "Material ID"
)

fig.tight_layout()

figure = (
    OUT_RESULTS
    / "figure_01_toe_erosion_bands.png"
)

fig.savefig(
    figure,
    dpi=220,
)

plt.close(
    fig
)


# ============================================================
# SUMMARY
# ============================================================

summary = (
    OUT_RESULTS
    / "phase05a_summary.txt"
)

lines = [
    "PHASE 05A EROSION-READY TOE MESH",
    "",
    (
        f"Toe x-coordinate [m]: "
        f"{TOE_X:.3f}"
    ),
    (
        f"Band width [m]: "
        f"{BAND_WIDTH:.3f}"
    ),
    (
        f"Number of bands: "
        f"{N_BANDS}"
    ),
    (
        f"Maximum coarse recession E [m]: "
        f"{N_BANDS*BAND_WIDTH:.3f}"
    ),
    (
        f"Nominal notch depth [m]: "
        f"{NOTCH_DEPTH:.3f}"
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
        f"Total candidate erosion cells: "
        f"{erosion_cells}"
    ),
    (
        f"Total candidate erosion area [m2]: "
        f"{erosion_area:.8f}"
    ),
    (
        "Erosion/domain area [%]: "
        f"{100*erosion_area/domain_area:.6f}"
    ),
    "",
    "EROSION-BAND TOPOLOGY: PASS",
    "",
    (
        "Interpretation: material IDs define "
        "sequential toe-recession increments."
    ),
    (
        "These bands represent process-informed "
        "geometric erosion, not wave-resolved "
        "sediment transport."
    ),
]

summary.write_text(
    "\n".join(lines) + "\n",
    encoding="utf-8",
)


print()
print(
    "PASS:",
    manifest,
)

print(
    "PASS:",
    cumulative_path,
)

print(
    "PASS:",
    figure,
)

print(
    "PASS:",
    summary,
)

print()
print(
    "EROSION-BAND TOPOLOGY: PASS"
)
