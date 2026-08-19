from pathlib import Path
import csv
import shutil
import xml.etree.ElementTree as ET

import numpy as np
import pyvista as pv


ROOT = Path.cwd()

SRC_MODEL = (
    ROOT
    / "model"
    / "phase05b_v2_moving_front"
)

SRC_PRJ = (
    SRC_MODEL
    / "reference_moving_front.prj"
)

MODEL = (
    ROOT
    / "model"
    / "phase05c_transition_refinement"
)

OUT = (
    ROOT
    / "results"
    / "phase05c_analysis"
)

MODEL.mkdir(
    parents=True,
    exist_ok=True,
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


if not SRC_PRJ.exists():
    raise SystemExit(
        f"FAIL: missing {SRC_PRJ}"
    )


# ============================================================
# COPY MESHES
# ============================================================

for name in [
    "slope_toe_notch_ready.vtu",
    "slope_left.vtu",
    "slope_right.vtu",
    "slope_top.vtu",
    "slope_bottom.vtu",
]:

    src = SRC_MODEL / name

    if not src.exists():

        raise SystemExit(
            f"FAIL: missing {src}"
        )

    shutil.copy2(
        src,
        MODEL / name,
    )


# ============================================================
# GEOMETRIC DEACTIVATION-EVENT AUDIT
#
# OGS moving plane uses element centre.
# Therefore each unique recession coordinate
# is a discrete topology-change event.
# ============================================================

mesh = pv.read(
    MODEL
    / "slope_toe_notch_ready.vtu"
)

material_ids = np.asarray(
    mesh.cell_data["MaterialIDs"],
    dtype=int,
)

centers = (
    mesh
    .cell_centers()
    .points
)

areas_mesh = mesh.compute_cell_sizes(
    length=False,
    area=True,
    volume=False,
)

areas = np.asarray(
    areas_mesh.cell_data["Area"],
    dtype=float,
)

TOE_X = 22.0

candidate = (
    material_ids > 0
)

recession = (
    TOE_X
    - centers[:, 0]
)

audit_mask = (
    candidate
    & (recession >= 0.35)
    & (recession <= 0.60)
)

event_values = np.unique(
    np.round(
        recession[audit_mask],
        8,
    )
)

event_rows = []

print(
    "========================================"
)

print(
    "DEACTIVATION EVENTS NEAR TRANSITION"
)

print(
    "========================================"
)

for E_event in event_values:

    mask = (
        audit_mask
        & np.isclose(
            recession,
            E_event,
            atol=1e-7,
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

    mids = sorted(
        set(
            material_ids[
                mask
            ].tolist()
        )
    )

    event_rows.append(
        {
            "E_event_m":
                float(E_event),

            "cell_count":
                count,

            "area_m2":
                area,

            "material_ids":
                " ".join(
                    str(x)
                    for x in mids
                ),
        }
    )

    print(
        f"E={E_event:.6f} m | "
        f"cells={count} | "
        f"area={area:.8f} m2 | "
        f"MaterialIDs={mids}"
    )


audit_csv = (
    OUT
    / "deactivation_events_035_060.csv"
)

with audit_csv.open(
    "w",
    newline="",
    encoding="utf-8",
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "E_event_m",
            "cell_count",
            "area_m2",
            "material_ids",
        ],
    )

    writer.writeheader()
    writer.writerows(
        event_rows
    )


# ============================================================
# BUILD TWO TEMPORAL-RESOLUTION CASES
#
# Physics and erosion law are identical.
# Only dt changes.
#
# Output every 0.5 s:
# dE between saved states = 0.01 m.
# ============================================================

cases = {
    "dt0p10": {
        "dt": 0.10,
        "output_each": 5,
    },
    "dt0p05": {
        "dt": 0.05,
        "output_each": 10,
    },
}


for case, cfg in cases.items():

    tree = ET.parse(
        SRC_PRJ
    )

    root = tree.getroot()

    # --------------------------------------------------------
    # Time stepping
    # --------------------------------------------------------

    ts = root.find(
        "./time_loop/processes/"
        "process/time_stepping"
    )

    if ts is None:

        raise SystemExit(
            "FAIL: time stepping missing"
        )

    ts.clear()

    ET.SubElement(
        ts,
        "type",
    ).text = "FixedTimeStepping"

    ET.SubElement(
        ts,
        "t_initial",
    ).text = "0"

    # E = 0.7 m at t = 45 s.
    # We expect current model to fail earlier,
    # but this leaves room if refinement converges.
    ET.SubElement(
        ts,
        "t_end",
    ).text = "45"

    timesteps = ET.SubElement(
        ts,
        "timesteps",
    )

    pair = ET.SubElement(
        timesteps,
        "pair",
    )

    n_steps = int(
        round(
            45.0
            / cfg["dt"]
        )
    )

    ET.SubElement(
        pair,
        "repeat",
    ).text = str(
        n_steps
    )

    ET.SubElement(
        pair,
        "delta_t",
    ).text = str(
        cfg["dt"]
    )

    # --------------------------------------------------------
    # Output every 0.5 s.
    # --------------------------------------------------------

    out_steps = root.find(
        "./time_loop/output/timesteps"
    )

    if out_steps is None:

        raise SystemExit(
            "FAIL: output timesteps missing"
        )

    out_steps.clear()

    pair = ET.SubElement(
        out_steps,
        "pair",
    )

    ET.SubElement(
        pair,
        "repeat",
    ).text = str(
        n_steps
    )

    ET.SubElement(
        pair,
        "each_steps",
    ).text = str(
        cfg["output_each"]
    )

    prefix = root.find(
        "./time_loop/output/prefix"
    )

    if prefix is None:

        raise SystemExit(
            "FAIL: output prefix missing"
        )

    prefix.text = (
        f"phase05c_{case}"
    )

    # --------------------------------------------------------
    # Newton
    # --------------------------------------------------------

    solver = root.find(
        "./nonlinear_solvers/"
        "nonlinear_solver"
    )

    if solver is None:

        raise SystemExit(
            "FAIL: nonlinear solver missing"
        )

    max_iter = solver.find(
        "max_iter"
    )

    if max_iter is None:

        max_iter = ET.SubElement(
            solver,
            "max_iter",
        )

    max_iter.text = "100"

    # --------------------------------------------------------
    # Save
    # --------------------------------------------------------

    ET.indent(
        tree,
        space="    ",
    )

    dst = (
        MODEL
        / f"{case}.prj"
    )

    tree.write(
        dst,
        encoding="UTF-8",
        xml_declaration=True,
    )

    ET.parse(
        dst
    )

    print()
    print(
        f"PASS BUILD: {case}"
    )

    print(
        f"dt = {cfg['dt']} s"
    )

    print(
        "saved-state spacing = 0.5 s"
    )

    print(
        "saved nominal dE = 0.01 m"
    )


print()
print(
    "PHASE 05C BUILD: PASS"
)

print(
    "Audit:",
    audit_csv,
)
