from pathlib import Path
import shutil
import xml.etree.ElementTree as ET


ROOT = Path.cwd()

SRC_PRJ = (
    ROOT
    / "model"
    / "phase04d_v2"
    / "screening.prj"
)

TOE_MODEL = (
    ROOT
    / "model"
    / "phase05a_v2_toe_notch"
)

MODEL = (
    ROOT
    / "model"
    / "phase05b_reference_continuation"
)

MODEL.mkdir(
    parents=True,
    exist_ok=True,
)

DST = (
    MODEL
    / "reference_continuation.prj"
)


# ============================================================
# COPY EROSION-READY INITIALIZED MESH
# ============================================================

for name in [
    "slope_toe_notch_ready.vtu",
    "slope_left.vtu",
    "slope_right.vtu",
    "slope_top.vtu",
    "slope_bottom.vtu",
]:

    src = TOE_MODEL / name

    if not src.exists():
        raise SystemExit(
            f"FAIL: missing {src}"
        )

    shutil.copy2(
        src,
        MODEL / name,
    )


# ============================================================
# LOAD VERIFIED MC-HM PROJECT
# ============================================================

tree = ET.parse(
    SRC_PRJ
)

root = tree.getroot()


# ============================================================
# BULK MESH
# ============================================================

meshes = root.find(
    "meshes"
)

if meshes is None:
    raise SystemExit(
        "FAIL: meshes missing"
    )

bulk_replaced = False

for mesh in meshes.findall(
    "mesh"
):

    text = (
        mesh.text.strip()
        if mesh.text
        else ""
    )

    if text in {
        "slope_initialized.vtu",
        "slope.vtu",
    }:

        mesh.text = (
            "slope_toe_notch_ready.vtu"
        )

        bulk_replaced = True

if not bulk_replaced:
    raise SystemExit(
        "FAIL: bulk mesh not replaced"
    )


# ============================================================
# SAME CONSTITUTIVE MODEL FOR IDs 0..5
# ============================================================

process = None

for candidate in root.findall(
    "./processes/process"
):

    if (
        candidate.findtext("type")
        == "RICHARDS_MECHANICS"
    ):
        process = candidate
        break

if process is None:
    raise SystemExit(
        "FAIL: RM process missing"
    )

cr = process.find(
    "constitutive_relation"
)

if cr is None:
    raise SystemExit(
        "FAIL: constitutive relation missing"
    )

cr.set(
    "id",
    "0,1,2,3,4,5",
)


medium = root.find(
    "./media/medium"
)

if medium is None:
    raise SystemExit(
        "FAIL: medium missing"
    )

medium.set(
    "id",
    "0,1,2,3,4,5",
)


# ============================================================
# SEQUENTIAL DEACTIVATION
#
# t=10 -> ID1 -> E=0.4
# t=20 -> ID2 -> E=0.8
# t=30 -> ID3 -> E=1.2
# t=40 -> ID4 -> E=1.6
# t=50 -> ID5 -> E=2.0
#
# End=60 s gives 10 s hold after final event.
#
# Numerical event time is NOT physical erosion duration.
# ============================================================

root_pvs = root.find(
    "process_variables"
)

if root_pvs is None:
    raise SystemExit(
        "FAIL: root process_variables missing"
    )


for pv in root_pvs.findall(
    "process_variable"
):

    name = pv.findtext(
        "name"
    )

    if name not in {
        "displacement",
        "pressure",
    }:
        continue

    old = pv.find(
        "deactivated_subdomains"
    )

    if old is not None:
        pv.remove(
            old
        )

    ds = ET.SubElement(
        pv,
        "deactivated_subdomains",
    )

    for material_id in range(
        1,
        6,
    ):

        start = (
            material_id
            * 10
        )

        sub = ET.SubElement(
            ds,
            "deactivated_subdomain",
        )

        interval = ET.SubElement(
            sub,
            "time_interval",
        )

        ET.SubElement(
            interval,
            "start",
        ).text = str(
            start
        )

        ET.SubElement(
            interval,
            "end",
        ).text = "60"

        ET.SubElement(
            sub,
            "material_ids",
        ).text = str(
            material_id
        )


# ============================================================
# TIME STEPPING
# ============================================================

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

ET.SubElement(
    ts,
    "t_end",
).text = "60"

steps = ET.SubElement(
    ts,
    "timesteps",
)

pair = ET.SubElement(
    steps,
    "pair",
)

ET.SubElement(
    pair,
    "repeat",
).text = "60"

ET.SubElement(
    pair,
    "delta_t",
).text = "1"


# ============================================================
# NEWTON
# ============================================================

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


# ============================================================
# OUTPUT EVERY STEP
# ============================================================

prefix = root.find(
    "./time_loop/output/prefix"
)

if prefix is None:
    raise SystemExit(
        "FAIL: prefix missing"
    )

prefix.text = (
    "phase05b_reference"
)


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
).text = "60"

ET.SubElement(
    pair,
    "each_steps",
).text = "1"


variables = root.find(
    "./time_loop/output/variables"
)

if variables is None:
    raise SystemExit(
        "FAIL: output variables missing"
    )

existing = {
    v.text
    for v in variables.findall(
        "variable"
    )
}

for name in [
    "displacement",
    "pressure",
    "saturation",
    "sigma",
    "EquivalentPlasticStrain",
]:

    if name not in existing:

        ET.SubElement(
            variables,
            "variable",
        ).text = name


# ============================================================
# SAVE
# ============================================================

ET.indent(
    tree,
    space="    ",
)

tree.write(
    DST,
    encoding="UTF-8",
    xml_declaration=True,
)

ET.parse(
    DST
)


print(
    "PHASE 05B REFERENCE CONTINUATION BUILD: PASS"
)

print()

for material_id in range(
    1,
    6,
):

    print(
        f"t={material_id*10:2d} s | "
        f"remove ID={material_id} | "
        f"E={material_id*0.4:.1f} m"
    )

print()
print(
    "No boundary_parameter is imposed "
    "on the new notch surfaces."
)
