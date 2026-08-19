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
    / "phase05b_v2_moving_front"
)

MODEL.mkdir(
    parents=True,
    exist_ok=True,
)

DST = (
    MODEL
    / "reference_moving_front.prj"
)


# ============================================================
# CHECK / COPY VERIFIED MESHES
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


if not SRC_PRJ.exists():

    raise SystemExit(
        f"FAIL: missing {SRC_PRJ}"
    )


# ============================================================
# LOAD VERIFIED INITIALIZED MC-HM PROJECT
# ============================================================

tree = ET.parse(
    SRC_PRJ
)

root = tree.getroot()


# ============================================================
# USE TOE-NOTCH MATERIAL-ID MESH
# ============================================================

meshes = root.find(
    "meshes"
)

if meshes is None:

    raise SystemExit(
        "FAIL: meshes block missing"
    )


replaced = False

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

        replaced = True


if not replaced:

    raise SystemExit(
        "FAIL: bulk mesh reference "
        "not replaced"
    )


# ============================================================
# RICHARDS_MECHANICS
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
        "FAIL: RICHARDS_MECHANICS "
        "process missing"
    )


# ============================================================
# SAME CONSTITUTIVE MODEL FOR IDs 0..5
# ============================================================

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
# MOVING EROSION FRONT
#
# OGS line-segment logic:
#
# start = x=22
# end   = x=20
#
# direction points inland.
#
# curve position:
# 0 m -> no recession
# 2 m -> full candidate notch recession.
#
# Selected IDs 1..5 merely restrict the
# moving plane to the predefined toe-notch zone.
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

    sub = ET.SubElement(
        ds,
        "deactivated_subdomain",
    )

    ET.SubElement(
        sub,
        "time_curve",
    ).text = "erosion_front_curve"

    line = ET.SubElement(
        sub,
        "line_segment",
    )

    ET.SubElement(
        line,
        "start",
    ).text = "22 2.6 0"

    ET.SubElement(
        line,
        "end",
    ).text = "20 2.6 0"

    ET.SubElement(
        sub,
        "material_ids",
    ).text = "1 2 3 4 5"


# ============================================================
# ROOT-LEVEL CURVE
#
# t = 0..10 s:
# intact hold.
#
# t = 10..110 s:
# E ramps continuously 0 -> 2 m.
#
# t = 110..130 s:
# final hold.
#
# THIS TIME IS NUMERICAL CONTINUATION TIME,
# NOT PHYSICAL COASTAL EROSION DURATION.
# ============================================================

curves = None

for child in list(
    root
):

    if child.tag == "curves":

        curves = child
        break


if curves is None:

    curves = ET.Element(
        "curves"
    )

    children = list(
        root
    )

    pv_index = None

    for i, child in enumerate(
        children
    ):

        if child.tag == "process_variables":

            pv_index = i
            break

    if pv_index is None:

        raise SystemExit(
            "FAIL: root process_variables "
            "not found"
        )

    root.insert(
        pv_index,
        curves,
    )


for curve in list(
    curves.findall("curve")
):

    if (
        curve.findtext("name")
        == "erosion_front_curve"
    ):

        curves.remove(
            curve
        )


curve = ET.SubElement(
    curves,
    "curve",
)

ET.SubElement(
    curve,
    "name",
).text = "erosion_front_curve"

ET.SubElement(
    curve,
    "coords",
).text = (
    "0 10 110 130"
)

ET.SubElement(
    curve,
    "values",
).text = (
    "0 0 2 2"
)


# ============================================================
# TIME STEPPING
#
# dt = 0.5 s
#
# During erosion:
#
# dE = 2 / 100 * 0.5
#    = 0.01 m numerical front movement
#      per timestep.
#
# The mesh itself remains discrete, so cells
# deactivate when their centres are crossed.
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
).text = "130"

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
).text = "260"

ET.SubElement(
    pair,
    "delta_t",
).text = "0.5"


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
# OUTPUT
#
# Every 10 steps = every 5 s.
# During erosion this corresponds to
# nominal dE = 0.1 m between saved states.
# ============================================================

prefix = root.find(
    "./time_loop/output/prefix"
)

if prefix is None:

    raise SystemExit(
        "FAIL: output prefix missing"
    )

prefix.text = (
    "phase05b_v2_reference"
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
).text = "260"

ET.SubElement(
    pair,
    "each_steps",
).text = "10"


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
# SAVE / VERIFY
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


check = ET.parse(
    DST
).getroot()


for pv in check.findall(
    "./process_variables/process_variable"
):

    name = pv.findtext(
        "name"
    )

    if name not in {
        "pressure",
        "displacement",
    }:

        continue

    curve_name = pv.findtext(
        "./deactivated_subdomains/"
        "deactivated_subdomain/"
        "time_curve"
    )

    ids = pv.findtext(
        "./deactivated_subdomains/"
        "deactivated_subdomain/"
        "material_ids"
    )

    if curve_name != "erosion_front_curve":

        raise SystemExit(
            f"FAIL: {name} moving front absent"
        )

    if ids != "1 2 3 4 5":

        raise SystemExit(
            f"FAIL: {name} material IDs wrong"
        )


print(
    "PHASE 05B-V2 MOVING FRONT BUILD: PASS"
)

print()

print(
    "Erosion front:"
)

print(
    "x = 22 m -> x = 20 m"
)

print(
    "E = 0 -> 2.0 m"
)

print(
    "Numerical ramp:"
)

print(
    "t = 10 -> 110 s"
)

print(
    "Numerical front increment per step:"
)

print(
    "dE = 0.01 m"
)

print()

print(
    "Selected toe-notch MaterialIDs:"
)

print(
    "1 2 3 4 5"
)

print()

print(
    "No hydraulic boundary_parameter "
    "on newly exposed notch."
)
