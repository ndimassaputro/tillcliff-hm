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

SRC_MESH_DIR = (
    ROOT
    / "model"
    / "phase05a_v2_toe_notch"
)

MODEL = (
    ROOT
    / "model"
    / "phase05b0_single_erosion"
)

MODEL.mkdir(
    parents=True,
    exist_ok=True,
)

DST = (
    MODEL
    / "reference_E0p4.prj"
)


# ============================================================
# INPUT CHECK
# ============================================================

if not SRC_PRJ.exists():
    raise SystemExit(
        f"FAIL: missing {SRC_PRJ}"
    )

bulk = (
    SRC_MESH_DIR
    / "slope_toe_notch_ready.vtu"
)

if not bulk.exists():
    raise SystemExit(
        f"FAIL: missing {bulk}"
    )


# ============================================================
# COPY MESHES
# ============================================================

shutil.copy2(
    bulk,
    MODEL
    / "slope_toe_notch_ready.vtu",
)

for name in [
    "slope_left.vtu",
    "slope_right.vtu",
    "slope_top.vtu",
    "slope_bottom.vtu",
]:

    src = (
        SRC_MESH_DIR
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
# LOAD VERIFIED ELASTOPLASTIC HM PROJECT
# ============================================================

tree = ET.parse(
    SRC_PRJ
)

root = tree.getroot()


# ============================================================
# USE EROSION-READY BULK MESH
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
        "was not replaced"
    )


# ============================================================
# PROCESS / CONSTITUTIVE RELATION
#
# MaterialIDs now range 0..5.
# All IDs use the same till constitutive model.
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


# ============================================================
# MEDIUM MAPPING
#
# Same hydraulic/material properties for all bands.
# Material IDs here are geometric labels for erosion,
# not different soil types.
# ============================================================

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
# ADD DEACTIVATION TO BOTH PRIMARY VARIABLES
#
# MaterialID 1:
# nominal E = 0.4 m.
#
# t < 10 s  : intact
# t >= 10 s : band 1 removed
#
# Numerical time here is an event-continuation
# coordinate, NOT physical erosion duration.
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

    deactivated = ET.SubElement(
        pv,
        "deactivated_subdomains",
    )

    subdomain = ET.SubElement(
        deactivated,
        "deactivated_subdomain",
    )

    interval = ET.SubElement(
        subdomain,
        "time_interval",
    )

    ET.SubElement(
        interval,
        "start",
    ).text = "10"

    ET.SubElement(
        interval,
        "end",
    ).text = "20"

    ET.SubElement(
        subdomain,
        "material_ids",
    ).text = "1"


# ============================================================
# TIME STEPPING
#
# 1-second steps guarantee outputs immediately
# before and after the erosion event.
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
).text = "20"

timesteps = ET.SubElement(
    ts,
    "timesteps",
)

pair = ET.SubElement(
    timesteps,
    "pair",
)

ET.SubElement(
    pair,
    "repeat",
).text = "20"

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

max_iter.text = "80"


# ============================================================
# OUTPUT EVERY STEP
# ============================================================

prefix = root.find(
    "./time_loop/output/prefix"
)

if prefix is None:
    raise SystemExit(
        "FAIL: output prefix missing"
    )

prefix.text = (
    "phase05b0_reference_E0p4"
)


output_timesteps = root.find(
    "./time_loop/output/timesteps"
)

if output_timesteps is not None:

    output_timesteps.clear()

    pair = ET.SubElement(
        output_timesteps,
        "pair",
    )

    ET.SubElement(
        pair,
        "repeat",
    ).text = "20"

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
# SAVE / STRUCTURAL VALIDATION
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

check = ET.parse(
    DST
)

check_root = check.getroot()


check_cr = check_root.find(
    "./processes/process/"
    "constitutive_relation"
)

if (
    check_cr is None
    or check_cr.attrib.get("id")
    != "0,1,2,3,4,5"
):

    raise SystemExit(
        "FAIL: constitutive MaterialID "
        "mapping incorrect"
    )


check_medium = check_root.find(
    "./media/medium"
)

if (
    check_medium is None
    or check_medium.attrib.get("id")
    != "0,1,2,3,4,5"
):

    raise SystemExit(
        "FAIL: medium MaterialID "
        "mapping incorrect"
    )


for pv in check_root.findall(
    "./process_variables/"
    "process_variable"
):

    name = pv.findtext(
        "name"
    )

    if name not in {
        "pressure",
        "displacement",
    }:
        continue

    ids = pv.findtext(
        "./deactivated_subdomains/"
        "deactivated_subdomain/"
        "material_ids"
    )

    if ids != "1":

        raise SystemExit(
            f"FAIL: {name} deactivation "
            "not installed"
        )


print(
    "PHASE 05B-0 BUILD: PASS"
)

print()
print(
    "Reference antecedent hydraulic state"
)

print(
    "Screening Mohr-Coulomb: "
    "c=25 kPa, phi=28 deg"
)

print(
    "MaterialID removed: 1"
)

print(
    "Nominal erosion E: 0.4 m"
)

print(
    "Removal time: t=10 s"
)

print()
print(
    "No hydraulic boundary_parameter "
    "added on new notch surface."
)

print(
    "This isolates the first-order "
    "geometric/mechanical erosion response."
)
