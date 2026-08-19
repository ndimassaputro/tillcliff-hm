from pathlib import Path
import shutil
import xml.etree.ElementTree as ET

import numpy as np
import pyvista as pv
import ogstools as ot


ROOT = Path.cwd()

BASE_PRJ = (
    ROOT
    / "model"
    / "phase04a_slope"
    / "gravity_baseline.prj"
)

BASE_MESH = (
    ROOT
    / "model"
    / "phase04a_slope"
    / "slope.vtu"
)

BASE_RESULTS = (
    ROOT
    / "results"
    / "phase04a_gravity"
)

MODEL = (
    ROOT
    / "model"
    / "phase04d_v2"
)

MODEL.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# CHECK INPUTS
# ============================================================

for path in [
    BASE_PRJ,
    BASE_MESH,
]:
    if not path.exists():
        raise SystemExit(
            f"FAIL: missing {path}"
        )

pvds = sorted(
    BASE_RESULTS.glob("*.pvd")
)

if not pvds:
    raise SystemExit(
        "FAIL: no Phase 04A PVD"
    )


# ============================================================
# LOAD FINAL VERIFIED PHASE 04A STATE
# ============================================================

series = ot.MeshSeries(
    str(pvds[0])
)

times = np.asarray(
    series.timevalues,
    dtype=float,
)

final = series.mesh(
    len(times) - 1
)

print(
    "Phase 04A final time [days]:",
    times[-1] / 86400.0,
)


# ============================================================
# CREATE INITIALIZED BULK MESH
# ============================================================

bulk = pv.read(
    BASE_MESH
)

if (
    bulk.n_points
    != final.n_points
):
    raise SystemExit(
        "FAIL: point-count mismatch"
    )

if (
    bulk.n_cells
    != final.n_cells
):
    raise SystemExit(
        "FAIL: cell-count mismatch"
    )


# ------------------------------------------------------------
# PRESSURE -> nodal restart field p0
# ------------------------------------------------------------

if "pressure" in final.point_data:

    p0 = np.asarray(
        final.point_data["pressure"],
        dtype=float,
    ).squeeze()

elif "pressure" in final.cell_data:

    temp = final.cell_data_to_point_data()

    p0 = np.asarray(
        temp.point_data["pressure"],
        dtype=float,
    ).squeeze()

else:

    raise SystemExit(
        "FAIL: pressure missing "
        "from Phase 04A final state"
    )

if len(p0) != bulk.n_points:
    raise SystemExit(
        "FAIL: p0 size mismatch"
    )

bulk.point_data[
    "p0"
] = p0


# ------------------------------------------------------------
# STRESS -> element restart field sigma0
# ------------------------------------------------------------

if "sigma" in final.cell_data:

    sigma0 = np.asarray(
        final.cell_data["sigma"],
        dtype=float,
    )

elif "sigma" in final.point_data:

    temp = final.point_data_to_cell_data()

    sigma0 = np.asarray(
        temp.cell_data["sigma"],
        dtype=float,
    )

else:

    raise SystemExit(
        "FAIL: sigma missing "
        "from Phase 04A final state"
    )

if sigma0.shape[0] != bulk.n_cells:
    raise SystemExit(
        "FAIL: sigma0 cell-count mismatch"
    )

bulk.cell_data[
    "sigma0"
] = sigma0


# ------------------------------------------------------------
# Diagnostics
# ------------------------------------------------------------

if not np.all(
    np.isfinite(p0)
):
    raise SystemExit(
        "FAIL: non-finite p0"
    )

if not np.all(
    np.isfinite(sigma0)
):
    raise SystemExit(
        "FAIL: non-finite sigma0"
    )

INIT_MESH = (
    MODEL
    / "slope_initialized.vtu"
)

bulk.save(
    INIT_MESH
)

print()
print("INITIAL STATE TRANSFER")
print(
    "pressure range [kPa]:",
    np.min(p0) / 1000.0,
    "to",
    np.max(p0) / 1000.0,
)

print(
    "sigma0 shape:",
    sigma0.shape,
)

print(
    "sigma0 range [kPa]:",
    np.min(sigma0) / 1000.0,
    "to",
    np.max(sigma0) / 1000.0,
)

print(
    "PASS:",
    INIT_MESH,
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
        / "phase04a_slope"
        / name
    )

    dst = MODEL / name

    shutil.copy2(
        src,
        dst,
    )


# ============================================================
# BUILD TWO MOHR-COULOMB CASES
# ============================================================

cases = {
    "strong": {
        "c": 80000,
        "phi": 35,
    },
    "screening": {
        "c": 25000,
        "phi": 28,
    },
}


def find_parameter(
    parameters,
    name,
):

    for p in parameters.findall(
        "parameter"
    ):

        if (
            p.findtext("name")
            == name
        ):
            return p

    return None


def add_constant(
    parameters,
    name,
    value,
):

    old = find_parameter(
        parameters,
        name,
    )

    if old is not None:
        parameters.remove(
            old
        )

    p = ET.SubElement(
        parameters,
        "parameter",
    )

    ET.SubElement(
        p,
        "name",
    ).text = name

    ET.SubElement(
        p,
        "type",
    ).text = "Constant"

    ET.SubElement(
        p,
        "value",
    ).text = str(value)


for case, strength in cases.items():

    tree = ET.parse(
        BASE_PRJ
    )

    root = tree.getroot()

    # --------------------------------------------------------
    # BULK MESH
    # --------------------------------------------------------

    meshes = root.find(
        "meshes"
    )

    if meshes is None:
        raise SystemExit(
            "FAIL: meshes missing"
        )

    for mesh_tag in meshes.findall(
        "mesh"
    ):

        if (
            mesh_tag.text.strip()
            == "slope.vtu"
        ):
            mesh_tag.text = (
                "slope_initialized.vtu"
            )


    # --------------------------------------------------------
    # PROCESS
    # --------------------------------------------------------

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


    # --------------------------------------------------------
    # CONSTITUTIVE MODEL
    # --------------------------------------------------------

    old_cr = process.find(
        "constitutive_relation"
    )

    if old_cr is None:
        raise SystemExit(
            "FAIL: constitutive relation missing"
        )

    children = list(
        process
    )

    index = children.index(
        old_cr
    )

    process.remove(
        old_cr
    )

    cr = ET.Element(
        "constitutive_relation"
    )

    ET.SubElement(
        cr,
        "type",
    ).text = "MFront"

    ET.SubElement(
        cr,
        "behaviour",
    ).text = (
        "MohrCoulombAbboSloan"
    )

    mp = ET.SubElement(
        cr,
        "material_properties",
    )

    mapping = [
        (
            "YoungModulus",
            "E",
        ),
        (
            "PoissonRatio",
            "nu",
        ),
        (
            "Cohesion",
            "Cohesion",
        ),
        (
            "FrictionAngle",
            "FrictionAngle",
        ),
        (
            "DilatancyAngle",
            "DilatancyAngle",
        ),
        (
            "TransitionAngle",
            "TransitionAngle",
        ),
        (
            "TensionCutOffParameter",
            "TensionCutOffParameter",
        ),
    ]

    for prop_name, param_name in mapping:

        ET.SubElement(
            mp,
            "material_property",
            {
                "name": prop_name,
                "parameter": param_name,
            },
        )

    process.insert(
        index,
        cr,
    )


    # --------------------------------------------------------
    # INITIAL STRESS
    # --------------------------------------------------------

    old_initial_stress = process.find(
        "initial_stress"
    )

    if old_initial_stress is not None:
        process.remove(
            old_initial_stress
        )

    specific_body = process.find(
        "specific_body_force"
    )

    if specific_body is None:
        raise SystemExit(
            "FAIL: specific body force missing"
        )

    body_index = list(
        process
    ).index(
        specific_body
    )

    initial_stress = ET.Element(
        "initial_stress"
    )

    initial_stress.text = "sigma0"

    process.insert(
        body_index + 1,
        initial_stress,
    )


    # --------------------------------------------------------
    # SECONDARY VARIABLES
    # --------------------------------------------------------

    secondary = process.find(
        "secondary_variables"
    )

    if secondary is None:
        secondary = ET.SubElement(
            process,
            "secondary_variables",
        )

    existing = {
        s.attrib.get("name")
        for s in secondary.findall(
            "secondary_variable"
        )
    }

    for name in [
        "ElasticStrain",
        "EquivalentPlasticStrain",
    ]:

        if name not in existing:

            ET.SubElement(
                secondary,
                "secondary_variable",
                {
                    "name": name,
                },
            )


    # --------------------------------------------------------
    # KEEP SOLID DENSITY PHYSICAL AND CONSTANT
    # --------------------------------------------------------

    solid_phase = None

    for phase in root.findall(
        "./media/medium/phases/phase"
    ):

        if (
            phase.findtext("type")
            == "Solid"
        ):

            solid_phase = phase
            break

    if solid_phase is None:
        raise SystemExit(
            "FAIL: solid phase missing"
        )

    for prop in solid_phase.findall(
        "./properties/property"
    ):

        if (
            prop.findtext("name")
            == "density"
        ):

            prop.clear()

            ET.SubElement(
                prop,
                "name",
            ).text = "density"

            ET.SubElement(
                prop,
                "type",
            ).text = "Constant"

            ET.SubElement(
                prop,
                "value",
            ).text = "2650"


    # --------------------------------------------------------
    # PARAMETERS
    # --------------------------------------------------------

    parameters = root.find(
        "parameters"
    )

    if parameters is None:
        raise SystemExit(
            "FAIL: parameters missing"
        )

    add_constant(
        parameters,
        "Cohesion",
        strength["c"],
    )

    add_constant(
        parameters,
        "FrictionAngle",
        strength["phi"],
    )

    add_constant(
        parameters,
        "DilatancyAngle",
        0,
    )

    add_constant(
        parameters,
        "TransitionAngle",
        27,
    )

    add_constant(
        parameters,
        "TensionCutOffParameter",
        10000,
    )


    # sigma0 as MeshElement
    old = find_parameter(
        parameters,
        "sigma0",
    )

    if old is not None:
        parameters.remove(
            old
        )

    p = ET.SubElement(
        parameters,
        "parameter",
    )

    ET.SubElement(
        p,
        "name",
    ).text = "sigma0"

    ET.SubElement(
        p,
        "type",
    ).text = "MeshElement"

    ET.SubElement(
        p,
        "field_name",
    ).text = "sigma0"


    # pressure_ic as MeshNode
    old = find_parameter(
        parameters,
        "pressure_ic_restart",
    )

    if old is not None:
        parameters.remove(
            old
        )

    p = ET.SubElement(
        parameters,
        "parameter",
    )

    ET.SubElement(
        p,
        "name",
    ).text = (
        "pressure_ic_restart"
    )

    ET.SubElement(
        p,
        "type",
    ).text = "MeshNode"

    ET.SubElement(
        p,
        "field_name",
    ).text = "p0"


    # --------------------------------------------------------
    # PROCESS VARIABLE INITIAL PRESSURE
    # --------------------------------------------------------

    root_pvs = root.find(
        "process_variables"
    )

    if root_pvs is None:
        raise SystemExit(
            "FAIL: root process variables missing"
        )

    pressure_pv = None

    for pv in root_pvs.findall(
        "process_variable"
    ):

        if (
            pv.findtext("name")
            == "pressure"
        ):

            pressure_pv = pv
            break

    if pressure_pv is None:
        raise SystemExit(
            "FAIL: pressure PV missing"
        )

    ic = pressure_pv.find(
        "initial_condition"
    )

    if ic is None:
        raise SystemExit(
            "FAIL: pressure IC missing"
        )

    ic.text = (
        "pressure_ic_restart"
    )


    # --------------------------------------------------------
    # ADAPTIVE TIME STEPPING
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
    ).text = (
        "IterationNumberBasedTimeStepping"
    )

    ET.SubElement(
        ts,
        "t_initial",
    ).text = "0"

    # Two-day hold is enough for this transfer test.
    ET.SubElement(
        ts,
        "t_end",
    ).text = "172800"

    ET.SubElement(
        ts,
        "initial_dt",
    ).text = "60"

    ET.SubElement(
        ts,
        "minimum_dt",
    ).text = "0.01"

    ET.SubElement(
        ts,
        "maximum_dt",
    ).text = "21600"

    ET.SubElement(
        ts,
        "number_iterations",
    ).text = (
        "1 4 8 12 20 40 70"
    )

    ET.SubElement(
        ts,
        "multiplier",
    ).text = (
        "1.5 1.2 1.0 0.8 0.5 0.25 0.1"
    )


    # --------------------------------------------------------
    # NEWTON
    # --------------------------------------------------------

    nonlinear = root.find(
        "./nonlinear_solvers/"
        "nonlinear_solver"
    )

    if nonlinear is None:
        raise SystemExit(
            "FAIL: nonlinear solver missing"
        )

    max_iter = nonlinear.find(
        "max_iter"
    )

    if max_iter is None:

        max_iter = ET.SubElement(
            nonlinear,
            "max_iter",
        )

    max_iter.text = "80"


    # --------------------------------------------------------
    # OUTPUT
    # --------------------------------------------------------

    prefix = root.find(
        "./time_loop/output/prefix"
    )

    prefix.text = (
        f"phase04d_v2_{case}"
    )

    variables = root.find(
        "./time_loop/output/variables"
    )

    existing_variables = {
        x.text
        for x in variables.findall(
            "variable"
        )
    }

    for name in [
        "ElasticStrain",
        "EquivalentPlasticStrain",
    ]:

        if name not in existing_variables:

            ET.SubElement(
                variables,
                "variable",
            ).text = name


    # --------------------------------------------------------
    # SAVE
    # --------------------------------------------------------

    ET.indent(
        tree,
        space="    ",
    )

    out = (
        MODEL
        / f"{case}.prj"
    )

    tree.write(
        out,
        encoding="UTF-8",
        xml_declaration=True,
    )

    ET.parse(
        out
    )

    print()
    print(
        f"PASS BUILD {case.upper()}"
    )

    print(
        f"c = {strength['c']/1000:.1f} kPa"
    )

    print(
        f"phi = {strength['phi']} deg"
    )


print()
print(
    "PHASE 04D-V2 INITIALIZED BUILD: PASS"
)
