from pathlib import Path
import shutil
import xml.etree.ElementTree as ET


SRC = Path(
    "model/phase04a_slope/gravity_baseline.prj"
)

MODEL = Path(
    "model/phase04d_mc_hm"
)

DST = MODEL / "mc_hm_baseline.prj"

MODEL.mkdir(
    parents=True,
    exist_ok=True,
)

if not SRC.exists():
    raise SystemExit(
        f"FAIL: missing {SRC}"
    )


# ============================================================
# COPY CLEAN PHASE 04A BASELINE
# ============================================================

shutil.copy2(
    SRC,
    DST,
)

tree = ET.parse(DST)
root = tree.getroot()


# ============================================================
# FIND RICHARDS_MECHANICS PROCESS
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
        "FAIL: RICHARDS_MECHANICS process not found"
    )


# ============================================================
# REPLACE LINEAR ELASTICITY WITH MFRONT MOHR-COULOMB
# ============================================================

old_cr = process.find(
    "constitutive_relation"
)

if old_cr is None:
    raise SystemExit(
        "FAIL: constitutive_relation missing"
    )

children = list(process)
cr_index = children.index(old_cr)

process.remove(old_cr)

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
).text = "MohrCoulombAbboSloan"

material_properties = ET.SubElement(
    cr,
    "material_properties",
)

properties = [
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

for name, parameter in properties:

    ET.SubElement(
        material_properties,
        "material_property",
        {
            "name": name,
            "parameter": parameter,
        },
    )

process.insert(
    cr_index,
    cr,
)


# ============================================================
# FORCE PROCESS-VARIABLE ORDER TO MATCH OFFICIAL RM EXAMPLES
# ============================================================

process_pv = process.find(
    "process_variables"
)

if process_pv is None:
    raise SystemExit(
        "FAIL: process process_variables missing"
    )

process_pv.clear()

ET.SubElement(
    process_pv,
    "pressure",
).text = "pressure"

ET.SubElement(
    process_pv,
    "displacement",
).text = "displacement"


# ============================================================
# ADD PLASTIC SECONDARY VARIABLES
# ============================================================

secondary = process.find(
    "secondary_variables"
)

if secondary is None:
    secondary = ET.SubElement(
        process,
        "secondary_variables",
    )

existing_secondary = {
    x.attrib.get("name")
    for x in secondary.findall(
        "secondary_variable"
    )
}

for name in [
    "ElasticStrain",
    "EquivalentPlasticStrain",
]:

    if name not in existing_secondary:

        ET.SubElement(
            secondary,
            "secondary_variable",
            {
                "name": name,
            },
        )


# ============================================================
# PARAMETERS
# ============================================================

parameters = root.find(
    "parameters"
)

if parameters is None:
    raise SystemExit(
        "FAIL: parameters block missing"
    )


def existing_parameter_names():

    return {
        p.findtext("name")
        for p in parameters.findall(
            "parameter"
        )
    }


def add_constant(
    name,
    value,
):

    if name in existing_parameter_names():
        return

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


# ------------------------------------------------------------
# IMPORTANT:
# Screening parameters only.
# NOT site-calibrated Danish till parameters.
# ------------------------------------------------------------

add_constant(
    "Cohesion",
    "25000",
)

add_constant(
    "FrictionAngle",
    "28",
)

add_constant(
    "DilatancyAngle",
    "0",
)

add_constant(
    "TransitionAngle",
    "27",
)

add_constant(
    "TensionCutOffParameter",
    "10000",
)

add_constant(
    "rho_s_value",
    "2650",
)


# ============================================================
# NUMERICAL SELF-WEIGHT RAMP
# ============================================================

# Replace solid density:
#
# Constant 2650
#
# with
#
# Parameter -> rho_s_loading
#
# rho_s_loading itself is curve-scaled from
# 0 to full density over first 2 days.

solid_phase = None

for phase in root.findall(
    "./media/medium/phases/phase"
):

    if phase.findtext("type") == "Solid":
        solid_phase = phase
        break

if solid_phase is None:
    raise SystemExit(
        "FAIL: Solid phase not found"
    )

solid_properties = solid_phase.find(
    "properties"
)

if solid_properties is None:
    raise SystemExit(
        "FAIL: Solid properties missing"
    )

density_property = None

for prop in solid_properties.findall(
    "property"
):

    if prop.findtext("name") == "density":
        density_property = prop
        break

if density_property is None:
    raise SystemExit(
        "FAIL: solid density property missing"
    )

density_property.clear()

ET.SubElement(
    density_property,
    "name",
).text = "density"

ET.SubElement(
    density_property,
    "type",
).text = "Parameter"

ET.SubElement(
    density_property,
    "parameter_name",
).text = "rho_s_loading"


# Add CurveScaled parameter.

if "rho_s_loading" not in existing_parameter_names():

    p = ET.SubElement(
        parameters,
        "parameter",
    )

    ET.SubElement(
        p,
        "name",
    ).text = "rho_s_loading"

    ET.SubElement(
        p,
        "type",
    ).text = "CurveScaled"

    ET.SubElement(
        p,
        "curve",
    ).text = "gravity_density_ramp"

    ET.SubElement(
        p,
        "parameter",
    ).text = "rho_s_value"


# ============================================================
# ROOT-LEVEL CURVE BLOCK
# ============================================================

curves = root.find(
    "curves"
)

if curves is None:

    curves = ET.Element(
        "curves"
    )

    root_children = list(root)

    pv_index = None

    for i, child in enumerate(
        root_children
    ):
        if child.tag == "process_variables":
            pv_index = i
            break

    if pv_index is None:
        raise SystemExit(
            "FAIL: root process_variables missing"
        )

    root.insert(
        pv_index,
        curves,
    )


# Remove duplicate curve if rebuilding.
for curve in list(
    curves.findall("curve")
):

    if (
        curve.findtext("name")
        == "gravity_density_ramp"
    ):
        curves.remove(curve)


curve = ET.SubElement(
    curves,
    "curve",
)

ET.SubElement(
    curve,
    "name",
).text = "gravity_density_ramp"

# seconds:
# 0 d -> no solid self-weight
# 2 d -> full self-weight
# 30 d -> hold full self-weight
ET.SubElement(
    curve,
    "coords",
).text = (
    "0 "
    "172800 "
    "2592000"
)

ET.SubElement(
    curve,
    "values",
).text = (
    "0 "
    "1 "
    "1"
)


# ============================================================
# TIME STEPPING
# ============================================================

time_stepping = root.find(
    "./time_loop/processes/process/"
    "time_stepping"
)

if time_stepping is None:
    raise SystemExit(
        "FAIL: time_stepping missing"
    )

time_stepping.clear()

ET.SubElement(
    time_stepping,
    "type",
).text = "FixedTimeStepping"

ET.SubElement(
    time_stepping,
    "t_initial",
).text = "0"

ET.SubElement(
    time_stepping,
    "t_end",
).text = "2592000"

timesteps = ET.SubElement(
    time_stepping,
    "timesteps",
)

# 2-day gravity loading:
# 8 increments of 6 h.
pair1 = ET.SubElement(
    timesteps,
    "pair",
)

ET.SubElement(
    pair1,
    "repeat",
).text = "8"

ET.SubElement(
    pair1,
    "delta_t",
).text = "21600"

# Then 28 days hold/equilibration.
pair2 = ET.SubElement(
    timesteps,
    "pair",
)

ET.SubElement(
    pair2,
    "repeat",
).text = "28"

ET.SubElement(
    pair2,
    "delta_t",
).text = "86400"


# ============================================================
# NEWTON ROBUSTNESS
# ============================================================

solver = root.find(
    "./nonlinear_solvers/nonlinear_solver"
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
# OUTPUT
# ============================================================

prefix = root.find(
    "./time_loop/output/prefix"
)

if prefix is None:
    raise SystemExit(
        "FAIL: output prefix missing"
    )

prefix.text = "phase04d_mc_hm"

output_variables = root.find(
    "./time_loop/output/variables"
)

if output_variables is None:
    raise SystemExit(
        "FAIL: output variables missing"
    )

existing_output = {
    x.text
    for x in output_variables.findall(
        "variable"
    )
}

for name in [
    "ElasticStrain",
    "EquivalentPlasticStrain",
]:

    if name not in existing_output:

        ET.SubElement(
            output_variables,
            "variable",
        ).text = name


# ============================================================
# FORMAT + SAVE
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

# Reparse for XML sanity.
ET.parse(DST)


# ============================================================
# STRUCTURAL CHECKS
# ============================================================

check = ET.parse(DST)
check_root = check.getroot()

behaviour = check_root.findtext(
    "./processes/process/"
    "constitutive_relation/behaviour"
)

curve_names = [
    c.findtext("name")
    for c in check_root.findall(
        "./curves/curve"
    )
]

secondary_names = [
    s.attrib.get("name")
    for s in check_root.findall(
        "./processes/process/"
        "secondary_variables/"
        "secondary_variable"
    )
]

if behaviour != "MohrCoulombAbboSloan":
    raise SystemExit(
        "FAIL: MFront behaviour not installed"
    )

if (
    "gravity_density_ramp"
    not in curve_names
):
    raise SystemExit(
        "FAIL: gravity curve not ROOT-level"
    )

if (
    "EquivalentPlasticStrain"
    not in secondary_names
):
    raise SystemExit(
        "FAIL: plastic-strain output absent"
    )


print(
    "PHASE 04D PROJECT BUILD: PASS"
)

print()
print(
    "Constitutive model: "
    "MFront / MohrCoulombAbboSloan"
)

print(
    "Screening cohesion: 25 kPa"
)

print(
    "Screening friction angle: 28 deg"
)

print(
    "Screening dilatancy angle: 0 deg"
)

print(
    "Transition angle: 27 deg"
)

print(
    "Tension cut-off parameter: 10 kPa"
)

print()
print(
    "Self-weight loading: "
    "0 -> full over first 2 days"
)

print(
    "Full-gravity hold: "
    "day 2 -> day 30"
)

print()
print(
    "IMPORTANT: strength values are provisional "
    "screening parameters, not Danish site calibration."
)
