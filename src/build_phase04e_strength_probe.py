from pathlib import Path
import csv
import math
import shutil
import xml.etree.ElementTree as ET


ROOT = Path.cwd()

SRC_MODEL = (
    ROOT
    / "model"
    / "phase04d_v2"
)

SRC_PRJ = (
    SRC_MODEL
    / "screening.prj"
)

MODEL = (
    ROOT
    / "model"
    / "phase04e_strength_probe"
)

MODEL.mkdir(
    parents=True,
    exist_ok=True,
)

if not SRC_PRJ.exists():
    raise SystemExit(
        f"FAIL: missing {SRC_PRJ}"
    )


# ============================================================
# BASELINE SCREENING STRENGTH
# ============================================================

C0 = 25000.0
PHI0 = 28.0

SRFS = [
    1.00,
    1.25,
    1.50,
    1.75,
    2.00,
    2.50,
    3.00,
    4.00,
    5.00,
]


def case_name(srf):
    return (
        f"srf_{srf:.2f}"
        .replace(".", "p")
    )


# ============================================================
# COPY INITIALIZED MESH + BOUNDARIES
# ============================================================

for name in [
    "slope_initialized.vtu",
    "slope_left.vtu",
    "slope_right.vtu",
    "slope_top.vtu",
    "slope_bottom.vtu",
]:

    src = SRC_MODEL / name
    dst = MODEL / name

    if not src.exists():
        raise SystemExit(
            f"FAIL: missing {src}"
        )

    shutil.copy2(
        src,
        dst,
    )


# ============================================================
# XML HELPERS
# ============================================================

def find_parameter(
    parameters,
    name,
):

    for p in parameters.findall(
        "parameter"
    ):

        if p.findtext("name") == name:
            return p

    return None


def remove_parameter(
    parameters,
    name,
):

    old = find_parameter(
        parameters,
        name,
    )

    if old is not None:
        parameters.remove(
            old
        )


def add_constant(
    parameters,
    name,
    value,
):

    remove_parameter(
        parameters,
        name,
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


def add_curve_scaled(
    parameters,
    name,
    curve_name,
    base_parameter,
):

    remove_parameter(
        parameters,
        name,
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
    ).text = "CurveScaled"

    ET.SubElement(
        p,
        "curve",
    ).text = curve_name

    ET.SubElement(
        p,
        "parameter",
    ).text = base_parameter


def get_root_curves(
    root,
):

    curves = None

    for child in list(root):

        if child.tag == "curves":
            curves = child
            break

    if curves is None:

        curves = ET.Element(
            "curves"
        )

        children = list(root)

        pv_index = None

        for i, child in enumerate(
            children
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

    return curves


# ============================================================
# BUILD CASES
# ============================================================

case_rows = []

for srf in SRFS:

    c_target = (
        C0 / srf
    )

    phi_target = math.degrees(
        math.atan(
            math.tan(
                math.radians(PHI0)
            )
            / srf
        )
    )

    c_ratio = (
        c_target / C0
    )

    phi_ratio = (
        phi_target / PHI0
    )

    case = case_name(
        srf
    )

    tree = ET.parse(
        SRC_PRJ
    )

    root = tree.getroot()

    parameters = root.find(
        "parameters"
    )

    if parameters is None:
        raise SystemExit(
            "FAIL: parameters block missing"
        )

    # --------------------------------------------------------
    # Strength parameters become time-ramped.
    # --------------------------------------------------------

    add_constant(
        parameters,
        "CohesionBase",
        C0,
    )

    add_constant(
        parameters,
        "FrictionAngleBase",
        PHI0,
    )

    add_curve_scaled(
        parameters,
        "Cohesion",
        "cohesion_strength_ramp",
        "CohesionBase",
    )

    add_curve_scaled(
        parameters,
        "FrictionAngle",
        "friction_strength_ramp",
        "FrictionAngleBase",
    )

    # Keep these fixed during this diagnostic.
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

    # --------------------------------------------------------
    # Root-level strength curves.
    #
    # 0-10 s:
    #     hold full baseline strength
    #
    # 10-1010 s:
    #     numerical continuation to target strength
    #
    # 1010-1210 s:
    #     hold target strength
    #
    # Time here is NUMERICAL continuation time,
    # not physical geotechnical time.
    # --------------------------------------------------------

    curves = get_root_curves(
        root
    )

    for curve in list(
        curves.findall("curve")
    ):

        name = curve.findtext(
            "name"
        )

        if name in {
            "cohesion_strength_ramp",
            "friction_strength_ramp",
        }:

            curves.remove(
                curve
            )

    cohesion_curve = ET.SubElement(
        curves,
        "curve",
    )

    ET.SubElement(
        cohesion_curve,
        "name",
    ).text = "cohesion_strength_ramp"

    ET.SubElement(
        cohesion_curve,
        "coords",
    ).text = (
        "0 10 1010 1210"
    )

    ET.SubElement(
        cohesion_curve,
        "values",
    ).text = (
        f"1 1 "
        f"{c_ratio:.12f} "
        f"{c_ratio:.12f}"
    )

    friction_curve = ET.SubElement(
        curves,
        "curve",
    )

    ET.SubElement(
        friction_curve,
        "name",
    ).text = "friction_strength_ramp"

    ET.SubElement(
        friction_curve,
        "coords",
    ).text = (
        "0 10 1010 1210"
    )

    ET.SubElement(
        friction_curve,
        "values",
    ).text = (
        f"1 1 "
        f"{phi_ratio:.12f} "
        f"{phi_ratio:.12f}"
    )

    # --------------------------------------------------------
    # Adaptive numerical stepping.
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

    ET.SubElement(
        ts,
        "t_end",
    ).text = "1210"

    ET.SubElement(
        ts,
        "initial_dt",
    ).text = "1"

    ET.SubElement(
        ts,
        "minimum_dt",
    ).text = "0.001"

    ET.SubElement(
        ts,
        "maximum_dt",
    ).text = "50"

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
    # Newton.
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
    # Output.
    # --------------------------------------------------------

    prefix = root.find(
        "./time_loop/output/prefix"
    )

    if prefix is None:
        raise SystemExit(
            "FAIL: output prefix missing"
        )

    prefix.text = (
        f"phase04e_{case}"
    )

    output_variables = root.find(
        "./time_loop/output/variables"
    )

    if output_variables is None:
        raise SystemExit(
            "FAIL: output variables missing"
        )

    existing = {
        x.text
        for x in output_variables.findall(
            "variable"
        )
    }

    for variable in [
        "EquivalentPlasticStrain",
        "ElasticStrain",
        "pressure",
        "saturation",
        "displacement",
        "sigma",
    ]:

        if variable not in existing:

            ET.SubElement(
                output_variables,
                "variable",
            ).text = variable

    # --------------------------------------------------------
    # Save.
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

    case_rows.append(
        {
            "case": case,
            "srf": srf,
            "cohesion_kpa":
                c_target / 1000.0,
            "friction_angle_deg":
                phi_target,
            "c_ratio":
                c_ratio,
            "phi_ratio":
                phi_ratio,
        }
    )

    print(
        f"{case:10s} | "
        f"SRF={srf:.2f} | "
        f"c={c_target/1000:.3f} kPa | "
        f"phi={phi_target:.3f} deg"
    )


# ============================================================
# WRITE CASE MANIFEST
# ============================================================

manifest = (
    MODEL
    / "cases.csv"
)

with manifest.open(
    "w",
    newline="",
    encoding="utf-8",
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "case",
            "srf",
            "cohesion_kpa",
            "friction_angle_deg",
            "c_ratio",
            "phi_ratio",
        ],
    )

    writer.writeheader()
    writer.writerows(
        case_rows
    )


print()
print(
    "PHASE 04E STRENGTH-PROBE BUILD: PASS"
)

print(
    "Manifest:",
    manifest,
)

print()
print(
    "NOTE: continuation time is numerical only."
)

print(
    "NOTE: this is a strength-proximity diagnostic, "
    "not a final factor-of-safety claim."
)
