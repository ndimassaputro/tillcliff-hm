from pathlib import Path
import shutil
import xml.etree.ElementTree as ET


ROOT = Path.cwd()

SOURCE_ROOT = (
    ROOT
    / "model"
    / "phase05e2_mc_antecedent"
)

MODEL_ROOT = (
    ROOT
    / "model"
    / "phase05e3_antecedent_erosion"
)

MODEL_ROOT.mkdir(
    parents=True,
    exist_ok=True,
)

STATES = [
    "dry",
    "reference",
    "wet",
]

E_MAX = 0.55

HOLD_END = 20.0

EROSION_RATE = 0.01

EROSION_END = (
    HOLD_END
    + E_MAX / EROSION_RATE
)


def ensure_curves(root):

    curves = root.find(
        "curves"
    )

    if curves is not None:
        return curves

    curves = ET.Element(
        "curves"
    )

    children = list(root)

    insert_index = None

    for i, child in enumerate(
        children
    ):

        if child.tag == "process_variables":

            insert_index = i
            break

    if insert_index is None:

        raise RuntimeError(
            "Cannot locate insertion point "
            "for curves block"
        )

    root.insert(
        insert_index,
        curves,
    )

    return curves


print(
    "========================================"
)

print(
    "PHASE 05E-3 BUILD"
)

print(
    "========================================"
)

print(
    f"E_MAX = {E_MAX:.3f} m"
)

print(
    f"EROSION END TIME = "
    f"{EROSION_END:.3f} s"
)


for state in STATES:

    print()
    print(
        f"=== {state.upper()} ==="
    )


    src_dir = (
        SOURCE_ROOT
        / state
    )

    dst_dir = (
        MODEL_ROOT
        / state
    )

    dst_dir.mkdir(
        parents=True,
        exist_ok=True,
    )


    src_prj = (
        src_dir
        / "mc_antecedent_hold.prj"
    )

    if not src_prj.exists():

        raise RuntimeError(
            f"Missing {src_prj}"
        )


    # ========================================================
    # COPY VERIFIED MESH + BOUNDARIES
    # ========================================================

    for name in [
        "bulk.vtu",
        "slope_left.vtu",
        "slope_right.vtu",
        "slope_top.vtu",
        "slope_bottom.vtu",
    ]:

        src = (
            src_dir
            / name
        )

        dst = (
            dst_dir
            / name
        )

        if not src.exists():

            raise RuntimeError(
                f"Missing {src}"
            )

        shutil.copy2(
            src,
            dst,
        )


    # ========================================================
    # LOAD PROJECT
    # ========================================================

    tree = ET.parse(
        src_prj
    )

    root = tree.getroot()


    # ========================================================
    # MOVING TOE-EROSION FRONT
    #
    # Same geometry + rate for every antecedent state.
    #
    # Numerical continuation time only.
    # NOT physical erosion duration.
    # ========================================================

    pvs = root.find(
        "process_variables"
    )

    if pvs is None:

        raise RuntimeError(
            "process_variables missing"
        )


    for pv in pvs.findall(
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


        old = pv.find(
            "deactivated_subdomains"
        )

        if old is not None:

            pv.remove(
                old
            )


        block = ET.SubElement(
            pv,
            "deactivated_subdomains",
        )

        sub = ET.SubElement(
            block,
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
        ).text = "21 2.6 0"


        ET.SubElement(
            sub,
            "material_ids",
        ).text = "1 2 3 4 5"


    # ========================================================
    # CURVE
    #
    # 0 -> 20 s:
    # intact equilibrium hold.
    #
    # 20 -> 75 s:
    # recession 0 -> 0.55 m.
    # ========================================================

    curves = ensure_curves(
        root
    )


    for curve in list(
        curves.findall(
            "curve"
        )
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
        f"0 "
        f"{HOLD_END:.8f} "
        f"{EROSION_END:.8f}"
    )

    ET.SubElement(
        curve,
        "values",
    ).text = (
        f"0 0 {E_MAX:.8f}"
    )


    # ========================================================
    # TIMESTEPPING
    #
    # dt = 0.1 s
    # dE per timestep = 0.001 m.
    # ========================================================

    ts = root.find(
        "./time_loop/processes/"
        "process/time_stepping"
    )

    if ts is None:

        raise RuntimeError(
            "time stepping missing"
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
    ).text = (
        f"{EROSION_END:.8f}"
    )


    steps = ET.SubElement(
        ts,
        "timesteps",
    )

    pair = ET.SubElement(
        steps,
        "pair",
    )

    repeat = int(
        round(
            EROSION_END
            / 0.1
        )
    )

    ET.SubElement(
        pair,
        "repeat",
    ).text = str(
        repeat
    )

    ET.SubElement(
        pair,
        "delta_t",
    ).text = "0.1"


    # ========================================================
    # NONLINEAR ITERATIONS
    # ========================================================

    solver = root.find(
        "./nonlinear_solvers/"
        "nonlinear_solver"
    )

    if solver is None:

        raise RuntimeError(
            "nonlinear solver missing"
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


    # ========================================================
    # OUTPUT EVERY 1 SECOND
    #
    # During erosion:
    # nominal saved dE = 0.01 m.
    # ========================================================

    prefix = root.find(
        "./time_loop/output/prefix"
    )

    if prefix is None:

        raise RuntimeError(
            "output prefix missing"
        )

    prefix.text = (
        f"phase05e3_{state}"
    )


    out_ts = root.find(
        "./time_loop/output/timesteps"
    )

    if out_ts is None:

        raise RuntimeError(
            "output timesteps missing"
        )

    out_ts.clear()

    pair = ET.SubElement(
        out_ts,
        "pair",
    )

    ET.SubElement(
        pair,
        "repeat",
    ).text = str(
        repeat
    )

    ET.SubElement(
        pair,
        "each_steps",
    ).text = "10"


    # ========================================================
    # ENSURE USEFUL OUTPUT
    # ========================================================

    variables = root.find(
        "./time_loop/output/variables"
    )

    if variables is None:

        raise RuntimeError(
            "output variables missing"
        )

    existing = {
        v.text
        for v in variables.findall(
            "variable"
        )
    }

    for name in [
        "pressure",
        "saturation",
        "displacement",
        "sigma",
        "EquivalentPlasticStrain",
    ]:

        if name not in existing:

            ET.SubElement(
                variables,
                "variable",
            ).text = name


    # ========================================================
    # SAVE
    # ========================================================

    ET.indent(
        tree,
        space="    ",
    )

    dst_prj = (
        dst_dir
        / "antecedent_erosion.prj"
    )

    tree.write(
        dst_prj,
        encoding="UTF-8",
        xml_declaration=True,
    )

    ET.parse(
        dst_prj
    )


    print(
        "PASS:",
        dst_prj
    )


print()
print(
    "PHASE 05E-3 BUILD COMPLETE"
)
