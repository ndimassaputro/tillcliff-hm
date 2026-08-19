from pathlib import Path
import shutil
import xml.etree.ElementTree as ET

import numpy as np
import pyvista as pv


ROOT = Path.cwd()

SOURCE = (
    ROOT
    / "model"
    / "phase05e2_mc_antecedent"
)

MODEL = (
    ROOT
    / "model"
    / "phase05e4_decomposition"
)

MODEL.mkdir(
    parents=True,
    exist_ok=True,
)


CASES = {
    "refP_refS": {
        "pressure_state": "reference",
        "stress_state": "reference",
    },
    "dryP_refS": {
        "pressure_state": "dry",
        "stress_state": "reference",
    },
    "wetP_refS": {
        "pressure_state": "wet",
        "stress_state": "reference",
    },
    "refP_dryS": {
        "pressure_state": "reference",
        "stress_state": "dry",
    },
    "refP_wetS": {
        "pressure_state": "reference",
        "stress_state": "wet",
    },
}


HOLD_END = 20.0
E_MAX = 0.25
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

    children = list(
        root
    )

    idx = next(
        (
            i
            for i, child in enumerate(
                children
            )
            if child.tag
            == "process_variables"
        ),
        None,
    )

    if idx is None:

        raise RuntimeError(
            "Cannot insert curves block"
        )


    root.insert(
        idx,
        curves,
    )

    return curves


print(
    "========================================"
)

print(
    "PHASE 05E-4 PRESSURE/STRESS DECOMPOSITION"
)

print(
    "========================================"
)


for case, config in CASES.items():

    p_state = config[
        "pressure_state"
    ]

    s_state = config[
        "stress_state"
    ]


    print()
    print(
        f"=== {case} ==="
    )

    print(
        f"pressure state = {p_state}"
    )

    print(
        f"stress state   = {s_state}"
    )


    # ========================================================
    # INPUT FILES
    # ========================================================

    p_mesh_path = (
        SOURCE
        / p_state
        / "bulk.vtu"
    )

    s_mesh_path = (
        SOURCE
        / s_state
        / "bulk.vtu"
    )

    project_path = (
        SOURCE
        / p_state
        / "mc_antecedent_hold.prj"
    )


    for path in [
        p_mesh_path,
        s_mesh_path,
        project_path,
    ]:

        if not path.exists():

            raise RuntimeError(
                f"Missing {path}"
            )


    # ========================================================
    # LOAD PRESSURE-STATE AND STRESS-STATE MESHES
    #
    # IMPORTANT:
    # pv here ALWAYS remains the PyVista module.
    # XML loop variables below use pv_node.
    # ========================================================

    p_mesh = pv.read(
        p_mesh_path
    )

    s_mesh = pv.read(
        s_mesh_path
    )


    # ========================================================
    # SAME GEOMETRY IS MANDATORY
    # ========================================================

    if (
        p_mesh.n_points
        != s_mesh.n_points
        or p_mesh.n_cells
        != s_mesh.n_cells
    ):

        raise RuntimeError(
            f"{case}: topology mismatch"
        )


    coord_diff = float(
        np.max(
            np.abs(
                np.asarray(
                    p_mesh.points,
                    dtype=float,
                )
                - np.asarray(
                    s_mesh.points,
                    dtype=float,
                )
            )
        )
    )


    if coord_diff > 1e-10:

        raise RuntimeError(
            f"{case}: geometry mismatch "
            f"{coord_diff}"
        )


    # ========================================================
    # REQUIRED FIELDS
    # ========================================================

    if "p0" not in p_mesh.point_data:

        raise RuntimeError(
            f"{case}: p0 missing "
            f"from pressure-state mesh"
        )


    if (
        "antecedent_saturation"
        not in p_mesh.point_data
    ):

        raise RuntimeError(
            f"{case}: "
            "antecedent_saturation missing "
            "from pressure-state mesh"
        )


    if "sigma0" not in s_mesh.cell_data:

        raise RuntimeError(
            f"{case}: sigma0 missing "
            f"from stress-state mesh"
        )


    if (
        "MaterialIDs"
        not in p_mesh.cell_data
    ):

        raise RuntimeError(
            f"{case}: MaterialIDs missing"
        )


    # ========================================================
    # BUILD FACTORIAL STATE
    #
    # Pressure branch supplies:
    #   p0
    #   antecedent saturation
    #   hydraulic boundary condition
    #
    # Stress branch supplies:
    #   sigma0
    #
    # Everything else follows the verified medium mesh.
    # ========================================================

    mesh = p_mesh.copy(
        deep=True
    )


    mesh.cell_data[
        "sigma0"
    ] = np.asarray(
        s_mesh.cell_data[
            "sigma0"
        ],
        dtype=float,
    ).copy()


    out_dir = (
        MODEL
        / case
    )

    out_dir.mkdir(
        parents=True,
        exist_ok=True,
    )


    mesh.save(
        out_dir
        / "bulk.vtu"
    )


    # ========================================================
    # COPY VERIFIED MEDIUM-MESH BOUNDARIES
    # ========================================================

    for name in [
        "slope_left.vtu",
        "slope_right.vtu",
        "slope_top.vtu",
        "slope_bottom.vtu",
    ]:

        src_boundary = (
            SOURCE
            / p_state
            / name
        )

        if not src_boundary.exists():

            raise RuntimeError(
                f"{case}: missing "
                f"{src_boundary}"
            )


        shutil.copy2(
            src_boundary,
            out_dir
            / name,
        )


    # ========================================================
    # PROJECT INHERITS HYDRAULIC BC
    # FROM PRESSURE BRANCH
    # ========================================================

    tree = ET.parse(
        project_path
    )

    root = tree.getroot()


    # ========================================================
    # REMOVE OLD DEACTIVATION BLOCKS
    # ========================================================

    for pv_node in root.findall(
        "./process_variables/"
        "process_variable"
    ):

        old = pv_node.find(
            "deactivated_subdomains"
        )

        if old is not None:

            pv_node.remove(
                old
            )


    # ========================================================
    # MOVING TOE-EROSION FRONT
    #
    # Same geometry for all factorial cases.
    # ========================================================

    for pv_node in root.findall(
        "./process_variables/"
        "process_variable"
    ):

        name = pv_node.findtext(
            "name"
        )

        if name not in {
            "pressure",
            "displacement",
        }:

            continue


        block = ET.SubElement(
            pv_node,
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
    # EROSION CURVE
    #
    # t = 0..20 s:
    # intact equilibrium hold.
    #
    # t = 20..45 s:
    # nominal E = 0..0.25 m.
    #
    # Numerical continuation time only.
    # ========================================================

    curves = ensure_curves(
        root
    )


    for old_curve in list(
        curves.findall(
            "curve"
        )
    ):

        if (
            old_curve.findtext(
                "name"
            )
            == "erosion_front_curve"
        ):

            curves.remove(
                old_curve
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
    # dE = 0.001 m per step during erosion.
    # ========================================================

    repeat = int(
        round(
            EROSION_END
            / 0.1
        )
    )


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
    # OUTPUT EVERY 1 s
    #
    # saved dE = 0.01 m during erosion.
    # ========================================================

    prefix = root.find(
        "./time_loop/output/prefix"
    )


    if prefix is None:

        raise RuntimeError(
            "output prefix missing"
        )


    prefix.text = (
        f"phase05e4_{case}"
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
    # ENSURE REQUIRED OUTPUT VARIABLES
    # ========================================================

    variables = root.find(
        "./time_loop/output/variables"
    )


    if variables is None:

        raise RuntimeError(
            "output variables missing"
        )


    existing = {
        variable.text
        for variable in variables.findall(
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
    # SAVE PROJECT
    # ========================================================

    ET.indent(
        tree,
        space="    ",
    )


    prj = (
        out_dir
        / "decomposition.prj"
    )


    tree.write(
        prj,
        encoding="UTF-8",
        xml_declaration=True,
    )


    # Parse once again as XML validation.
    ET.parse(
        prj
    )


    # ========================================================
    # BUILD DIAGNOSTICS
    # ========================================================

    p0 = np.asarray(
        mesh.point_data[
            "p0"
        ],
        dtype=float,
    )


    sigma0 = np.asarray(
        mesh.cell_data[
            "sigma0"
        ],
        dtype=float,
    )


    saturation0 = np.asarray(
        mesh.point_data[
            "antecedent_saturation"
        ],
        dtype=float,
    )


    print(
        "p0 mean [kPa] = "
        f"{np.mean(p0)/1000:.6f}"
    )


    print(
        "saturation mean = "
        f"{np.mean(saturation0):.8f}"
    )


    print(
        "sigma0 mean [Pa] = "
        f"{np.mean(sigma0):.6f}"
    )


    print(
        "coordinate consistency [m] = "
        f"{coord_diff:.3e}"
    )


    print(
        "PASS:",
        prj,
    )


print()
print(
    "========================================"
)

print(
    "PHASE 05E-4 BUILD COMPLETE"
)

print(
    "========================================"
)
