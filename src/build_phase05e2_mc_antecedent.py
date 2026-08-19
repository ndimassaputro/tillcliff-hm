from pathlib import Path
import re
import shutil
import xml.etree.ElementTree as ET

import numpy as np
import pyvista as pv


ROOT = Path.cwd()

SOURCE_RESULTS = (
    ROOT
    / "results"
    / "phase04b_states"
)

TARGET_TEMPLATE = (
    ROOT
    / "model"
    / "phase05d1_restart_mesh"
    / "medium"
)

TEMPLATE_PRJ = (
    TARGET_TEMPLATE
    / "restart_check.prj"
)

TEMPLATE_BULK = (
    TARGET_TEMPLATE
    / "bulk.vtu"
)

MODEL_ROOT = (
    ROOT
    / "model"
    / "phase05e2_mc_antecedent"
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


# ============================================================
# HELPERS
# ============================================================

def time_from_name(path):

    m = re.search(
        r"_t_([0-9eE+.\-]+)\.vtu$",
        path.name,
    )

    if not m:
        return -np.inf

    return float(
        m.group(1)
    )


def final_vtu(state):

    folder = (
        SOURCE_RESULTS
        / state
    )

    files = list(
        folder.glob("*.vtu")
    )

    if not files:
        raise RuntimeError(
            f"No VTU for {state}"
        )

    return max(
        files,
        key=time_from_name,
    )


def structured_point_grid(
    mesh,
    field_name,
):

    xyz = np.asarray(
        mesh.points,
        dtype=float,
    )

    values = np.asarray(
        mesh.point_data[field_name],
        dtype=float,
    )

    x_values = np.unique(
        np.round(
            xyz[:, 0],
            10,
        )
    )

    x_values.sort()

    columns = []

    top = []

    for x in x_values:

        ids = np.where(
            np.isclose(
                xyz[:, 0],
                x,
                atol=1e-9,
            )
        )[0]

        ids = ids[
            np.argsort(
                xyz[
                    ids,
                    1,
                ]
            )
        ]

        columns.append(
            ids
        )

        top.append(
            float(
                np.max(
                    xyz[
                        ids,
                        1,
                    ]
                )
            )
        )

    counts = [
        len(ids)
        for ids in columns
    ]

    if len(
        set(counts)
    ) != 1:

        raise RuntimeError(
            "Source is not a uniform "
            "structured-column mesh"
        )

    top = np.asarray(
        top,
        dtype=float,
    )

    eta = (
        xyz[
            columns[0],
            1,
        ]
        / top[0]
    )

    if values.ndim == 1:

        grid = np.empty(
            (
                len(x_values),
                len(eta),
            ),
            dtype=float,
        )

    else:

        grid = np.empty(
            (
                len(x_values),
                len(eta),
                values.shape[1],
            ),
            dtype=float,
        )

    for i, ids in enumerate(
        columns
    ):

        grid[i] = values[
            ids
        ]

    return (
        x_values,
        eta,
        top,
        grid,
    )


def bilinear(
    x_grid,
    eta_grid,
    values,
    x,
    eta,
):

    i = int(
        np.searchsorted(
            x_grid,
            x,
            side="right",
        )
        - 1
    )

    j = int(
        np.searchsorted(
            eta_grid,
            eta,
            side="right",
        )
        - 1
    )

    i = int(
        np.clip(
            i,
            0,
            len(x_grid) - 2,
        )
    )

    j = int(
        np.clip(
            j,
            0,
            len(eta_grid) - 2,
        )
    )

    x0 = x_grid[i]
    x1 = x_grid[i + 1]

    e0 = eta_grid[j]
    e1 = eta_grid[j + 1]

    wx = (
        (x - x0)
        / (x1 - x0)
        if x1 != x0
        else 0.0
    )

    we = (
        (eta - e0)
        / (e1 - e0)
        if e1 != e0
        else 0.0
    )

    wx = float(
        np.clip(
            wx,
            0.0,
            1.0,
        )
    )

    we = float(
        np.clip(
            we,
            0.0,
            1.0,
        )
    )

    v00 = values[
        i,
        j,
    ]

    v10 = values[
        i + 1,
        j,
    ]

    v01 = values[
        i,
        j + 1,
    ]

    v11 = values[
        i + 1,
        j + 1,
    ]

    return (
        (1 - wx)
        * (1 - we)
        * v00
        + wx
        * (1 - we)
        * v10
        + (1 - wx)
        * we
        * v01
        + wx
        * we
        * v11
    )


def replace_or_add_constant_parameter(
    root,
    name,
    value,
):

    parameters = root.find(
        "parameters"
    )

    if parameters is None:

        raise RuntimeError(
            "parameters block missing"
        )

    for p in list(
        parameters.findall(
            "parameter"
        )
    ):

        if (
            p.findtext("name")
            == name
        ):

            parameters.remove(
                p
            )

    parameter = ET.SubElement(
        parameters,
        "parameter",
    )

    ET.SubElement(
        parameter,
        "name",
    ).text = name

    ET.SubElement(
        parameter,
        "type",
    ).text = "Constant"

    ET.SubElement(
        parameter,
        "value",
    ).text = (
        f"{value:.16g}"
    )


def set_top_pressure_bc(
    root,
    parameter_name,
):

    pressure_pv = None

    for pv in root.findall(
        "./process_variables/"
        "process_variable"
    ):

        if (
            pv.findtext("name")
            == "pressure"
        ):

            pressure_pv = pv
            break

    if pressure_pv is None:

        raise RuntimeError(
            "pressure process variable missing"
        )

    candidates = []

    for bc in pressure_pv.findall(
        "./boundary_conditions/"
        "boundary_condition"
    ):

        mesh_name = (
            bc.findtext("mesh")
            or ""
        )

        bc_type = (
            bc.findtext("type")
            or ""
        )

        if (
            bc_type == "Dirichlet"
            and "top"
            in mesh_name.lower()
        ):

            candidates.append(
                bc
            )

    if len(candidates) != 1:

        raise RuntimeError(
            "Expected exactly one top "
            "pressure Dirichlet BC, got "
            f"{len(candidates)}"
        )

    bc = candidates[0]

    parameter = bc.find(
        "parameter"
    )

    if parameter is None:

        parameter = ET.SubElement(
            bc,
            "parameter",
        )

    parameter.text = (
        parameter_name
    )


# ============================================================
# CHECK TEMPLATE
# ============================================================

for path in [
    TEMPLATE_PRJ,
    TEMPLATE_BULK,
]:

    if not path.exists():

        raise RuntimeError(
            f"Missing {path}"
        )


target_template = pv.read(
    TEMPLATE_BULK
)


# ============================================================
# BUILD STATES
# ============================================================

print(
    "========================================"
)

print(
    "PHASE 05E-2 MC ANTECEDENT BUILD"
)

print(
    "========================================"
)


for state in STATES:

    print()
    print(
        f"=== {state.upper()} ==="
    )

    source_path = final_vtu(
        state
    )

    source = pv.read(
        source_path
    )

    for field in [
        "pressure",
        "saturation",
        "sigma",
    ]:

        if field not in source.point_data:

            raise RuntimeError(
                f"{state}: missing {field}"
            )


    # --------------------------------------------------------
    # Structured source fields.
    # --------------------------------------------------------

    (
        x_source,
        eta_source,
        top_source,
        p_grid,
    ) = structured_point_grid(
        source,
        "pressure",
    )

    (
        x_sigma,
        eta_sigma,
        top_sigma,
        sigma_grid,
    ) = structured_point_grid(
        source,
        "sigma",
    )

    (
        x_sat,
        eta_sat,
        top_sat,
        sat_grid,
    ) = structured_point_grid(
        source,
        "saturation",
    )


    if not np.allclose(
        x_source,
        x_sigma,
    ):

        raise RuntimeError(
            f"{state}: x-grid mismatch"
        )

    if not np.allclose(
        eta_source,
        eta_sigma,
    ):

        raise RuntimeError(
            f"{state}: eta-grid mismatch"
        )


    # --------------------------------------------------------
    # Target mesh copy.
    # --------------------------------------------------------

    target = target_template.copy(
        deep=True
    )

    target_points = np.asarray(
        target.points,
        dtype=float,
    )

    target_centers = (
        target
        .cell_centers()
        .points
    )


    # --------------------------------------------------------
    # Target surface height by x.
    # --------------------------------------------------------

    target_x_unique = np.unique(
        np.round(
            target_points[:, 0],
            10,
        )
    )

    target_top = {}

    for xv in target_x_unique:

        ids = np.where(
            np.isclose(
                target_points[:, 0],
                xv,
                atol=1e-9,
            )
        )[0]

        target_top[
            float(xv)
        ] = float(
            np.max(
                target_points[
                    ids,
                    1,
                ]
            )
        )


    def top_at_x(x):

        return float(
            np.interp(
                x,
                x_source,
                top_source,
            )
        )


    # --------------------------------------------------------
    # p0 + diagnostic saturation at target nodes.
    # --------------------------------------------------------

    p0 = np.empty(
        target.n_points,
        dtype=float,
    )

    sat0 = np.empty(
        target.n_points,
        dtype=float,
    )


    for pid, xyz in enumerate(
        target_points
    ):

        x = float(
            xyz[0]
        )

        h = top_at_x(
            x
        )

        eta = (
            float(
                xyz[1]
            )
            / h
            if h > 0
            else 0.0
        )

        p0[pid] = bilinear(
            x_source,
            eta_source,
            p_grid,
            x,
            eta,
        )

        sat0[pid] = bilinear(
            x_sat,
            eta_sat,
            sat_grid,
            x,
            eta,
        )


    # --------------------------------------------------------
    # sigma0 at target element centres.
    #
    # Source sigma is nodal/extrapolated output.
    # We project it bilinearly to each target
    # cell centre to obtain MeshElement sigma0.
    # --------------------------------------------------------

    sigma0 = np.empty(
        (
            target.n_cells,
            sigma_grid.shape[2],
        ),
        dtype=float,
    )


    for cid, xyz in enumerate(
        target_centers
    ):

        x = float(
            xyz[0]
        )

        h = top_at_x(
            x
        )

        eta = (
            float(
                xyz[1]
            )
            / h
            if h > 0
            else 0.0
        )

        sigma0[
            cid,
            :,
        ] = bilinear(
            x_sigma,
            eta_sigma,
            sigma_grid,
            x,
            eta,
        )


    target.point_data[
        "p0"
    ] = p0

    target.point_data[
        "antecedent_saturation"
    ] = sat0

    target.cell_data[
        "sigma0"
    ] = sigma0


    # --------------------------------------------------------
    # Determine terminal surface pressure.
    # --------------------------------------------------------

    source_xyz = np.asarray(
        source.points,
        dtype=float,
    )

    top_mask = np.zeros(
        source.n_points,
        dtype=bool,
    )

    for i, x in enumerate(
        x_source
    ):

        ids = np.where(
            np.isclose(
                source_xyz[:, 0],
                x,
                atol=1e-9,
            )
        )[0]

        top_id = ids[
            np.argmax(
                source_xyz[
                    ids,
                    1,
                ]
            )
        ]

        top_mask[
            top_id
        ] = True


    top_pressures = np.asarray(
        source.point_data[
            "pressure"
        ],
        dtype=float,
    )[
        top_mask
    ]


    surface_pressure = float(
        np.median(
            top_pressures
        )
    )

    surface_spread = float(
        np.max(
            top_pressures
        )
        - np.min(
            top_pressures
        )
    )


    # --------------------------------------------------------
    # Save state mesh.
    # --------------------------------------------------------

    state_dir = (
        MODEL_ROOT
        / state
    )

    state_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    bulk_path = (
        state_dir
        / "bulk.vtu"
    )

    target.save(
        bulk_path
    )


    # --------------------------------------------------------
    # Copy verified medium boundaries.
    # --------------------------------------------------------

    for name in [
        "slope_left.vtu",
        "slope_right.vtu",
        "slope_top.vtu",
        "slope_bottom.vtu",
    ]:

        shutil.copy2(
            TARGET_TEMPLATE
            / name,
            state_dir
            / name,
        )


    # --------------------------------------------------------
    # Build MC intact-hold project.
    # --------------------------------------------------------

    tree = ET.parse(
        TEMPLATE_PRJ
    )

    root = tree.getroot()


    # Remove any erosion blocks defensively.
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


    replace_or_add_constant_parameter(
        root,
        "antecedent_surface_pressure",
        surface_pressure,
    )

    set_top_pressure_bc(
        root,
        "antecedent_surface_pressure",
    )


    # --------------------------------------------------------
    # 20 s intact MC equilibrium hold.
    # --------------------------------------------------------

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
    ).text = "20"

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
    ).text = "20"

    ET.SubElement(
        pair,
        "delta_t",
    ).text = "1"


    prefix = root.find(
        "./time_loop/output/prefix"
    )

    if prefix is None:

        raise RuntimeError(
            "output prefix missing"
        )

    prefix.text = (
        f"phase05e2_{state}"
    )


    ET.indent(
        tree,
        space="    ",
    )

    prj_path = (
        state_dir
        / "mc_antecedent_hold.prj"
    )

    tree.write(
        prj_path,
        encoding="UTF-8",
        xml_declaration=True,
    )

    ET.parse(
        prj_path
    )


    # --------------------------------------------------------
    # Diagnostics.
    # --------------------------------------------------------

    print(
        "Source final VTU:"
    )

    print(
        source_path
    )

    print(
        f"Source final time [s]: "
        f"{time_from_name(source_path):.0f}"
    )

    print(
        f"Target points/cells: "
        f"{target.n_points}/"
        f"{target.n_cells}"
    )

    print(
        "Transferred pressure [kPa]: "
        f"{np.min(p0)/1000:.6f} "
        f"to "
        f"{np.max(p0)/1000:.6f}"
    )

    print(
        "Transferred mean pressure [kPa]: "
        f"{np.mean(p0)/1000:.6f}"
    )

    print(
        "Transferred mean saturation: "
        f"{np.mean(sat0):.8f}"
    )

    print(
        "Surface pressure BC [kPa]: "
        f"{surface_pressure/1000:.6f}"
    )

    print(
        "Source top-pressure spread [Pa]: "
        f"{surface_spread:.6e}"
    )

    print(
        "sigma0 range [Pa]: "
        f"{np.min(sigma0):.6e} "
        f"to "
        f"{np.max(sigma0):.6e}"
    )

    print(
        "PASS BUILD:",
        prj_path,
    )


print()
print(
    "PHASE 05E-2 BUILD COMPLETE"
)
