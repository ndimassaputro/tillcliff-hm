from pathlib import Path
import csv
import xml.etree.ElementTree as ET


ROOT = Path.cwd()

SRC_ROOT = (
    ROOT
    / "model"
    / "phase05d1_restart_mesh"
)

MODEL = (
    ROOT
    / "model"
    / "phase05d2_mesh_erosion"
)

MODEL.mkdir(
    parents=True,
    exist_ok=True,
)

CASES = [
    "coarse",
    "medium",
    "fine",
]


def ensure_curves(root):

    for child in list(root):

        if child.tag == "curves":
            return child

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

        raise SystemExit(
            "FAIL: root process_variables missing"
        )

    root.insert(
        insert_index,
        curves,
    )

    return curves


rows = []


for case in CASES:

    print()
    print(
        "========================================"
    )

    print(
        f"BUILD {case.upper()}"
    )

    print(
        "========================================"
    )


    src_dir = (
        SRC_ROOT
        / case
    )

    src_prj = (
        src_dir
        / "restart_check.prj"
    )

    if not src_prj.exists():

        raise SystemExit(
            f"FAIL: missing {src_prj}"
        )


    dst_dir = (
        MODEL
        / case
    )

    dst_dir.mkdir(
        parents=True,
        exist_ok=True,
    )


    # ========================================================
    # Reuse verified mesh + boundary files from Phase 05D-1.
    # ========================================================

    for name in [
        "bulk.vtu",
        "slope_left.vtu",
        "slope_right.vtu",
        "slope_top.vtu",
        "slope_bottom.vtu",
        "erosion_events.csv",
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

            raise SystemExit(
                f"FAIL: missing {src}"
            )

        dst.write_bytes(
            src.read_bytes()
        )


    # ========================================================
    # Project.
    # ========================================================

    tree = ET.parse(
        src_prj
    )

    root = tree.getroot()


    # ========================================================
    # Install identical moving erosion front on
    # displacement AND pressure equations.
    #
    # No pressure boundary is imposed on the new surface.
    # ========================================================

    pvs = root.find(
        "process_variables"
    )

    if pvs is None:

        raise SystemExit(
            "FAIL: process_variables missing"
        )


    for pv in pvs.findall(
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

        # Horizontal front advancing inland.
        #
        # Start: toe x = 22 m
        # End  : x = 21.2 m
        #
        # Total probe recession = 0.8 m.
        ET.SubElement(
            line,
            "start",
        ).text = "22 2.6 0"

        ET.SubElement(
            line,
            "end",
        ).text = "21.2 2.6 0"


        ET.SubElement(
            sub,
            "material_ids",
        ).text = "1 2 3 4 5"


    # ========================================================
    # Erosion curve.
    #
    # t = 0..20 s:
    # intact equilibration hold.
    #
    # t = 20..100 s:
    # E = 0 -> 0.8 m.
    #
    # dE/dt = 0.01 m/s.
    #
    # t = 100..110 s:
    # hold final geometry.
    #
    # Time is numerical continuation time,
    # NOT physical coastal erosion duration.
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
    ).text = "0 20 100 110"

    ET.SubElement(
        curve,
        "values",
    ).text = "0 0 0.8 0.8"


    # ========================================================
    # Time stepping.
    #
    # Same dt for ALL meshes.
    #
    # dt = 0.1 s
    # numerical front movement / timestep:
    #
    # dE = 0.001 m
    #
    # Much smaller than:
    # coarse dx = 0.20
    # medium dx = 0.10
    # fine dx = 0.05
    # ========================================================

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
    ).text = "110"


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
    ).text = "1100"

    ET.SubElement(
        pair,
        "delta_t",
    ).text = "0.1"


    # ========================================================
    # Newton.
    # ========================================================

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


    # ========================================================
    # Output every 1 s.
    #
    # During erosion:
    # saved dE = 0.01 m.
    # ========================================================

    prefix = root.find(
        "./time_loop/output/prefix"
    )

    if prefix is None:

        raise SystemExit(
            "FAIL: output prefix missing"
        )


    prefix.text = (
        f"phase05d2_{case}"
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
    ).text = "1100"

    ET.SubElement(
        pair,
        "each_steps",
    ).text = "10"


    # ========================================================
    # Save.
    # ========================================================

    ET.indent(
        tree,
        space="    ",
    )


    dst_prj = (
        dst_dir
        / "erosion_probe.prj"
    )


    tree.write(
        dst_prj,
        encoding="UTF-8",
        xml_declaration=True,
    )


    ET.parse(
        dst_prj
    )


    # ========================================================
    # Read mesh event information.
    # ========================================================

    event_rows = []

    with (
        dst_dir
        / "erosion_events.csv"
    ).open(
        "r",
        encoding="utf-8",
    ) as f:

        reader = csv.DictReader(
            f
        )

        for row in reader:

            event_rows.append(
                {
                    "E_m":
                        float(
                            row["E_m"]
                        ),

                    "cell_count":
                        int(
                            row[
                                "cell_count"
                            ]
                        ),

                    "area_m2":
                        float(
                            row[
                                "area_m2"
                            ]
                        ),
                }
            )


    nearest = min(
        event_rows,
        key=lambda r:
            abs(
                r["E_m"]
                - 0.5
            ),
    )


    rows.append(
        {
            "case":
                case,

            "event_near_0p5_E":
                nearest[
                    "E_m"
                ],

            "event_near_0p5_cells":
                nearest[
                    "cell_count"
                ],

            "event_near_0p5_area":
                nearest[
                    "area_m2"
                ],
        }
    )


    print(
        "Project:",
        dst_prj
    )

    print(
        "Event nearest 0.5 m:"
    )

    print(
        f"  E = "
        f"{nearest['E_m']:.6f} m"
    )

    print(
        f"  cells = "
        f"{nearest['cell_count']}"
    )

    print(
        f"  area = "
        f"{nearest['area_m2']:.10f} m2"
    )


# ============================================================
# Manifest.
# ============================================================

manifest = (
    MODEL
    / "mesh_cases.csv"
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
            "event_near_0p5_E",
            "event_near_0p5_cells",
            "event_near_0p5_area",
        ],
    )

    writer.writeheader()
    writer.writerows(
        rows
    )


print()
print(
    "PHASE 05D-2 BUILD: PASS"
)

print(
    "Manifest:",
    manifest,
)
