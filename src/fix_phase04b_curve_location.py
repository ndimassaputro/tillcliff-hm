from pathlib import Path
import xml.etree.ElementTree as ET


MODEL = Path("model/phase04b_states")

states = {
    "dry": 1.20,
    "reference": 1.00,
    "wet": 0.30,
}


def indent(elem, level=0):
    """
    Pretty-print XML without requiring manual edits.
    """
    space = "\n" + level * "    "

    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = space + "    "

        for child in elem:
            indent(child, level + 1)

        if not elem[-1].tail or not elem[-1].tail.strip():
            elem[-1].tail = space

    if level and (
        not elem.tail
        or not elem.tail.strip()
    ):
        elem.tail = space


for state, multiplier in states.items():

    path = MODEL / f"{state}.prj"

    if not path.exists():
        raise SystemExit(
            f"FAIL: missing {path}"
        )

    tree = ET.parse(path)
    root = tree.getroot()

    if root.tag != "OpenGeoSysProject":
        raise SystemExit(
            f"FAIL: unexpected root tag in {path}"
        )

    # ========================================================
    # 1. REMOVE EVERY EXISTING antecedent curve block,
    #    including the wrongly nested one.
    # ========================================================

    removed = 0

    for parent in root.iter():
        for child in list(parent):

            if child.tag != "curves":
                continue

            names = [
                c.findtext("name")
                for c in child.findall("curve")
            ]

            if "antecedent_state_curve" in names:
                parent.remove(child)
                removed += 1

    print(
        f"{state}: removed misplaced curve blocks = {removed}"
    )

    # ========================================================
    # 2. CREATE CURVE AS DIRECT CHILD OF PROJECT ROOT
    # ========================================================

    curves = ET.Element("curves")

    curve = ET.SubElement(
        curves,
        "curve",
    )

    ET.SubElement(
        curve,
        "name",
    ).text = "antecedent_state_curve"

    ET.SubElement(
        curve,
        "coords",
    ).text = (
        "0 "
        "2592000 "
        "5184000 "
        "10368000"
    )

    ET.SubElement(
        curve,
        "values",
    ).text = (
        f"1.0 "
        f"1.0 "
        f"{multiplier:.8f} "
        f"{multiplier:.8f}"
    )

    # Put <curves> immediately before the ROOT-level
    # <process_variables>.
    root_children = list(root)

    root_pv_index = None

    for i, child in enumerate(root_children):
        if child.tag == "process_variables":
            root_pv_index = i
            break

    if root_pv_index is None:
        raise SystemExit(
            f"FAIL: root-level process_variables "
            f"not found in {state}"
        )

    root.insert(
        root_pv_index,
        curves,
    )

    indent(root)

    tree.write(
        path,
        encoding="UTF-8",
        xml_declaration=True,
    )

    # ========================================================
    # 3. HARD STRUCTURE VALIDATION
    # ========================================================

    check_tree = ET.parse(path)
    check_root = check_tree.getroot()

    direct_curves = [
        c
        for c in list(check_root)
        if c.tag == "curves"
    ]

    matching = []

    for block in direct_curves:
        for c in block.findall("curve"):

            if (
                c.findtext("name")
                == "antecedent_state_curve"
            ):
                matching.append(c)

    if len(matching) != 1:
        raise SystemExit(
            f"FAIL {state}: expected exactly "
            f"one ROOT-level antecedent curve, "
            f"found {len(matching)}"
        )

    curve = matching[0]

    coords = curve.findtext(
        "coords"
    ).split()

    values = curve.findtext(
        "values"
    ).split()

    if len(coords) != 4:
        raise SystemExit(
            f"FAIL {state}: bad coords"
        )

    if len(values) != 4:
        raise SystemExit(
            f"FAIL {state}: bad values"
        )

    print(
        f"PASS {state:9s}: "
        f"ROOT-level curve installed, "
        f"target suction = "
        f"{100 * multiplier:.1f} kPa"
    )


print()
print(
    "PHASE 04B CURVE-LOCATION FIX: PASS"
)
