from pathlib import Path

import numpy as np
import pyvista as pv


ROOT = Path.cwd()

MODEL_ROOT = (
    ROOT
    / "model"
    / "phase05d1_restart_mesh"
)

CASES = [
    "coarse",
    "medium",
    "fine",
]


required_point = {
    "bulk_node_ids",
}

required_cell = {
    "bulk_element_ids",
    "bulk_face_ids",
}


for case in CASES:

    folder = (
        MODEL_ROOT
        / case
    )

    path = (
        folder
        / "boundary_all.vtu"
    )

    if not path.exists():

        raise SystemExit(
            f"FAIL: missing {path}"
        )


    boundary = pv.read(
        path
    )

    centres = (
        boundary
        .cell_centers()
        .points
    )

    xmin = float(
        np.min(
            boundary.points[:, 0]
        )
    )

    xmax = float(
        np.max(
            boundary.points[:, 0]
        )
    )

    ymin = float(
        np.min(
            boundary.points[:, 1]
        )
    )

    tol = 1e-8


    left = np.isclose(
        centres[:, 0],
        xmin,
        atol=tol,
    )

    right = np.isclose(
        centres[:, 0],
        xmax,
        atol=tol,
    )

    bottom = np.isclose(
        centres[:, 1],
        ymin,
        atol=tol,
    )

    top = ~(
        left
        | right
        | bottom
    )


    groups = {
        "slope_left.vtu":
            left,

        "slope_right.vtu":
            right,

        "slope_bottom.vtu":
            bottom,

        "slope_top.vtu":
            top,
    }


    print()
    print(
        "========================================"
    )

    print(
        f"BOUNDARY SPLIT {case.upper()}"
    )

    print(
        "========================================"
    )


    for filename, mask in (
        groups.items()
    ):

        ids = np.where(
            mask
        )[0]

        if len(ids) == 0:

            raise SystemExit(
                f"FAIL: empty {filename}"
            )


        sub = boundary.extract_cells(
            ids
        )


        missing_point = (
            required_point
            - set(
                sub.point_data.keys()
            )
        )

        missing_cell = (
            required_cell
            - set(
                sub.cell_data.keys()
            )
        )


        if missing_point:

            raise SystemExit(
                f"FAIL: {filename} missing "
                f"point arrays {missing_point}"
            )

        if missing_cell:

            raise SystemExit(
                f"FAIL: {filename} missing "
                f"cell arrays {missing_cell}"
            )


        out = (
            folder
            / filename
        )

        sub.save(
            out
        )


        print(
            f"{filename}: "
            f"points={sub.n_points}, "
            f"cells={sub.n_cells}"
        )


print()
print(
    "PHASE 05D-1 BOUNDARY SPLIT: PASS"
)
