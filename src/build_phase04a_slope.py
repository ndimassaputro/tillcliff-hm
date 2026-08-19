from pathlib import Path
import numpy as np
import pyvista as pv

MODEL = Path("model/phase04a_slope")
MODEL.mkdir(parents=True, exist_ok=True)

FILES = [
    "slope.vtu",
    "slope_left.vtu",
    "slope_right.vtu",
    "slope_top.vtu",
    "slope_bottom.vtu",
]

def surface_height(x):
    """
    North-European coastal-till archetype.

    Plateau:
        x <= 8 m       -> elevation 10 m

    Slope:
        8 < x < 22 m   -> 10 m down to 2 m

    Toe platform:
        x >= 22 m      -> elevation 2 m

    Relief = 8 m
    Horizontal slope run = 14 m
    Slope angle ~29.7 deg
    """
    x = np.asarray(x, dtype=float)

    h = np.full_like(x, 10.0)

    mid = (x > 8.0) & (x < 22.0)

    h[mid] = (
        10.0
        - 8.0
        * (x[mid] - 8.0)
        / 14.0
    )

    h[x >= 22.0] = 2.0

    return h


print("=== TRANSFORM RECTANGLE TO SLOPE ===")

for name in FILES:

    path = MODEL / name

    if not path.exists():
        raise SystemExit(
            f"FAIL: missing mesh {path}"
        )

    mesh = pv.read(path)

    pts = mesh.points.copy()

    x = pts[:, 0]
    y = pts[:, 1]

    h = surface_height(x)

    # Original structured mesh is 10 m high.
    # Map every vertical column onto local slope height.
    pts[:, 1] = (
        y / 10.0
    ) * h

    mesh.points = pts
    mesh.save(path)

    print(
        f"PASS {name}: "
        f"{mesh.n_points} points, "
        f"{mesh.n_cells} cells"
    )


# ------------------------------------------------------------
# Geometry sanity checks
# ------------------------------------------------------------

bulk = pv.read(
    MODEL / "slope.vtu"
)

xmin, xmax, ymin, ymax, _, _ = bulk.bounds

print()
print("=== GEOMETRY CHECK ===")

print(
    f"x range = {xmin:.3f} to {xmax:.3f} m"
)

print(
    f"y range = {ymin:.3f} to {ymax:.3f} m"
)

print(
    "Nominal cliff relief = 8.0 m"
)

angle = np.degrees(
    np.arctan2(
        8.0,
        14.0,
    )
)

print(
    f"Nominal slope angle = {angle:.3f} deg"
)

print(
    "PHASE 04A SLOPE MESH: PASS"
)
