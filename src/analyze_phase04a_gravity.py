from pathlib import Path

import numpy as np
import ogstools as ot


ROOT = Path.cwd()

RUN = (
    ROOT
    / "results"
    / "phase04a_gravity"
)

pvds = sorted(
    RUN.glob("*.pvd")
)

if not pvds:
    raise SystemExit(
        "FAIL: no phase04a PVD found"
    )

series = ot.MeshSeries(
    str(pvds[0])
)

times = np.asarray(
    series.timevalues
)

mesh = series.mesh(
    len(times) - 1
)


def get(mesh, name):

    if name in mesh.point_data:
        return np.asarray(
            mesh.point_data[name]
        )

    if name in mesh.cell_data:
        return np.asarray(
            mesh.cell_data[name]
        )

    raise KeyError(
        f"{name} not found"
    )


disp = get(
    mesh,
    "displacement",
)

sat = get(
    mesh,
    "saturation",
).squeeze()

por = get(
    mesh,
    "porosity",
).squeeze()

pressure = get(
    mesh,
    "pressure",
).squeeze()


if disp.ndim == 1:
    disp_mag = np.abs(
        disp
    )
else:
    disp_mag = np.linalg.norm(
        disp,
        axis=1,
    )


summary = Path(
    "results/phase04a_gravity/"
    "phase04a_summary.txt"
)

lines = [
    "PHASE 04A 2D SLOPE GRAVITY BASELINE",
    "",
    f"Timesteps: {len(times)}",
    (
        "Final time [days]: "
        f"{times[-1]/86400:.6f}"
    ),
    "",
    (
        "Maximum displacement [mm]: "
        f"{np.nanmax(disp_mag)*1000:.8f}"
    ),
    (
        "Pressure range [kPa]: "
        f"{np.nanmin(pressure)/1000:.6f} "
        f"to "
        f"{np.nanmax(pressure)/1000:.6f}"
    ),
    (
        "Saturation range [-]: "
        f"{np.nanmin(sat):.8f} "
        f"to "
        f"{np.nanmax(sat):.8f}"
    ),
    (
        "Porosity range [-]: "
        f"{np.nanmin(por):.8f} "
        f"to "
        f"{np.nanmax(por):.8f}"
    ),
    "",
]


finite = (
    np.all(
        np.isfinite(disp_mag)
    )
    and np.all(
        np.isfinite(sat)
    )
    and np.all(
        np.isfinite(por)
    )
    and np.all(
        np.isfinite(pressure)
    )
)


if not finite:
    status = "FAIL: non-finite solution"

elif np.nanmin(sat) < -1e-8:
    status = "FAIL: negative saturation"

elif np.nanmax(sat) > 1.000001:
    status = "FAIL: saturation > 1"

else:
    status = (
        "PHASE 04A GRAVITY BASELINE: PASS"
    )


lines.append(
    status
)

summary.write_text(
    "\n".join(lines) + "\n",
    encoding="utf-8",
)

print(
    "\n".join(lines)
)
