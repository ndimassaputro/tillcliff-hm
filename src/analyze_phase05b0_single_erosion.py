from pathlib import Path

import numpy as np
import ogstools as ot


ROOT = Path.cwd()

RUN = (
    ROOT
    / "results"
    / "phase05b0_single_erosion"
)

pvds = sorted(
    RUN.glob("*.pvd")
)

if not pvds:
    raise SystemExit(
        "FAIL: no PVD output"
    )


series = ot.MeshSeries(
    str(pvds[0])
)

times = np.asarray(
    series.timevalues,
    dtype=float,
)


def index_at_time(
    target,
):

    return int(
        np.argmin(
            np.abs(
                times
                - target
            )
        )
    )


i_pre = index_at_time(
    9.0
)

i_event = index_at_time(
    10.0
)

i_post = index_at_time(
    20.0
)


pre = series.mesh(
    i_pre
)

event = series.mesh(
    i_event
)

post = series.mesh(
    i_post
)


def field(
    mesh,
    name,
):

    if name in mesh.point_data:

        return np.asarray(
            mesh.point_data[name],
            dtype=float,
        )

    if name in mesh.cell_data:

        return np.asarray(
            mesh.cell_data[name],
            dtype=float,
        )

    raise KeyError(
        f"{name!r} missing"
    )


def max_displacement(
    mesh,
):

    u = field(
        mesh,
        "displacement",
    )

    mag = np.linalg.norm(
        u,
        axis=1,
    )

    return float(
        np.nanmax(
            mag
        )
    )


def plastic_metrics(
    mesh,
):

    epsp = field(
        mesh,
        "EquivalentPlasticStrain",
    ).squeeze()

    finite = epsp[
        np.isfinite(
            epsp
        )
    ]

    if finite.size == 0:

        return (
            np.nan,
            np.nan,
            np.nan,
        )

    return (
        float(
            np.nanmax(
                finite
            )
        ),

        float(
            np.mean(
                finite > 1e-6
            )
        ),

        float(
            np.mean(
                finite > 1e-4
            )
        ),
    )


u_pre = max_displacement(
    pre
)

u_event = max_displacement(
    event
)

u_post = max_displacement(
    post
)


epsp_pre, f6_pre, f4_pre = (
    plastic_metrics(
        pre
    )
)

epsp_post, f6_post, f4_post = (
    plastic_metrics(
        post
    )
)


# ============================================================
# HYDRAULIC CHECK
# ============================================================

pressure_pre = field(
    pre,
    "pressure",
).squeeze()

pressure_post = field(
    post,
    "pressure",
).squeeze()


# ============================================================
# PRINT
# ============================================================

print(
    "========================================"
)

print(
    "PHASE 05B-0 SINGLE TOE-EROSION TEST"
)

print(
    "========================================"
)

print()

print(
    f"Output times: {len(times)}"
)

print(
    f"Pre-event time: "
    f"{times[i_pre]:.6f} s"
)

print(
    f"Event time: "
    f"{times[i_event]:.6f} s"
)

print(
    f"Post-event time: "
    f"{times[i_post]:.6f} s"
)

print()

print(
    "=== MESH ==="
)

print(
    f"Pre-event cells : {pre.n_cells}"
)

print(
    f"Event cells     : {event.n_cells}"
)

print(
    f"Post-event cells: {post.n_cells}"
)

print()

print(
    "=== DEFORMATION ==="
)

print(
    "Pre-event max |u| [mm]: "
    f"{u_pre*1000:.10f}"
)

print(
    "At-event max |u| [mm]: "
    f"{u_event*1000:.10f}"
)

print(
    "Post-event max |u| [mm]: "
    f"{u_post*1000:.10f}"
)

print(
    "Change post-pre [mm]: "
    f"{(u_post-u_pre)*1000:+.10f}"
)

print()

print(
    "=== PLASTICITY ==="
)

print(
    "Pre epsp_max: "
    f"{epsp_pre:.10e}"
)

print(
    "Post epsp_max: "
    f"{epsp_post:.10e}"
)

print(
    "Post fraction epsp > 1e-6: "
    f"{100*f6_post:.6f}%"
)

print(
    "Post fraction epsp > 1e-4: "
    f"{100*f4_post:.6f}%"
)

print()

print(
    "=== HYDRAULIC CHECK ==="
)

print(
    "Pre pressure [kPa]: "
    f"{np.nanmin(pressure_pre)/1000:.6f} "
    f"to "
    f"{np.nanmax(pressure_pre)/1000:.6f}"
)

print(
    "Post pressure [kPa]: "
    f"{np.nanmin(pressure_post)/1000:.6f} "
    f"to "
    f"{np.nanmax(pressure_post)/1000:.6f}"
)


# ============================================================
# CLASSIFY
# ============================================================

finite_ok = all(
    [
        np.isfinite(
            u_pre
        ),
        np.isfinite(
            u_post
        ),
        np.isfinite(
            epsp_post
        ),
    ]
)


if not finite_ok:

    status = (
        "FAIL: NON-FINITE RESPONSE"
    )

elif epsp_post <= 1e-10:

    status = (
        "PASS — ELASTIC EROSION RESPONSE"
    )

elif f4_post <= 0.01:

    status = (
        "PASS — LOCAL YIELD ONSET"
    )

elif f4_post <= 0.20:

    status = (
        "PASS — LOCALIZED PLASTIC RESPONSE"
    )

else:

    status = (
        "PASS — WIDESPREAD PLASTIC RESPONSE"
    )


print()
print(
    "EROSION RESPONSE STATUS:"
)

print(
    status
)


summary = (
    RUN
    / "phase05b0_summary.txt"
)

summary.write_text(
    "\n".join(
        [
            (
                "PHASE 05B-0 SINGLE "
                "TOE-EROSION TEST"
            ),
            "",
            "Antecedent state: REFERENCE",
            "Nominal E [m]: 0.4",
            "Removed MaterialID: 1",
            "",
            (
                "Pre-event cells: "
                f"{pre.n_cells}"
            ),
            (
                "Post-event cells: "
                f"{post.n_cells}"
            ),
            "",
            (
                "Pre max displacement [mm]: "
                f"{u_pre*1000:.10f}"
            ),
            (
                "Post max displacement [mm]: "
                f"{u_post*1000:.10f}"
            ),
            (
                "Delta max displacement [mm]: "
                f"{(u_post-u_pre)*1000:+.10f}"
            ),
            "",
            (
                "Pre epsp_max: "
                f"{epsp_pre:.12e}"
            ),
            (
                "Post epsp_max: "
                f"{epsp_post:.12e}"
            ),
            (
                "Post frac epsp > 1e-6 [%]: "
                f"{100*f6_post:.8f}"
            ),
            (
                "Post frac epsp > 1e-4 [%]: "
                f"{100*f4_post:.8f}"
            ),
            "",
            (
                "STATUS: "
                f"{status}"
            ),
            "",
            (
                "NOTE: no hydraulic boundary "
                "condition was imposed on the "
                "new erosion surface."
            ),
        ]
    )
    + "\n",
    encoding="utf-8",
)

print()
print(
    "PASS:",
    summary,
)
