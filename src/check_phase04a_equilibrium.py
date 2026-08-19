from pathlib import Path

import numpy as np
import ogstools as ot


ROOT = Path.cwd()
RUN = ROOT / "results" / "phase04a_gravity"

pvds = sorted(RUN.glob("*.pvd"))

if not pvds:
    raise SystemExit("FAIL: no Phase 04A PVD found")

series = ot.MeshSeries(str(pvds[0]))

times_s = np.asarray(
    series.timevalues,
    dtype=float,
)

times_d = times_s / 86400.0


def displacement(mesh):
    if "displacement" not in mesh.point_data:
        raise KeyError(
            "displacement not found in point data"
        )

    return np.asarray(
        mesh.point_data["displacement"],
        dtype=float,
    )


def max_displacement(mesh):
    u = displacement(mesh)

    mag = np.linalg.norm(
        u,
        axis=1,
    )

    idx = int(
        np.nanargmax(mag)
    )

    return (
        float(mag[idx]),
        mesh.points[idx],
        u[idx],
    )


print("========================================")
print("PHASE 04A EQUILIBRIUM CHECK")
print("========================================")

print()
print("Final output times [days]:")

for t in times_d[-8:]:
    print(f"{t:.6f}")


# ------------------------------------------------------------
# Maximum displacement at final state
# ------------------------------------------------------------

final_mesh = series.mesh(
    len(times_d) - 1
)

umax, xyz, uvec = max_displacement(
    final_mesh
)

print()
print("=== FINAL MAXIMUM DISPLACEMENT ===")

print(
    f"|u|max = {umax * 1000:.8f} mm"
)

print(
    "location [m] = "
    f"({xyz[0]:.6f}, {xyz[1]:.6f})"
)

print(
    "components [mm] = "
    f"({uvec[0]*1000:.8f}, "
    f"{uvec[1]*1000:.8f})"
)


# ------------------------------------------------------------
# Incremental displacement between late outputs
# ------------------------------------------------------------

print()
print("=== LATE-TIME INCREMENTS ===")

records = []

start = max(
    1,
    len(times_d) - 8,
)

for i in range(
    start,
    len(times_d),
):

    mesh_a = series.mesh(i - 1)
    mesh_b = series.mesh(i)

    ua = displacement(mesh_a)
    ub = displacement(mesh_b)

    du = ub - ua

    du_mag = np.linalg.norm(
        du,
        axis=1,
    )

    max_du = float(
        np.nanmax(du_mag)
    )

    rms_du = float(
        np.sqrt(
            np.nanmean(
                du_mag ** 2
            )
        )
    )

    dt_days = (
        times_d[i]
        - times_d[i - 1]
    )

    records.append(
        (
            times_d[i],
            dt_days,
            max_du,
            rms_du,
        )
    )

    print(
        f"t={times_d[i]:.3f} d | "
        f"dt={dt_days:.3f} d | "
        f"max Δu={max_du*1000:.8f} mm | "
        f"RMS Δu={rms_du*1000:.8f} mm"
    )


# ------------------------------------------------------------
# Compare final increment with total displacement
# ------------------------------------------------------------

last_max_du = records[-1][2]

ratio = (
    last_max_du / umax
    if umax > 0
    else np.nan
)

print()
print("=== EQUILIBRIUM METRIC ===")

print(
    "Final-step max Δu / total max |u| = "
    f"{ratio:.10e}"
)

# Screening criterion only:
# final daily increment < 0.1% of accumulated displacement.
threshold = 1.0e-3

if np.isfinite(ratio) and ratio < threshold:

    result = "EQUILIBRIUM SCREENING: PASS"

    interpretation = (
        "Late-time displacement increment is small "
        "relative to accumulated gravity deformation."
    )

else:

    result = "EQUILIBRIUM SCREENING: REVIEW"

    interpretation = (
        "Gravity state is still evolving materially; "
        "extend equilibration before seasonal branching."
    )


print()
print(result)
print(interpretation)


summary = RUN / "phase04a_equilibrium_summary.txt"

summary.write_text(
    "\n".join(
        [
            "PHASE 04A EQUILIBRIUM CHECK",
            "",
            f"Final time [days]: {times_d[-1]:.8f}",
            f"Total max displacement [mm]: {umax*1000:.10f}",
            (
                "Final max displacement location [m]: "
                f"{xyz[0]:.8f}, {xyz[1]:.8f}"
            ),
            (
                "Final displacement components [mm]: "
                f"{uvec[0]*1000:.10f}, "
                f"{uvec[1]*1000:.10f}"
            ),
            (
                "Final-step max incremental displacement [mm]: "
                f"{last_max_du*1000:.12f}"
            ),
            (
                "Final-step increment / total displacement: "
                f"{ratio:.12e}"
            ),
            "",
            result,
            interpretation,
        ]
    )
    + "\n",
    encoding="utf-8",
)

print()
print("PASS:", summary)
