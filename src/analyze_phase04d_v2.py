from pathlib import Path

import numpy as np
import ogstools as ot


ROOT = Path.cwd()

RUN_ROOT = (
    ROOT
    / "results"
    / "phase04d_v2"
)


def load_case(
    case,
):

    folder = (
        RUN_ROOT
        / case
    )

    pvds = sorted(
        folder.glob("*.pvd")
    )

    if not pvds:
        return None

    series = ot.MeshSeries(
        str(pvds[0])
    )

    times = np.asarray(
        series.timevalues,
        dtype=float,
    )

    mesh = series.mesh(
        len(times) - 1
    )

    return (
        series,
        times,
        mesh,
    )


def get_field(
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
        f"{name} missing"
    )


print(
    "========================================"
)

print(
    "PHASE 04D-V2 INITIALIZED MC-HM"
)

print(
    "========================================"
)


for case in [
    "strong",
    "screening",
]:

    result = load_case(
        case
    )

    print()
    print(
        f"=== {case.upper()} ==="
    )

    if result is None:

        print(
            "NO VALID PVD"
        )

        continue

    series, times, mesh = result

    u = get_field(
        mesh,
        "displacement",
    )

    epsp = get_field(
        mesh,
        "EquivalentPlasticStrain",
    ).squeeze()

    pressure = get_field(
        mesh,
        "pressure",
    ).squeeze()

    saturation = get_field(
        mesh,
        "saturation",
    ).squeeze()

    umag = np.linalg.norm(
        u,
        axis=1,
    )

    umax = float(
        np.nanmax(
            umag
        )
    )

    epsp_max = float(
        np.nanmax(
            epsp
        )
    )

    epsp_fraction_1e4 = float(
        np.mean(
            epsp > 1e-4
        )
    )

    print(
        "Final time [days]: "
        f"{times[-1]/86400:.8f}"
    )

    print(
        "Max incremental displacement [mm]: "
        f"{umax*1000:.10f}"
    )

    print(
        "Max EquivalentPlasticStrain: "
        f"{epsp_max:.10e}"
    )

    print(
        "Fraction eps_p > 1e-4: "
        f"{100*epsp_fraction_1e4:.6f}%"
    )

    print(
        "Pressure range [kPa]: "
        f"{np.nanmin(pressure)/1000:.6f} "
        f"to "
        f"{np.nanmax(pressure)/1000:.6f}"
    )

    print(
        "Saturation range [-]: "
        f"{np.nanmin(saturation):.8f} "
        f"to "
        f"{np.nanmax(saturation):.8f}"
    )


print()
print(
    "PHASE 04D-V2 ANALYSIS COMPLETE"
)
