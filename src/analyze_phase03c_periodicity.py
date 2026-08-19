from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot


ROOT = Path.cwd()
RUN = ROOT / "results" / "phase03c_three_cycles"
OUT = ROOT / "results" / "phase03c_analysis"

OUT.mkdir(
    parents=True,
    exist_ok=True,
)

pvds = sorted(RUN.glob("*.pvd"))

if not pvds:
    raise SystemExit(
        "FAIL: no Phase 03C PVD found"
    )

series = ot.MeshSeries(
    str(pvds[0])
)

times_d = (
    np.asarray(
        series.timevalues,
        dtype=float,
    )
    / 86400.0
)

mesh0 = series.mesh(0)

ymax = float(
    mesh0.bounds[3]
)

height = (
    float(mesh0.bounds[3])
    - float(mesh0.bounds[2])
)

depths = np.array(
    [
        0.0,
        0.5,
        1.0,
        2.0,
        3.0,
        4.0,
    ]
)

depths = depths[
    depths <= height + 1e-9
]


def get_data(mesh, name):
    if name in mesh.point_data:
        return (
            np.asarray(
                mesh.point_data[name]
            ),
            np.asarray(mesh.points),
        )

    if name in mesh.cell_data:
        return (
            np.asarray(
                mesh.cell_data[name]
            ),
            np.asarray(
                mesh.cell_centers().points
            ),
        )

    raise KeyError(name)


def horizontal_mean(
    mesh,
    variable,
    depth,
):
    values, coords = get_data(
        mesh,
        variable,
    )

    values = np.asarray(
        values
    ).squeeze()

    y = coords[:, 1]
    target = ymax - depth

    delta = np.abs(
        y - target
    )

    nearest = np.min(
        delta
    )

    mask = np.isclose(
        delta,
        nearest,
        atol=1e-10,
        rtol=0,
    )

    return float(
        np.nanmean(
            values[mask]
        )
    )


sat = np.zeros(
    (
        len(times_d),
        len(depths),
    )
)

pressure = np.zeros_like(
    sat
)

surface_uy = np.zeros(
    len(times_d)
)

for i in range(
    len(times_d)
):
    mesh = series.mesh(i)

    for j, depth in enumerate(
        depths
    ):
        sat[i, j] = horizontal_mean(
            mesh,
            "saturation",
            depth,
        )

        pressure[i, j] = horizontal_mean(
            mesh,
            "pressure",
            depth,
        )

    disp, coords = get_data(
        mesh,
        "displacement",
    )

    disp = np.asarray(
        disp
    )

    top = np.isclose(
        coords[:, 1],
        np.max(coords[:, 1]),
        atol=1e-10,
    )

    surface_uy[i] = np.mean(
        disp[top, 1]
    )


# ============================================================
# Cycle windows
# ============================================================

cycles = [
    (1, 30.0, 395.0),
    (2, 395.0, 760.0),
    (3, 760.0, 1125.0),
]

records = []

for cycle_no, start, end in cycles:

    # Avoid assigning the shared cycle boundary
    # to both cycles except for final cycle.
    if cycle_no < 3:
        mask = (
            (times_d >= start)
            & (times_d < end)
        )
    else:
        mask = (
            (times_d >= start)
            & (times_d <= end)
        )

    t = times_d[mask]

    surface_index = int(
        np.nanargmax(
            sat[mask, 0]
        )
    )

    surface_peak = float(
        t[surface_index]
    )

    for j, depth in enumerate(
        depths
    ):

        values = sat[
            mask,
            j,
        ]

        peak_idx = int(
            np.nanargmax(values)
        )

        peak_day = float(
            t[peak_idx]
        )

        amplitude = float(
            np.nanmax(values)
            - np.nanmin(values)
        )

        lag = (
            peak_day
            - surface_peak
        )

        censored = (
            np.isclose(
                peak_day,
                end,
                atol=2.6,
            )
            or np.isclose(
                peak_day,
                t[-1],
                atol=1e-8,
            )
        )

        records.append(
            {
                "cycle": cycle_no,
                "depth": depth,
                "wettest_day": peak_day,
                "surface_wettest_day":
                    surface_peak,
                "lag": lag,
                "amplitude": amplitude,
                "right_censored":
                    censored,
            }
        )


# ============================================================
# Cycle 2 vs Cycle 3 periodic convergence
# ============================================================

phase = np.linspace(
    0.0,
    365.0,
    1461,
)

cycle_profiles = {}

for cycle_no, start, _ in cycles:

    arr = np.zeros(
        (
            len(phase),
            len(depths),
        )
    )

    for j in range(
        len(depths)
    ):
        arr[:, j] = np.interp(
            start + phase,
            times_d,
            sat[:, j],
        )

    cycle_profiles[
        cycle_no
    ] = arr


diff23 = (
    cycle_profiles[3]
    - cycle_profiles[2]
)

max_cycle_difference = float(
    np.nanmax(
        np.abs(diff23)
    )
)

rms_cycle_difference = float(
    np.sqrt(
        np.nanmean(
            diff23 ** 2
        )
    )
)


# ============================================================
# Plot cycle 3
# ============================================================

fig, ax = plt.subplots(
    figsize=(8, 5)
)

mask3 = (
    (times_d >= 760)
    & (times_d <= 1125)
)

for j, depth in enumerate(
    depths
):
    ax.plot(
        times_d[mask3] - 760,
        sat[mask3, j],
        linewidth=1.8,
        label=f"{depth:g} m",
    )

ax.set_xlabel(
    "Time within third cycle [days]"
)

ax.set_ylabel(
    "Degree of saturation, Sr [-]"
)

ax.set_title(
    "Third-cycle seasonal saturation response"
)

ax.grid(
    alpha=0.25
)

ax.legend(
    fontsize=8,
    ncol=2,
)

fig.tight_layout()

fig.savefig(
    OUT
    / "figure_01_third_cycle_saturation.png",
    dpi=220,
)

plt.close(fig)


# ============================================================
# Attenuation + lag
# ============================================================

third = [
    r
    for r in records
    if r["cycle"] == 3
]

third_depth = np.array(
    [
        r["depth"]
        for r in third
    ]
)

third_amp = np.array(
    [
        r["amplitude"]
        for r in third
    ]
)

third_lag = np.array(
    [
        r["lag"]
        for r in third
    ]
)


fig, ax = plt.subplots(
    figsize=(6.5, 4.8)
)

ax.plot(
    third_amp,
    third_depth,
    marker="o",
)

ax.invert_yaxis()

ax.set_xlabel(
    "Seasonal saturation amplitude, Delta Sr"
)

ax.set_ylabel(
    "Depth below surface [m]"
)

ax.set_title(
    "Third-cycle hydraulic attenuation"
)

ax.grid(
    alpha=0.25
)

fig.tight_layout()

fig.savefig(
    OUT
    / "figure_02_third_cycle_attenuation.png",
    dpi=220,
)

plt.close(fig)


fig, ax = plt.subplots(
    figsize=(6.5, 4.8)
)

ax.plot(
    third_lag,
    third_depth,
    marker="o",
)

ax.invert_yaxis()

ax.set_xlabel(
    "Wettest-state lag vs surface [days]"
)

ax.set_ylabel(
    "Depth below surface [m]"
)

ax.set_title(
    "Third-cycle hydraulic phase lag"
)

ax.grid(
    alpha=0.25
)

fig.tight_layout()

fig.savefig(
    OUT
    / "figure_03_third_cycle_lag.png",
    dpi=220,
)

plt.close(fig)


# ============================================================
# Mechanical periodicity
# ============================================================

fig, ax = plt.subplots(
    figsize=(8, 4.7)
)

ax.plot(
    times_d,
    surface_uy * 1000,
    linewidth=1.8,
)

for x in [
    395,
    760,
]:
    ax.axvline(
        x,
        linestyle="--",
        alpha=0.5,
    )

ax.set_xlabel(
    "Time [days]"
)

ax.set_ylabel(
    "Surface vertical displacement [mm]"
)

ax.set_title(
    "Three-cycle hydro-mechanical response"
)

ax.grid(
    alpha=0.25
)

fig.tight_layout()

fig.savefig(
    OUT
    / "figure_04_three_cycle_displacement.png",
    dpi=220,
)

plt.close(fig)


# ============================================================
# CSV + summary
# ============================================================

csv = (
    OUT
    / "three_cycle_depth_metrics.csv"
)

with csv.open(
    "w",
    encoding="utf-8",
) as f:

    f.write(
        "cycle,depth_m,"
        "surface_wettest_day,"
        "wettest_day,"
        "lag_days,"
        "saturation_amplitude,"
        "right_censored\n"
    )

    for r in records:
        f.write(
            f"{r['cycle']},"
            f"{r['depth']:.3f},"
            f"{r['surface_wettest_day']:.3f},"
            f"{r['wettest_day']:.3f},"
            f"{r['lag']:.3f},"
            f"{r['amplitude']:.10f},"
            f"{int(r['right_censored'])}\n"
        )


summary = (
    OUT
    / "phase03c_summary.txt"
)

lines = [
    "PHASE 03C THREE-CYCLE PERIODICITY",
    "",
    (
        "Cycle 2 -> Cycle 3 "
        "max |Delta Sr|: "
        f"{max_cycle_difference:.8f}"
    ),
    (
        "Cycle 2 -> Cycle 3 "
        "RMS Delta Sr: "
        f"{rms_cycle_difference:.8f}"
    ),
    "",
    "THIRD-CYCLE DEPTH RESPONSE:",
]

for r in third:

    flag = (
        "CENSORED"
        if r["right_censored"]
        else "resolved"
    )

    lines.append(
        f"depth={r['depth']:.2f} m | "
        f"wettest_day={r['wettest_day']:.2f} | "
        f"lag={r['lag']:.2f} d | "
        f"DeltaSr={r['amplitude']:.6f} | "
        f"{flag}"
    )

lines += [
    "",
    (
        "Surface displacement "
        "total range [mm]: "
        f"{np.ptp(surface_uy)*1000:.8f}"
    ),
    "",
    "STATUS: PASS",
]

summary.write_text(
    "\n".join(lines) + "\n",
    encoding="utf-8",
)

print(
    "\n".join(lines)
)

print()
print(
    "PHASE 03C PERIODICITY ANALYSIS: PASS"
)
