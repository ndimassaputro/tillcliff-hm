from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot


ROOT = Path.cwd()
PVD_DIR = ROOT / "results" / "phase03_seasonal_column"
OUT = ROOT / "results" / "phase03_analysis"
OUT.mkdir(parents=True, exist_ok=True)

pvd_files = sorted(PVD_DIR.glob("*.pvd"))
if not pvd_files:
    raise SystemExit("FAIL: no Phase 03 PVD found")

pvd = pvd_files[0]
series = ot.MeshSeries(str(pvd))

times_s = np.asarray(series.timevalues, dtype=float)
times_d = times_s / 86400.0

print("=== PHASE 03B SEASONAL RESPONSE ANALYSIS ===")
print("PVD:", pvd)
print("Timesteps:", len(times_d))
print(
    f"Time range: {times_d.min():.2f}"
    f" to {times_d.max():.2f} days"
)


def get_data(mesh, name):
    if name in mesh.point_data:
        return (
            np.asarray(mesh.point_data[name]),
            np.asarray(mesh.points),
            "point",
        )

    if name in mesh.cell_data:
        centers = mesh.cell_centers().points
        return (
            np.asarray(mesh.cell_data[name]),
            np.asarray(centers),
            "cell",
        )

    raise KeyError(
        f"{name!r} not found.\n"
        f"Point data: {list(mesh.point_data.keys())}\n"
        f"Cell data: {list(mesh.cell_data.keys())}"
    )


mesh0 = series.mesh(0)

print()
print("=== AVAILABLE ARRAYS ===")
print("Point:", list(mesh0.point_data.keys()))
print("Cell :", list(mesh0.cell_data.keys()))

ymin = float(mesh0.bounds[2])
ymax = float(mesh0.bounds[3])
height = ymax - ymin

print()
print("=== GEOMETRY ===")
print(f"Column ymin   = {ymin:.6g} m")
print(f"Column ymax   = {ymax:.6g} m")
print(f"Column height = {height:.6g} m")

sensor_depths = np.array(
    [0.0, 0.5, 1.0, 2.0, 3.0, 4.0],
    dtype=float,
)

sensor_depths = sensor_depths[
    sensor_depths <= height + 1e-9
]


def horizontal_mean(mesh, variable, depth):
    data, coords, _ = get_data(mesh, variable)

    data = np.asarray(data).squeeze()
    y = coords[:, 1]

    target_y = ymax - depth

    dy = np.abs(y - target_y)
    nearest = np.min(dy)

    mask = np.isclose(
        dy,
        nearest,
        atol=1e-10,
        rtol=0,
    )

    values = data[mask]

    if values.ndim > 1:
        return np.nanmean(values, axis=0)

    return float(np.nanmean(values))


pressure = np.zeros(
    (len(times_d), len(sensor_depths))
)

saturation = np.zeros_like(pressure)

surface_uy = np.zeros(len(times_d))
mean_porosity = np.zeros(len(times_d))

for i in range(len(times_d)):
    mesh = series.mesh(i)

    for j, depth in enumerate(sensor_depths):
        pressure[i, j] = horizontal_mean(
            mesh,
            "pressure",
            depth,
        )

        saturation[i, j] = horizontal_mean(
            mesh,
            "saturation",
            depth,
        )

    disp, disp_coords, _ = get_data(
        mesh,
        "displacement",
    )

    disp = np.asarray(disp)
    y = disp_coords[:, 1]

    top_mask = np.isclose(
        y,
        np.max(y),
        atol=1e-10,
        rtol=0,
    )

    if disp.ndim == 2:
        surface_uy[i] = np.mean(
            disp[top_mask, 1]
        )
    else:
        surface_uy[i] = np.mean(
            disp[top_mask]
        )

    porosity, _, _ = get_data(
        mesh,
        "porosity",
    )

    mean_porosity[i] = np.nanmean(
        np.asarray(porosity)
    )


# ============================================================
# 1. PRESSURE / SUCTION TIME SERIES
# ============================================================

fig, ax = plt.subplots(figsize=(8.0, 5.0))

for j, depth in enumerate(sensor_depths):
    suction_kpa = -pressure[:, j] / 1000.0

    ax.plot(
        times_d,
        suction_kpa,
        linewidth=1.8,
        label=f"{depth:g} m depth",
    )

ax.set_xlabel("Time [days]")
ax.set_ylabel("Matric suction [kPa]")
ax.set_title(
    "Seasonal suction propagation through clayey till"
)
ax.grid(alpha=0.25)
ax.legend(fontsize=8, ncol=2)

fig.tight_layout()
fig.savefig(
    OUT / "figure_01_suction_timeseries.png",
    dpi=220,
)
plt.close(fig)


# ============================================================
# 2. SATURATION TIME SERIES
# ============================================================

fig, ax = plt.subplots(figsize=(8.0, 5.0))

for j, depth in enumerate(sensor_depths):
    ax.plot(
        times_d,
        saturation[:, j],
        linewidth=1.8,
        label=f"{depth:g} m depth",
    )

ax.set_xlabel("Time [days]")
ax.set_ylabel("Degree of saturation, Sr [-]")
ax.set_title(
    "Seasonal saturation response and attenuation"
)
ax.grid(alpha=0.25)
ax.legend(fontsize=8, ncol=2)

fig.tight_layout()
fig.savefig(
    OUT / "figure_02_saturation_timeseries.png",
    dpi=220,
)
plt.close(fig)


# ============================================================
# 3. DEPTH PROFILES AT KEY TIMES
# ============================================================

profile_days = np.array(
    [30, 90, 180, 270, 360, 395],
    dtype=float,
)

fig, ax = plt.subplots(figsize=(7.0, 5.3))

for target_day in profile_days:
    idx = int(
        np.argmin(np.abs(times_d - target_day))
    )

    mesh = series.mesh(idx)

    sat, coords, _ = get_data(
        mesh,
        "saturation",
    )

    sat = np.asarray(sat).squeeze()
    depth = ymax - coords[:, 1]

    order = np.argsort(depth)

    # Collapse duplicate horizontal locations.
    rounded_depth = np.round(
        depth[order],
        decimals=8,
    )

    unique_depths = np.unique(rounded_depth)

    profile_sat = []
    profile_z = []

    for d in unique_depths:
        mask = np.isclose(
            rounded_depth,
            d,
            atol=1e-10,
        )

        profile_z.append(d)
        profile_sat.append(
            np.nanmean(sat[order][mask])
        )

    ax.plot(
        profile_sat,
        profile_z,
        linewidth=1.8,
        label=f"{times_d[idx]:.0f} d",
    )

ax.invert_yaxis()
ax.set_xlabel("Degree of saturation, Sr [-]")
ax.set_ylabel("Depth below surface [m]")
ax.set_title(
    "Seasonal saturation profiles"
)
ax.grid(alpha=0.25)
ax.legend(fontsize=8)

fig.tight_layout()
fig.savefig(
    OUT / "figure_03_saturation_depth_profiles.png",
    dpi=220,
)
plt.close(fig)


# ============================================================
# 4. MECHANICAL RESPONSE
# ============================================================

fig, ax = plt.subplots(figsize=(8.0, 4.8))

ax.plot(
    times_d,
    surface_uy * 1000.0,
    linewidth=2,
)

ax.set_xlabel("Time [days]")
ax.set_ylabel("Surface vertical displacement [mm]")
ax.set_title(
    "Hydro-mechanical surface response"
)
ax.grid(alpha=0.25)

fig.tight_layout()
fig.savefig(
    OUT / "figure_04_surface_displacement.png",
    dpi=220,
)
plt.close(fig)


fig, ax = plt.subplots(figsize=(8.0, 4.8))

ax.plot(
    times_d,
    mean_porosity,
    linewidth=2,
)

ax.set_xlabel("Time [days]")
ax.set_ylabel("Spatial mean porosity [-]")
ax.set_title(
    "Deformation-dependent porosity evolution"
)
ax.grid(alpha=0.25)

fig.tight_layout()
fig.savefig(
    OUT / "figure_05_mean_porosity.png",
    dpi=220,
)
plt.close(fig)


# ============================================================
# 5. PHASE LAG + ATTENUATION
# ============================================================

analysis_mask = times_d >= 30.0

peak_days = []
sat_amplitude = []

for j, depth in enumerate(sensor_depths):
    s = saturation[analysis_mask, j]
    t = times_d[analysis_mask]

    max_idx = int(np.nanargmax(s))

    peak_days.append(
        float(t[max_idx])
    )

    sat_amplitude.append(
        float(np.nanmax(s) - np.nanmin(s))
    )

peak_days = np.asarray(peak_days)
sat_amplitude = np.asarray(sat_amplitude)

surface_peak_day = peak_days[0]
phase_lag = peak_days - surface_peak_day


# ============================================================
# 6. WRITE CSV
# ============================================================

csv = OUT / "seasonal_depth_metrics.csv"

with csv.open("w", encoding="utf-8") as f:
    f.write(
        "depth_m,"
        "peak_saturation_day,"
        "lag_vs_surface_days,"
        "saturation_amplitude\n"
    )

    for d, tp, lag, amp in zip(
        sensor_depths,
        peak_days,
        phase_lag,
        sat_amplitude,
    ):
        f.write(
            f"{d:.3f},"
            f"{tp:.6f},"
            f"{lag:.6f},"
            f"{amp:.10f}\n"
        )


# ============================================================
# 7. SUMMARY
# ============================================================

summary = OUT / "phase03_summary.txt"

lines = [
    "PHASE 03B SEASONAL HM RESPONSE",
    "",
    f"Column height: {height:.3f} m",
    f"Simulation duration: {times_d[-1]:.2f} days",
    "",
    "DEPTH RESPONSE:",
]

for d, tp, lag, amp in zip(
    sensor_depths,
    peak_days,
    phase_lag,
    sat_amplitude,
):
    lines.append(
        f"depth={d:.2f} m | "
        f"wettest_day={tp:.2f} | "
        f"lag={lag:.2f} d | "
        f"DeltaSr={amp:.6f}"
    )

lines += [
    "",
    (
        "Surface vertical displacement range [mm]: "
        f"{np.ptp(surface_uy) * 1000.0:.8f}"
    ),
    (
        "Initial mean porosity: "
        f"{mean_porosity[0]:.10f}"
    ),
    (
        "Minimum mean porosity: "
        f"{np.min(mean_porosity):.10f}"
    ),
    (
        "Maximum mean porosity: "
        f"{np.max(mean_porosity):.10f}"
    ),
    (
        "Net mean porosity change: "
        f"{mean_porosity[-1]-mean_porosity[0]:.10e}"
    ),
    "",
    "STATUS: PASS",
]

summary.write_text(
    "\n".join(lines) + "\n",
    encoding="utf-8",
)

print()
print("=== DEPTH RESPONSE ===")

for line in lines[6:6 + len(sensor_depths)]:
    print(line)

print()
print("=== MECHANICAL RESPONSE ===")
print(
    "Surface displacement range [mm]:",
    np.ptp(surface_uy) * 1000.0,
)
print(
    "Mean porosity range:",
    np.min(mean_porosity),
    "to",
    np.max(mean_porosity),
)

print()
print("PASS:", csv)
print("PASS:", summary)

print()
print("========================================")
print("PHASE 03B SEASONAL RESPONSE: PASS")
print("========================================")
