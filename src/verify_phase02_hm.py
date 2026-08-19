from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot


ROOT = Path.cwd()

LIA_DIR = ROOT / "results" / "phase02_liakopoulos_hm"
POR_DIR = ROOT / "results" / "phase02_deformation_porosity"
OUT = ROOT / "results" / "phase02_verification"

OUT.mkdir(parents=True, exist_ok=True)


# ============================================================
# Helpers
# ============================================================

def find_pvd(folder: Path) -> Path:
    files = sorted(folder.glob("*.pvd"))

    if not files:
        raise FileNotFoundError(
            f"No PVD file found in {folder}"
        )

    print(f"PVD selected: {files[0]}")
    return files[0]


def array_names(mesh):
    return {
        "point": list(mesh.point_data.keys()),
        "cell": list(mesh.cell_data.keys()),
    }


def get_array(mesh, name):
    if name in mesh.point_data:
        return np.asarray(mesh.point_data[name]), "point"

    if name in mesh.cell_data:
        return np.asarray(mesh.cell_data[name]), "cell"

    raise KeyError(
        f"{name!r} not found.\n"
        f"Point arrays: {list(mesh.point_data.keys())}\n"
        f"Cell arrays: {list(mesh.cell_data.keys())}"
    )


def dominant_axis(mesh):
    bounds = np.array(mesh.bounds).reshape(3, 2)
    spans = bounds[:, 1] - bounds[:, 0]
    axis = int(np.argmax(spans))

    names = ["x", "y", "z"]

    print(
        "Dominant geometry axis:",
        names[axis],
        "span =",
        spans[axis],
    )

    return axis


def line_endpoints(mesh, axis):
    bounds = np.array(mesh.bounds).reshape(3, 2)
    centre = bounds.mean(axis=1)

    p0 = centre.copy()
    p1 = centre.copy()

    p0[axis] = bounds[axis, 0]
    p1[axis] = bounds[axis, 1]

    return p0, p1


def selected_indices(n, count=5):
    if n <= count:
        return np.arange(n)

    return np.unique(
        np.linspace(
            0,
            n - 1,
            count,
            dtype=int,
        )
    )


def sample_profile(mesh, array_name, axis):
    p0, p1 = line_endpoints(mesh, axis)

    sampled = mesh.sample_over_line(
        p0,
        p1,
        resolution=250,
    )

    if array_name not in sampled.point_data:
        raise KeyError(
            f"{array_name!r} unavailable after line sampling.\n"
            f"Available: {list(sampled.point_data.keys())}"
        )

    coord = sampled.points[:, axis]
    values = np.asarray(
        sampled.point_data[array_name]
    )

    order = np.argsort(coord)

    return coord[order], values[order]


# ============================================================
# 1. LIAKOPOULOS COUPLED HM
# ============================================================

print()
print("========================================")
print("LIAKOPOULOS COUPLED HM")
print("========================================")

lia_pvd = find_pvd(LIA_DIR)
lia = ot.MeshSeries(str(lia_pvd))

lia_times = np.asarray(lia.timevalues)

print("Timesteps:", len(lia_times))
print("Time range:", lia_times[0], "to", lia_times[-1])

lia_mesh0 = lia.mesh(0)

print("Arrays:")
print(array_names(lia_mesh0))

axis = dominant_axis(lia_mesh0)
axis_name = ["x", "y", "z"][axis]

if "saturation" not in lia_mesh0.point_data:
    raise SystemExit(
        "FAIL: saturation is not point data in Liakopoulos result."
    )

if "displacement" not in lia_mesh0.point_data:
    raise SystemExit(
        "FAIL: displacement is not point data in Liakopoulos result."
    )

ids = selected_indices(len(lia_times), count=5)


# ------------------------------------------------------------
# Saturation profiles
# ------------------------------------------------------------

fig, ax = plt.subplots(figsize=(7.0, 5.0))

sat_global_min = np.inf
sat_global_max = -np.inf

for idx in ids:
    mesh = lia.mesh(int(idx))

    coord, sat = sample_profile(
        mesh,
        "saturation",
        axis,
    )

    sat = np.asarray(sat).squeeze()

    sat_global_min = min(
        sat_global_min,
        float(np.nanmin(sat)),
    )

    sat_global_max = max(
        sat_global_max,
        float(np.nanmax(sat)),
    )

    ax.plot(
        sat,
        coord,
        linewidth=2,
        label=f"t = {lia_times[idx]:.3g} s",
    )

ax.set_xlabel("Saturation, Sr [-]")
ax.set_ylabel(f"{axis_name}-coordinate [m]")
ax.set_title(
    "Liakopoulos benchmark — saturation profile"
)
ax.grid(alpha=0.25)
ax.legend(fontsize=8)
fig.tight_layout()

sat_png = OUT / "liakopoulos_saturation_profiles.png"
fig.savefig(sat_png, dpi=220)
plt.close(fig)

print("PASS:", sat_png)


# ------------------------------------------------------------
# Vertical / dominant-axis displacement profiles
# ------------------------------------------------------------

fig, ax = plt.subplots(figsize=(7.0, 5.0))

disp_global_min = np.inf
disp_global_max = -np.inf

for idx in ids:
    mesh = lia.mesh(int(idx))

    coord, disp = sample_profile(
        mesh,
        "displacement",
        axis,
    )

    disp = np.asarray(disp)

    if disp.ndim == 1:
        component = disp
    else:
        component = disp[:, axis]

    disp_global_min = min(
        disp_global_min,
        float(np.nanmin(component)),
    )

    disp_global_max = max(
        disp_global_max,
        float(np.nanmax(component)),
    )

    ax.plot(
        component * 1000.0,
        coord,
        linewidth=2,
        label=f"t = {lia_times[idx]:.3g} s",
    )

ax.set_xlabel(
    f"Displacement u_{axis_name} [mm]"
)
ax.set_ylabel(f"{axis_name}-coordinate [m]")
ax.set_title(
    "Liakopoulos benchmark — displacement profile"
)
ax.grid(alpha=0.25)
ax.legend(fontsize=8)
fig.tight_layout()

disp_png = OUT / "liakopoulos_displacement_profiles.png"
fig.savefig(disp_png, dpi=220)
plt.close(fig)

print("PASS:", disp_png)


# ============================================================
# 2. DEFORMATION-DEPENDENT POROSITY
# ============================================================

print()
print("========================================")
print("DEFORMATION-DEPENDENT POROSITY")
print("========================================")

por_pvd = find_pvd(POR_DIR)
por = ot.MeshSeries(str(por_pvd))

por_times = np.asarray(por.timevalues)

print("Timesteps:", len(por_times))
print("Time range:", por_times[0], "to", por_times[-1])

por_mesh0 = por.mesh(0)

print("Arrays:")
print(array_names(por_mesh0))


# ------------------------------------------------------------
# Porosity evolution
# ------------------------------------------------------------

mean_porosity = []
min_porosity = []
max_porosity = []

max_displacement = []

for i in range(len(por_times)):
    mesh = por.mesh(i)

    phi, phi_location = get_array(
        mesh,
        "porosity",
    )

    phi = np.asarray(phi).squeeze()

    mean_porosity.append(
        float(np.nanmean(phi))
    )
    min_porosity.append(
        float(np.nanmin(phi))
    )
    max_porosity.append(
        float(np.nanmax(phi))
    )

    disp, _ = get_array(
        mesh,
        "displacement",
    )

    disp = np.asarray(disp)

    if disp.ndim == 1:
        mag = np.abs(disp)
    else:
        mag = np.linalg.norm(
            disp,
            axis=1,
        )

    max_displacement.append(
        float(np.nanmax(mag))
    )

print("Porosity location:", phi_location)

mean_porosity = np.asarray(mean_porosity)
min_porosity = np.asarray(min_porosity)
max_porosity = np.asarray(max_porosity)

max_displacement = np.asarray(
    max_displacement
)


fig, ax = plt.subplots(figsize=(7.0, 4.8))

ax.plot(
    por_times,
    mean_porosity,
    marker="o",
    label="spatial mean",
)

ax.fill_between(
    por_times,
    min_porosity,
    max_porosity,
    alpha=0.2,
    label="spatial min–max",
)

ax.set_xlabel("Time [s]")
ax.set_ylabel("Porosity, n [-]")
ax.set_title(
    "Deformation-dependent porosity benchmark"
)
ax.grid(alpha=0.25)
ax.legend()

fig.tight_layout()

por_png = OUT / "porosity_evolution.png"
fig.savefig(por_png, dpi=220)
plt.close(fig)

print("PASS:", por_png)


fig, ax = plt.subplots(figsize=(7.0, 4.8))

ax.plot(
    por_times,
    max_displacement * 1000.0,
    marker="o",
)

ax.set_xlabel("Time [s]")
ax.set_ylabel("Maximum displacement [mm]")
ax.set_title(
    "Mechanical response during porosity benchmark"
)
ax.grid(alpha=0.25)

fig.tight_layout()

por_disp_png = OUT / "porosity_benchmark_displacement.png"
fig.savefig(
    por_disp_png,
    dpi=220,
)
plt.close(fig)

print("PASS:", por_disp_png)


# ============================================================
# 3. NUMERICAL SANITY CHECKS
# ============================================================

print()
print("========================================")
print("NUMERICAL SANITY CHECKS")
print("========================================")

print(
    "Liakopoulos saturation range:",
    sat_global_min,
    "to",
    sat_global_max,
)

print(
    "Liakopoulos displacement range [m]:",
    disp_global_min,
    "to",
    disp_global_max,
)

print(
    "Porosity initial mean:",
    mean_porosity[0],
)

print(
    "Porosity final mean:",
    mean_porosity[-1],
)

print(
    "Porosity mean change:",
    mean_porosity[-1] - mean_porosity[0],
)

print(
    "Maximum displacement reached [m]:",
    np.nanmax(max_displacement),
)


if sat_global_min < -1e-8:
    raise SystemExit(
        "FAIL: negative saturation detected"
    )

if sat_global_max > 1.0001:
    raise SystemExit(
        "FAIL: saturation above physical upper bound"
    )

if not np.all(
    np.isfinite(mean_porosity)
):
    raise SystemExit(
        "FAIL: non-finite porosity detected"
    )

if np.ptp(mean_porosity) == 0:
    raise SystemExit(
        "FAIL: porosity did not evolve"
    )


# ============================================================
# 4. SUMMARY FILE
# ============================================================

summary = OUT / "verification_summary.txt"

summary.write_text(
    "\n".join(
        [
            "PHASE 02 COUPLED HM VERIFICATION",
            "",
            f"Liakopoulos PVD: {lia_pvd}",
            f"Liakopoulos timesteps: {len(lia_times)}",
            (
                "Saturation range: "
                f"{sat_global_min:.8g} "
                f"to {sat_global_max:.8g}"
            ),
            (
                "Displacement range [m]: "
                f"{disp_global_min:.8g} "
                f"to {disp_global_max:.8g}"
            ),
            "",
            f"Porosity PVD: {por_pvd}",
            f"Porosity timesteps: {len(por_times)}",
            (
                "Initial mean porosity: "
                f"{mean_porosity[0]:.10g}"
            ),
            (
                "Final mean porosity: "
                f"{mean_porosity[-1]:.10g}"
            ),
            (
                "Mean porosity change: "
                f"{mean_porosity[-1]-mean_porosity[0]:.10g}"
            ),
            (
                "Maximum displacement [m]: "
                f"{np.nanmax(max_displacement):.10g}"
            ),
            "",
            "STATUS: PASS",
        ]
    )
    + "\n",
    encoding="utf-8",
)

print("PASS:", summary)

print()
print("========================================")
print("PHASE 02 HM FIELD VERIFICATION: PASS")
print("========================================")
