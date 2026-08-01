#!/usr/bin/env python3
"""Regenerate the README gallery from the repository's simulation data."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.animation as animation
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
LABS = ROOT / "labs"
QUANTUM_DATA = ROOT / "labs" / "12" / "data"
ASSETS = ROOT / "docs" / "assets"
MAX_FRAMES = 60
FPS = 15
LOG_PROBABILITY_DECADES = 6


def frame_indices(length: int, maximum: int = MAX_FRAMES) -> np.ndarray:
    """Return evenly spaced, unique frame indices."""
    return np.unique(np.linspace(0, length - 1, min(length, maximum), dtype=int))


def log_probability_scale(probabilities: np.ndarray) -> tuple[LogNorm, float]:
    """Return a fixed log scale and floor for nonnegative probability data."""
    if not np.all(np.isfinite(probabilities)):
        raise ValueError("probability data contain a non-finite value")
    if np.any(probabilities < 0.0):
        raise ValueError("probability data contain a negative value")

    maximum = float(probabilities.max())
    if maximum <= 0.0:
        raise ValueError("probability data contain no positive values")

    floor = maximum * 10.0 ** -LOG_PROBABILITY_DECADES
    return LogNorm(vmin=floor, vmax=maximum, clip=True), floor


def metadata_float(path: Path, key: str) -> float:
    """Read a ``# key value`` floating-point entry from a data file."""
    prefix = f"# {key} "
    with path.open(encoding="utf-8") as source:
        for line in source:
            if line.startswith(prefix):
                return float(line[len(prefix):].strip())
    raise ValueError(f"{path} has no {prefix.strip()!r} metadata entry")


def save_animation(
    figure: plt.Figure,
    update,
    frames,
    filename: str,
    *,
    fps: int = FPS,
) -> None:
    """Render one compact, GitHub-friendly GIF."""
    movie = animation.FuncAnimation(
        figure,
        update,
        frames=frames,
        interval=1000 / fps,
        blit=False,
        repeat=True,
        cache_frame_data=False,
    )
    output = ASSETS / filename
    movie.save(output, writer=animation.PillowWriter(fps=fps), dpi=82)
    plt.close(figure)
    print(f"generated {output.relative_to(ROOT)}")


def probability_chain(
    data_file: Path,
    title: str,
    filename: str,
    *,
    barrier_site: int,
    barrier_height: float,
    t_max: float,
) -> None:
    probabilities = np.loadtxt(data_file, ndmin=2)
    indices = frame_indices(len(probabilities))
    x = np.arange(probabilities.shape[1])
    times = np.linspace(0.0, t_max, len(probabilities))
    if not 0 <= barrier_site < len(x):
        raise ValueError(f"barrier site {barrier_site} is outside the chain")

    figure, axis = plt.subplots(figsize=(6.4, 3.8))
    axis.axvspan(
        barrier_site - 0.45,
        barrier_site + 0.45,
        color="#f97316",
        alpha=0.18,
        zorder=0,
    )
    axis.axvline(
        barrier_site,
        color="#c2410c",
        linestyle="--",
        linewidth=1.2,
        zorder=1,
    )
    axis.text(
        barrier_site,
        0.98,
        rf"barrier $\Delta E=+{barrier_height:g}$",
        transform=axis.get_xaxis_transform(),
        color="#9a3412",
        fontsize=8.5,
        ha="center",
        va="top",
        bbox={
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": 0.8,
            "pad": 2,
        },
    )
    (curve,) = axis.plot(
        x,
        probabilities[0],
        color="#4f46e5",
        linewidth=2.2,
        marker="o",
        markersize=3.2,
        zorder=3,
    )
    axis.fill_between(x, 0, probabilities[0], color="#a5b4fc", alpha=0.35)
    fill = [axis.collections[-1]]
    axis.set(
        xlim=(-0.5, len(x) - 0.5),
        ylim=(0, max(1.0, float(probabilities.max()) * 1.05)),
        xlabel="Lattice site (zero-based)",
        ylabel=r"Probability $|\psi|^2$",
        title=title,
    )
    time_text = axis.text(
        0.98,
        0.92,
        "",
        ha="right",
        va="top",
        transform=axis.transAxes,
        fontsize=8.5,
        color="#334155",
        bbox={
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": 0.8,
            "pad": 2,
        },
    )

    def update(frame: int):
        fill[0].remove()
        curve.set_ydata(probabilities[frame])
        fill[0] = axis.fill_between(
            x, 0, probabilities[frame], color="#a5b4fc", alpha=0.35
        )
        time_text.set_text(
            f"t = {times[frame]:.2f}\n"
            rf"$\sum |\psi|^2$ = {probabilities[frame].sum():.6f}"
        )
        return curve, fill[0], time_text

    figure.tight_layout()
    save_animation(figure, update, indices, filename)


def square_lattice() -> None:
    probabilities = np.loadtxt(QUANTUM_DATA / "timeevolution2d.txt", ndmin=2)
    side = math.isqrt(probabilities.shape[1])
    grids = probabilities.reshape(-1, side, side)
    indices = frame_indices(len(grids))
    times = np.linspace(0.0, 10.0, len(grids))
    color_scale, probability_floor = log_probability_scale(probabilities)

    figure, axis = plt.subplots(figsize=(5.2, 4.6))
    image = axis.imshow(
        np.maximum(grids[0], probability_floor),
        origin="lower",
        interpolation="nearest",
        cmap="magma",
        norm=color_scale,
    )
    figure.colorbar(
        image,
        ax=axis,
        extend="min",
        label=r"$|\psi(x,y)|^2$ (log scale)",
    )
    axis.set(xlabel="x site", ylabel="y site")
    axis.set_xticks(np.arange(-0.5, side, 1), minor=True)
    axis.set_yticks(np.arange(-0.5, side, 1), minor=True)
    axis.grid(which="minor", color="white", linewidth=0.3, alpha=0.2)
    axis.tick_params(which="minor", bottom=False, left=False)
    title = axis.set_title("Square-lattice quantum walk · log probability")

    def update(frame: int):
        image.set_data(np.maximum(grids[frame], probability_floor))
        title.set_text(f"Square-lattice quantum walk · t = {times[frame]:.2f}")
        return image, title

    figure.tight_layout()
    save_animation(figure, update, indices, "square-lattice.gif")


def honeycomb_geometry(side: int) -> tuple[np.ndarray, list[tuple]]:
    coordinates = []
    for y in range(side):
        for x in range(side):
            display_x = x + x // 2
            if y % 2:
                display_x += -0.5 if x % 2 == 0 else 0.5
            coordinates.append((display_x, math.sqrt(3) * y / 2))

    coordinates_array = np.asarray(coordinates)
    bonds = []
    for y in range(side):
        for x in range(side):
            site = x + side * y
            if y + 1 < side:
                bonds.append((coordinates_array[site], coordinates_array[site + side]))
            if x + 1 < side and (x + y) % 2 == 0:
                bonds.append((coordinates_array[site], coordinates_array[site + 1]))
    return coordinates_array, bonds


def honeycomb_lattice(
    data_file: Path,
    title_text: str,
    filename: str,
    *,
    duration: float,
    endpoints: tuple[tuple[int, int], tuple[int, int]] | None = None,
) -> None:
    probabilities = np.loadtxt(data_file, ndmin=2)
    side = math.isqrt(probabilities.shape[1])
    coordinates, bonds = honeycomb_geometry(side)
    indices = frame_indices(len(probabilities))
    times = np.linspace(0.0, duration, len(probabilities))
    color_scale, probability_floor = log_probability_scale(probabilities)

    figure, axis = plt.subplots(figsize=(5.5, 4.8))
    axis.add_collection(LineCollection(bonds, colors="#cbd5e1", linewidths=1.5))
    sites = axis.scatter(
        coordinates[:, 0],
        coordinates[:, 1],
        c=np.maximum(probabilities[0], probability_floor),
        cmap="viridis",
        norm=color_scale,
        s=150,
        edgecolor="#1e293b",
        linewidth=0.5,
        zorder=2,
    )
    figure.colorbar(
        sites,
        ax=axis,
        extend="min",
        label=r"$|\psi|^2$ (log scale)",
    )
    axis.set_aspect("equal")
    axis.set_xticks([])
    axis.set_yticks([])

    if endpoints is not None:
        endpoint_handles = []
        for label, (x, y), color in zip(
            ("source", "target"),
            endpoints,
            ("#0284c7", "#ea580c"),
        ):
            if not (0 <= x < side and 0 <= y < side):
                raise ValueError(f"{label} site ({x}, {y}) is outside the lattice")
            endpoint = coordinates[x + side * y]
            endpoint_handles.append(
                axis.scatter(
                    [endpoint[0]],
                    [endpoint[1]],
                    s=245,
                    facecolors="none",
                    edgecolors=color,
                    linewidths=1.8,
                    label=label,
                    zorder=3,
                )
            )
        axis.legend(
            handles=endpoint_handles,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.025),
            ncol=2,
            frameon=False,
            fontsize=8,
            handletextpad=0.35,
            columnspacing=1.3,
        )

    title = axis.set_title(f"{title_text} · log probability")

    def update(frame: int):
        sites.set_array(np.maximum(probabilities[frame], probability_floor))
        title.set_text(f"{title_text} · t = {times[frame]:.2f}")
        return sites, title

    figure.tight_layout()
    save_animation(figure, update, indices, filename)


def kepler_comparison() -> None:
    first = np.loadtxt(LABS / "06" / "data" / "kepler1")
    second = np.loadtxt(LABS / "06" / "data" / "kepler2")
    count = min(len(first), len(second))
    indices = frame_indices(count)

    figure, axis = plt.subplots(figsize=(5.2, 4.8))
    axis.set(
        xlim=(-1.15, 1.15),
        ylim=(-1.15, 1.15),
        xlabel="x (AU)",
        ylabel="y (AU)",
        title="Kepler orbit: integrator stability",
    )
    axis.set_aspect("equal")
    axis.scatter(
        [0],
        [0],
        color="#f59e0b",
        edgecolor="#b45309",
        linewidth=0.7,
        marker="*",
        s=180,
        label="Sun",
        zorder=4,
    )
    (path_one,) = axis.plot(
        [],
        [],
        color="#2563eb",
        linewidth=2.0,
        label="Symplectic Euler (drift–kick)",
    )
    (path_two,) = axis.plot(
        [],
        [],
        color="#dc2626",
        linewidth=1.8,
        linestyle="--",
        label="Forward Euler",
    )
    planet_one = axis.scatter([], [], color="#2563eb", s=38, zorder=5)
    planet_two = axis.scatter(
        [],
        [],
        facecolor="white",
        edgecolor="#dc2626",
        linewidth=1.4,
        s=42,
        zorder=5,
    )
    axis.legend(loc="upper right", fontsize=8, framealpha=0.92)
    time_text = axis.text(
        0.03,
        0.95,
        "",
        va="top",
        transform=axis.transAxes,
        bbox={"boxstyle": "round,pad=0.25", "fc": "white", "alpha": 0.88},
    )

    x_min = np.minimum.accumulate(
        np.minimum(np.minimum(first[:, 1], second[:, 1]), 0.0)
    )
    x_max = np.maximum.accumulate(
        np.maximum(np.maximum(first[:, 1], second[:, 1]), 0.0)
    )
    y_min = np.minimum.accumulate(
        np.minimum(np.minimum(first[:, 2], second[:, 2]), 0.0)
    )
    y_max = np.maximum.accumulate(
        np.maximum(np.maximum(first[:, 2], second[:, 2]), 0.0)
    )

    def update(frame: int):
        path_one.set_data(first[: frame + 1, 1], first[: frame + 1, 2])
        path_two.set_data(second[: frame + 1, 1], second[: frame + 1, 2])
        planet_one.set_offsets([first[frame, 1:3]])
        planet_two.set_offsets([second[frame, 1:3]])
        separation = np.linalg.norm(first[frame, 1:3] - second[frame, 1:3])
        time_text.set_text(
            f"t = {first[frame, 0]:.3f} yr\nmethod gap = {separation:.2e} AU"
        )

        center_x = 0.5 * (x_min[frame] + x_max[frame])
        center_y = 0.5 * (y_min[frame] + y_max[frame])
        span = max(
            x_max[frame] - x_min[frame],
            y_max[frame] - y_min[frame],
            2.0,
        )
        half_limit = 0.575 * span
        axis.set_xlim(center_x - half_limit, center_x + half_limit)
        axis.set_ylim(center_y - half_limit, center_y + half_limit)
        return path_one, path_two, planet_one, planet_two, time_text

    figure.tight_layout()
    save_animation(figure, update, indices, "kepler-comparison.gif")


def three_body() -> None:
    data = np.loadtxt(LABS / "07" / "data" / "3body")
    indices = frame_indices(len(data), 72)
    positions = [data[:, 1:3], data[:, 3:5], data[:, 5:7]]
    colors = ["#f59e0b", "#2563eb", "#dc2626"]

    limit = 1.1 * max(np.max(np.abs(data[:, 1:7])), 1.0)
    figure, axis = plt.subplots(figsize=(5.2, 4.8))
    axis.set(
        xlim=(-limit, limit),
        ylim=(-limit, limit),
        xlabel="x (AU)",
        ylabel="y (AU)",
        title="Perturbed equal-mass figure-eight",
    )
    axis.set_aspect("equal")
    paths = [
        axis.plot([], [], color=color, linewidth=1.4, alpha=0.5)[0]
        for color in colors
    ]
    bodies = [
        axis.scatter(
            [],
            [],
            color=color,
            edgecolor="white",
            linewidth=0.7,
            s=60,
            label=f"Body {label}",
            zorder=4,
        )
        for color, label in zip(colors, "ABC")
    ]
    barycenter = axis.scatter(
        [0],
        [0],
        color="#334155",
        marker="+",
        linewidth=1.2,
        s=52,
        label="Barycenter",
        zorder=3,
    )
    axis.legend(loc="lower right", ncol=2, fontsize=8, framealpha=0.9)
    time_text = axis.text(
        0.03,
        0.95,
        "",
        va="top",
        transform=axis.transAxes,
        bbox={"boxstyle": "round,pad=0.22", "fc": "white", "alpha": 0.85},
    )

    def update(frame: int):
        start = max(0, frame - 10000)
        for path, body, position in zip(paths, bodies, positions):
            trail = position[start : frame + 1 : 40]
            path.set_data(trail[:, 0], trail[:, 1])
            body.set_offsets([position[frame]])
        time_text.set_text(f"t = {data[frame, 0]:.2f} yr")
        return *paths, *bodies, barycenter, time_text

    figure.tight_layout()
    save_animation(figure, update, indices, "three-body.gif")


def diffusion_comparison() -> None:
    explicit = np.loadtxt(LABS / "08" / "data" / "diffusion.txt", ndmin=2)
    alternate = np.loadtxt(
        LABS / "08" / "data" / "diffusion2.txt", ndmin=2
    )
    count = min(len(explicit), len(alternate))
    indices = frame_indices(count)
    x = np.linspace(0, 1, explicit.shape[1])
    maximum = max(float(explicit.max()), float(alternate.max()))

    figure, axis = plt.subplots(figsize=(6.4, 3.8))
    (line_one,) = axis.plot(
        x,
        explicit[0],
        color="#2563eb",
        linewidth=2.2,
        label="Fixed-value endpoints",
    )
    (line_two,) = axis.plot(
        x,
        alternate[0],
        color="#dc2626",
        linewidth=2.0,
        linestyle="--",
        label="Zero-flux endpoints",
    )
    difference_fill = [
        axis.fill_between(
            x,
            explicit[0],
            alternate[0],
            color="#94a3b8",
            alpha=0.2,
            linewidth=0,
        )
    ]
    axis.set(
        xlim=(0, 1),
        ylim=(0, maximum * 1.05),
        xlabel="Normalized position",
        ylabel=r"Density $\rho(x)$",
        title="FTCS diffusion: boundary-condition comparison",
    )
    axis.legend(loc="upper left", fontsize=8, framealpha=0.92)
    time_text = axis.text(
        0.98,
        0.94,
        "",
        ha="right",
        va="top",
        transform=axis.transAxes,
        bbox={"boxstyle": "round,pad=0.22", "fc": "white", "alpha": 0.85},
    )

    def update(frame: int):
        difference_fill[0].remove()
        line_one.set_ydata(explicit[frame])
        line_two.set_ydata(alternate[frame])
        difference_fill[0] = axis.fill_between(
            x,
            explicit[frame],
            alternate[frame],
            color="#94a3b8",
            alpha=0.2,
            linewidth=0,
        )
        time_text.set_text(f"saved snapshot {frame + 1}/{count}")
        return line_one, line_two, difference_fill[0], time_text

    figure.tight_layout()
    save_animation(figure, update, indices, "diffusion.gif")


def molecular_dynamics() -> None:
    data = np.loadtxt(
        LABS / "11" / "data" / "positions.csv",
        delimiter=",",
        skiprows=1,
    )
    steps = data[:, 0]
    positions = data[:, 1:].reshape(len(data), -1, 2)
    indices = frame_indices(len(data), 72)
    lower = positions.min(axis=(0, 1))
    upper = positions.max(axis=(0, 1))

    step_deltas = np.diff(steps)
    step_deltas[step_deltas == 0] = 1
    velocities = np.empty_like(positions)
    velocities[1:] = np.diff(positions, axis=0) / step_deltas[:, None, None]
    velocities[0] = velocities[1]
    speeds = np.linalg.norm(velocities, axis=2)
    smoothed_speeds = np.empty_like(speeds)
    smoothed_speeds[0] = speeds[0]
    for frame in range(1, len(steps)):
        smoothed_speeds[frame] = (
            0.15 * speeds[frame] + 0.85 * smoothed_speeds[frame - 1]
        )

    positive_speeds = smoothed_speeds[smoothed_speeds > 0]
    speed_floor = max(float(positive_speeds.min()), 1e-12)
    log_speeds = np.log10(np.maximum(smoothed_speeds, speed_floor))
    color_min, color_max = np.quantile(log_speeds, [0.02, 0.98])

    figure, axis = plt.subplots(figsize=(5.2, 4.8))
    bonds = LineCollection(
        [], colors="#38bdf8", linewidths=0.55, alpha=0.35, zorder=1
    )
    axis.add_collection(bonds)
    points = axis.scatter(
        positions[0, :, 0],
        positions[0, :, 1],
        s=16,
        c=log_speeds[0],
        cmap="plasma",
        vmin=color_min,
        vmax=color_max,
        alpha=0.9,
        zorder=2,
    )
    figure.colorbar(
        points,
        ax=axis,
        pad=0.025,
        fraction=0.05,
        label=r"$\log_{10}$(speed per saved step)",
    )
    axis.plot(
        [lower[0], upper[0], upper[0], lower[0], lower[0]],
        [lower[1], lower[1], upper[1], upper[1], lower[1]],
        color="#334155",
        linewidth=1.2,
        zorder=3,
    )
    axis.set(
        xlim=(lower[0] - 0.14, upper[0] + 0.14),
        ylim=(lower[1] - 0.14, upper[1] + 0.14),
        xlabel="x",
        ylabel="y",
        title="Lennard–Jones cooling and crystallization",
    )
    axis.set_aspect("equal")
    axis.grid(False)
    step_text = axis.text(
        0.03,
        0.96,
        "",
        va="top",
        transform=axis.transAxes,
        bbox={"boxstyle": "round,pad=0.22", "fc": "white", "alpha": 0.82},
    )

    def lattice_bonds(frame: int) -> tuple[list[list[np.ndarray]], int]:
        cool_indices = np.flatnonzero(smoothed_speeds[frame] < 2.0e-4)
        if len(cool_indices) < 2:
            return [], len(cool_indices)

        cool_positions = positions[frame, cool_indices]
        offsets = cool_positions[:, None, :] - cool_positions[None, :, :]
        distances = np.linalg.norm(offsets, axis=2)
        np.fill_diagonal(distances, np.inf)
        neighbor_count = min(6, len(cool_indices) - 1)
        nearest = np.argpartition(
            distances, neighbor_count - 1, axis=1
        )[:, :neighbor_count]
        pairs = {
            tuple(sorted((index, int(neighbor))))
            for index, neighbors in enumerate(nearest)
            for neighbor in neighbors
            if distances[index, neighbor] <= 0.27
        }
        return [
            [cool_positions[first], cool_positions[second]]
            for first, second in pairs
        ], len(cool_indices)

    def update(frame: int):
        points.set_offsets(positions[frame])
        points.set_array(log_speeds[frame])
        segments, cool_count = lattice_bonds(frame)
        bonds.set_segments(segments)
        step_text.set_text(
            f"step {int(steps[frame])} · cool atoms {cool_count}/{len(positions[frame])}"
        )
        return bonds, points, step_text

    figure.tight_layout()
    save_animation(figure, update, indices, "molecular-dynamics.gif")


def projectile_motion() -> None:
    data = np.loadtxt(LABS / "05" / "data" / "projmotion")
    final_above_ground = np.flatnonzero(data[:, 2] >= 0)[-1]
    if final_above_ground + 1 < len(data):
        above = data[final_above_ground]
        below = data[final_above_ground + 1]
        fraction = above[2] / (above[2] - below[2])
        landing = above + fraction * (below - above)
        data = np.vstack((data[: final_above_ground + 1], landing))
    else:
        data = data[: final_above_ground + 1]
    indices = frame_indices(len(data))
    x, y = data[:, 1], data[:, 2]
    apex = int(np.argmax(y))
    vertical_margin = 0.08 * y.max()

    figure, axis = plt.subplots(figsize=(6.2, 3.8))
    axis.plot(x, y, color="#94a3b8", linewidth=1.2, linestyle=":")
    axis.axhline(0, color="#475569", linewidth=1.2)
    axis.fill_between(
        [x.min(), x.max() * 1.03],
        -vertical_margin,
        0,
        color="#d6d3d1",
        alpha=0.55,
        linewidth=0,
    )
    axis.scatter([x[apex]], [y[apex]], color="#475569", marker="x", s=36)
    axis.annotate(
        "apex",
        (x[apex], y[apex]),
        xytext=(6, -14),
        textcoords="offset points",
        fontsize=8,
        color="#475569",
    )
    (trajectory,) = axis.plot([], [], color="#2563eb", linewidth=2.4)
    projectile = axis.scatter(
        [],
        [],
        color="#ef4444",
        edgecolor="white",
        linewidth=0.7,
        s=58,
        zorder=4,
    )
    axis.set(
        xlim=(x.min(), x.max() * 1.03),
        ylim=(-vertical_margin, y.max() * 1.08),
        xlabel="Horizontal position x",
        ylabel="Height y",
        title="Ballistic trajectory under uniform gravity",
    )
    time_text = axis.text(
        0.03,
        0.93,
        "",
        va="top",
        transform=axis.transAxes,
        bbox={"boxstyle": "round,pad=0.22", "fc": "white", "alpha": 0.85},
    )

    def update(frame: int):
        trajectory.set_data(x[: frame + 1], y[: frame + 1])
        projectile.set_offsets([[x[frame], y[frame]]])
        time_text.set_text(f"t = {data[frame, 0]:.2f}")
        return trajectory, projectile, time_text

    figure.tight_layout()
    save_animation(figure, update, indices, "projectile-motion.gif")


def mass_spring() -> None:
    data = np.loadtxt(LABS / "04" / "data" / "mass_spring")
    data = data[np.concatenate(([True], np.diff(data[:, 0]) > 0))]
    indices = frame_indices(len(data))
    time, position = data[:, 0], data[:, 1]

    figure, axis = plt.subplots(figsize=(6.2, 3.8))
    axis.axhline(0, color="#64748b", linestyle="--", linewidth=1.1, alpha=0.8)
    axis.plot(time, position, color="#94a3b8", linewidth=1.0, alpha=0.45)
    (history,) = axis.plot([], [], color="#059669", linewidth=2.4)
    time_guide = axis.axvline(
        time[0], color="#f59e0b", linewidth=1.0, alpha=0.55
    )
    marker = axis.scatter(
        [],
        [],
        color="#f59e0b",
        edgecolor="white",
        linewidth=0.7,
        s=58,
        zorder=4,
    )
    axis.set(
        xlim=(time.min(), time.max()),
        ylim=(position.min() * 1.08, position.max() * 1.08),
        xlabel="Time t",
        ylabel="Displacement x",
        title="Mass–spring oscillator",
    )
    time_text = axis.text(
        0.98,
        0.93,
        "",
        ha="right",
        va="top",
        transform=axis.transAxes,
        bbox={"boxstyle": "round,pad=0.22", "fc": "white", "alpha": 0.85},
    )

    def update(frame: int):
        history.set_data(time[: frame + 1], position[: frame + 1])
        time_guide.set_xdata([time[frame], time[frame]])
        marker.set_offsets([[time[frame], position[frame]]])
        time_text.set_text(f"t = {time[frame]:.3f}")
        return history, time_guide, marker, time_text

    figure.tight_layout()
    save_animation(figure, update, indices, "mass-spring.gif")


def density_of_states() -> None:
    values = [
        np.loadtxt(QUANTUM_DATA / "eigenvalues.dat"),
        np.loadtxt(QUANTUM_DATA / "eigenvalues2.dat"),
        np.loadtxt(QUANTUM_DATA / "eigenvalues2d.dat"),
    ]
    labels = ["1D uniform", "1D alternating", "2D square"]
    colors = ["#2563eb", "#7c3aed", "#059669"]

    figure, axis = plt.subplots(figsize=(6.4, 3.8))
    for eigenvalues, label, color in zip(values, labels, colors):
        axis.hist(
            eigenvalues,
            bins=48,
            density=True,
            histtype="step",
            linewidth=1.8,
            label=label,
            color=color,
        )
    axis.set(
        xlabel="Energy",
        ylabel="Normalized density",
        title="Tight-binding density of states",
    )
    axis.legend()
    figure.tight_layout()
    output = ASSETS / "density-of-states.png"
    figure.savefig(output, dpi=160)
    plt.close(figure)
    print(f"generated {output.relative_to(ROOT)}")


def main() -> None:
    ASSETS.mkdir(parents=True, exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")

    probability_chain(
        QUANTUM_DATA / "timeevolution.txt",
        "Barrier scattering on a 1D chain",
        "quantum-chain.gif",
        barrier_site=10,
        barrier_height=2.0,
        t_max=25.1,
    )
    probability_chain(
        QUANTUM_DATA / "annealed_barrier_evolution.txt",
        "Annealed end-to-end state transfer",
        "annealed-barrier.gif",
        barrier_site=10,
        barrier_height=2.0,
        t_max=8.0 * math.pi,
    )
    square_lattice()
    honeycomb_lattice(
        QUANTUM_DATA / "timehex.txt",
        "Honeycomb-lattice quantum walk",
        "honeycomb-lattice.gif",
        duration=2.0 * math.pi,
    )
    honeycomb_lattice(
        QUANTUM_DATA / "annealed_evolution.txt",
        "Annealed honeycomb state transfer",
        "annealed-honeycomb.gif",
        duration=metadata_float(QUANTUM_DATA / "annealed_couplings.dat", "time"),
        endpoints=((0, 0), (5, 0)),
    )
    kepler_comparison()
    three_body()
    diffusion_comparison()
    molecular_dynamics()
    projectile_motion()
    mass_spring()
    density_of_states()


if __name__ == "__main__":
    main()
