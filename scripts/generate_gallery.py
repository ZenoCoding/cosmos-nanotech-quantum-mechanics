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
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
LABS = ROOT / "labs"
QUANTUM_DATA = ROOT / "labs" / "12" / "data"
ASSETS = ROOT / "docs" / "assets"
MAX_FRAMES = 60
FPS = 15


def frame_indices(length: int, maximum: int = MAX_FRAMES) -> np.ndarray:
    """Return evenly spaced, unique frame indices."""
    return np.unique(np.linspace(0, length - 1, min(length, maximum), dtype=int))


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


def probability_chain(data_file: Path, title: str, filename: str) -> None:
    probabilities = np.loadtxt(data_file, ndmin=2)
    indices = frame_indices(len(probabilities))
    x = np.arange(1, probabilities.shape[1] + 1)

    figure, axis = plt.subplots(figsize=(6.4, 3.8))
    (curve,) = axis.plot(x, probabilities[0], color="#5b5bd6", linewidth=2.2)
    axis.fill_between(x, 0, probabilities[0], color="#a5b4fc", alpha=0.35)
    fill = [axis.collections[-1]]
    axis.set(
        xlim=(1, len(x)),
        ylim=(0, max(1.0, float(probabilities.max()) * 1.05)),
        xlabel="Lattice site",
        ylabel=r"Probability $|\psi|^2$",
        title=title,
    )
    time_text = axis.text(0.98, 0.92, "", ha="right", transform=axis.transAxes)

    def update(frame: int):
        fill[0].remove()
        curve.set_ydata(probabilities[frame])
        fill[0] = axis.fill_between(
            x, 0, probabilities[frame], color="#a5b4fc", alpha=0.35
        )
        time_text.set_text(f"frame {frame + 1}/{len(probabilities)}")
        return curve, fill[0], time_text

    figure.tight_layout()
    save_animation(figure, update, indices, filename)


def square_lattice() -> None:
    probabilities = np.loadtxt(QUANTUM_DATA / "timeevolution2d.txt", ndmin=2)
    side = math.isqrt(probabilities.shape[1])
    grids = probabilities.reshape(-1, side, side)
    indices = frame_indices(len(grids))

    figure, axis = plt.subplots(figsize=(5.2, 4.6))
    image = axis.imshow(
        grids[0],
        origin="lower",
        cmap="magma",
        vmin=0,
        vmax=max(float(probabilities.max()), 1e-12),
    )
    figure.colorbar(image, ax=axis, label=r"$|\psi(x,y)|^2$")
    axis.set(xlabel="x site", ylabel="y site")
    title = axis.set_title("Square-lattice quantum walk")

    def update(frame: int):
        image.set_data(grids[frame])
        title.set_text(f"Square-lattice quantum walk · frame {frame + 1}")
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


def honeycomb_lattice(data_file: Path, title_text: str, filename: str) -> None:
    probabilities = np.loadtxt(data_file, ndmin=2)
    side = math.isqrt(probabilities.shape[1])
    coordinates, bonds = honeycomb_geometry(side)
    indices = frame_indices(len(probabilities))

    figure, axis = plt.subplots(figsize=(5.5, 4.8))
    axis.add_collection(LineCollection(bonds, colors="#cbd5e1", linewidths=1.5))
    sites = axis.scatter(
        coordinates[:, 0],
        coordinates[:, 1],
        c=probabilities[0],
        cmap="viridis",
        vmin=0,
        vmax=max(float(probabilities.max()), 1e-12),
        s=150,
        edgecolor="#1e293b",
        linewidth=0.5,
        zorder=2,
    )
    figure.colorbar(sites, ax=axis, label=r"$|\psi|^2$")
    axis.set_aspect("equal")
    axis.set_xticks([])
    axis.set_yticks([])
    title = axis.set_title(title_text)

    def update(frame: int):
        sites.set_array(probabilities[frame])
        title.set_text(f"{title_text} · frame {frame + 1}")
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
        title="Kepler integrator comparison",
    )
    axis.set_aspect("equal")
    axis.scatter([0], [0], color="#f59e0b", s=110, label="Sun")
    (path_one,) = axis.plot(
        [], [], color="#2563eb", label="Symplectic Euler (drift–kick)"
    )
    (path_two,) = axis.plot(
        [], [], color="#dc2626", label="Forward Euler"
    )
    planet_one = axis.scatter([], [], color="#2563eb", s=35)
    planet_two = axis.scatter([], [], color="#dc2626", s=35)
    axis.legend(loc="upper right")
    time_text = axis.text(0.03, 0.94, "", transform=axis.transAxes)

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
        time_text.set_text(f"t = {first[frame, 0]:.3f} · Δr = {separation:.2e}")

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
        xlabel="x",
        ylabel="y",
        title="Three-body dynamics",
    )
    axis.set_aspect("equal")
    paths = [axis.plot([], [], color=color, alpha=0.65)[0] for color in colors]
    bodies = [axis.scatter([], [], color=color, s=55) for color in colors]
    time_text = axis.text(0.03, 0.94, "", transform=axis.transAxes)

    def update(frame: int):
        start = max(0, frame - 10000)
        for path, body, position in zip(paths, bodies, positions):
            trail = position[start : frame + 1 : 40]
            path.set_data(trail[:, 0], trail[:, 1])
            body.set_offsets([position[frame]])
        time_text.set_text(f"t = {data[frame, 0]:.2f}")
        return *paths, *bodies, time_text

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
        label="FTCS · fixed-value boundaries",
    )
    (line_two,) = axis.plot(
        x,
        alternate[0],
        color="#dc2626",
        label="FTCS · zero-flux boundaries",
    )
    axis.set(
        xlim=(0, 1),
        ylim=(0, maximum * 1.05),
        xlabel="Normalized position",
        ylabel="Concentration",
        title="Diffusion equation",
    )
    axis.legend()
    time_text = axis.text(0.98, 0.92, "", ha="right", transform=axis.transAxes)

    def update(frame: int):
        line_one.set_ydata(explicit[frame])
        line_two.set_ydata(alternate[frame])
        time_text.set_text(f"step {frame}")
        return line_one, line_two, time_text

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

    figure, axis = plt.subplots(figsize=(5.2, 4.8))
    points = axis.scatter(
        positions[0, :, 0],
        positions[0, :, 1],
        s=16,
        color="#7c3aed",
        alpha=0.8,
    )
    axis.set(
        xlim=(lower[0] - 0.25, upper[0] + 0.25),
        ylim=(lower[1] - 0.25, upper[1] + 0.25),
        xlabel="x",
        ylabel="y",
        title="Lennard–Jones molecular dynamics",
    )
    axis.set_aspect("equal")
    step_text = axis.text(0.03, 0.94, "", transform=axis.transAxes)

    def update(frame: int):
        points.set_offsets(positions[frame])
        step_text.set_text(f"step {int(steps[frame])}")
        return points, step_text

    figure.tight_layout()
    save_animation(figure, update, indices, "molecular-dynamics.gif")


def projectile_motion() -> None:
    data = np.loadtxt(LABS / "05" / "data" / "projmotion")
    indices = frame_indices(len(data))
    x, y = data[:, 1], data[:, 2]

    figure, axis = plt.subplots(figsize=(6.2, 3.8))
    (trajectory,) = axis.plot([], [], color="#2563eb", linewidth=2)
    projectile = axis.scatter([], [], color="#ef4444", s=55)
    axis.set(
        xlim=(x.min(), x.max() * 1.03),
        ylim=(0, y.max() * 1.08),
        xlabel="Horizontal position",
        ylabel="Height",
        title="Projectile motion",
    )
    time_text = axis.text(0.03, 0.92, "", transform=axis.transAxes)

    def update(frame: int):
        trajectory.set_data(x[: frame + 1], y[: frame + 1])
        projectile.set_offsets([[x[frame], y[frame]]])
        time_text.set_text(f"t = {data[frame, 0]:.2f}")
        return trajectory, projectile, time_text

    figure.tight_layout()
    save_animation(figure, update, indices, "projectile-motion.gif")


def mass_spring() -> None:
    data = np.loadtxt(LABS / "04" / "data" / "mass_spring")
    indices = frame_indices(len(data))
    time, position = data[:, 0], data[:, 1]

    figure, axis = plt.subplots(figsize=(6.2, 3.8))
    (history,) = axis.plot([], [], color="#059669", linewidth=2)
    marker = axis.scatter([], [], color="#f59e0b", s=55)
    axis.set(
        xlim=(time.min(), time.max()),
        ylim=(position.min() * 1.08, position.max() * 1.08),
        xlabel="Time",
        ylabel="Displacement",
        title="Mass–spring oscillator",
    )

    def update(frame: int):
        history.set_data(time[: frame + 1], position[: frame + 1])
        marker.set_offsets([[time[frame], position[frame]]])
        return history, marker

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
    )
    probability_chain(
        QUANTUM_DATA / "annealed_barrier_evolution.txt",
        "Annealed end-to-end state transfer",
        "annealed-barrier.gif",
    )
    square_lattice()
    honeycomb_lattice(
        QUANTUM_DATA / "timehex.txt",
        "Honeycomb-lattice quantum walk",
        "honeycomb-lattice.gif",
    )
    honeycomb_lattice(
        QUANTUM_DATA / "annealed_evolution.txt",
        "Annealed honeycomb state transfer",
        "annealed-honeycomb.gif",
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
