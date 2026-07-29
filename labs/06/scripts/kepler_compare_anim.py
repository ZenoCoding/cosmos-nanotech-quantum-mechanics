from pathlib import Path

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np


DATA_DIR = Path(__file__).resolve().parents[1] / "data"
kepler1 = np.loadtxt(DATA_DIR / "kepler1")
kepler2 = np.loadtxt(DATA_DIR / "kepler2")

frame_count = min(len(kepler1), len(kepler2))
kepler1 = kepler1[:frame_count]
kepler2 = kepler2[:frame_count]

t = kepler1[:, 0]
x1, y1 = kepler1[:, 1], kepler1[:, 2]
x2, y2 = kepler2[:, 1], kepler2[:, 2]

fig, ax = plt.subplots()
ax.set_xlim(-1.15, 1.15)
ax.set_ylim(-1.15, 1.15)
ax.set_aspect("equal")
ax.set_xlabel("x (AU)")
ax.set_ylabel("y (AU)")
ax.set_title("Kepler Integration Order Comparison")

path1, = ax.plot(
    [],
    [],
    color="royalblue",
    linewidth=1.5,
    label="Kepler 1: move, new force, accelerate",
)
path2, = ax.plot(
    [],
    [],
    color="firebrick",
    linewidth=1.5,
    label="Kepler 2: old force, move, accelerate",
)
sun = ax.scatter([0.0], [0.0], color="orange", s=100, label="Sun")
planet1 = ax.scatter([], [], color="royalblue", s=40)
planet2 = ax.scatter([], [], color="firebrick", s=40)
time_text = ax.text(0.02, 0.96, "", transform=ax.transAxes)
separation_text = ax.text(0.02, 0.91, "", transform=ax.transAxes)
ax.legend(loc="upper right")

# Precompute growing bounds so each frame can rescale without scanning all
# previous points. Including zero ensures that the Sun always remains visible.
x_min = np.minimum.accumulate(np.minimum(np.minimum(x1, x2), 0.0))
x_max = np.maximum.accumulate(np.maximum(np.maximum(x1, x2), 0.0))
y_min = np.minimum.accumulate(np.minimum(np.minimum(y1, y2), 0.0))
y_max = np.maximum.accumulate(np.maximum(np.maximum(y1, y2), 0.0))


def update(frame):
    path1.set_data(x1[: frame + 1], y1[: frame + 1])
    path2.set_data(x2[: frame + 1], y2[: frame + 1])
    planet1.set_offsets([[x1[frame], y1[frame]]])
    planet2.set_offsets([[x2[frame], y2[frame]]])

    separation = np.hypot(x1[frame] - x2[frame], y1[frame] - y2[frame])
    time_text.set_text(f"Time: {t[frame]:.3f} years")
    separation_text.set_text(f"Path difference: {separation:.6f} AU")

    # Use the same span on both axes so orbital shapes are not distorted.
    center_x = 0.5 * (x_min[frame] + x_max[frame])
    center_y = 0.5 * (y_min[frame] + y_max[frame])
    span = max(x_max[frame] - x_min[frame], y_max[frame] - y_min[frame], 2.0)
    half_limit = 0.575 * span
    ax.set_xlim(center_x - half_limit, center_x + half_limit)
    ax.set_ylim(center_y - half_limit, center_y + half_limit)

    return path1, path2, sun, planet1, planet2, time_text, separation_text


ani = animation.FuncAnimation(
    fig,
    update,
    frames=frame_count,
    interval=16,
    blit=False,
    repeat=True,
    cache_frame_data=False,
)

plt.show()
