import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.collections import LineCollection
import numpy as np
from pathlib import Path

# Lattice-highlighting controls. Speeds are measured in position units per
# saved simulation step because positions.csv does not contain velocities.
SPEED_THRESHOLD = 2.0e-4
SPEED_SMOOTHING = 0.15
NEAREST_NEIGHBORS = 6
MAX_BOND_LENGTH = 0.27
BOND_COLOR = "#4a90e2"
BOND_ALPHA = 0.28

# Load CSV data
df = pd.read_csv(Path(__file__).resolve().parents[1] / "data" / "positions.csv")
steps = df["step"].values

# Extract 2D positions
pos_cols = [col for col in df.columns if col != "step"]
positions = df[pos_cols].values.reshape(len(df), -1, 2)  # shape: (num_steps, num_particles, 2)
num_particles = positions.shape[1]

# Estimate and smooth speeds from the saved positions. Smoothing keeps bonds
# from flickering when an atom briefly vibrates above the threshold.
step_deltas = np.diff(steps)
step_deltas[step_deltas == 0] = 1
speeds = np.empty(positions.shape[:2])
speeds[1:] = (
    np.linalg.norm(np.diff(positions, axis=0), axis=2) / step_deltas[:, None]
)
speeds[0] = speeds[1]
smoothed_speeds = np.empty_like(speeds)
smoothed_speeds[0] = speeds[0]
for frame in range(1, len(steps)):
    smoothed_speeds[frame] = (
        SPEED_SMOOTHING * speeds[frame]
        + (1.0 - SPEED_SMOOTHING) * smoothed_speeds[frame - 1]
    )

# Determine bounds
x_min, x_max = positions[:, :, 0].min(), positions[:, :, 0].max()
y_min, y_max = positions[:, :, 1].min(), positions[:, :, 1].max()
margin = 0.5

# Setup plot
fig, ax = plt.subplots()
bonds = LineCollection(
    [], colors=BOND_COLOR, linewidths=0.7, alpha=BOND_ALPHA, zorder=1
)
ax.add_collection(bonds)
sc = ax.scatter(
    positions[0, :, 0], positions[0, :, 1], s=40, c="#9a05fc", zorder=2
)

ax.set_xlim(x_min - margin, x_max + margin)
ax.set_ylim(y_min - margin, y_max + margin)
ax.set_xlabel("x")
ax.set_ylabel("y")
title = ax.set_title("2D Molecular Dynamics with Lennard-Jones")


def lattice_bonds(frame):
    """Return bonds between cool atoms that are nearby lattice neighbors."""
    cool_indices = np.flatnonzero(smoothed_speeds[frame] < SPEED_THRESHOLD)
    if len(cool_indices) < 2:
        return [], len(cool_indices)

    cool_positions = positions[frame, cool_indices]
    offsets = cool_positions[:, None, :] - cool_positions[None, :, :]
    distances = np.linalg.norm(offsets, axis=2)
    np.fill_diagonal(distances, np.inf)

    k = min(NEAREST_NEIGHBORS, len(cool_indices) - 1)
    nearest = np.argpartition(distances, k - 1, axis=1)[:, :k]

    # A set removes duplicate i->j and j->i connections.
    pairs = set()
    for i, neighbors in enumerate(nearest):
        for j in neighbors:
            if distances[i, j] <= MAX_BOND_LENGTH:
                pairs.add(tuple(sorted((i, int(j)))))

    segments = [[cool_positions[i], cool_positions[j]] for i, j in pairs]
    return segments, len(cool_indices)

# Init function
def init():
    sc.set_offsets(positions[0])
    segments, cool_count = lattice_bonds(0)
    bonds.set_segments(segments)
    title.set_text(f"Step {steps[0]}  |  cool atoms: {cool_count}/{num_particles}")
    return bonds, sc, title

# Frame update
def update(frame):
    sc.set_offsets(positions[frame])
    segments, cool_count = lattice_bonds(frame)
    bonds.set_segments(segments)
    title.set_text(
        f"Step {steps[frame]}  |  cool atoms: {cool_count}/{num_particles}"
    )
    return bonds, sc, title

# Animate
skip = 1
ani = animation.FuncAnimation(
    fig,
    update,
    frames=range(0, len(steps), skip),
    init_func=init,
    blit=True,
    interval=16,
    repeat=False
)

plt.tight_layout()
plt.show()
