import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from pathlib import Path

# Load CSV data
df = pd.read_csv(Path(__file__).resolve().parents[1] / "data" / "positions.csv")
steps = df["step"].values

# Extract 2D positions
pos_cols = [col for col in df.columns if col != "step"]
positions = df[pos_cols].values.reshape(len(df), -1, 2)  # shape: (num_steps, num_particles, 2)
num_particles = positions.shape[1]

# Determine bounds
x_min, x_max = positions[:, :, 0].min(), positions[:, :, 0].max()
y_min, y_max = positions[:, :, 1].min(), positions[:, :, 1].max()
margin = 0.5

# Setup plot
fig, ax = plt.subplots()
sc = ax.scatter(positions[0, :, 0], positions[0, :, 1], s=40, c='#9a05fc')

ax.set_xlim(x_min - margin, x_max + margin)
ax.set_ylim(y_min - margin, y_max + margin)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("2D Molecular Dynamics with Lennard-Jones")

# Init function
def init():
    sc.set_offsets(positions[0])
    return sc,

# Frame update
def update(frame):
    sc.set_offsets(positions[frame])
    ax.set_title(f"Step {steps[frame]}")
    return sc,

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
