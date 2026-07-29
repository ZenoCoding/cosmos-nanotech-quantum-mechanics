from pathlib import Path

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np


DATA_FILE = Path(__file__).resolve().parents[1] / "data" / "3body"
data = np.loadtxt(DATA_FILE)
# Keep every tenth simulation step for smooth motion without excessive frames.
data = data[::10]
t = data[:, 0]
sun_x, sun_y = data[:, 1], data[:, 2]
earth_x, earth_y = data[:, 3], data[:, 4]
jupiter_x, jupiter_y = data[:, 5], data[:, 6]

fig, ax = plt.subplots()

# Equal, centered limits keep the orbit's shape accurate on screen.
limit = 1.15 * max(
    1.0,
    np.max(np.abs(sun_x)), np.max(np.abs(sun_y)),
    np.max(np.abs(earth_x)), np.max(np.abs(earth_y)),
    np.max(np.abs(jupiter_x)), np.max(np.abs(jupiter_y)),
)
ax.set_xlim(-limit, limit)
ax.set_ylim(-limit, limit)
ax.set_aspect("equal")
ax.set_xlabel("x (AU)")
ax.set_ylabel("y (AU)")
ax.set_title("Perturbed Three-Body Orbit")

# Line trails are substantially faster to update than growing scatter plots.
sun_orbit, = ax.plot([], [], color="orange", linewidth=1.0, alpha=0.55)
earth_orbit, = ax.plot([], [], color="royalblue", linewidth=1.0, alpha=0.55)
jupiter_orbit, = ax.plot([], [], color="firebrick", linewidth=1.0, alpha=0.55)
sun = ax.scatter([], [], color="orange", s=60, label="Body A")
earth = ax.scatter([], [], color="royalblue", s=60, label="Body B")
jupiter = ax.scatter([], [], color="firebrick", s=60, label="Body C")
time_text = ax.text(0.02, 0.96, "", transform=ax.transAxes)
ax.legend(loc="upper right")

# Show a moving trail instead of redrawing the entire history each frame.
TRAIL_LENGTH = 500


def update(frame):
    start = max(0, frame - TRAIL_LENGTH)
    sun_orbit.set_data(sun_x[start : frame + 1], sun_y[start : frame + 1])
    earth_orbit.set_data(earth_x[start : frame + 1], earth_y[start : frame + 1])
    jupiter_orbit.set_data(
        jupiter_x[start : frame + 1], jupiter_y[start : frame + 1]
    )
    sun.set_offsets([[sun_x[frame], sun_y[frame]]])
    earth.set_offsets([[earth_x[frame], earth_y[frame]]])
    jupiter.set_offsets([[jupiter_x[frame], jupiter_y[frame]]])
    time_text.set_text(f"Time: {t[frame]:.2f} years")
    return sun_orbit, earth_orbit, jupiter_orbit, sun, earth, jupiter, time_text


ani = animation.FuncAnimation(
    fig,
    update,
    frames=len(t),
    interval=16,
    blit=True,
    repeat=True,
    cache_frame_data=False,
)

plt.show()
