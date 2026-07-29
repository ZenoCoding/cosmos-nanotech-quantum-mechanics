from pathlib import Path

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np


DATA_FILE = Path(__file__).resolve().parents[1] / "data" / "kepler2"
data = np.loadtxt(DATA_FILE)
t = data[:, 0]
x = data[:, 1]
y = data[:, 2]

fig, ax = plt.subplots()
limit = 1.15 * max(1.0, np.max(np.abs(x)), np.max(np.abs(y)))
ax.set_xlim(-limit, limit)
ax.set_ylim(-limit, limit)
ax.set_aspect("equal")
ax.set_xlabel("x (AU)")
ax.set_ylabel("y (AU)")
ax.set_title("Kepler Orbit")

orbit, = ax.plot([], [], color="royalblue", linewidth=1.2, alpha=0.65)
sun = ax.scatter([0.0], [0.0], color="orange", s=100, label="Sun")
planet = ax.scatter([], [], color="royalblue", s=40, label="Planet")
time_text = ax.text(0.02, 0.96, "", transform=ax.transAxes)
ax.legend(loc="upper right")


def update(frame):
    orbit.set_data(x[: frame + 1], y[: frame + 1])
    planet.set_offsets([[x[frame], y[frame]]])
    time_text.set_text(f"Time: {t[frame]:.3f} years")
    return orbit, sun, planet, time_text


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
