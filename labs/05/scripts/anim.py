import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pathlib import Path

data = np.loadtxt(Path(__file__).resolve().parents[1] / "data" / "projmotion")

t = data[:, 0]
x = data[:, 1]
y = data[:, 2]

fig, ax = plt.subplots()

ax.set_xlim(x.min(), x.max())
ax.set_ylim(y.min(), y.max())
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title('Particle Trajectory in 2D')

scat = ax.scatter([], [], color='blue', s=10)
time_text = ax.text(
    0.05, 0.95, "", transform=ax.transAxes, fontsize=14, verticalalignment='top'
)

#update plot each frame
def update(frame):
    points = np.column_stack((x[:frame+1], y[:frame+1]))
    time_text.set_text("Time:" + str(t[frame]))
    scat.set_offsets(points)
    return scat, time_text

ani = animation.FuncAnimation(fig, update, frames=len(x), interval=30, blit=True)

plt.show()
