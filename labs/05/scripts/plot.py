import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

data = np.loadtxt(Path(__file__).resolve().parents[1] / "data" / "projmotion")

t = data[:,0]
y = data[:,1]
# e = data[:,2]

plt.xlim(t.min(), t.max())
plt.ylim(y.min(), y.max()*(1+0.1))

plt.scatter(t, y, s=1)
plt.tick_params(direction='in', top=True, right=True)
plt.xlabel("Time t")
plt.ylabel("Height y")
plt.title("Particle Height vs Time")

plt.show()

# plt.scatter(t, e, s=1)
# plt.tick_params(direction='in', top=True, right=True)
# plt.xlabel("Time t")
# plt.ylabel("Energy E")
# plt.title("Particle Energy vs Time")

# plt.show()
