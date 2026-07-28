"""Animate one-dimensional probability data produced by the C evolvers."""

import argparse
import math

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Animate probability on a one-dimensional chain."
    )
    parser.add_argument("input", help="probability data file")
    parser.add_argument(
        "--tmax",
        type=float,
        default=8.0 * math.pi,
        help="final simulation time (default: 8*pi)",
    )
    parser.add_argument(
        "--interval",
        type=int,
        default=None,
        help="milliseconds between frames (default: about 1 second total)",
    )
    parser.add_argument(
        "--save",
        metavar="FILE",
        help="save the animation instead of opening a window",
    )
    return parser.parse_args()


def initial(curve, timer):
    curve.set_data([], [])
    timer.set_text("")
    return curve, timer

def f(frame, ts, xs, fxs, curve, timer):
    curve.set_data(xs, fxs[frame, :])
    timer.set_text(f"t = {ts[frame]:.3f}")
    return curve, timer

args = parse_arguments()
if args.tmax < 0:
    raise ValueError("--tmax must be nonnegative")
if args.interval is not None and args.interval < 1:
    raise ValueError("--interval must be at least 1 millisecond")

fxs = np.loadtxt(args.input, ndmin=2)
nts, nxs = fxs.shape
interval = args.interval if args.interval is not None else max(1, 1000 / nts)

fig, ax = plt.subplots(figsize=(8,6), dpi=200)

ax.set_ylim(0,1)
ax.set_xlim(1, nxs)
ax.set_title("Quantum State Transfer")

xs = np.linspace(1, nxs, nxs)
ts = np.linspace(0, args.tmax, nts)

ax.set_xlabel(r"$x$", fontsize=10)
ax.set_ylabel(r"$|\psi(x)|^2$", fontsize=10)

ax.tick_params(axis="both", which="major", labelsize=10, length=6, width=1)

#ax.grid(True)

curve, = ax.plot([], [], lw=1.5, color="blue", alpha=0.7)
timer = ax.text(x=0.98, y=0.95, s="", transform=ax.transAxes, fontsize=10, ha="right")

animation = ani.FuncAnimation(fig,
                                               f,
                                               frames=nts,
                                               init_func=lambda: initial(curve, timer),
                                               fargs=(ts, xs, fxs, curve, timer),
                                               interval=interval,
                                               blit=True)

fig.tight_layout()
if args.save:
    animation.save(args.save)
    print(f"Saved animation to {args.save}")
else:
    plt.show()
