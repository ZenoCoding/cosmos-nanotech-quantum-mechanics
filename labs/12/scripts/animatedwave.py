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
    parser.add_argument(
        "--barrier-site",
        type=int,
        default=None,
        help="draw an on-site barrier at this zero-based lattice site",
    )
    parser.add_argument(
        "--barrier-height",
        type=float,
        default=2.0,
        help="on-site energy shift at the barrier (default: 2)",
    )
    return parser.parse_args()


def initial(curve, timer, norm_text):
    curve.set_data([], [])
    timer.set_text("")
    norm_text.set_text("")
    return curve, timer, norm_text


def f(frame, ts, xs, fxs, curve, timer, norm_text):
    curve.set_data(xs, fxs[frame, :])
    timer.set_text(f"t = {ts[frame]:.3f}")
    norm_text.set_text(f"total probability = {fxs[frame, :].sum():.6f}")
    return curve, timer, norm_text


args = parse_arguments()
if args.tmax < 0:
    raise ValueError("--tmax must be nonnegative")
if args.interval is not None and args.interval < 1:
    raise ValueError("--interval must be at least 1 millisecond")
if args.barrier_height < 0:
    raise ValueError("--barrier-height must be nonnegative")

fxs = np.loadtxt(args.input, ndmin=2)
nts, nxs = fxs.shape
interval = args.interval if args.interval is not None else max(1, 1000 / nts)
barrier_site = args.barrier_site
if barrier_site is not None and not 0 <= barrier_site < nxs:
    raise ValueError(f"--barrier-site must be between 0 and {nxs - 1}")

fig, ax = plt.subplots(figsize=(8, 6), dpi=200)

ax.set_ylim(0, max(1.0, float(fxs.max()) * 1.05))
ax.set_xlim(-0.5, nxs - 0.5)
ax.set_title("Quantum evolution on a one-dimensional chain")

xs = np.arange(nxs)
ts = np.linspace(0, args.tmax, nts)

ax.set_xlabel("Lattice site (zero-based)", fontsize=10)
ax.set_ylabel(r"$|\psi(x)|^2$", fontsize=10)

ax.tick_params(axis="both", which="major", labelsize=10, length=6, width=1)
ax.grid(axis="y", alpha=0.25)

if barrier_site is not None:
    ax.axvspan(
        barrier_site - 0.45,
        barrier_site + 0.45,
        color="#f97316",
        alpha=0.18,
        zorder=0,
    )
    ax.axvline(barrier_site, color="#c2410c", linestyle="--", linewidth=1.2)
    ax.text(
        barrier_site,
        0.98,
        rf"barrier $\Delta E=+{args.barrier_height:g}$",
        transform=ax.get_xaxis_transform(),
        color="#9a3412",
        fontsize=9,
        ha="center",
        va="top",
        bbox={
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": 0.8,
            "pad": 2,
        },
    )

curve, = ax.plot([], [], lw=2.0, color="#4f46e5", alpha=0.9)
timer = ax.text(
    x=0.98,
    y=0.95,
    s="",
    transform=ax.transAxes,
    fontsize=10,
    ha="right",
)
norm_text = ax.text(
    x=0.98,
    y=0.89,
    s="",
    transform=ax.transAxes,
    fontsize=8,
    color="#475569",
    ha="right",
)

animation = ani.FuncAnimation(
    fig,
    f,
    frames=nts,
    init_func=lambda: initial(curve, timer, norm_text),
    fargs=(ts, xs, fxs, curve, timer, norm_text),
    interval=interval,
    blit=True,
)

fig.tight_layout()
if args.save:
    animation.save(args.save)
    print(f"Saved animation to {args.save}")
else:
    plt.show()
