#
# Animate probability data produced by timeevolve2d.c.
#
# Run like:
#
#     python3 animatedwave2d.py timeevolution2d.txt
#
# Each input row is one time slice. The number of columns must be a perfect
# square because the values are reshaped into an L-by-L lattice.
#
import argparse
import math

import matplotlib.animation as animation
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np


DEFAULT_LOG_DECADES = 6


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Animate probability on a square lattice."
    )
    parser.add_argument(
        "input",
        nargs="?",
        default="timeevolution2d.txt",
        help="time-evolution data file (default: timeevolution2d.txt)",
    )
    parser.add_argument(
        "--tmax",
        type=float,
        default=10.0,
        help="final simulation time shown in the title (default: 10)",
    )
    parser.add_argument(
        "--interval",
        type=int,
        default=100,
        help="milliseconds between frames (default: 100)",
    )
    parser.add_argument(
        "--scale",
        choices=("log", "linear"),
        default="log",
        help="probability color scale (default: log)",
    )
    parser.add_argument(
        "--log-floor",
        type=float,
        metavar="PROBABILITY",
        help=(
            "absolute probability shown as the bottom of the log scale "
            f"(default: {DEFAULT_LOG_DECADES} decades below the global maximum)"
        ),
    )
    parser.add_argument(
        "--save",
        metavar="FILE",
        help="save the animation instead of opening a window",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    probabilities = np.loadtxt(args.input, ndmin=2)

    frame_count, site_count = probabilities.shape
    side_length = math.isqrt(site_count)
    if side_length * side_length != site_count:
        raise ValueError(
            f"{args.input} has {site_count} columns; "
            "a square lattice requires a perfect square"
        )
    if args.tmax < 0:
        raise ValueError("--tmax must be nonnegative")
    if args.interval < 1:
        raise ValueError("--interval must be at least 1 millisecond")
    if not np.all(np.isfinite(probabilities)):
        raise ValueError(f"{args.input} contains a non-finite probability")
    if np.any(probabilities < 0.0):
        raise ValueError(f"{args.input} contains a negative probability")

    maximum_probability = float(probabilities.max())
    if maximum_probability <= 0.0:
        raise ValueError(f"{args.input} contains no positive probabilities")

    log_floor = None
    if args.scale == "log":
        log_floor = (
            maximum_probability * 10.0 ** -DEFAULT_LOG_DECADES
            if args.log_floor is None
            else args.log_floor
        )
        if not 0.0 < log_floor < maximum_probability:
            raise ValueError(
                "--log-floor must be positive and below the global maximum "
                f"probability ({maximum_probability:.6g})"
            )

    lattice_data = probabilities.reshape(
        frame_count, side_length, side_length
    )
    times = np.linspace(0.0, args.tmax, frame_count)

    fig, ax = plt.subplots(figsize=(7, 6), dpi=150)
    initial_data = (
        np.maximum(lattice_data[0], log_floor)
        if log_floor is not None
        else lattice_data[0]
    )
    color_scale = (
        {"norm": LogNorm(log_floor, maximum_probability, clip=True)}
        if log_floor is not None
        else {"vmin": 0.0, "vmax": maximum_probability}
    )
    image = ax.imshow(
        initial_data,
        origin="lower",
        interpolation="nearest",
        cmap="magma",
        extent=(-0.5, side_length - 0.5, -0.5, side_length - 0.5),
        **color_scale,
    )
    colorbar = fig.colorbar(
        image,
        ax=ax,
        extend="min" if log_floor is not None else "neither",
    )
    colorbar.set_label(
        r"$|\psi(x,y)|^2$"
        + (" (log scale)" if log_floor is not None else "")
    )

    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    ax.set_xticks(range(side_length))
    ax.set_yticks(range(side_length))
    ax.set_xticks(np.arange(-0.5, side_length, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, side_length, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.35, alpha=0.22)
    ax.tick_params(which="minor", bottom=False, left=False)
    title = ax.set_title(f"Square-lattice evolution, t = {times[0]:.3f}")

    if log_floor is not None:
        fig.text(
            0.5,
            0.015,
            f"Probabilities at or below {log_floor:.1e} are shown at the color floor.",
            ha="center",
            fontsize=8,
            color="0.35",
        )

    def update(frame):
        frame_data = lattice_data[frame]
        if log_floor is not None:
            frame_data = np.maximum(frame_data, log_floor)
        image.set_data(frame_data)
        title.set_text(f"Square-lattice evolution, t = {times[frame]:.3f}")
        return image, title

    movie = animation.FuncAnimation(
        fig,
        update,
        frames=frame_count,
        interval=args.interval,
        blit=False,
    )

    fig.tight_layout(rect=(0.0, 0.035 if log_floor is not None else 0.0, 1.0, 1.0))
    if args.save:
        movie.save(args.save)
        print(f"Saved animation to {args.save}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
