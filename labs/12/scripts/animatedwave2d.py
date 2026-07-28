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
import numpy as np


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

    lattice_data = probabilities.reshape(
        frame_count, side_length, side_length
    )
    times = np.linspace(0.0, args.tmax, frame_count)

    fig, ax = plt.subplots(figsize=(7, 6), dpi=150)
    image = ax.imshow(
        lattice_data[0],
        origin="lower",
        interpolation="nearest",
        cmap="viridis",
        vmin=0.0,
        vmax=probabilities.max(),
        extent=(-0.5, side_length - 0.5, -0.5, side_length - 0.5),
    )
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label(r"$|\psi(x,y)|^2$")

    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    ax.set_xticks(range(side_length))
    ax.set_yticks(range(side_length))
    title = ax.set_title(f"Square-lattice evolution, t = {times[0]:.3f}")

    def update(frame):
        image.set_data(lattice_data[frame])
        title.set_text(f"Square-lattice evolution, t = {times[frame]:.3f}")
        return image, title

    movie = animation.FuncAnimation(
        fig,
        update,
        frames=frame_count,
        interval=args.interval,
        blit=False,
    )

    fig.tight_layout()
    if args.save:
        movie.save(args.save)
        print(f"Saved animation to {args.save}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
