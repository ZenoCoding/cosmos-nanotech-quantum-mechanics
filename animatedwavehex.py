#
# Animate probability data produced by timehex.c.
#
# Run interactively:
#
#     python3 animatedwavehex.py
#
# Or save the animation:
#
#     python3 animatedwavehex.py --save hex.gif
#
# Each input row is one time slice. Sites are ordered as in timehex.c:
#
#     (0,0), (1,0), ..., (L-1,0), (0,1), ..., (L-1,L-1).
#
import argparse
import math

import matplotlib.animation as animation
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Animate probability on a honeycomb lattice."
    )
    parser.add_argument(
        "input",
        nargs="?",
        default="timehex.txt",
        help="time-evolution data file (default: timehex.txt)",
    )
    parser.add_argument(
        "--tmax",
        type=float,
        default=8.0 * math.pi,
        help="final simulation time shown in the title (default: 8*pi)",
    )
    parser.add_argument(
        "--interval",
        type=int,
        default=100,
        help="milliseconds between frames (default: 100)",
    )
    parser.add_argument(
        "--fixed-scale",
        action="store_true",
        help="keep one color scale for all frames instead of rescaling each frame",
    )
    parser.add_argument(
        "--save",
        metavar="FILE",
        help="save the animation instead of opening a window",
    )
    return parser.parse_args()


def lattice_geometry(side_length):
    """Return coordinates and bonds drawn as a regular honeycomb."""
    coordinates = []
    row_spacing = math.sqrt(3.0) / 2.0

    for y in range(side_length):
        for x in range(side_length):
            # Along even rows, bonded pairs lie at (0, 1), (3, 4), ...
            # Odd rows shift alternate sites outward. Consequently, the two
            # vertical brick-wall edges become the sloping sides of a regular
            # hexagon while every nearest-neighbor bond retains unit length.
            display_x = x + x // 2
            if y % 2 == 1:
                display_x += -0.5 if x % 2 == 0 else 0.5
            coordinates.append((display_x, row_spacing * y))

    coordinates = np.asarray(coordinates, dtype=float)
    bonds = []

    for y in range(side_length):
        for x in range(side_length):
            site = x + side_length * y

            # Every site is connected vertically.
            if y + 1 < side_length:
                above = site + side_length
                bonds.append((coordinates[site], coordinates[above]))

            # Horizontal bonds alternate from one row to the next.
            if x + 1 < side_length and (x + y) % 2 == 0:
                right = site + 1
                bonds.append((coordinates[site], coordinates[right]))

    return coordinates, bonds


def main():
    args = parse_arguments()
    probabilities = np.loadtxt(args.input, ndmin=2)

    frame_count, site_count = probabilities.shape
    side_length = math.isqrt(site_count)
    if side_length * side_length != site_count:
        raise ValueError(
            f"{args.input} has {site_count} columns; "
            "the number of sites must be a perfect square"
        )
    if args.tmax < 0:
        raise ValueError("--tmax must be nonnegative")
    if args.interval < 1:
        raise ValueError("--interval must be at least 1 millisecond")

    coordinates, bonds = lattice_geometry(side_length)
    times = np.linspace(0.0, args.tmax, frame_count)
    maximum_probability = max(float(probabilities.max()), 1.0e-15)
    initial_maximum = max(float(probabilities[0].max()), 1.0e-15)

    fig, ax = plt.subplots(figsize=(7, 7), dpi=150)
    ax.add_collection(
        LineCollection(bonds, colors="0.72", linewidths=1.5, zorder=1)
    )
    sites = ax.scatter(
        coordinates[:, 0],
        coordinates[:, 1],
        c=probabilities[0],
        cmap="viridis",
        vmin=0.0,
        vmax=maximum_probability if args.fixed_scale else initial_maximum,
        s=170,
        edgecolors="0.15",
        linewidths=0.5,
        zorder=2,
    )
    colorbar = fig.colorbar(sites, ax=ax)
    colorbar.set_label(r"$|\psi(x,y)|^2$")

    margin = 0.7
    ax.set_xlim(coordinates[:, 0].min() - margin,
                coordinates[:, 0].max() + margin)
    ax.set_ylim(coordinates[:, 1].min() - margin,
                coordinates[:, 1].max() + margin)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    ax.set_xticks([])
    ax.set_yticks([])
    scale_name = "fixed scale" if args.fixed_scale else "adaptive scale"
    title = ax.set_title(
        "Honeycomb-lattice evolution\n"
        f"t = {times[0]:.3f}, max probability = {initial_maximum:.3f} "
        f"({scale_name})"
    )

    def update(frame):
        frame_probabilities = probabilities[frame]
        frame_maximum = max(float(frame_probabilities.max()), 1.0e-15)
        sites.set_array(frame_probabilities)
        if not args.fixed_scale:
            sites.set_clim(0.0, frame_maximum)
            colorbar.update_normal(sites)
        title.set_text(
            "Honeycomb-lattice evolution\n"
            f"t = {times[frame]:.3f}, max probability = "
            f"{frame_maximum:.3f} ({scale_name})"
        )
        return sites, title

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
