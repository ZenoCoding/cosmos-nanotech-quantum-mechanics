"""Display a density-of-states histogram from a file of eigenvalues."""

import argparse

import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot a histogram of eigenvalues."
    )
    parser.add_argument(
        "filename",
        help="text file containing one real eigenvalue per line",
    )
    parser.add_argument(
        "bin_width",
        type=float,
        help="energy width of each histogram bin (for example, 0.1)",
    )
    parser.add_argument(
        "--save",
        metavar="IMAGE",
        help="save the plot to an image instead of opening the viewer",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if args.bin_width <= 0.0:
        raise SystemExit("bin_width must be greater than zero")

    eigenvalues = np.atleast_1d(np.loadtxt(args.filename, dtype=float))
    eigenvalues = eigenvalues[np.isfinite(eigenvalues)]
    if eigenvalues.size == 0:
        raise SystemExit("the input file contains no finite eigenvalues")

    minimum = eigenvalues.min()
    maximum = eigenvalues.max()
    if minimum == maximum:
        edges = np.array([
            minimum - args.bin_width / 2.0,
            maximum + args.bin_width / 2.0,
        ])
    else:
        start = np.floor(minimum / args.bin_width) * args.bin_width
        stop = np.ceil(maximum / args.bin_width) * args.bin_width
        edges = np.arange(start, stop + args.bin_width, args.bin_width)

    fig, axis = plt.subplots()
    axis.hist(eigenvalues, bins=edges, edgecolor="black")
    axis.tick_params(direction="out")
    axis.set_xlabel("Energy, E")
    axis.set_ylabel("Number of states")
    axis.set_title("Density of States")
    fig.tight_layout()

    if args.save:
        fig.savefig(args.save, dpi=160)
        print(f"Saved histogram to {args.save}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
