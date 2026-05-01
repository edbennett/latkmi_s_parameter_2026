#!/usr/bin/env python3

"""
Plot comparing the time moment method (coloured points)
with the VP-momentum method (grey points)
"""

from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np

from ..io import read_numpy
from ..plot import Props, save_or_show
from ..utils import nested_get


def get_args():
    parser = ArgumentParser()
    parser.add_argument("primary_input_files", nargs="+", metavar="input_file")
    parser.add_argument(
        "--secondary_input_file",
        action="append",
        dest="secondary_input_files",
        metavar="input_file",
    )
    parser.add_argument("--plot_styles", default="styles/prd.mplstyle")
    parser.add_argument("--output_file", default=None)
    return parser.parse_args()


def plot_single_series(ax, data, secondary_data, key, props):
    """
    Plot the given data on the given ax.
    Access the relevant quantity from each datum using key.
    Use the given Props instance to decide the colour, label, and marker,
    unless secondary_data = True,
    in which case the colour is always grey and the label not set.
    (The marker is still taken from props.)
    Shift secondary_data horizontally by a small amount to allow easier comparison.
    """
    sizes = reversed(sorted(set(datum["Nx"] for datum in data)))
    for size in sizes:
        subset = [datum for datum in data if datum["Nx"] == size]
        mass_offset = 0.001 if secondary_data else 0
        masses = [datum["mass"] + mass_offset for datum in subset]

        results = np.array([nested_get(datum, key) for datum in subset])
        values = results[:, 0]
        errors = (results[:, 1:] ** 2).sum(axis=1) ** 0.5

        colour, marker, label = props.get(subset[0])
        ax.errorbar(
            masses,
            values,
            errors,
            color="grey" if secondary_data else colour,
            marker=marker,
            label=None if secondary_data else label,
            linestyle="none",
        )


def plot(primary, secondary):
    fig, ax = plt.subplots()
    props = Props(length_only=True)

    plot_single_series(ax, primary, False, "S_infinite_t", props)
    plot_single_series(
        ax, secondary, True, ["pade_fit_result", "Conserved", "S"], props
    )

    ax.set_xlabel("$am_f$")
    ax.set_ylabel(r"$S|_{\textnormal{1-doublet}}$")

    ax.set_ylim(0.2, 0.3)
    ax.set_xlim(0, None)

    ax.legend(loc="best")
    return fig


def main():
    args = get_args()
    plt.style.use(args.plot_styles)
    primary = [read_numpy(input_file) for input_file in args.primary_input_files]
    secondary = [read_numpy(input_file) for input_file in args.secondary_input_files]

    save_or_show(plot(primary, secondary), args.output_file)


if __name__ == "__main__":
    main()
