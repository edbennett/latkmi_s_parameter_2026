#!/usr/bin/env python3

from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np

from ..plot import save_or_show
from ..io import read_numpy


colours_markers = {
    18: ("C5", "p"),
    24: ("C4", "D"),
    30: ("C3", "v"),
    36: ("C2", "^"),
    42: ("C1", "o"),
    48: ("C0", "s"),
}


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_files", nargs="+", metavar="input_file")
    parser.add_argument("--plot_styles", default="styles/prd.mplstyle")
    parser.add_argument("--output_file", default=None)
    return parser.parse_args()


def plot_single_series(
    ax,
    data,
    offset,
    lattice_size,
    upper_bound,
    colour,
    fillstyle,
    label_suffix,
    datum_getter,
):
    default_colour, marker = colours_markers[lattice_size]
    volume_data = [
        datum
        for datum in data
        if datum["Nx"] == lattice_size and datum["upper_bound"] == upper_bound
    ]
    if not volume_data:
        return

    assert all(
        datum["Ny"] == lattice_size and datum["Nz"] == lattice_size
        for datum in volume_data
    )
    masses = [datum["mass"] + offset for datum in volume_data]
    value, error, _ = map(
        np.array, zip(*[datum_getter(datum) for datum in volume_data])
    )
    ax.errorbar(
        masses,
        value,
        yerr=error,
        marker=marker,
        linestyle="none",
        label=None if label_suffix is None else f"$L = {lattice_size}${label_suffix}",
        fillstyle=fillstyle,
        color=colour or default_colour,
    )


def get_S(datum):
    return datum["pade_fit_result"]["Conserved"]["S"]


def get_suffix(upper_bound_key):
    if upper_bound_key.isdigit():
        return f", $q^2 < {upper_bound_key}$"
    return rf", $q^2 < q^2_{{\mathrm{{{upper_bound_key}}}}}$"


def plot(data, datum_getter, ylabel, ylim, yscale):
    fig, ax = plt.subplots()
    upper_bounds = sorted(set(datum["upper_bound"] for datum in data))
    lattice_sizes = reversed(sorted(set(datum["Nx"] for datum in data)))
    if len(upper_bounds) <= 2:
        suffixes = ["", None]
    else:
        suffixes = [get_suffix(upper_bound_key) for upper_bound_key in upper_bounds]

    for lattice_size in lattice_sizes:
        for upper_bound, colour, fillstyle, offset, label_suffix in zip(
            upper_bounds,
            [None, None, "grey"],
            ["none", "full", "full"],
            [0, 0.001, -0.001],
            suffixes,
        ):
            plot_single_series(
                ax,
                data,
                offset,
                lattice_size,
                upper_bound,
                colour,
                fillstyle,
                label_suffix,
                datum_getter,
            )

    ax.set_xlabel("$am_f$")
    ax.set_ylabel(ylabel)
    ax.set_yscale(yscale)
    ax.legend(loc="best")

    if ax.get_xlim()[0] < 0.03:
        ax.set_xlim(0, None)
    if len(upper_bounds) > 1:
        ax.set_ylim(ylim)
    else:
        ax.set_ylim(0, 0.35)

    return fig


def main(datum_getter, ylabel, ylim, yscale="linear"):
    args = get_args()
    plt.style.use(args.plot_styles)
    data = [read_numpy(filename) for filename in args.input_files]

    fig = plot(data, datum_getter, ylabel, ylim, yscale)
    save_or_show(fig, args.output_file)


if __name__ == "__main__":
    main(get_S, r"$S|_{\textnormal{1-doublet}}$", (0.2, 0.3))
