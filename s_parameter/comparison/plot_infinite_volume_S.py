#!/usr/bin/env python3

from argparse import ArgumentParser

import matplotlib.pyplot as plt

from ..io import read_numpy
from ..plot import save_or_show
from ..utils import nested_get


def get_args():
    parser = ArgumentParser()
    parser.add_argument(
        "--time_moment_infinite_volume", metavar="input_file", nargs="+"
    )
    parser.add_argument("--time_moment_finite_volume", metavar="input_file", nargs="+")
    parser.add_argument("--vp_momentum_finite_volume", metavar="input_file", nargs="+")
    parser.add_argument("--plot_styles", default="styles/prd.mplstyle")
    parser.add_argument("--output_file", default=None)
    return parser.parse_args()


def plot(
    time_moment_infinite_volume, time_moment_finite_volume, vp_momentum_finite_volume
):
    fig, ax = plt.subplots()

    ax.set_xlabel("$am_f$")
    ax.set_ylabel(r"$S|_{\textnormal{\scriptsize{1-doublet}}}$")

    for data, value_key, label, marker, offset in (
        (
            time_moment_infinite_volume,
            "S_infinite_volume",
            r"Time-moment $L=\infty$",
            "s",
            0.0,
        ),
        (
            time_moment_finite_volume,
            "S_infinite_t",
            r"Time-moment finite volume",
            "o",
            -0.0005,
        ),
        (
            vp_momentum_finite_volume,
            ("pade_fit_result", "Conserved", "S"),
            "VP-momentum finite volume",
            "^",
            0.0005,
        ),
    ):
        ax.errorbar(
            [datum["mass"] + offset for datum in data],
            [nested_get(datum, value_key)[0] for datum in data],
            yerr=[nested_get(datum, value_key)[1] for datum in data],
            marker=marker,
            label=label,
            linestyle="none",
        )

    ax.set_xlim(0, None)
    ax.set_ylim(0.2, 0.32)
    ax.legend(loc="best")

    return fig


def main():
    args = get_args()
    plt.style.use(args.plot_styles)

    time_moment_infinite_volume = [
        read_numpy(input_file) for input_file in args.time_moment_infinite_volume
    ]
    time_moment_finite_volume = [
        read_numpy(input_file) for input_file in args.time_moment_finite_volume
    ]
    vp_momentum_finite_volume = [
        read_numpy(input_file) for input_file in args.vp_momentum_finite_volume
    ]

    fig = plot(
        time_moment_infinite_volume,
        time_moment_finite_volume,
        vp_momentum_finite_volume,
    )
    save_or_show(fig, args.output_file)


if __name__ == "__main__":
    main()
