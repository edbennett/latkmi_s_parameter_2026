#!/usr/bin/env python3

from argparse import ArgumentParser

import matplotlib.pyplot as plt

from ..io import read_numpy
from ..plot import save_or_show
from ..stats import add_quadrature


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_file", nargs="+", metavar="input_file")
    parser.add_argument("--plot_styles", default="styles/prd.mplstyle")
    parser.add_argument("--output_file", default=None)
    return parser.parse_args()


def symmetrise_y_axis(ax):
    initial_ymin, initial_ymax = ax.get_ylim()
    ymax = max(abs(initial_ymin), abs(initial_ymax))
    ax.set_ylim(-ymax, ymax)


def plot(data):
    fig, ax = plt.subplots()

    ax.set_xlabel("$am_f$")
    ax.set_ylabel(r"$\overline{S}_{\mathrm{IR}}$")

    masses = [datum["mass"] for datum in data]
    sign_names = {"+": "plus", "-": "minus"}
    for channel, parity, marker in [
        ("V", "+", "s"),
        ("A", "-", "o"),
        ("A", "+", "^"),
        ("V", "-", "v"),
    ]:
        values_errors = [datum[f"S_{channel}_{sign_names[parity]}"] for datum in data]
        values = [contributions[0] for contributions in values_errors]
        combined_errors = [add_quadrature(*errors) for _, *errors in values_errors]
        label = (
            rf"${{\overline{{S}}_\mathrm{{IR}}^{{\mathrm{{{channel}}}}}}}^{{{parity}}}$"
        )
        ax.errorbar(
            masses,
            values,
            yerr=combined_errors,
            label=label,
            linestyle="none",
            marker=marker,
        )

    ax.set_xlim(0, None)
    ax.axhline(0, color="black")
    ax.legend(loc="best", ncol=4)

    symmetrise_y_axis(ax)
    return fig


def main():
    args = get_args()
    plt.style.use(args.plot_styles)
    data = [read_numpy(filename) for filename in args.input_file]

    fig = plot(data)
    save_or_show(fig, args.output_file)


if __name__ == "__main__":
    main()
