#!/usr/bin/env python3

from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np

from ..io import read_numpy
from ..plot import save_or_show


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_file", metavar="input_file")
    parser.add_argument("--plot_styles", default="styles/prd.mplstyle")
    parser.add_argument("--output_file", default=None)
    return parser.parse_args()


def plot(data):
    fig, ax = plt.subplots(layout="constrained", figsize=(3.4, 3))

    ax.set_yscale("log")
    ax.set_xlabel("$t$")
    ax.set_ylabel(r"$-\hat{G}$")

    for channel, marker in [("V", "o"), ("A", "^"), ("V-A", "s")]:
        value, error = data["Conserved"][f"{channel}_renormalised"]
        ax.errorbar(
            np.arange(data["Nt"] // 2 + 1),
            -value,  # Correlators come in negative; negate them to show on a log scale
            error,
            marker=marker,
            linestyle="none",
            label=rf"$\hat{{G}}^{{{channel}}}$",
        )

    ax.legend(loc="best")
    return fig


def main():
    args = get_args()
    plt.style.use(args.plot_styles)
    data = read_numpy(args.input_file)

    fig = plot(data)
    save_or_show(fig, args.output_file)


if __name__ == "__main__":
    main()
