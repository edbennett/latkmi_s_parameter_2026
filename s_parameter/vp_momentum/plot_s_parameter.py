#!/usr/bin/env python3

from argparse import ArgumentParser

import matplotlib.pyplot as plt

from ..plot import save_or_show
from .plot_vpf import read


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_files", nargs="+", metavar="input_file")
    parser.add_argument("--plot_styles", default="styles/prd.mplstyle")
    parser.add_argument("--output_file", default=None)
    return parser.parse_args()


def plot(data):
    fig, ax = plt.subplots(layout="constrained", figsize=(3.4, 3))

    for lattice_size, marker in zip(
        sorted(set(datum["Nx"] for datum in data)), "os^vDp"
    ):
        volume_data = [datum for datum in data if datum["Nx"] == lattice_size]
        assert all(
            datum["Ny"] == lattice_size and datum["Nz"] == lattice_size
            for datum in volume_data
        )
        masses = [datum["mass"] for datum in volume_data]
        s_parameter, s_parameter_error = zip(
            *[datum["pade_fit_result"]["Conserved"]["S"] for datum in volume_data]
        )
        ax.errorbar(
            masses,
            s_parameter,
            yerr=s_parameter_error,
            marker=marker,
            linestyle="none",
            label=f"$L = {lattice_size}$",
        )

    ax.set_xlabel("$m_f$")
    ax.set_ylabel(r"$S|_{\textnormal{1-doublet}}$")
    ax.legend(loc="best")

    ax.set_xlim(0, None)
    ax.set_ylim(0, None)

    return fig


def main():
    args = get_args()
    plt.style.use(args.plot_styles)
    data = [read(filename) for filename in args.input_files]

    fig = plot(data)
    save_or_show(fig, args.output_file)


if __name__ == "__main__":
    main()
