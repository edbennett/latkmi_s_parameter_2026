#!/usr/bin/env python3

from argparse import ArgumentParser

import matplotlib as mpl
import matplotlib.pyplot as plt

from ..plot import save_or_show
from .plot_vpf import read, plot_fit_result


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_files", nargs="+", metavar="input_file")
    parser.add_argument(
        "--fit_result", action="append", dest="fit_results", default=None
    )
    parser.add_argument("--plot_styles", default="styles/prd.mplstyle")
    parser.add_argument("--output_file", default=None)
    parser.add_argument("--q_squared_upper_bound", type=float, default=None)
    return parser.parse_args()


class Props:
    _markers = ["o", "s", "^", "v", "<", "H", ">", "+", "x", "D", "h", "p", "8", "*"]
    _colours = mpl.color_sequences["tab10"] + ["black", "seagreen", "violet"]

    def __init__(self):
        self._props = {}

    def get(self, datum):
        lattice_size = datum["Nx"]
        assert lattice_size == datum["Ny"] and lattice_size == datum["Nz"]
        mass = datum["mass"]
        if (mass, lattice_size) not in self._props:
            self._props[mass, lattice_size] = (
                self._colours.pop(0),
                self._markers.pop(0),
                f"$am_f={mass}$, $L={lattice_size}$",
            )

        return self._props[mass, lattice_size]


def plot_single_datum(ax, datum, props):
    colour, marker, label = props.get(datum)
    vpf = datum["renormalised_vpf"]["Conserved"]
    ax.errorbar(
        datum["momentum_squared"],
        vpf[0],
        yerr=vpf[1],
        linestyle="none",
        marker=marker,
        label=label,
        color=colour,
    )


def plot_single_fit(ax, datum, props):
    colour, _, _ = props.get(datum)
    plot_fit_result(ax, datum, colour=colour, target_currents=["Conserved"])


def plot(data=[], fit_results=[], q_squared_upper_bound=None):
    fig, ax = plt.subplots(layout="constrained", figsize=(3.4, 3))
    props = Props()

    for fit_result in fit_results:
        plot_single_fit(ax, fit_result, props)

    for datum in data:
        plot_single_datum(ax, datum, props)

    ax.set_ylabel(r"$\Pi(q^2)$")
    ax.set_xlabel("$q^2$")
    ax.set_ylim(None, 0)
    ax.set_xlim(0, q_squared_upper_bound)
    ax.legend(loc="best")

    return fig


def main():
    args = get_args()
    plt.style.use(args.plot_styles)
    data = [read(input_file) for input_file in args.input_files]
    fit_results = [read(fit_result) for fit_result in args.fit_results]

    fig = plot(
        data,
        fit_results,
        args.q_squared_upper_bound,
    )
    save_or_show(fig, args.output_file)


if __name__ == "__main__":
    main()
