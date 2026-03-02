#!/usr/bin/env python3

from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np

from ..io import read_numpy_optional
from ..plot import Props, save_or_show
from .pade import get_momentum_filter
from .plot_vpf import read, plot_fit_result


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_files", nargs="+", metavar="input_file")
    parser.add_argument(
        "--fit_result", action="append", dest="fit_results", default=None
    )
    parser.add_argument("--Z_A", default=None)
    parser.add_argument("--plot_styles", default="styles/prd.mplstyle")
    parser.add_argument("--output_file", default=None)
    parser.add_argument("--q_squared_upper_bound", type=float, default=None)
    return parser.parse_args()


def plot_single_datum(ax, datum, props, Z_A):
    colour, marker, label = props.get(datum)
    momentum_filter = get_momentum_filter(
        datum["reordered_momentum"],
        [datum[key] for key in ["Nx", "Ny", "Nz", "Nt"]],
    )
    vpf = datum["vpf"]["Conserved"]
    Z_A_value, Z_A_error = Z_A["central"]["Z_A_0"] if Z_A else (1, 0)
    vpf_value = vpf[0] * Z_A_value
    vpf_error = (
        np.abs(vpf_value)
        * ((vpf[1] / vpf[0]) ** 2 + (Z_A_error / Z_A_value) ** 2) ** 0.5
    )
    ax.errorbar(
        datum["momentum_squared"][momentum_filter],
        vpf_value[momentum_filter],
        yerr=vpf_error[momentum_filter],
        linestyle="none",
        marker=marker,
        label=label,
        color=colour,
    )


def plot_single_fit(ax, datum, props, Z_A):
    colour, _, _ = props.get(datum)
    plot_fit_result(ax, datum, colour=colour, target_currents=["Conserved"], Z_A=Z_A)


def plot(data=[], fit_results=[], Z_A=None, q_squared_upper_bound=None):
    fig, ax = plt.subplots(layout="constrained", figsize=(3.4, 3))
    props = Props()

    for fit_result in fit_results:
        plot_single_fit(ax, fit_result, props, Z_A["central"]["Z_A_0"][0])

    for datum in data:
        plot_single_datum(ax, datum, props, Z_A)

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
    Z_A = read_numpy_optional(args.Z_A)

    fig = plot(
        data,
        fit_results,
        Z_A,
        args.q_squared_upper_bound,
    )
    save_or_show(fig, args.output_file)


if __name__ == "__main__":
    main()
