#!/usr/bin/env python3

from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np

from ..io import read_numpy, read_numpy_optional
from ..plot import save_or_show
from .pade import pade, get_momentum_filter


colours = {"Conserved": "C0", "OneLink": "C1"}


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_file", metavar="input_file")
    parser.add_argument("--fit_result", default=None)
    parser.add_argument("--plot_styles", default="styles/prd.mplstyle")
    parser.add_argument("--output_file", default=None)
    parser.add_argument("--renormalised", action="store_true")
    parser.add_argument("--q_squared_upper_bound", type=float, default=None)
    return parser.parse_args()


def plot_fit_result(
    ax, fit_results, colour=None, target_currents=["Conserved", "OneLink"], Z_A=1
):
    xmin, xmax = ax.get_xlim()
    x_range = np.linspace(
        xmin,
        min(xmax, fit_results["pade_fit_result"]["max_momentum_squared"] * 1.05),
        1000,
    )
    for current, result in fit_results["pade_fit_result"].items():
        if current not in target_currents:
            continue
        ax.plot(
            x_range,
            pade(x_range, *result["params"][0]) * Z_A,
            color=colour if colour else colours[current],
        )
    ax.set_xlim(xmin, xmax)


def plot(data, renormalised=False, q_squared_upper_bound=None, fit_result=None):
    fig, ax = plt.subplots()

    key = {True: "renormalised_vpf", False: "vpf"}[renormalised]
    title = {
        True: "Renormalized one-link operator",
        False: "Bare one-link operator",
    }[renormalised]
    labels = {
        True: {
            "Conserved": r"$Z_A \times \textnormal{Conserved--OneLink}$",
            "OneLink": r"$Z_A^2 \times \textnormal{OneLink--OneLink}$",
        },
        False: {
            "Conserved": "$Conserved--OneLink",
            "OneLink": "$OneLink--OneLink",
        },
    }
    momentum_filter = get_momentum_filter(
        data["reordered_momentum"],
        [data[key] for key in ["Nx", "Ny", "Nz", "Nt"]],
    )

    for current, marker in [("Conserved", "o"), ("OneLink", "s")]:
        vpf = data[key][current]
        ax.errorbar(
            data["momentum_squared"][momentum_filter],
            vpf[0][momentum_filter],
            yerr=vpf[1][momentum_filter],
            linestyle="none",
            marker=marker,
            label=labels[renormalised][current],
            color=colours[current],
        )

    assert data["Nx"] == data["Ny"] == data["Nz"]
    ax.set_title(rf"$m_f = {data['mass']}$, ${data['Nx']}^3 \times {data['Nt']}$")
    ax.set_ylabel(r"$\Pi^{V-A}(q^2)$")
    ax.set_xlabel("$q^2$")
    ax.set_xlim(0, q_squared_upper_bound)

    # Ensure that full fit result fits on plot vertically,
    # but scan full horizontal range,
    # so do this after x-limit but before y-limit
    if fit_result:
        plot_fit_result(ax, fit_result)

    ax.set_ylim(None, 0)

    fig.suptitle(title)

    return fig


def main():
    args = get_args()
    plt.style.use(args.plot_styles)
    data = read_numpy(args.input_file)
    fit_result = read_numpy_optional(args.fit_result)

    fig = plot(data, args.renormalised, args.q_squared_upper_bound, fit_result)
    save_or_show(fig, args.output_file)


if __name__ == "__main__":
    main()
