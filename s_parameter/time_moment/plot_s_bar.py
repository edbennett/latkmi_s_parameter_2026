#!/usr/bin/env python3

#!/usr/bin/env python3

from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np

from ..io import read_numpy
from ..plot import save_or_show


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_file", metavar="input_file")
    parser.add_argument("--fit_result", required=True)
    parser.add_argument("--plot_styles", default="styles/prd.mplstyle")
    parser.add_argument("--output_file", default=None)
    return parser.parse_args()


def plot_single_fit(ax, data, label, colour):
    time = data["times"]
    s_eff, s_eff_error = data["S_effective"]
    ax.fill_between(
        time,
        s_eff + s_eff_error,
        s_eff - s_eff_error,
        label=label,
        color=colour,
        alpha=0.5,
    )


def plot_single_data(ax, data, label, colour):
    s_eff, s_eff_error = data
    time = np.arange(len(s_eff))
    ax.errorbar(
        time,
        s_eff,
        yerr=s_eff_error,
        linestyle="none",
        label=label,
        marker="o",
        color=colour,
    )


def plot(data, fit_result):
    fig, ax = plt.subplots()

    ax.set_xlabel(r"$\overline{t}$")
    ax.set_ylabel(r"$\overline{S}(\overline{t})$")

    plot_single_data(ax, data["Conserved"]["S_parameter_eff"], "Data", "C0")
    plot_single_fit(ax, fit_result["S_correlator_interpolation"], "Fit", "C1")
    plot_single_fit(
        ax, fit_result["S_extrapolation_large_t"], r"$T \rightarrow \infty$", "C2"
    )

    ax.set_xlim(0, data["Nt"])
    ax.set_ylim(0, 0.3)
    ax.legend(loc="best")
    ax.grid()

    return fig


def main():
    args = get_args()
    plt.style.use(args.plot_styles)
    data = read_numpy(args.input_file)
    fit_result = read_numpy(args.fit_result)

    fig = plot(data, fit_result)
    save_or_show(fig, args.output_file)


if __name__ == "__main__":
    main()
