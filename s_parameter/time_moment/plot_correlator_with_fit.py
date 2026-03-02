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


def plot_fit(ax, fit_result, colour="C0", label="Fit"):
    single_result = fit_result["S_correlator_interpolation"]
    correlator_fit, correlator_fit_error = single_result["correlator"]
    ax.fill_between(
        single_result["times"],
        correlator_fit - correlator_fit_error,
        correlator_fit + correlator_fit_error,
        color=colour,
        label=label,
        alpha=0.5,
    )


def plot(data, fit_result):
    fig, ax = plt.subplots(layout="constrained", figsize=(3.4, 3))

    ax.set_yscale("log")
    ax.set_xlabel("$t$")
    ax.set_ylabel(r"$-\hat{G}^{V-A}(t)$")

    value, error = data["Conserved"]["V-A_renormalised"]
    ax.errorbar(
        np.arange(data["Nt"] // 2 + 1),
        -value,  # Correlators come in negative; negate them to show on a log scale
        error,
        marker="s",
        linestyle="none",
        label="Data",
        color="C2",
    )

    num_timeslices = data["Conserved"]["V-A_renormalised"].shape[-1] - 1
    max_timeslice = 5 * num_timeslices / 4
    # plot_fit(ax, fit_result, "C1", r"$T\rightarrow\infty$", max_timeslice=max_timeslice)
    plot_fit(ax, fit_result, "C0", "Fit")
    ax.set_xlim(None, max_timeslice)

    ax.legend(loc="best")
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
