#!/usr/bin/env python3

import matplotlib.pyplot as plt

from ..plot import comparison_plot_main, add_qcd_value, plot_new_series, iterate_lengths


ylabels = {
    "frho-fpi": r"$F_{\rho} / F_{\pi}$",
    "frho-fa1": r"$F_{\rho} / F_{a_1}$",
    "wsr-i-normalised": (
        r"$(F_{\rho}^2 - F_{a_1}^2 - F_{\pi}^2)"
        r"/ (F_{\rho}^2 + F_{a_1}^2 + F_{\pi}^2)$"
    ),
    "wsr-ii-normalised": (
        r"$(F_{\rho}^2 M_{\rho}^2 - F_{a_1}^2 M_{a_1}^2)"
        r"/ (F_{\rho}^2 M_{\rho}^2 + F_{a_1}^2 M_{a_1}^2)$"
    ),
}
ylims = {
    "frho-fpi": (0, 2),
    "frho-fa1": (0, 2),
    "wsr-i-normalised": (-1, 1),
    "wsr-ii-normalised": (-1, 1),
}


def plot(data, plot_type):
    fig, ax = plt.subplots()

    ax.set_xlabel(r"$am_{f}$")
    ax.set_ylabel(ylabels[plot_type])

    for subset, (colour, marker, label) in iterate_lengths(data):
        plot_new_series(ax, subset, plot_type, colour, marker, label)

    add_qcd_value(ax, plot_type)

    ax.set_ylim(*ylims[plot_type])
    ax.legend(loc="best")

    return fig


if __name__ == "__main__":
    comparison_plot_main(plot, select_plot=True)
