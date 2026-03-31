#!/usr/bin/env python3

import matplotlib.pyplot as plt

from ..plot import comparison_plot_main, plot_new_series, iterate_lengths


def plot(data):
    fig, ax = plt.subplots()

    ax.set_xlabel(r"$am_{f}$")
    ax.set_ylabel(r"$g_{\rho\pi\pi}^{\mathrm{KSRF}}$")

    for subset, (colour, marker, label) in iterate_lengths(data):
        plot_new_series(ax, subset, "ksrf-i", colour, marker, label)
        plot_new_series(ax, subset, "ksrf-ii", "gray", marker, label)

    ax.set_xlim(0, None)
    ax.set_ylim(4, 7)

    return fig


if __name__ == "__main__":
    comparison_plot_main(plot)
