#!/usr/bin/env python3

import matplotlib.pyplot as plt

from ..plot import comparison_plot_main, add_qcd_value, plot_new_series, iterate_lengths


def plot(data):
    fig, ax = plt.subplots()

    ax.set_xlabel(r"$am_{f}$")
    ax.set_ylabel(r"$F_{\rho} / F_{a_1}$")

    for subset, (colour, marker, label) in iterate_lengths(data):
        plot_new_series(ax, subset, "frho-fa1", colour, marker, label)

    add_qcd_value(ax, "rho_decay_const", "a_1_decay_const")

    ax.set_ylim(0, 2)

    ax.legend(loc="best")
    return fig


if __name__ == "__main__":
    comparison_plot_main(plot, old_data=False)
