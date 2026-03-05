#!/usr/bin/env python3

import matplotlib.pyplot as plt

from ..plot import Props, comparison_plot_main, add_qcd_value

from .plot_a_1_rho_ratio import get_ratio


def add_data(ax, data, colour, marker, label):
    masses, ratios, errors = get_ratio(data, "rho_decay_const", "a_1_decay_const")
    ax.errorbar(
        masses,
        ratios,
        yerr=errors,
        marker=marker,
        color=colour,
        linestyle="none",
        label=label,
    )


def plot(data):
    fig, ax = plt.subplots()

    ax.set_xlabel(r"$am_{f}$")
    ax.set_ylabel(r"$F_{\rho} / F_{a_1}$")

    lengths = reversed(sorted(set(datum["Nx"] for datum in data)))
    props = Props(length_only=True)

    for length in lengths:
        subset = [datum for datum in data if datum["Nx"] == length]
        colour, marker, label = props.get(subset[0])
        add_data(ax, subset, colour, marker, label)

    add_qcd_value(ax, "rho_decay_const", "a_1_decay_const")

    ax.set_ylim(0, 2)

    ax.legend(loc="best")
    return fig


if __name__ == "__main__":
    comparison_plot_main(plot, old_data=False)
