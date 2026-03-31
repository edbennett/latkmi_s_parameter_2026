#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd

from ..plot import comparison_plot_main, iterate_lengths, plot_new_series, add_qcd_value


def get_unique_masses(data):
    return sorted(set(datum["mass"] for datum in data))


def product_error_contribution(data, *keys):
    return (
        sum([(data[f"error_{key}"] / data[f"value_{key}"]) ** 2 for key in keys]).sum()
        ** 0.5
    )


def add_old_data(ax, target_ensembles, source_data, marker):
    matching_data = []
    (spatial_size,) = list(set(datum["Nx"] for datum in target_ensembles))
    (temporal_size,) = list(set(datum["Nt"] for datum in target_ensembles))
    masses = get_unique_masses(target_ensembles)

    for mass in masses:
        result = source_data.query(
            "Nf == 8 & beta == 3.8 & "
            f"L == {spatial_size} & T == {temporal_size} & mf == {mass}"
        )
        if len(result.index) > 1:
            raise ValueError("Multiple ensembles found.")
        matching_data.append(result)

    data = pd.concat(matching_data)
    ratio = data["value_ma1"] / data["value_mrho"]
    ax.errorbar(
        data.mf + 0.001,
        ratio,
        yerr=ratio * product_error_contribution(data, "ma1", "mrho"),
        marker=marker,
        color="gray",
        linestyle="none",
    )


def plot(new_data, old_data):
    fig, ax = plt.subplots()

    ax.set_xlabel(r"$am_{f}$")
    ax.set_ylabel(r"$M_{a_1} / M_{\mathrm{\rho}}$")

    for subset, (colour, marker, label) in iterate_lengths(new_data):
        plot_new_series(ax, subset, "ma1-mrho", colour, marker, label)
        add_old_data(ax, subset, old_data, marker)

    add_qcd_value(ax, "a_1_mass", "rho_mass")

    ax.set_ylim(0, 2)

    ax.legend(loc="best")
    return fig


if __name__ == "__main__":
    comparison_plot_main(plot, old_data=True)
