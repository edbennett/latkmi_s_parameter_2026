#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd

from ..io import get_samples
from ..stats import jackknife_mean_variance, product_error_contribution
from ..plot import Props, comparison_plot_main, add_qcd_value


def get_unique_masses(data):
    return sorted(set(datum["mass"] for datum in data))


def get_ratio(data, numerator="a_1_mass", denominator="rho_mass"):
    result = []

    masses = get_unique_masses(data)
    for mass in masses:
        numerator_masses = get_samples(data, mass, numerator)
        denominator_masses = get_samples(data, mass, denominator)
        ratio_samples = numerator_masses / denominator_masses
        result.append(jackknife_mean_variance(ratio_samples))

    return masses, *zip(*result)


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


def add_new_data(ax, data, colour, marker, label):
    masses, ratios, errors = get_ratio(data)
    ax.errorbar(
        masses,
        ratios,
        yerr=errors,
        marker=marker,
        color=colour,
        linestyle="none",
        label=label,
    )


def plot(new_data, old_data):
    fig, ax = plt.subplots()

    ax.set_xlabel(r"$am_{f}$")
    ax.set_ylabel(r"$M_{a_1} / M_{\mathrm{\rho}}$")

    lengths = reversed(sorted(set(datum["Nx"] for datum in new_data)))
    props = Props(length_only=True)

    for length in lengths:
        subset = [datum for datum in new_data if datum["Nx"] == length]
        colour, marker, label = props.get(subset[0])
        add_new_data(ax, subset, colour, marker, label)
        add_old_data(ax, subset, old_data, marker)

    add_qcd_value(ax, "a_1_mass", "rho_mass")

    ax.set_ylim(0, 2)

    ax.legend(loc="best")
    return fig


if __name__ == "__main__":
    comparison_plot_main(plot)
