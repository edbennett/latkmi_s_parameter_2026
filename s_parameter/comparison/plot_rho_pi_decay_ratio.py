#!/usr/bin/env python3

import matplotlib.pyplot as plt

from ..stats import jackknife_mean_variance, generate_jackknife
from ..plot import Props, comparison_plot_main, add_qcd_value


def compute_ratio_with_old(datum, new_key, old_data, old_key):
    result = old_data.query(
        "Nf == 8 & beta == 3.8 & "
        f"L == {datum['Nx']} & T == {datum['Nt']} & mf == {datum['mass']}"
    )
    assert datum["Nx"] == datum["Ny"] and datum["Nx"] == datum["Nz"]
    if len(result.index) > 1:
        raise ValueError("Multiple ensembles found.")
    if len(result.index) == 0:
        raise ValueError("No old data found")
    denominator_value = result[f"value_{old_key}"].iloc[0]
    denominator_error = result[f"error_{old_key}"].iloc[0]

    numerator_samples = datum["fit_result_samples"][new_key]
    denominator_samples = generate_jackknife(
        denominator_value,
        denominator_error,
        datum,
        len(numerator_samples),
    )
    ratio_samples = numerator_samples / denominator_samples
    return jackknife_mean_variance(ratio_samples)


def plot_ratio_with_old(ax, props, new_data, new_key, old_data, old_key):
    colour, marker, label = props.get(new_data[0])
    masses = [datum["mass"] for datum in new_data]
    values, errors = zip(
        *[
            compute_ratio_with_old(datum, new_key, old_data, old_key)
            for datum in new_data
        ]
    )
    ax.errorbar(
        masses,
        values,
        yerr=errors,
        color=colour,
        marker=marker,
        label=label,
        linestyle="none",
    )


def plot(new_data, old_data):
    fig, ax = plt.subplots()

    ax.set_xlabel(r"$am_{f}$")
    ax.set_ylabel(r"$F_{\rho} / F_{\pi}$")

    lengths = reversed(sorted(set(datum["Nx"] for datum in new_data)))
    props = Props(length_only=True)
    for length in lengths:
        subset = [datum for datum in new_data if datum["Nx"] == length]
        plot_ratio_with_old(ax, props, subset, "rho_decay_const", old_data, "fpi")

    add_qcd_value(ax, "rho_decay_const", "pi_decay_const")

    ax.set_ylim(0, 2)

    return fig


if __name__ == "__main__":
    comparison_plot_main(plot, old_data=True)
