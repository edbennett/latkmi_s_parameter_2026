#!/usr/bin/env python3

import matplotlib.pyplot as plt

from ..stats import jackknife_mean_variance
from ..plot import Props, comparison_plot_main

from .plot_rho_pi_decay_ratio import (
    plot_ratio_with_old,
    plot_new_old_derived,
    get_old_samples,
)


def get_ksrf_i(new_datum, old_data):
    rho_mass_samples = new_datum["fit_result_samples"]["rho_mass"]
    rho_decay_const_samples = new_datum["fit_result_samples"]["rho_decay_const"]
    pi_decay_const_samples = get_old_samples(new_datum, old_data, "fpi")
    ratio_samples = (
        rho_mass_samples
        * rho_decay_const_samples
        / (2**0.5 * pi_decay_const_samples**2)
    )
    return jackknife_mean_variance(ratio_samples)


def plot_ksrf_i(ax, props, new_data, old_data):
    plot_new_old_derived(ax, props, get_ksrf_i, new_data, old_data)


def plot_ksrf_ii(ax, props, new_data, old_data, **kwargs):
    plot_ratio_with_old(ax, props, new_data, "rho_mass", old_data, "fpi", **kwargs)


def plot(new_data, old_data):
    fig, ax = plt.subplots()

    ax.set_xlabel(r"$am_{f}$")
    ax.set_ylabel(r"$g_{\rho\pi\pi}^{\mathrm{KSRF}}$")

    lengths = reversed(sorted(set(datum["Nx"] for datum in new_data)))
    props = Props(length_only=True)
    for length in lengths:
        subset = [datum for datum in new_data if datum["Nx"] == length]
        plot_ksrf_i(ax, props, subset, old_data)
        plot_ksrf_ii(ax, props, subset, old_data, colour="gray", label=None)

    ax.set_xlim(0, None)
    ax.set_ylim(4, 7)

    return fig


if __name__ == "__main__":
    comparison_plot_main(plot, old_data=True)
