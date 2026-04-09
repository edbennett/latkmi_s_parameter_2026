#!/usr/bin/env python3

from argparse import ArgumentParser

import matplotlib.pyplot as plt

from ..io import read_numpy
from ..plot import iterate_attribute, save_or_show, Props
from ..stats import jackknife_mean_variance


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_files", metavar="input_file", nargs="+")
    parser.add_argument("--fit_result", required=True)
    parser.add_argument("--fit_form", required=True)
    parser.add_argument("--plot_styles", default="styles/prd.mplstyle")
    parser.add_argument("--output_file", default=None)
    return parser.parse_args()


def add_fit_band(ax, fit_result, fit_form, props):
    const_coefficient_value, const_coefficient_error = fit_result["C"]
    L_M_pi = fit_form["masses"] * fit_form["length"]
    delta_fv_S = fit_form["delta_fv_S"]
    colour, _ = props.burn()
    ax.fill_between(
        L_M_pi,
        delta_fv_S * (const_coefficient_value - const_coefficient_error),
        delta_fv_S * (const_coefficient_value + const_coefficient_error),
        color=colour,
        alpha=0.3,
        label=r"Fit: $-\Delta^{\mathrm{FV}} S(LM_\pi)$",
    )


def add_data(ax, data, fit_result, props):
    for subset, (colour, marker, label) in iterate_attribute(data, "mass", props):
        mass = subset[0]["mass"]
        mass_index = list(fit_result["masses"]).index(mass)
        infinite_volume_S_value, infinite_volume_S_error = fit_result[
            "S_infinite_volume"
        ][mass_index]
        L_M_pi = [
            datum["Nx"] * jackknife_mean_variance(datum["pi_mass_samples"])[0]
            for datum in subset
        ]
        difference_value = [
            infinite_volume_S_value - datum["S_infinite_t"][0] for datum in subset
        ]
        difference_error = [
            (infinite_volume_S_error**2 + datum["S_infinite_t"][1] ** 2)
            ** 0.5  # TODO systematic
            for datum in subset
        ]
        ax.errorbar(
            L_M_pi,
            difference_value,
            difference_error,
            linestyle="none",
            color=colour,
            marker=marker,
            label=label,
        )


def plot(data, fit_result, fit_form):
    fig, ax = plt.subplots()

    ax.set_xlim(5, 12)
    ax.set_ylim(-0.002, 0.04)

    ax.set_xlabel(r"$LM_\pi$")
    ax.set_ylabel(r"$S_\infty(M_\pi)-S(L, M_\pi)$")

    props = Props(mass_only=True)
    add_fit_band(ax, fit_result["fit_result"], fit_form, props)
    add_data(ax, data, fit_result["fit_result"], props)

    ax.legend(loc="best")
    return fig


def main():
    args = get_args()
    plt.style.use(args.plot_styles)

    data = [read_numpy(input_file) for input_file in args.input_files]
    fit_result = read_numpy(args.fit_result)
    fit_form = read_numpy(args.fit_form)

    fig = plot(data, fit_result, fit_form)
    save_or_show(fig, args.output_file)


if __name__ == "__main__":
    main()
