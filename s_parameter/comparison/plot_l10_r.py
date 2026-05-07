#!/usr/bin/env python3

from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..io import read_numpy
from ..plot import save_or_show
from ..stats import add_quadrature


def get_args():
    parser = ArgumentParser()
    parser.add_argument("latkmi_input_files", metavar="input_file", nargs="+")
    parser.add_argument("--rbc_ukqcd_data", required=True)
    parser.add_argument("--jlqcd_data", required=True)
    parser.add_argument("--plot_styles", default="styles/prd.mplstyle")
    parser.add_argument("--output_file", default=None)
    return parser.parse_args()


def add_band(ax, data, colour, label):
    # If the mass is finite, we should plot points rather than a band
    (m_f,) = data.m_f
    assert m_f == 0.0

    (value,) = data.value_l10_r
    (statistical_uncertainty,) = data.error_l10_r

    if "systematic_upper_l10_r" in data.columns:
        (systematic_upper,) = data.systematic_upper_l10_r
        (systematic_lower,) = data.systematic_lower_l10_r
    elif "systematic_l10_r" in data.columns:
        (systematic_upper,) = data.systematic_l10_r
        systematic_lower = systematic_upper
    else:
        systematic_upper = 0
        systematic_lower = 0

    total_upper_uncertainty = add_quadrature(systematic_upper, statistical_uncertainty)
    total_lower_uncertainty = add_quadrature(systematic_lower, statistical_uncertainty)

    ax.axhline(value, color=colour)
    ax.axhspan(
        value - total_lower_uncertainty,
        value + total_upper_uncertainty,
        color=colour,
        label=label,
        alpha=0.4,
    )


def add_data(ax, data, key, offset=0, **props):
    x_values, *x_errors = zip(*[datum["mpi-mrho"] for datum in data])
    y_values, *y_errors = zip(*[datum[key] for datum in data])
    ax.errorbar(
        [value + offset for value in x_values],
        y_values,
        xerr=np.linalg.norm(x_errors),
        yerr=np.linalg.norm(y_errors),
        linestyle="none",
        **props,
    )


def plot(latkmi_data, rbc_ukqcd_data, jlqcd_data):
    fig, ax = plt.subplots()

    ax.set_xlabel(r"$M_\pi / M_\rho$")
    ax.set_ylabel(r"$L_{10}^r (M_\rho)$")

    add_band(ax, rbc_ukqcd_data, "C0", "RBC/UKQCD 2010 ($N_f=2+1$)")
    add_band(ax, jlqcd_data, "C1", "JLQCD 2008 ($N_f=2$)")

    add_data(
        ax,
        latkmi_data,
        "l10-r",
        color="C2",
        label=r"LatKMI w. $S_{\mathrm{TM},\infty}$ ($N_f=8$)",
        marker="s",
        fillstyle="full",
    )
    add_data(
        ax,
        latkmi_data,
        "dmo-l10-r",
        color="gray",
        label=r"LatKMI w. $S_{\mathrm{DMO}}$ ($N_f=8$)",
        marker="o",
        fillstyle="none",
        offset=2e-3,
    )

    ax.set_xlim(0.6, 0.72)
    ax.set_ylim(-8e-3, 0)

    ax.legend(loc="best")

    return fig


def main():
    args = get_args()
    plt.style.use(args.plot_styles)

    latkmi_data = [read_numpy(input_file) for input_file in args.latkmi_input_files]
    rbc_ukqcd_data = pd.read_csv(args.rbc_ukqcd_data, comment="#")
    jlqcd_data = pd.read_csv(args.jlqcd_data, comment="#")

    fig = plot(latkmi_data, rbc_ukqcd_data, jlqcd_data)
    save_or_show(fig, args.output_file)


if __name__ == "__main__":
    main()
