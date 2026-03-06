#!/usr/bin/env python3

from argparse import ArgumentParser

import matplotlib.pyplot as plt
import pandas as pd

from ..io import read_numpy
from ..stats import jackknife_mean_variance
from ..plot import save_or_show

from .plot_rho_pi_decay_ratio import compute_ratio_with_old, get_old_samples


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_files", metavar="input_file", nargs="+")
    parser.add_argument("--previous_data", required=True)
    parser.add_argument("--lsd_data", required=True)
    parser.add_argument("--plot_styles", default="styles/prd.mplstyle")
    parser.add_argument("--output_file", default=None)
    return parser.parse_args()


def get_mpi_s8t0(new_data, old_data):
    result = []
    for datum in new_data:
        mpi_samples = get_old_samples(datum, old_data, "mpi")
        t0_samples = get_old_samples(datum, old_data, "t0c")
        result.append(jackknife_mean_variance(mpi_samples * (8 * t0_samples) ** 0.5))
    return zip(*result)


def plot(new_data, old_data, lsd_data):
    fig, ax = plt.subplots()

    ax.set_xlabel(r"$am_{f}$")
    ax.set_ylabel(r"$g_{\rho\pi\pi}^{\textnormal{\scriptsize{KSRF-II}}}$")

    latkmi_mpi_s8t0_values, latkmi_mpi_s8t0_errors = get_mpi_s8t0(new_data, old_data)
    latkmi_ksrf_values, latkmi_ksrf_errors = zip(
        *[
            compute_ratio_with_old(datum, "rho_mass", old_data, "fpi")
            for datum in new_data
        ]
    )
    ax.errorbar(
        latkmi_mpi_s8t0_values,
        latkmi_ksrf_values,
        xerr=latkmi_mpi_s8t0_errors,
        yerr=latkmi_ksrf_errors,
        linestyle="none",
        marker="s",
        label="LatKMI",
    )

    ax.errorbar(
        lsd_data.value_mpi_s8t0,
        lsd_data.value_g_rho_pi_pi,
        xerr=lsd_data.error_mpi_s8t0,
        yerr=lsd_data.error_g_rho_pi_pi,
        linestyle="none",
        marker="o",
        label="LSD Phys. Rev. D99, 014509",
    )
    ax.set_xlim(0, None)
    ax.set_ylim(4, 7)

    ax.legend(loc="best")

    return fig


def main():
    args = get_args()
    plt.style.use(args.plot_styles)

    new_data = [read_numpy(input_file) for input_file in args.input_files]

    old_data = pd.read_csv(args.previous_data)
    lsd_data = pd.read_csv(args.lsd_data)

    fig = plot(new_data, old_data, lsd_data)
    save_or_show(fig, args.output_file)


if __name__ == "__main__":
    main()
