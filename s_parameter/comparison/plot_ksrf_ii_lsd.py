#!/usr/bin/env python3

from argparse import ArgumentParser

import matplotlib.pyplot as plt
import pandas as pd

from ..io import read_numpy
from ..plot import save_or_show


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_files", metavar="input_file", nargs="+")
    parser.add_argument("--lsd_data", required=True)
    parser.add_argument("--plot_styles", default="styles/prd.mplstyle")
    parser.add_argument("--output_file", default=None)
    return parser.parse_args()


def plot(latkmi_data, lsd_data):
    fig, ax = plt.subplots()

    ax.set_xlabel(r"$M_\rho \sqrt{8t_0}$")
    ax.set_ylabel(r"$g_{\rho\pi\pi}^{\textnormal{\scriptsize{KSRF-II}}}$")

    latkmi_mpi_s8t0_values, latkmi_mpi_s8t0_errors = zip(
        *[datum["mrho_s8t0"] for datum in latkmi_data]
    )
    latkmi_ksrf_values, latkmi_ksrf_errors = zip(
        *[datum["ksrf-ii"] for datum in latkmi_data]
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
        lsd_data.value_mrho_s8t0,
        lsd_data.value_g_rho_pi_pi,
        xerr=lsd_data.error_mrho_s8t0,
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

    latkmi_data = [read_numpy(input_file) for input_file in args.input_files]
    lsd_data = pd.read_csv(args.lsd_data)

    fig = plot(latkmi_data, lsd_data)
    save_or_show(fig, args.output_file)


if __name__ == "__main__":
    main()
