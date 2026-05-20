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
    parser.add_argument("input_files", metavar="input_file", nargs="+")
    parser.add_argument("--lsd_data", required=True)
    parser.add_argument("--plot_styles", default="styles/prd.mplstyle")
    parser.add_argument("--output_file", default=None)
    return parser.parse_args()


def plot(latkmi_data, lsd_data):
    fig, ax = plt.subplots()

    ax.set_xlabel(r"$M_\rho \sqrt{8t_0}$")
    ax.set_ylabel(r"$g_{\rho\pi\pi}^{\textnormal{\scriptsize{KSRF-II}}}$")

    latkmi_mpi_s8t0_values, *latkmi_mpi_s8t0_errors = map(
        np.array,
        zip(*[datum["mrho_s8t0"] for datum in latkmi_data]),
    )
    latkmi_ksrf_values, *latkmi_ksrf_errors = map(
        np.array,
        zip(*[datum["ksrf-ii"] for datum in latkmi_data]),
    )
    for xerrors, yerrors, label in [
        (latkmi_mpi_s8t0_errors[0], latkmi_ksrf_errors[0], "LatKMI"),
        (
            add_quadrature(*latkmi_mpi_s8t0_errors),
            add_quadrature(*latkmi_ksrf_errors),
            None,
        ),
    ]:
        ax.errorbar(
            latkmi_mpi_s8t0_values,
            latkmi_ksrf_values,
            xerr=xerrors,
            yerr=yerrors,
            linestyle="none",
            marker="s",
            label=label,
            color="C0",
        )

    ax.errorbar(
        lsd_data.value_mrho_s8t0,
        lsd_data.value_g_rho_pi_pi,
        xerr=lsd_data.error_mrho_s8t0,
        yerr=lsd_data.error_g_rho_pi_pi,
        linestyle="none",
        marker="o",
        label="LSD Phys. Rev. D99, 014509",
        color="C1",
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
