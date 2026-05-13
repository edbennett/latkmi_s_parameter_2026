#!/usr/bin/env python3

from argparse import ArgumentParser

import matplotlib.pyplot as plt
import pandas as pd

from ..io import read_numpy
from ..plot import save_or_show


def get_args():
    parser = ArgumentParser()
    parser.add_argument(
        "--time_moment_infinite_volume", metavar="input_file", nargs="+"
    )
    parser.add_argument("--chiral_fit_result", required=True)
    parser.add_argument("--spectrum_data", required=True)
    parser.add_argument("--lsd_data", required=True)
    parser.add_argument("--plot_styles", default="styles/prd.mplstyle")
    parser.add_argument("--output_file", default=None)
    return parser.parse_args()


def get_rho_mass(chiral_data):
    linear_result = chiral_data.query(
        "Nf == 8 & state == 'mrho' & fit_form == 'Linear'"
    )
    quadratic_result = chiral_data.query(
        "Nf == 8 & state == 'mrho' & fit_form == 'Quadratic'"
    )

    assert len(linear_result) == 1
    assert len(quadratic_result) == 1
    result = quadratic_result.iloc[0]
    params = list(map(float, result["params"].split(";")))
    errors = list(map(float, result["errors"].split(";")))

    linear_params = list(map(float, linear_result.iloc[0]["params"].split(";")))
    systematic = abs(linear_params[0] - params[0])
    return params[0], errors[0], systematic


def get_pi_mass(spectrum_data, datum):
    results = spectrum_data.query(f"Nf == {datum['Nf']} & mf == {datum['mass']}")
    assert len(results) >= 1
    result = results.sort_values(by="L", ascending=False).iloc[0]
    return result.value_mpi, result.error_mpi


def get_mpi_over_mrho(chiral_data, spectrum_data, datum):
    value_rho_mass_chiral, error_rho_mass_chiral, _ = get_rho_mass(chiral_data)
    value_pi_mass, error_pi_mass = get_pi_mass(spectrum_data, datum)
    value_ratio = value_pi_mass / value_rho_mass_chiral
    error_ratio = (
        value_ratio
        * (
            (error_pi_mass / value_pi_mass) ** 2
            + (error_rho_mass_chiral / value_rho_mass_chiral) ** 2
        )
        ** 0.5
    )
    return value_ratio, error_ratio


def add_latkmi_data(ax, infinite_volume_data, chiral_data, spectrum_data):
    data_to_plot = []
    for infinite_volume_datum in infinite_volume_data:
        data_to_plot.append(
            (
                *get_mpi_over_mrho(chiral_data, spectrum_data, infinite_volume_datum),
                *infinite_volume_datum["S_infinite_volume"],
            )
        )
    value_mpi_over_mrho, error_mpi_over_mrho, value_S, error_S = zip(*data_to_plot)
    ax.errorbar(
        value_mpi_over_mrho,
        value_S,
        xerr=error_mpi_over_mrho,
        yerr=error_S,
        color="C0",
        marker="s",
        linestyle="none",
        label=r"LatKMI $L=\infty$",
    )


def add_lsd_data(ax, data):
    value_mpi_over_mrho = data["value_pi_mass_over_chiral_rho_mass"]
    error_mpi_over_mrho = data["error_pi_mass_over_chiral_rho_mass"]

    ax.errorbar(
        value_mpi_over_mrho,
        data["value_S_infinite_volume"],
        xerr=error_mpi_over_mrho,
        yerr=data["error_S_infinite_volume"],
        color="C1",
        marker="x",
        linestyle="none",
        label=r"LSD $L=\infty$",
    )
    ax.errorbar(
        value_mpi_over_mrho + 0.02,
        data["value_S_lattice"],
        xerr=error_mpi_over_mrho,
        yerr=data["error_S_lattice"],
        color="silver",
        marker=",",
        linestyle="none",
        label="LSD $L=32$",
    )


def plot(time_moment_data, chiral_data, spectrum_data, lsd_data):
    fig, ax = plt.subplots()

    ax.set_xlabel(r"$M_\pi / M_\rho^{\chi\mathrm{lim}}$")
    ax.set_ylabel(r"$S|_{\textnormal{\scriptsize{1-doublet}}}$")

    add_latkmi_data(ax, time_moment_data, chiral_data, spectrum_data)
    add_lsd_data(ax, lsd_data)

    ax.set_xlim(0, None)
    ax.set_ylim(0.2, 0.3)

    ax.legend(loc="best")
    return fig


def main():
    args = get_args()
    plt.style.use(args.plot_styles)

    time_moment_infinite_volume = [
        read_numpy(input_file) for input_file in args.time_moment_infinite_volume
    ]
    chiral_data = pd.read_csv(args.chiral_fit_result, comment="#")
    spectrum_data = pd.read_csv(args.spectrum_data, comment="#")
    lsd_data = pd.read_csv(args.lsd_data, comment="#")

    fig = plot(time_moment_infinite_volume, chiral_data, spectrum_data, lsd_data)
    save_or_show(fig, args.output_file)


if __name__ == "__main__":
    main()
