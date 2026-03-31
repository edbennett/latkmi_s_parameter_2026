#!/usr/bin/env python3

from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..io import read_numpy, collate_ensembles
from ..plot import iterate_attribute, save_or_show, Props


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_files", metavar="input_file", nargs="+")
    parser.add_argument("--spectrum_data", required=True)
    parser.add_argument("--lsd_data", default=None)
    parser.add_argument("--plot_styles", default="styles/prd.mplstyle")
    parser.add_argument("--output_file", default=None)
    parser.add_argument("--group_by", required=True, choices=["mass", "length"])
    return parser.parse_args()


def get_only(data, key):
    (result,) = set(data[key])
    return result


def add_lsd_data(ax, data, props):
    mock_datum = {
        "Nx": get_only(data, "Nx"),
        "Ny": get_only(data, "Ny"),
        "Nz": get_only(data, "Nz"),
        "Nt": get_only(data, "Nt"),
        "mass": None,
    }
    colour, marker, label = props.get(mock_datum)
    ax.errorbar(
        data.value_pi_mass * data.Nx,
        data.value_S_lattice,
        xerr=data.error_pi_mass * data.Nx,
        yerr=data.error_S_lattice,
        color=colour,
        marker=marker,
        label=f"LSD DWF {label}",
    )


def get_mpi_L(datum, spectrum_data):
    L = datum["Nx"]
    assert L == datum["Ny"] and L == datum["Nz"]

    subset = spectrum_data.query(f"Nf == 8 & beta == 3.8 & mf == {datum['mass']}")
    if len(subset) == 0:
        return np.nan, np.nan

    largest_volume_datum = subset.sort_values(by="L").iloc[-1]
    return largest_volume_datum.value_mpi * L, largest_volume_datum.error_mpi * L


def plot(latkmi_data, spectrum_data, group_by, lsd_data=None):
    fig, ax = plt.subplots()

    ax.set_xlabel(r"$M_\pi L$")
    ax.set_ylabel(r"$4\pi\Pi^\prime(0)$")
    props = Props(**{f"{group_by}_only": True})

    sorted_data = sorted(
        collate_ensembles(latkmi_data),
        key=lambda datum: tuple(get_mpi_L(datum, spectrum_data)),
    )

    for subset, (colour, marker, label) in iterate_attribute(
        sorted_data, group_by, props
    ):
        mpi_L, mpi_L_error = zip(*[get_mpi_L(datum, spectrum_data) for datum in subset])
        S_param, S_param_error = zip(
            *[datum["pade_fit_result"]["Conserved"]["S"] for datum in subset]
        )
        ax.errorbar(
            mpi_L,
            S_param,
            xerr=mpi_L_error,
            yerr=S_param_error,
            marker=marker,
            color=colour,
            label=label,
            dashes=(1, 3),
        )

    if lsd_data is not None:
        add_lsd_data(ax, lsd_data, props)

    ax.set_xlim(4.5, 12)
    ax.set_ylim(0.15, 0.32)

    ax.legend(loc="best", ncols=2)

    return fig


def main():
    args = get_args()
    plt.style.use(args.plot_styles)

    latkmi_data = [read_numpy(input_file) for input_file in args.input_files]
    spectrum_data = pd.read_csv(args.spectrum_data)
    lsd_data = pd.read_csv(args.lsd_data) if args.lsd_data else None

    fig = plot(latkmi_data, spectrum_data, args.group_by, lsd_data)
    save_or_show(fig, args.output_file)


if __name__ == "__main__":
    main()
