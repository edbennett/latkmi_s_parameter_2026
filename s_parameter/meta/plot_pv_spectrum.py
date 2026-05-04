#!/usr/bin/env python3

from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..plot import save_or_show, Props


def get_args():
    parser = ArgumentParser()
    parser.add_argument("--input_metadata", required=True)
    parser.add_argument("--input_spectrum", required=True)
    parser.add_argument("--plot_styles", default="styles/prd.mplstyle")
    parser.add_argument("--output_file", default=None)
    return parser.parse_args()


def add_single_size(ax, lattice_size, metadata, data, offset, props):
    colour, marker = props.burn()
    ax.errorbar(
        [np.nan],
        [np.nan],
        yerr=[np.nan],
        marker=marker,
        color=colour,
        linestyle="none",
        label=f"$L = {lattice_size}$",
    )

    subset_metadata = metadata[metadata.Nx == lattice_size]
    subset_data = []

    for metadatum in subset_metadata.to_dict(orient="records"):
        assert metadatum["Ny"] == lattice_size
        assert metadatum["Nz"] == lattice_size

        matched_data = data.query(
            f"Nf == {metadatum['Nf']} & mf == {metadatum['mf']} & "
            f"T == {metadatum['Nt']} & L == {metadatum['Nx']}"
        )
        assert len(matched_data) == 1
        subset_data.append(matched_data.iloc[0])

    subset_df = pd.DataFrame(subset_data)

    for channel in "rho", "a1":
        ax.errorbar(
            subset_df.mf + offset * 0.0005,
            subset_df[f"value_m{channel}"],
            subset_df[f"error_m{channel}"],
            color=colour,
            marker=marker,
            linestyle="--",
        )


def plot(metadata, spectrum):
    fig, ax = plt.subplots()
    props = Props()

    ax.set_xlabel("$am_f$")
    ax.set_ylabel("$aM_X$")
    offsets = {18: -1, 24: 1, 30: 0, 36: -1, 42: 0, 48: 0}

    for lattice_size in sorted(set(metadata.Nx), reverse=True):
        add_single_size(
            ax, lattice_size, metadata, spectrum, offsets[lattice_size], props
        )

    ax.set_xlim(0, None)
    ax.set_ylim(0, None)

    ax.text(0.015, 0.5, r"$aM_{a_1(\mathrm{PV})}$")
    ax.text(0.03, 0.35, r"$aM_{\rho(\mathrm{PV})}$")

    ax.legend(loc="best")
    return fig


def main():
    args = get_args()
    plt.style.use(args.plot_styles)
    metadata = pd.read_csv(args.input_metadata, comment="#")
    spectrum = pd.read_csv(args.input_spectrum, comment="#")
    save_or_show(plot(metadata[metadata.plot_pv_masses], spectrum), args.output_file)


if __name__ == "__main__":
    main()
