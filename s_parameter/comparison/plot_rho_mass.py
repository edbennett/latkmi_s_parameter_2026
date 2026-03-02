#!/usr/bin/env python3

from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..io import read_numpy
from ..plot import Props, save_or_show


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_files", metavar="input_file", nargs="+")
    parser.add_argument("--previous_data", required=True)
    parser.add_argument("--plot_styles", default="styles/prd.mplstyle")
    parser.add_argument("--output_file", default=None)
    return parser.parse_args()


def add_old_data(ax, target_ensembles, source_data, marker):
    matching_data = []
    for ensemble in target_ensembles:
        result = source_data.query(
            "Nf == 8 & beta == 3.8 & "
            f"L == {ensemble['Nx']} & T == {ensemble['Nt']} & mf == {ensemble['mass']}"
        )
        if len(result.index) > 1:
            raise ValueError("Multiple ensembles found.")
        matching_data.append(result)

    data = pd.concat(matching_data)
    ax.errorbar(
        data.mf + 0.001,
        data.value_mrho,
        data.error_mrho,
        marker=marker,
        color="grey",
        linestyle="none",
    )


def add_new_data(ax, data, colour, marker, label):
    masses = np.array([datum["mass"] for datum in data])
    rho = np.array([datum["fit_result"]["rho_mass"] for datum in data])
    ax.errorbar(
        masses,
        rho[:, 0],
        yerr=rho[:, 1],
        marker=marker,
        color=colour,
        label=label,
        linestyle="none",
    )
    ax.errorbar(
        masses,
        rho[:, 0],
        yerr=(rho[:, 1:] ** 2).sum() ** 0.5,
        color=colour,
        linestyle="none",
    )


def plot(new_data, old_data):
    fig, ax = plt.subplots()

    ax.set_xlabel(r"$am_{f}$")
    ax.set_ylabel(r"$aM_{\mathrm{\rho}}$")

    props = Props(length_only=True)
    lengths = reversed(sorted(set(datum["Nx"] for datum in new_data)))
    for length in lengths:
        subset = [datum for datum in new_data if datum["Nx"] == length]
        colour, marker, label = props.get(subset[0])
        add_new_data(ax, subset, colour, marker, label)
        add_old_data(ax, subset, old_data, marker)

    ax.set_xlim(0, None)
    ax.set_ylim(0.2, 0.5)

    ax.legend(loc="best")
    return fig


def main():
    args = get_args()
    plt.style.use(args.plot_styles)
    new_data = [read_numpy(input_file) for input_file in args.input_files]
    old_data = pd.read_csv(args.previous_data)

    fig = plot(new_data, old_data)
    save_or_show(fig, args.output_file)


if __name__ == "__main__":
    main()
