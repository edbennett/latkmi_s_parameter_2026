#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

import matplotlib.pyplot as plt
import numpy as np

from ..define import define
from ..io import read_Z_A
from ..plot import save_or_show


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_files", nargs="+", metavar="input_file")
    parser.add_argument("--plot_styles", default="styles/prd.mplstyle")
    parser.add_argument("--output_file", default=None)
    parser.add_argument("--Z_A_lower_bound", type=float, default=None)
    parser.add_argument("--output_definitions", type=FileType("w"), default=None)
    return parser.parse_args()


def plot_plateau_lines(ax, ensemble, colour):
    error = ensemble["Z_A_error"]
    for offset, dashes in [(-error, (3, 2)), (0, (None, None)), (error, (3, 2))]:
        y_position = ensemble["Z_A"] + offset
        ax.plot(
            [ensemble["plateau_start"] - 0.5, ensemble["plateau_end"] + 0.5],
            [y_position, y_position],
            dashes=dashes,
            color=colour,
            linewidth=0.7,
        )


def plot(data, Z_A_lower_bound=None):
    fig, ax = plt.subplots()
    markers = "ovsD^"
    for colour_index, (marker, ensemble) in enumerate(
        zip(markers, sorted(data, key=lambda v: v["mass"]))
    ):
        colour = f"C{colour_index}"
        lattice_size = ensemble["Nx"]
        assert lattice_size == ensemble["Ny"] == ensemble["Nz"]
        ax.errorbar(
            np.arange(2, len(ensemble["Z_A_eff"]) - 2),
            ensemble["Z_A_eff"][2:-2],
            yerr=ensemble["Z_A_eff_error"][2:-2],
            marker=marker,
            color=colour,
            linestyle="none",
            label=f"$m_f = {ensemble['mass']}$, $L={lattice_size}$",
        )
        plot_plateau_lines(ax, ensemble, colour)

    ax.set_xlabel("$t / a$")
    ax.set_ylabel("$Z_A$")
    ax.set_ylim(Z_A_lower_bound, None)

    ax.legend(loc="best")

    return fig


def read(filenames):
    return [read_Z_A(filename) for filename in filenames]


def get_definitions(data):
    components = [
        f"({datum['Nx']}, {datum['mass']})"
        for datum in sorted(data, key=lambda d: d["Nx"], reverse=True)
    ]
    content = ", ".join(components)
    return define("Z_A_Ensembles", f"\\left\\{{{content}\\right\\}}")


def main():
    args = get_args()
    plt.style.use(args.plot_styles)
    data = read(args.input_files)
    save_or_show(plot(data, args.Z_A_lower_bound), args.output_file)
    if args.output_definitions:
        print(get_definitions(data), file=args.output_definitions)


if __name__ == "__main__":
    main()
