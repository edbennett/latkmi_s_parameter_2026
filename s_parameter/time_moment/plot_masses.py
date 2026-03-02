#!/usr/bin/env python3

from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np

from ..io import read_numpy
from ..plot import save_or_show


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_files", metavar="input_file", nargs="+")
    parser.add_argument("--plot_styles", default="styles/prd.mplstyle")
    parser.add_argument("--output_file", default=None)
    return parser.parse_args()


def plot(data):
    fig, axes = plt.subplots(ncols=2, figsize=(7, 3))
    for observable, label, ax in zip(["mass", "decay_const"], ["M", "F"], axes):
        ax.set_xlabel("$am_f$")
        ax.set_ylabel(f"$a{label}_X$")
        for channel, channel_label, colour, marker in [
            ("rho", r"\rho", "C0", "o"),
            ("a_1", "a_1", "C1", "s"),
        ]:
            channel_data = [
                datum for datum in data if f"{channel}_mass" in datum["fit_result"]
            ]
            masses = np.array([datum["mass"] for datum in channel_data])
            values, errors, systematics = map(
                np.array,
                zip(
                    *[
                        datum["fit_result"][f"{channel}_{observable}"]
                        for datum in channel_data
                    ]
                ),
            )
            combined_errors = (errors**2 + systematics**2) ** 0.5
            valid_data = abs(errors / values) < 0.2
            ax.errorbar(
                masses[valid_data],
                values[valid_data],
                yerr=errors[valid_data],
                label=f"${channel_label}$",
                linestyle="none",
                marker=marker,
                color=colour,
            )
            ax.errorbar(
                masses[valid_data],
                values[valid_data],
                yerr=combined_errors[valid_data],
                linestyle="none",
                color=colour,
            )
        ax.set_xlim(0, None)
        ax.set_ylim(0, None)
        ax.legend(loc="best")

    return fig


def main():
    args = get_args()
    plt.style.use(args.plot_styles)
    data = [read_numpy(input_file) for input_file in args.input_files]

    fig = plot(data)
    save_or_show(fig, args.output_file)


if __name__ == "__main__":
    main()
