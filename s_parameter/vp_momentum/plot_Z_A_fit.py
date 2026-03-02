#!/usr/bin/env python3

from argparse import ArgumentParser
import json

import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat

from ..fits import fit_form_linear, fit_form_quadratic
from ..io import read_Z_A
from ..plot import save_or_show


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_files", nargs="+", metavar="input_file")
    parser.add_argument(
        "--extra_input_file",
        dest="extra_input_files",
        action="append",
        metavar="input_file",
    )
    parser.add_argument("--plot_styles", default="styles/prd.mplstyle")
    parser.add_argument("--output_file", default=None)
    parser.add_argument("--fit_result", required=True)
    return parser.parse_args()


def read(filenames):
    return [read_Z_A(filename) for filename in filenames]


def plot_Z_A_data(ax, data, marker, label_offset):
    mf_values = [ensemble["mass"] for ensemble in data]
    Z_A_values = [ensemble["Z_A"] for ensemble in data]
    Z_A_errors = [ensemble["Z_A_error"] for ensemble in data]

    errorbar = ax.errorbar(
        mf_values,
        Z_A_values,
        yerr=Z_A_errors,
        linestyle="none",
        marker=marker,
    )

    for ensemble in data:
        assert ensemble["Nx"] == ensemble["Ny"] == ensemble["Nz"]
        ax.annotate(
            f"$L = {ensemble['Nx']}$",
            (ensemble["mass"], ensemble["Z_A"]),
            xytext=(0, label_offset),
            textcoords="offset points",
            color=errorbar[0]._color,
            horizontalalignment="right" if label_offset < 0 else "left",
        )


def plot_fit_data(ax, fit_results):
    xmin, xmax = ax.get_xlim()
    xmax *= 1.15
    for tag, fit_form, label in [
        ("linear", fit_form_linear, "$Z_A + cm_f$"),
        ("quadratic", fit_form_quadratic, "$Z_A + cm_f^2$"),
    ]:
        fit_result = fit_results[tag]
        mf_range = np.linspace(xmin, fit_result["mf_max"] * 1.2, 1000)
        ax.plot(
            mf_range,
            fit_form(mf_range, fit_result["Z_A_0"][0], fit_result["c"][0]),
            label=(
                rf"{label}, $\chi^2/\mathrm{{dof}} = "
                f"{fit_result['chisquare']:.1f}/{fit_result['dof']}$"
            ),
        )

    ax.set_xlim(xmin, xmax)


def add_arrow(ax, fit_result):
    result_ufloat = ufloat(*fit_result)
    ax.annotate(
        f"$Z_A = {result_ufloat:.01uSL}$",
        (0, fit_result[0]),
        xytext=(16, 12),
        textcoords="offset points",
        arrowprops={
            "facecolor": "black",
            "width": 0.2,
            "headwidth": 3,
            "headlength": 4,
        },
    )


def plot(data, extra_data, fit_result):
    fig, ax = plt.subplots()
    plot_Z_A_data(ax, data, "o", 6)
    plot_Z_A_data(ax, extra_data, "s", -12)
    ax.set_xlim(0, None)
    plot_fit_data(ax, fit_result)
    add_arrow(ax, fit_result["central"]["Z_A_0"])

    ax.set_xlabel("$am_f$")
    ax.set_ylabel("$Z_A$")
    ax.legend(loc="best")

    return fig


def main():
    args = get_args()
    plt.style.use(args.plot_styles)

    Z_A_data = read(args.input_files)
    extra_Z_A_data = read(args.extra_input_files)
    with open(args.fit_result) as file_object:
        Z_A_fit = json.load(file_object)

    save_or_show(plot(Z_A_data, extra_Z_A_data, Z_A_fit), args.output_file)


if __name__ == "__main__":
    main()
