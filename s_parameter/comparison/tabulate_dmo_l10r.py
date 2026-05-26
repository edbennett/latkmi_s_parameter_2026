#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

from format_multiple_errors import format_multiple_errors
import numpy as np

from .tabulate_sum_rules import get_with_attribute
from ..io import read_numpy


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_files", metavar="input_file", nargs="+")
    parser.add_argument("--output_file", type=FileType("w"), default="-")
    return parser.parse_args()


def get_row(data, spatial_size, temporal_size, mass):
    time_moment_datum = get_with_attribute(
        data, spatial_size, temporal_size, mass, "S_infinite_volume"
    )
    sum_rule_datum = get_with_attribute(
        data,
        spatial_size,
        temporal_size,
        mass,
        "l10-r",
    )

    def formatter(value_with_errors, scale_factor=1):
        if value_with_errors is None:
            return "---"
        if np.isnan(value_with_errors[-1]):
            slug = "(-)"
            value_with_errors = value_with_errors[:-1]
        else:
            slug = ""

        return "${}{}$".format(
            format_multiple_errors(
                *(value_with_errors * scale_factor),
                abbreviate=True,
                length_control="central",
                significant_figures=3,
            ),
            slug,
        )

    time_moment_S = formatter(time_moment_datum["S_infinite_volume"])
    time_moment_l10_r = formatter(sum_rule_datum.get("l10-r"), 1000)
    dmo_S = formatter(sum_rule_datum.get("dmo"))
    dmo_l10_r = formatter(sum_rule_datum.get("dmo-l10-r"), 1000)

    row_data = [
        spatial_size,
        temporal_size,
        mass,
        time_moment_S,
        time_moment_l10_r,
        dmo_S,
        dmo_l10_r,
    ]
    return " & ".join(map(str, row_data)) + r" \\"


def tabulate(data):
    ensembles = sorted(
        sorted(
            set(
                [
                    (datum["Nx"], datum["Nt"], datum["mass"])
                    for datum in data
                    if "Nx" in datum
                ]
            ),
            key=lambda x: x[0],  # Reverse sort by volume first
            reverse=True,
        ),
        key=lambda x: x[2],  # Then sort by mass
    )
    header = "\n".join(
        [
            r"\begin{tabular}{lllllll}",
            r"\toprule",
            " & ".join(
                [
                    "$L$",
                    "$T$",
                    "$am_f$",
                    r"$S_{\mathrm{TM},\infty}$",
                    r"$L_{10}^r(M_\rho)\cdot{10}^3$",
                    r"$S_{\mathrm{DMO}}$",
                    r"$L_{10}^r(M_\rho)|_{\mathrm{DMO}}\cdot{10}^3$",
                ]
            )
            + r" \\",
            r"\midrule",
        ]
    )
    footer = "\n".join([r"\bottomrule", r"\end{tabular}"])

    content = [get_row(data, *ensemble) for ensemble in ensembles]
    return "\n".join([header, *content, footer])


def main():
    args = get_args()
    data = [read_numpy(input_file) for input_file in args.input_files]
    print(tabulate(data), file=args.output_file)


if __name__ == "__main__":
    main()
