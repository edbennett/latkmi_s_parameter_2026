#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

from format_multiple_errors import format_multiple_errors

from ..io import read_numpy


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_files", metavar="input_file", nargs="+")
    parser.add_argument("--output_file", type=FileType("w"), default="-")
    return parser.parse_args()


def get_with_attribute(data, spatial_size, temporal_size, mass, attribute):
    subset = [
        datum
        for datum in data
        if ("Nx" not in datum or datum["Nx"] == spatial_size)
        and ("Ny" not in datum or datum["Ny"] == spatial_size)
        and ("Nz" not in datum or datum["Nz"] == spatial_size)
        and ("Nt" not in datum or datum["Nt"] == temporal_size)
        and datum["mass"] == mass
        and attribute in datum
    ]
    if len(subset) == 0:
        raise ValueError("Datum not found")
    if len(subset) > 1:
        raise ValueError("Multiple results found")
    return subset[0]


def get_row(datum):
    spatial_size = datum["Nx"]
    assert spatial_size == datum["Ny"] and spatial_size == datum["Nz"]
    temporal_size = datum["Nt"]
    mass = datum["mass"]

    row_data = [spatial_size, temporal_size, mass] + [
        format_multiple_errors(*datum[key], abbreviate=True, significant_figures=2)
        for key in ["ksrf-i", "ksrf-ii", "wsr-i", "wsr-ii"]
    ]
    return " & ".join(map(str, row_data)) + r" \\"


def tabulate(data):
    header = "\n".join(
        [
            r"\begin{tabular}{lllllll}",
            r"\toprule",
            " & ".join(
                [
                    "$L$",
                    "$T$",
                    "$am_f$",
                    r"$g_{\rho\pi\pi}^{\textnormal{\scriptsize{KSRF-I}}}$",
                    r"$g_{\rho\pi\pi}^{\textnormal{\scriptsize{KSRF-II}}}$",
                    "WSR-I",
                    "WSR-II",
                ]
            )
            + r" \\",
            r"\midrule",
        ]
    )
    footer = "\n".join([r"\bottomrule", r"\end{tabular}"])

    content = [get_row(datum) for datum in sorted(data, key=lambda x: x["mass"])]
    return "\n".join([header, *content, footer])


def main():
    args = get_args()
    data = [read_numpy(input_file) for input_file in args.input_files]
    print(tabulate(data), file=args.output_file)


if __name__ == "__main__":
    main()
