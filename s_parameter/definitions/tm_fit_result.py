#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

from format_multiple_errors import format_multiple_errors

from ..define import define_many
from ..io import read_numpy


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_data")
    parser.add_argument("--prefix", default="")
    parser.add_argument("--output_definitions", type=FileType("w"), default="-")
    return parser.parse_args()


def get_definitions(fit_result, prefix):
    plateau_start = fit_result["min_timeslice"]
    plateau_end = fit_result["max_timeslice"] - 1

    return define_many(
        (f"{prefix}_Plateau_Start", plateau_start),
        (f"{prefix}_Plateau_End", plateau_end),
        (f"{prefix}_Plateau_Range", f"[{plateau_start}, {plateau_end}]"),
        (
            f"{prefix}_Chisquare_Per_Dof",
            f"{fit_result['chisquare']:.02g}/{fit_result['dof']}",
        ),
        (
            f"{prefix}_S_Parameter",
            format_multiple_errors(
                *fit_result["S_infinite_t"],
                significant_figures=1,
                latex=True,
                abbreviate=True,
            ),
        ),
    )


def main():
    args = get_args()
    fit_result = read_numpy(args.input_data)
    print(
        get_definitions(fit_result, args.prefix),
        file=args.output_definitions,
    )


if __name__ == "__main__":
    main()
