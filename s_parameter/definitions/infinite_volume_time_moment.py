#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

from ..define import define_many
from ..io import read_numpy


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_data")
    parser.add_argument("--output_definitions", type=FileType("w"), default="-")
    return parser.parse_args()


def get_definitions(fit_result):
    chisquare = fit_result["chisquare"]
    dof = fit_result["dof"]
    return define_many(
        ("Time_Moment_Infinite_Volume_C", tuple(fit_result["C"])),
        ("Time_Moment_Infinite_Volume_Chisquare_Dof", f"{chisquare:.02g}/{dof}"),
    )


def main():
    args = get_args()
    fit_result = read_numpy(args.input_data)
    print(get_definitions(fit_result["fit_result"]), file=args.output_definitions)


if __name__ == "__main__":
    main()
