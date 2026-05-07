#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

from ..define import define_many
from ..io import read_numpy


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_data", nargs="+")
    parser.add_argument("--output_definitions", type=FileType("w"), default="-")
    return parser.parse_args()


def get_definitions(data):
    chisquare_values = [
        datum["pade_fit_result"]["Conserved"]["chisquare"][0]
        / datum["pade_fit_result"]["Conserved"]["dof"]
        for datum in data
    ]
    second_largest_chisquare = sorted(chisquare_values)[-2]
    largest_chisquare = max(chisquare_values)

    return define_many(
        ("Pade_Fit_Largest_Chisquare", f"{largest_chisquare:.01g}"),
        ("Pade_Fit_Second_Largest_Chisquare", f"{second_largest_chisquare:.01g}"),
    )


def main():
    args = get_args()
    data = [read_numpy(input_file) for input_file in args.input_data]
    print(get_definitions(data), file=args.output_definitions)


if __name__ == "__main__":
    main()
