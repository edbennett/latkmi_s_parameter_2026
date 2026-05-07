#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

from uncertainties import ufloat

from ..define import define_many
from ..io import read_numpy


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_data")
    parser.add_argument("--output_definitions", type=FileType("w"), default="-")
    return parser.parse_args()


def get_definitions(Z_A):
    Z_A_ufloat = ufloat(Z_A[0], Z_A[2])
    Z_A_statistical_percentage = Z_A[1] / Z_A[0] * 100
    Z_A_systematic_percentage = Z_A[2] / Z_A[0] * 100
    return define_many(
        ("Z_A_Value", f"{Z_A_ufloat:.01uSL}"),
        ("Z_A_Statistical_Percentage", f"{Z_A_statistical_percentage:.01g}\\%"),
        ("Z_A_Systematic_Percentage", f"{Z_A_systematic_percentage:.01g}\\%"),
    )


def main():
    args = get_args()
    Z_A = read_numpy(args.input_data)
    print(get_definitions(Z_A["central"]["Z_A_0"]), file=args.output_definitions)


if __name__ == "__main__":
    main()
