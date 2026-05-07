#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

from uncertainties import ufloat

from ..define import define
from ..io import read_numpy


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_data")
    parser.add_argument("--output_definitions", type=FileType("w"), default="-")
    return parser.parse_args()


def get_definitions(datum):
    S_ufloat = ufloat(*datum["pade_fit_result"]["Conserved"]["S"])
    return define(f"S_Parameter_VP_{datum['mass']}", f"{S_ufloat:.01uSL}")


def main():
    args = get_args()
    datum = read_numpy(args.input_data)
    print(get_definitions(datum), file=args.output_definitions)


if __name__ == "__main__":
    main()
