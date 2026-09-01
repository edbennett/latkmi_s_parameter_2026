#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

from ..define import define
from ..io import read_numpy


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_data")
    parser.add_argument("--output_definitions", type=FileType("w"), default="-")
    return parser.parse_args()


def get_definitions(data):
    assert len(data["systematics"] == 1)
    systematic = data["systematics"][0]["systematic"]
    return define(
        f"S_Parameter_Infinite_Volume_Time_Moment_Nf{data['Nf']}_mf{data['mass']}",
        tuple(list(data["S_infinite_volume"]) + [systematic]),
        length_control="largest",
        significant_figures=2,
    )


def main():
    args = get_args()
    infinite_volume_S = read_numpy(args.input_data)
    print(get_definitions(infinite_volume_S), file=args.output_definitions)


if __name__ == "__main__":
    main()
