#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

from ..define import define, name_ensemble
from ..io import read_numpy


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_data")
    parser.add_argument("--output_definitions", type=FileType("w"), default="-")
    return parser.parse_args()


def get_definitions(sum_rules):
    return define(f"l10_r_{name_ensemble(sum_rules)}", tuple(sum_rules["l10-r"]))


def main():
    args = get_args()
    sum_rules = read_numpy(args.input_data)
    print(get_definitions(sum_rules), file=args.output_definitions)


if __name__ == "__main__":
    main()
