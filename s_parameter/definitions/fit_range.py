#!/usr/bin/env python3

from argparse import ArgumentParser, FileType


from ..define import define_many


def get_args():
    parser = ArgumentParser()
    parser.add_argument("plateau_start", type=int)
    parser.add_argument("plateau_end", type=int)
    parser.add_argument("--prefix", default="")
    parser.add_argument("--output_definitions", type=FileType("w"), default="-")
    return parser.parse_args()


def get_definitions(plateau_start, plateau_end, prefix):
    return define_many(
        (f"{prefix}_Plateau_Start", plateau_start),
        (f"{prefix}_Plateau_End", plateau_end),
        (f"{prefix}_Plateau_Range", f"[{plateau_start}, {plateau_end}]"),
    )


def main():
    args = get_args()
    print(
        get_definitions(args.plateau_start, args.plateau_end, args.prefix),
        file=args.output_definitions,
    )


if __name__ == "__main__":
    main()
