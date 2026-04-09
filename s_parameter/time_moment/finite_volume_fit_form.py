#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

import numpy as np

from ..io import dump_numpy
from .finite_volume import delta_fv_S


def get_args():
    parser = ArgumentParser()
    parser.add_argument("--output_file", default="-", type=FileType("w"))
    return parser.parse_args()


def main():
    args = get_args()
    masses = np.linspace(0.2, 1, 81)
    length = 16  # Dummy length, to give m_pi L in (0, 16]
    result = delta_fv_S(length, masses, outer_projected_decay_time=42)

    dump_numpy(
        {
            "length": length,
            "masses": masses,
            "delta_fv_S": result,
        },
        args.output_file,
    )


if __name__ == "__main__":
    main()
