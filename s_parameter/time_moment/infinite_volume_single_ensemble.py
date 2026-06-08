#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

from ..io import read_numpy, dump_numpy
from ..stats import generate_jackknife, jackknife_mean_variance


METADATA_KEYS = ["mass", "Nt", "Nx", "Ny", "Nz", "Nf", "bin_size"]


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_files", metavar="INPUT_FILE", nargs="+")
    parser.add_argument("--fit_result", required=True)
    parser.add_argument("--output_file", default="-", type=FileType("w"))
    return parser.parse_args()


def compute_systematics(datum, infinite_volume_S):
    result = {key: datum[key] for key in ["Nx", "Ny", "Nz", "Nt"]}
    finite_volume_S = datum["S_infinite_t_samples"]
    result["systematic"] = abs(finite_volume_S[0] - infinite_volume_S[0])
    return result


def get_infinite_volume_S(data, result):
    assert len(data) > 0
    (mass,) = set(datum["mass"] for datum in data)
    (Nf,) = set(datum["Nf"] for datum in data)

    if len(data) > 1:
        # We already fitted these data; don't need to re-compute
        mass_index = list(result["fit_result"]["masses"]).index(mass)
        infinite_volume_S = result["fit_result"]["S_infinite_volume"][mass_index]
        return {
            "S_infinite_volume": infinite_volume_S,
            "method": "fit_result",
            "mass": mass,
            "Nf": Nf,
            "source_metadata": result["source_metadata"],
            "systematics": [
                compute_systematics(datum, infinite_volume_S) for datum in data
            ],
        }

    # We need to perform an extrapolation from our one point
    (datum,) = data
    finite_volume_S = datum["S_infinite_t_samples"]
    finite_volume_factor = datum["delta_fv_S_samples"]
    const_coefficient = generate_jackknife(
        *result["fit_result"]["C"], datum, finite_volume_S.shape
    )
    # Apply Eq. (53) in reverse using fit result to get S_\infty
    infinite_volume_S_samples = (
        finite_volume_S - const_coefficient * finite_volume_factor
    )
    infinite_volume_S = jackknife_mean_variance(infinite_volume_S_samples)
    return {
        "S_infinite_volume": infinite_volume_S,
        "S_infinite_volume_samples": infinite_volume_S_samples,
        "method": "extrapolation",
        "mass": mass,
        "Nf": Nf,
        "source_metadata": [{key: datum[key] for key in METADATA_KEYS}],
        "systematics": [compute_systematics(datum, infinite_volume_S)],
    }


def main():
    args = get_args()
    data = [read_numpy(input_file) for input_file in args.input_files]
    result = read_numpy(args.fit_result)

    dump_numpy(get_infinite_volume_S(data, result), args.output_file)


if __name__ == "__main__":
    main()
