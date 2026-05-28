#!/usr/bin/env python3

from argparse import ArgumentParser, FileType
from functools import partial
from multiprocessing import Pool

import pandas as pd

from ..io import read_numpy
from ..stats import generate_jackknife, jackknife_mean_variance, add_quadrature
from ..time_moment.finite_volume import delta_fv_S


JACKKNIFE_SAMPLE_SIZE = 50


def _generate_jackknife(datum, channel):
    return generate_jackknife(
        datum[f"value_{channel}"],
        datum[f"error_{channel}"],
        datum,
        JACKKNIFE_SAMPLE_SIZE,
    )


def get_args():
    parser = ArgumentParser()
    parser.add_argument("--s_parameter_data", required=True)
    parser.add_argument("--chiral_data", required=True)
    parser.add_argument("--fit_result", required=True)
    parser.add_argument("--output_file", default="-", type=FileType("w"))
    return parser.parse_args()


def extrapolate_single(datum, fit_result):
    assert datum["Nx"] == datum["Ny"]
    assert datum["Nx"] == datum["Nz"]
    # Add metadata for compatibility with generate_jackknife_samples
    datum["mass"] = datum["m_f"]
    datum["bin_size"] = 1
    pi_mass = _generate_jackknife(datum, "pi_mass")
    finite_volume_factor = delta_fv_S(datum["Nx"], pi_mass)
    finite_volume_S = _generate_jackknife(datum, "S_lattice")
    const_coefficient = generate_jackknife(
        *fit_result["fit_result"]["C"], datum, JACKKNIFE_SAMPLE_SIZE
    )
    return jackknife_mean_variance(
        finite_volume_S - const_coefficient * finite_volume_factor
    )


def extrapolate(s_parameter_data, fit_result):
    with Pool() as pool:
        results = pool.map(
            partial(extrapolate_single, fit_result=fit_result),
            s_parameter_data.to_dict(orient="records"),
        )
    result_df = s_parameter_data.copy()
    result_df["value_S_infinite_volume"], result_df["error_S_infinite_volume"] = zip(
        *results
    )
    return result_df


def add_pi_over_chiral_rho_mass(input_data, chiral_data):
    chiral_result = chiral_data.query("Nf == 8")
    assert len(chiral_result) == 1
    chiral_datum = chiral_result.iloc[0]
    value_chiral_rho_mass = chiral_datum.value_rho_mass
    error_chiral_rho_mass = chiral_datum.error_rho_mass

    data = input_data.copy()

    value_ratio = data["value_pi_mass"] / value_chiral_rho_mass
    error_ratio = value_ratio * add_quadrature(
        (data["value_pi_mass"], data["error_pi_mass"]),
        (value_chiral_rho_mass, error_chiral_rho_mass),
    )
    data["value_pi_mass_over_chiral_rho_mass"] = value_ratio
    data["error_pi_mass_over_chiral_rho_mass"] = error_ratio
    return data


def process(s_parameter_data, chiral_data, fit_result):
    new_data = add_pi_over_chiral_rho_mass(s_parameter_data, chiral_data)
    extrapolated_data = extrapolate(new_data, fit_result)
    return extrapolated_data


def main():
    args = get_args()

    s_parameter_data = pd.read_csv(args.s_parameter_data, comment="#")
    chiral_data = pd.read_csv(args.chiral_data, comment="#")
    fit_result = read_numpy(args.fit_result)
    process(s_parameter_data, chiral_data, fit_result).to_csv(
        args.output_file, index=False
    )


if __name__ == "__main__":
    main()
