#!/usr/bin/env python3

from argparse import ArgumentParser, FileType
import logging

import numpy as np
import pandas as pd

from ..io import read_numpy, dump_numpy
from ..stats import jackknife_mean_error, generate_jackknife


METADATA_KEYS = ["mass", "Nf", "Nt", "Nx", "Ny", "Nz", "bin_size"]


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_files", metavar="input_file", nargs="+")
    parser.add_argument("--previous_data", required=True)
    parser.add_argument("--skip_rule", dest="skip_rules", action="append", default=[])
    parser.add_argument("--output_file", default="-", type=FileType("w"))
    return parser.parse_args()


def get_common_shape(all_samples):
    (result,) = set(samples.shape for samples in all_samples.values())
    return result


def get_systematic(candidates, match_result):
    (matched_candidate,) = [
        candidate
        for candidate in candidates
        if all(candidate[key] == match_result[key] for key in ["Nx", "Ny", "Nz", "Nt"])
    ]
    return matched_candidate["systematic"]


def combine_samples(data):
    """
    Combine results of multiple different fits into a single dict for ease of processing.
    """
    if not data:
        return {}

    to_remove = []
    result = data[0]
    if "fit_result" not in result:
        result["fit_result"] = {}
    if "fit_result_samples" not in result:
        result["fit_result_samples"] = {}

    for datum in data[1:]:
        for key in METADATA_KEYS:
            if key in datum and key not in result:
                result[key] = datum[key]
            if key in datum and result[key] != datum[key]:
                raise ValueError(f"Data mismatch: {result[key]} != {datum[key]}")

        if "S_infinite_volume" in datum:
            value, error = datum["S_infinite_volume"]
            if "S_infinite_volume_samples" in datum:
                result["fit_result_samples"]["S_infinite_volume"] = datum[
                    "S_infinite_volume_samples"
                ]
            systematic = get_systematic(datum["systematics"], result)
            result["fit_result"]["S_infinite_volume"] = np.array(
                [value, error, systematic]
            )
            continue

        for key in datum["fit_result"]:
            if key in result["fit_result"]:
                to_remove.append(key)
                del result["fit_result"][key]
                del result["fit_result_samples"][key]
                del result["fit_result_systematic_samples"][key]

            if key in to_remove:
                continue

            result["fit_result"][key] = datum["fit_result"][key]
            result["fit_result_samples"][key] = datum["fit_result_samples"][key]
            systematic_samples = datum["fit_result_systematic_samples"][key]
            result["fit_result_systematic_samples"][key] = systematic_samples

    for key, value in result["fit_result_systematic_samples"].items():
        if np.isnan(value).all():
            result["fit_result_systematic_samples"][key] = result["fit_result"][key][0]

    for key, samples in result["fit_result_systematic_samples"].items():
        # Randomise the direction of the variation,
        # to capture the effect of systematics in unknown directions
        if not isinstance(samples, np.ndarray):
            continue
        num_samples = len(samples)
        result["fit_result_systematic_samples"][key] = (-1) ** np.arange(
            num_samples
        ) * samples / num_samples**0.5 + result["fit_result"][key][0]

    if "S_infinite_volume" in result["fit_result"]:
        shape = get_common_shape(result["fit_result_samples"])
        value, error, systematic = result["fit_result"]["S_infinite_volume"]
        if "S_infinite_volume" not in result["fit_result_samples"]:
            result["fit_result_samples"]["S_infinite_volume"] = generate_jackknife(
                value, error, result, shape
            )
        result["fit_result_systematic_samples"]["S_infinite_volume"] = (
            generate_jackknife(value, systematic, result, shape)
        )

    return result


def get_old_data(datum, old_data, key):
    result = old_data.query(
        "Nf == 8 & beta == 3.8 & "
        f"L == {datum['Nx']} & T == {datum['Nt']} & mf == {datum['mass']}"
    )
    assert datum["Nx"] == datum["Ny"] and datum["Nx"] == datum["Nz"]
    if len(result.index) > 1:
        raise ValueError("Multiple ensembles found.")
    if len(result.index) == 0:
        return None, None
    value = result[f"value_{key}"].iloc[0]
    error = result[f"error_{key}"].iloc[0]
    return value, error


def add_generated_samples(new_datum, old_data):
    """
    Take the preexisting data from old_data for the ensemble described in new_datum,
    generate jackknife samples for relevant observables,
    and add these to new_datum.
    """
    num_samples = len(next(iter(new_datum["fit_result_samples"].values())))

    for old_key, new_key in [
        ("fpi", "pi_decay_const"),
        ("mpi", "pi_mass"),
        ("t0c", "t0"),
        ("mrho", "rho_pv_mass"),
    ]:
        value, error = get_old_data(new_datum, old_data, old_key)
        if value:
            samples = generate_jackknife(value, error, new_datum, num_samples)
        else:
            samples, value, error = np.array([np.nan]), np.nan, np.nan

        new_datum["fit_result_samples"][new_key] = samples
        new_datum["fit_result"][new_key] = [value, error]
        new_datum["fit_result_systematic_samples"][new_key] = value


def wsr_i(samples):
    """
    First Weinberg sum rule; see Eq. (45) of the paper
    """
    f_rho = samples["rho_decay_const"]
    f_pi = samples["pi_decay_const"]
    f_a_1 = samples["a_1_decay_const"]

    return f_rho**2 - f_a_1**2 - f_pi**2


def wsr_i_normalised(samples):
    """
    First Weinberg sum rule with normalisation; see Fig. 12 of the paper.
    """
    f_rho = samples["rho_decay_const"]
    f_pi = samples["pi_decay_const"]
    f_a_1 = samples["a_1_decay_const"]

    return wsr_i(samples) / (f_rho**2 + f_a_1**2 + f_pi**2)


def wsr_ii(samples):
    """
    Second Weinberg sum rule with normalisation; see Eq. (46) and Fig. 12 of the paper
    """
    f_m_rho_squared = (samples["rho_decay_const"] * samples["rho_mass"]) ** 2
    f_m_a_1_squared = (samples["a_1_decay_const"] * samples["a_1_mass"]) ** 2

    return f_m_rho_squared - f_m_a_1_squared


def wsr_ii_normalised(samples):
    """
    Second Weinberg sum rule with normalisation; see Fig. 12 of the paper.
    """
    f_m_rho_squared = (samples["rho_decay_const"] * samples["rho_mass"]) ** 2
    f_m_a_1_squared = (samples["a_1_decay_const"] * samples["a_1_mass"]) ** 2

    return wsr_ii(samples) / (f_m_rho_squared + f_m_a_1_squared)


def ksrf_i(samples):
    """
    First Kawarabayashi-Suzuki-Riazuddin-Fayyazuddin relation;
    see Eq. (35) of the paper.
    """
    return (
        samples["rho_mass"]
        * samples["rho_decay_const"]
        / (2**0.5 * samples["pi_decay_const"] ** 2)
    )


def ksrf_ii(samples):
    """
    Second Kawarabayashi-Suzuki-Riazuddin-Fayyazuddin relation;
    see Eq. (36) of the paper.
    """
    return samples["rho_mass"] / samples["pi_decay_const"]


def dmo(samples):
    """
    Das-Mathur-Okubo sum rule; see Eq. (44) of the paper.
    """
    f_rho = samples["rho_decay_const"]
    f_a_1 = samples["a_1_decay_const"]
    m_rho = samples["rho_mass"]
    m_a_1 = samples["a_1_mass"]
    return 2 * np.pi * (f_rho**2 / m_rho**2 - f_a_1**2 / m_a_1**2)


def l10_r(samples, S, rho_type):
    """
    The low energy constant $L_{10}^r$; see Eq. (55) of the paper.
    """
    rho_keys = {
        "pv": "rho_pv_mass",
        "vt": "rho_mass",
    }
    Nf = 8  # Number of flavours
    m_pi = samples["pi_mass"]
    m_rho = samples[rho_keys[rho_type]]
    return -S / (16 * np.pi) - 1 / (192 * np.pi**2) * (Nf / 2) * (
        np.log(m_pi**2 / m_rho**2) + 1
    )


rules = {
    "frho-fpi": lambda s: s["rho_decay_const"]
    / s["pi_decay_const"],  # TODO check normalisation
    "frho-fa1": lambda s: s["rho_decay_const"] / s["a_1_decay_const"],
    "mpi-mrho": lambda s: s["pi_mass"] / s["rho_mass"],
    "ma1-mrho": lambda s: s["a_1_mass"] / s["rho_mass"],
    "mrho_s8t0": lambda s: s["rho_mass"] * (8 * s["t0"]) ** 0.5,
    "mpi_L": lambda s: s["pi_mass"] * s["Nx"],
    "wsr-i": wsr_i,
    "wsr-ii": wsr_ii,
    "wsr-i-normalised": wsr_i_normalised,
    "wsr-ii-normalised": wsr_ii_normalised,
    "ksrf-i": ksrf_i,
    "ksrf-ii": ksrf_ii,
    "dmo": dmo,
    "l10-r": lambda s: l10_r(s, s["S_infinite_volume"], "pv"),
    "dmo-l10-r": lambda s: l10_r(s, dmo(s), "vt"),
}


def compute_sum_rules(data, prefixes_to_skip=[]):
    result = {}
    for name, func in rules.items():
        if any(name.startswith(prefix) for prefix in prefixes_to_skip):
            continue
        try:
            central, statistical = jackknife_mean_error(
                func({**data, **data["fit_result_samples"]})
            )
            systematic_samples = func({**data, **data["fit_result_systematic_samples"]})
            if isinstance(systematic_samples, float):
                systematic = np.nan
            else:
                _, systematic = jackknife_mean_error(systematic_samples)
            result[name] = [central, statistical, systematic]
        except KeyError as key:
            message = f"Key {key} not found in data. Skipping computation of {name}."
            logging.warning(message)

    return result


def main():
    args = get_args()

    new_data = combine_samples(
        [read_numpy(input_file) for input_file in args.input_files]
    )
    old_data = pd.read_csv(args.previous_data)
    add_generated_samples(new_data, old_data)

    result = compute_sum_rules(new_data, args.skip_rules)
    dump_numpy(
        {
            **result,
            **{key: new_data[key] for key in METADATA_KEYS},
        },
        args.output_file,
    )


if __name__ == "__main__":
    main()
