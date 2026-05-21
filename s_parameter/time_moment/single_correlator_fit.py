#!/usr/bin/env python3

from argparse import ArgumentParser, FileType
from functools import partial
from itertools import product

import numpy as np
from scipy.optimize import curve_fit

from ..io import read_numpy, dump_numpy
from ..stats import (
    jackknife_mean_variance,
    jackknife_systematic_error,
    jackknife_systematic_intermediary,
    sample_systematics,
)


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_file", metavar="input_file")
    parser.add_argument("--min_timeslice", type=int, default=0)
    parser.add_argument("--max_timeslice", type=int, default=None)
    parser.add_argument("--channel", required=True, choices=["V", "A"])
    parser.add_argument("--output_file", type=FileType("w"), default="-")
    return parser.parse_args()


def fit_form(time, mass_main, amplitude_main, mass_osc, amplitude_osc, max_time):
    return amplitude_main * (
        np.exp(-mass_main * time) + np.exp(-mass_main * (max_time - time))
    ) + amplitude_osc * (-1) ** time * (
        np.exp(-mass_osc * time) + np.exp(-mass_osc * (max_time - time))
    )


def decay_const(mass, amplitude):
    num_staggered_tastes = 4
    return np.abs((2.0 * amplitude / mass) * (1 / num_staggered_tastes)) ** 0.5


def fit_single(samples, min_timeslice, max_timeslice, starting_guess=None):
    _, data_uncertainty = jackknife_mean_variance(samples)

    range_fit_form = partial(fit_form, max_time=(samples[0].shape[-1] - 1) * 2)
    target_range = np.arange(min_timeslice, max_timeslice)
    target_slice = slice(min_timeslice, max_timeslice)

    if starting_guess is None:
        starting_guess, _ = curve_fit(
            range_fit_form,
            target_range,
            samples.mean(axis=0)[target_slice],
            sigma=data_uncertainty[target_slice],
            p0=[1.0, 1.0, 1.0, -1.0],
        )

    results = []
    for sample in samples:
        (mass_main, amplitude_main, mass_osc, amplitude_osc), _ = curve_fit(
            range_fit_form,
            target_range,
            sample[target_slice],
            sigma=data_uncertainty[target_slice],
            p0=starting_guess,
        )
        decay_const_main = decay_const(mass_main, amplitude_main)
        decay_const_osc = decay_const(mass_osc, amplitude_osc)
        results.append(
            [
                mass_main,
                amplitude_main,
                decay_const_main,
                mass_osc,
                amplitude_osc,
                decay_const_osc,
            ]
        )
    return results


def fit_systematic(samples, result, min_timeslice, max_timeslice):
    starting_guess = [result[i] for i in [0, 1, 3, 4]]

    def fit_time_range(min_timeslice, max_timeslice):
        return fit_single(
            samples, min_timeslice, max_timeslice, starting_guess=starting_guess
        )

    return jackknife_systematic_intermediary(
        sample_systematics(fit_time_range, min_timeslice, max_timeslice), result
    )


def fit(full_data, min_timeslice, max_timeslice, channel):
    # Negate samples to have positive data to fit
    data = -full_data["Conserved"][f"{channel}_renormalised_samples"]

    fit_samples = fit_single(data, min_timeslice, max_timeslice)
    values, errors = jackknife_mean_variance(fit_samples)
    systematic_samples = fit_systematic(data, values, min_timeslice, max_timeslice)
    systematic_errors = jackknife_systematic_error(systematic_samples, values)

    states = {
        "V": {"main": "rho", "osc": "V_osc"},
        "A": {"main": "a_1", "osc": "A_osc"},
    }
    return {
        key: {
            f"{states[channel][name]}_{observable}": data
            for (name, observable), data in zip(
                product(["main", "osc"], ["mass", "amplitude", "decay_const"]),
                zip(*values),
            )
        }
        for key, values in [
            ("fit_result_samples", fit_samples),
            ("fit_result_systematic_samples", systematic_samples),
            ("fit_result", [values, errors, systematic_errors]),
        ]
    }


def main():
    args = get_args()
    data = read_numpy(args.input_file)

    result = fit(data, args.min_timeslice, args.max_timeslice, args.channel)
    dump_numpy(
        {
            **result,
            **{
                key: data[key]
                for key in ["mass", "Nf", "Nt", "Nx", "Ny", "Nz", "bin_size"]
            },
            "channel": args.channel,
            "min_timeslice": args.min_timeslice,
            "max_timeslice": args.max_timeslice,
        },
        args.output_file,
    )


if __name__ == "__main__":
    main()
