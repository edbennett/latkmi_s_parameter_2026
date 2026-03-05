#!/usr/bin/env python3

from argparse import ArgumentParser, FileType
from functools import partial

import numpy as np
from scipy.optimize import curve_fit

from ..io import read_numpy, dump_numpy
from ..stats import (
    jackknife_mean_variance,
    jackknife_systematic_error,
    sample_systematics,
)


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_file", metavar="input_file")
    parser.add_argument("--input_m_rho", required=True)
    parser.add_argument("--input_m_a_1", required=True)
    parser.add_argument("--min_timeslice", type=int, default=0)
    parser.add_argument("--max_timeslice", type=int, default=None)
    parser.add_argument("--output_file", type=FileType("w"), default="-")
    return parser.parse_args()


def fit_form(time, C_A_plus, C_A_minus, C_V_plus, C_V_minus, m_a_1, m_rho, max_time):
    """
    Fit form described in Eq. (26) of the paper.
    `time` is expected to be an integer.
    The four `C_` parameters are the fit parameters.
    `m_a_1` and `m_rho` are the masses of the a1 and rho meson,
    and should be preset using `partial` or similar.
    `max_time` is the value $T$ from the paper,
    and should be preset using `partial` or similar.
    """
    return (C_V_plus - (-1) ** time * C_A_plus) * (
        np.exp(-m_rho * time)
        + (0 if max_time is None else np.exp(-m_rho * (max_time - time)))
    ) - (C_A_minus - (-1) ** time * C_V_minus) * (
        np.exp(-m_a_1 * time)
        + (0 if max_time is None else np.exp(-m_a_1 * (max_time - time)))
    )


def correlator_infinite_t_I(mass, t0):
    """
    Extrapolate the forward correlator,
    implementing  Eq. (29) of the paper.

    t0: start of fit region
    mass: mass of state whose correlator we're extrapolating
    """
    exp_mass = np.exp(mass)
    return (
        np.exp(-mass * (t0 - 1))
        / (exp_mass - 1) ** 3
        * (t0**2 * (exp_mass - 1) ** 2 + 2 * t0 * (exp_mass - 1) + (exp_mass + 1))
    )


def correlator_infinite_t_J(t0, mass):
    """
    Extrapolate the forward correlator,
    implementing  Eq. (30) of the paper.

    t0: start of fit region
    mass: mass of state whose correlator we're extrapolating
    """
    exp_mass = np.exp(mass)
    return (
        (-1) ** t0
        * np.exp(-mass * (t0 - 1))
        / (exp_mass + 1) ** 3
        * (t0**2 * (exp_mass + 1) ** 2 - 2 * t0 * (exp_mass + 1) - (exp_mass - 1))
    )


def extrapolate_S_infinite_t(fit_samples, m_rho, m_a_1, t_0, S_t_0):
    """
    Take a jackknife sample set of fit results,
    and apply Eq. (27) of the paper to obtain
    the extrapolated S parameter at infinite T.
    fit_samples: The fit results for the four parameters C_{+,-}^{V,A}
    m_rho, m_a_1: Samples of the masses of the rho and a_1 states
    t_0: The start of the fit interval
    S_t_0: The cumulative summed value of S at the start of the fit interval
    """

    # Take transpose to match shape of fit_samples: (samples, parameters)
    factors = np.array(
        [
            correlator_infinite_t_I(t_0, m_rho),
            -correlator_infinite_t_I(t_0, m_a_1),
            -correlator_infinite_t_J(t_0, m_rho),
            correlator_infinite_t_I(t_0, m_a_1),
        ]
    ).T

    # Eq. (28) uses 4pi * (-1/2!); this is simplified to -2pi here.
    contributions = -2 * np.pi * (fit_samples * factors)
    result = {
        f"S_{name}": (value, error)
        for name, value, error in zip(
            ["V_plus", "A_minus", "A_plus", "V_minus"],
            *jackknife_mean_variance(contributions),
        )
    }
    S_samples = S_t_0 + contributions.sum(axis=1)
    result["S_infinite_t_samples"] = S_samples
    result["S_infinite_t"] = jackknife_mean_variance(S_samples)

    return result


def extrapolate_S_finite_t(fit_samples, m_rho, m_a_1, t_0, t_max, S_t_0):
    r"""
    Take a jackknife sample set of fit results,
    and apply Eqs. (27)--(30) of the paper without the T->\infty limit
    to obtain the extrapolated effective S parameter as a function of T.
    """
    times = np.arange(t_0, t_max + 1)
    masses = np.array([m_rho, m_a_1, m_rho, m_a_1])

    # Desired axis ordering:
    # bootstrap sample, contribution (rho/a1; A/V), time
    factors = np.cumsum(
        np.array([1, 1, -1, -1])[:, np.newaxis] ** times
        * times**2
        * np.exp(-masses.T[:, :, np.newaxis] * times),
        axis=2,
    )
    S_eff_samples = S_t_0[:, np.newaxis] - 2 * np.pi * (
        fit_samples[..., np.newaxis] * factors
    ).sum(axis=1)
    return [times, *jackknife_mean_variance(S_eff_samples)]


def interpolate_correlator_and_S(fit_samples, m_rho, m_a_1, t_0, t_1, S_t_0):
    """
    Take a jackknife sample set of fit results,
    and apply Eq. (22) of the paper
    to obtain the interpolated effective S parameter as a function of T.
    """
    times = np.arange(t_0, t_1)[:, np.newaxis]
    correlator = fit_form(times, *(np.array(fit_samples).T), m_a_1, m_rho, t_1 - 1).T
    S_eff_samples = (S_t_0 + np.cumsum(-2 * np.pi * times**2 * correlator.T, axis=0)).T
    return {
        "times": times[:, 0],
        "correlator": jackknife_mean_variance(correlator),
        "S_effective": jackknife_mean_variance(S_eff_samples),
    }


def fit_single(full_data, m_rho, m_a_1, min_timeslice, max_timeslice):
    # Negate samples to have positive data to fit
    data = -full_data["Conserved"]["V-A_renormalised_samples"]

    assert data.shape[:-1] == m_rho.shape
    assert data.shape[:-1] == m_a_1.shape
    _, data_uncertainty = jackknife_mean_variance(data)

    fit_samples = []
    for sample, m_rho_sample, m_a_1_sample in zip(data, m_rho, m_a_1):
        result, _ = curve_fit(
            partial(
                fit_form,
                m_a_1=m_a_1_sample,
                m_rho=m_rho_sample,
                max_time=2 * (data.shape[-1] - 1),
            ),
            np.arange(min_timeslice, max_timeslice),
            sample[min_timeslice:max_timeslice],
            sigma=data_uncertainty[min_timeslice:max_timeslice],
        )
        fit_samples.append(result)

    values, errors = jackknife_mean_variance(fit_samples)
    C_A_plus, C_A_minus, C_V_plus, C_V_minus = zip(values, errors)

    start_S_samples = full_data["Conserved"]["S_parameter_eff_samples"][
        :, min_timeslice
    ]
    extrapolated_S = extrapolate_S_infinite_t(
        fit_samples,
        m_rho,
        m_a_1,
        min_timeslice,
        start_S_samples,
    )
    return {
        "fit_result_samples": np.array(fit_samples),
        "C_V_plus": C_V_plus,
        "C_A_plus": C_A_plus,
        "C_V_minus": C_V_minus,
        "C_A_minus": C_A_minus,
        **extrapolated_S,
    }


def add_systematic(result, all_samples, key):
    value, error = result[key]
    samples = [sample[key][0] for sample in all_samples]
    systematic = jackknife_systematic_error(samples, value)
    return value, error, systematic


def fit(full_data, m_rho, m_a_1, min_timeslice, max_timeslice):
    def fit_time_range(min_timeslice, max_timeslice):
        return fit_single(full_data, m_rho, m_a_1, min_timeslice, max_timeslice)

    full_range_result = fit_time_range(min_timeslice, max_timeslice)
    systematic_samples = sample_systematics(
        fit_time_range, min_timeslice, max_timeslice
    )

    for key in ["C_V_plus", "C_A_plus", "C_V_minus", "C_A_minus", "S_infinite_t"]:
        full_range_result[key] = add_systematic(
            full_range_result, systematic_samples, key
        )

    fit_samples = full_range_result["fit_result_samples"]
    start_S_samples = full_data["Conserved"]["S_parameter_eff_samples"][
        :, min_timeslice
    ]
    full_range_result["S_extrapolation_large_t"] = extrapolate_S_finite_t(
        fit_samples, m_rho, m_a_1, min_timeslice, full_data["Nt"], start_S_samples
    )
    full_range_result["S_correlator_interpolation"] = interpolate_correlator_and_S(
        fit_samples, m_rho, m_a_1, min_timeslice, max_timeslice, start_S_samples
    )

    return full_range_result


def get_max_timeslice(data, max_timeslice):
    if max_timeslice is not None:
        return max_timeslice + 1

    return data["Conserved"]["V-A_renormalised_samples"].shape[-1]


def main():
    args = get_args()
    data = read_numpy(args.input_file)

    m_rho_data = read_numpy(args.input_m_rho)
    m_a_1_data = read_numpy(args.input_m_a_1)
    max_timeslice = get_max_timeslice(data, args.max_timeslice)
    result = fit(
        data,
        m_rho_data["fit_result_samples"]["rho_mass"],
        m_a_1_data["fit_result_samples"]["a_1_mass"],
        args.min_timeslice,
        max_timeslice,
    )
    dump_numpy(
        {
            **result,
            **{key: data[key] for key in ["mass", "Nt", "Nx", "Ny", "Nz", "bin_size"]},
            "m_rho": m_rho_data["fit_result"]["rho_mass"],
            "m_a_1": m_a_1_data["fit_result"]["a_1_mass"],
            "min_timeslice": args.min_timeslice,
            "max_timeslice": max_timeslice,
        },
        args.output_file,
    )


if __name__ == "__main__":
    main()
