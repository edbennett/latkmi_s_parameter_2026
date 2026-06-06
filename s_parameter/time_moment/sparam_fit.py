#!/usr/bin/env python3

from argparse import ArgumentParser, FileType
from functools import partial

import numpy as np
from scipy.optimize import curve_fit

from .sparam_tm import compute_S_parameter_contribution
from ..io import read_numpy, dump_numpy
from ..stats import (
    jackknife_mean_variance,
    jackknife_statistical_intermediary,
    jackknife_systematic_intermediary,
    sample_systematics,
)


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_file", metavar="input_file")
    parser.add_argument("--input_mass_samples", required=True)
    parser.add_argument("--min_timeslice", type=int, default=0)
    parser.add_argument("--max_timeslice", type=int, default=None)
    parser.add_argument("--m_rho_vt_samples", default=None)
    parser.add_argument("--m_a_1_samples", default=None)
    parser.add_argument("--output_file", type=FileType("w"), default="-")
    return parser.parse_args()


def fit_form(time, C_A_plus, C_A_minus, C_V_plus, C_V_minus, m_rho, m_a_1, max_time):
    """
    Fit form described in Eq. (26) of the paper.
    `time` is expected to be an integer.
    The four `C_` parameters are the fit parameters.
    `m_rho` and `m_a_1` are the masses of the rho and a_1 meson,
    and should be preset using `partial` or similar for fitting.
    `max_time` is the value $T$ from the paper,
    and should be preset using `partial` or similar for fitting.
    `max_time = None` corresponds to the infinite $T$ limit.
    """
    # C_A_plus = np.clip(C_A_plus, 0, 1)
    # C_A_minus = np.clip(C_A_minus, -1, 0)
    # C_V_plus = np.clip(C_V_plus, -1, 0)
    # C_V_minus = np.clip(C_V_minus, 0, 1)
    return (C_V_plus - (-1) ** time * C_A_plus) * (
        np.exp(-m_rho * time)
        + (0 if max_time is None else np.exp(-m_rho * (max_time - time)))
    ) - (C_A_minus - (-1) ** time * C_V_minus) * (
        np.exp(-m_a_1 * time)
        + (0 if max_time is None else np.exp(-m_a_1 * (max_time - time)))
    )


def correlator_infinite_t_I(t0, mass):
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
    and apply Eq. (27)/(28) of the paper to obtain
    the extrapolated S parameter at infinite T.
    fit_samples: The fit results for the four parameters C_{+,-}^{V,A}
    m_rho, m_a_1: Samples of the masses of the rho and a_1 states
    t_0: The start of the fit interval
    S_t_0: The cumulative summed value of S at the start of the fit interval
    """

    # Take transpose to match shape of fit_samples: (samples, parameters)
    factors = np.array(
        [
            -correlator_infinite_t_J(t_0, m_rho),
            -correlator_infinite_t_I(t_0, m_a_1),
            correlator_infinite_t_I(t_0, m_rho),
            correlator_infinite_t_J(t_0, m_a_1),
        ]
    ).T
    # Eq. (28) uses 4pi * (-1/2!); this is simplified to -2pi here.
    contributions = -np.pi / 2 * (fit_samples * factors)
    result = {
        name: contributions[..., idx]
        for idx, name in enumerate(["S_A_plus", "S_A_minus", "S_V_plus", "S_V_minus"])
    }
    S_samples = S_t_0 + contributions.sum(axis=1)
    result["S_infinite_t"] = S_samples

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
    factors = (
        np.cumsum(
            np.array([-1, 1, 1, -1])[:, np.newaxis] ** times
            * times**2
            * np.exp(-masses.T[:, :, np.newaxis] * times),
            axis=2,
        )
        * np.array([-1, -1, 1, 1])[:, np.newaxis]
    )
    S_eff_samples = S_t_0[:, np.newaxis] - np.pi / 2 * (
        fit_samples[..., np.newaxis] * factors
    ).sum(axis=1)
    return {
        "times": times,
        "S_effective": jackknife_mean_variance(S_eff_samples),
    }


def interpolate_correlator_and_S(fit_samples, m_rho, m_a_1, t_0, t_1, max_time, S_t_0):
    """
    Take a jackknife sample set of fit results,
    and apply Eq. (22) of the paper
    to obtain the interpolated effective S parameter as a function of T.
    """
    times = np.arange(t_0, t_1)[:, np.newaxis]
    correlator = -fit_form(times, *(np.array(fit_samples).T), m_rho, m_a_1, max_time).T

    S_eff_samples = (
        S_t_0 - np.cumsum(compute_S_parameter_contribution(correlator, t_0), axis=1).T
    ).T
    return {
        "times": times[:, 0],
        "correlator": jackknife_mean_variance(correlator),
        "S_effective": jackknife_mean_variance(S_eff_samples),
    }


def fit_single(full_data, m_rho, m_a_1, min_timeslice, max_timeslice):
    data = full_data["Conserved"]["V-A_renormalised_samples"]

    assert data.shape[:-1] == m_rho.shape
    assert data.shape[:-1] == m_a_1.shape
    data_values, data_uncertainty = jackknife_mean_variance(data)

    fit_samples = []
    chisquares = []
    times = np.arange(min_timeslice, max_timeslice)
    for sample, m_rho_sample, m_a_1_sample in zip(data, m_rho, m_a_1):
        result, _, info, _, _ = curve_fit(
            partial(
                fit_form,
                m_a_1=m_a_1_sample,
                m_rho=m_rho_sample,
                max_time=2 * (data.shape[-1] - 1),
            ),
            times,
            sample[min_timeslice:max_timeslice],
            sigma=data_uncertainty[min_timeslice:max_timeslice],
            bounds=([0, -1, -1, 0], [1, 0, 0, 1]),
            full_output=True,
        )
        fit_samples.append(result)
        chisquares.append((info["fvec"] ** 2).sum())

    start_S_samples = full_data["Conserved"]["S_parameter_eff_samples"][
        :, min_timeslice - 1
    ]
    result_samples = extrapolate_S_infinite_t(
        fit_samples,
        m_rho,
        m_a_1,
        min_timeslice,
        start_S_samples,
    )
    model_result_samples = extrapolate_S_infinite_t(fit_samples, m_rho, m_a_1, 1, 0)

    fit_samples_array = np.array(fit_samples)

    result = {
        "fit_result": jackknife_mean_variance(fit_samples),
        "fit_result_samples": result_samples,
        "chisquare": jackknife_mean_variance(chisquares),
        "dof": max_timeslice - min_timeslice + 1 - len(result),
        "S_infinite_t": jackknife_mean_variance(result_samples["S_infinite_t"]),
        "bare_result_samples": fit_samples_array,
        "min_timeslice": min_timeslice,
        "max_timeslice": max_timeslice,
    }

    for idx, name in enumerate(["A_plus", "A_minus", "V_plus", "V_minus"]):
        samples = fit_samples_array[..., idx]
        result_samples[f"C_{name}"] = samples
        result[f"C_{name}"] = jackknife_mean_variance(samples)
        result[f"S_{name}"] = jackknife_mean_variance(result_samples[f"S_{name}"])
        result_samples[f"model_S_{name}"] = model_result_samples[f"S_{name}"]
        result[f"model_S_{name}"] = jackknife_mean_variance(
            model_result_samples[f"S_{name}"]
        )

    return result


def get_chisquare(full_data, fit_samples, m_rho, m_a_1, min_timeslice, max_timeslice):
    data = full_data["Conserved"]["V-A_renormalised_samples"]
    data_values, data_uncertainty = jackknife_mean_variance(data)
    fit_values, _ = jackknife_mean_variance(fit_samples)
    times = np.arange(min_timeslice, max_timeslice)

    func_at_fit_result = fit_form(
        times, *fit_values, m_rho.mean(), m_a_1.mean(), 2 * (data.shape[-1] - 1)
    )
    return (
        (
            (func_at_fit_result - data_values[min_timeslice:max_timeslice])
            / data_uncertainty[min_timeslice:max_timeslice]
        )
        ** 2
    ).sum()


def add_systematic(result, all_samples, key):
    value, error = result[key]
    samples = np.array([sample["fit_result_samples"][key] for sample in all_samples])
    systematic = jackknife_systematic_intermediary(samples).mean()
    return value, error, systematic


def fit(full_data, m_rho, m_a_1, min_timeslice, max_timeslice):
    def fit_time_range(min_timeslice, max_timeslice):
        return fit_single(full_data, m_rho, m_a_1, min_timeslice, max_timeslice)

    result = {
        "fit_result_samples": {},
        "fit_result_systematic_samples": {},
    }
    all_samples = sample_systematics(
        fit_time_range,
        min_timeslice,
        max_timeslice,
        min_end_timeslice=full_data["Nt"] // 2 - 2,
    )

    for key in [
        *[
            f"{obs}_{channel}_{sign}"
            for obs in ["model_S", "S", "C"]
            for channel in ["V", "A"]
            for sign in ["plus", "minus"]
        ],
        "S_infinite_t",
    ]:
        obs_samples = np.array(
            [sample["fit_result_samples"][key] for sample in all_samples]
        )
        statistical_samples = jackknife_statistical_intermediary(obs_samples, "flat")
        result["fit_result_samples"][key] = statistical_samples

        value, error = jackknife_mean_variance(statistical_samples)

        systematic_samples = jackknife_systematic_intermediary(obs_samples)
        result["fit_result_systematic_samples"][key] = systematic_samples

        systematic_error = systematic_samples.mean()
        result[key] = (value, error, systematic_error)

    result["S_infinite_t_systematic_samples"] = jackknife_systematic_intermediary(
        np.array(
            [sample["fit_result_samples"]["S_infinite_t"] for sample in all_samples]
        )
    )
    fit_samples = jackknife_statistical_intermediary(
        np.array([sample["bare_result_samples"] for sample in all_samples]), "flat"
    )
    result["chisquare"] = get_chisquare(
        full_data, fit_samples, m_rho, m_a_1, min_timeslice, max_timeslice
    )
    result["dof"] = max_timeslice - min_timeslice + 1 - 4
    start_S_samples = full_data["Conserved"]["S_parameter_eff_samples"][
        :, min_timeslice - 1
    ]
    result["bare_result_samples"] = fit_samples
    result["S_extrapolation_large_t"] = extrapolate_S_finite_t(
        fit_samples, m_rho, m_a_1, min_timeslice, full_data["Nt"], start_S_samples
    )
    result["S_correlator_interpolation"] = interpolate_correlator_and_S(
        fit_samples,
        m_rho,
        m_a_1,
        min_timeslice,
        max_timeslice,
        full_data["Nt"],
        start_S_samples,
    )
    result["S_correlator_extrapolation"] = interpolate_correlator_and_S(
        fit_samples,
        m_rho,
        m_a_1,
        min_timeslice,
        full_data["Nt"] + 1,
        None,
        start_S_samples,
    )

    return result


def get_max_timeslice(data, max_timeslice):
    if max_timeslice is not None:
        return max_timeslice + 1

    return data["Conserved"]["V-A_renormalised_samples"].shape[-1]


def main():
    args = get_args()
    data = read_numpy(args.input_file)

    mass_data = read_numpy(args.input_mass_samples)
    m_rho_data = mass_data["rho_mass"]
    m_a_1_data = mass_data["a_1_mass"]
    if args.m_rho_vt_samples:
        m_rho_data = read_numpy(args.m_rho_vt_samples)["fit_result_samples"]["rho_mass"]
    if args.m_a_1_samples:
        m_a_1_data = read_numpy(args.m_a_1_samples)["fit_result_samples"]["a_1_mass"]

    max_timeslice = get_max_timeslice(data, args.max_timeslice)
    result = fit(
        data,
        m_rho_data,
        m_a_1_data,
        args.min_timeslice,
        max_timeslice,
    )
    dump_numpy(
        {
            **result,
            **{
                key: data[key]
                for key in ["mass", "Nt", "Nx", "Ny", "Nz", "Nf", "bin_size"]
            },
            "m_rho": m_rho_data,
            "m_a_1": m_a_1_data,
            "min_timeslice": args.min_timeslice,
            "max_timeslice": max_timeslice,
        },
        args.output_file,
    )


if __name__ == "__main__":
    main()
