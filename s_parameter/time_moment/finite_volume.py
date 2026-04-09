#!/usr/bin/env python3

from argparse import ArgumentParser, FileType
from functools import cache

import numpy as np
import pandas as pd
from scipy.integrate import quad_vec

from ..io import read_numpy, dump_numpy
from ..stats import generate_jackknife, jackknife_mean_variance


METADATA_KEYS = ["mass", "Nt", "Nx", "Ny", "Nz", "bin_size"]

taste_multiplicities = {
    None: 1,
    "45": 1,
    "i5": 3,
    "i4": 3,
    "ij": 3,  # 12, 13, 23 - others cancel or are duplicates
    "4": 1,
    "i": 3,
    "id": 1,
}
total_taste_multiplicity = sum(taste_multiplicities.values())


@cache
def wrapping_multiplicities():
    max_num_wrappings_squared = 100
    max_num_wrappings = int(max_num_wrappings_squared**0.5) + 1

    result = {
        num_wrappings_squared: 0
        for num_wrappings_squared in range(max_num_wrappings_squared + 1)
    }

    def add_single_contribution(direction_count, x_count, y_count=0, z_count=0):
        num_wrappings_squared = x_count**2 + y_count**2 + z_count**2
        if num_wrappings_squared <= max_num_wrappings_squared:
            result[num_wrappings_squared] += direction_count

    for x_count in range(1, max_num_wrappings):
        add_single_contribution(6, x_count)  # x only
        add_single_contribution(12, x_count, x_count)  # Equal x and y
        for z_count in range(1, max_num_wrappings):
            add_single_contribution(8, x_count, x_count, z_count)

        for y_count in range(x_count + 1, max_num_wrappings):
            add_single_contribution(24, x_count, y_count)  # Different x and y
            for z_count in range(1, max_num_wrappings):
                add_single_contribution(16, x_count, y_count, z_count)

    return result


class EnsembleManager:
    def __init__(self, datum, spectrum):
        self._datum = datum
        self._samples = {}

        result = spectrum.query(f"Nf == 8 & beta == 3.8 & mf == {datum['mass']}")
        if len(result) != 1:
            raise ValueError("Unique ensemble not found.")
        self._spectrum = result.iloc[0]

    @property
    def L(self):
        L = self._datum["Nx"]
        assert self._datum["Ny"] == L
        assert self._datum["Nz"] == L
        return L

    def get_samples(self, channel):
        if channel not in self._samples:
            value = self._spectrum[f"value_{channel}"]
            error = self._spectrum[f"error_{channel}"]
            self._samples[channel] = generate_jackknife(
                value, error, self._datum, self._datum["fit_result_samples"].shape[0]
            )

        return self._samples[channel]


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_file")
    parser.add_argument("--previous_data", required=True)
    parser.add_argument("--output_file", default="-", type=FileType("w"))
    return parser.parse_args()


def integrand(momentum, time, length, pi_mass, num_wrappings):
    """
    Inner portion of Eq. (49)
    """
    # Rightmost portion, setting \hat{M}_{\pi,\xi} = 1 as discussed in text,
    # cancelling the L, and multiplying the M_\pi through
    momentum_factor = momentum**2 + pi_mass**2
    summand = np.exp(-2 * time * (momentum_factor) ** 0.5) / momentum_factor

    return (
        momentum**3
        / (2 * np.pi**2)
        * np.sin(momentum * num_wrappings * length)
        / (num_wrappings * length)
        * total_taste_multiplicity
        * summand
    )


def half_infinite_sum(data):
    """
    Take data defined on the range [0, \\infty) computed with some arbitrary upper bound.
    Compute their sum,
    verifying with some basic checks that the summand has correctly decayed to zero.
    """
    result = np.sum(data, axis=0)

    # Verify there are sufficient data to have a tail at all
    count = len(data)
    assert count >= 40

    # Make it easier to split the data
    half_count = int(count // 2)
    quarter_count = int(count // 4)
    abs_data = np.abs(data)

    # This could be made more robust
    # by also imposing a lower bound relative to the start of the data,
    # but makes it fragile when things get very small in an inner infinite sum
    # For now, assert based on having looked at the intermediate data
    # that this threshold is sufficient
    zero_threshold = max(1e-15, abs_data[0].mean() / 1e8)

    # Ensure that the tail is decaying or zero
    try:
        assert abs_data[:half_count].sum() > abs_data[half_count:].sum() * 4
        assert abs_data[:-quarter_count].sum() > abs_data[-quarter_count:].sum() * 8
        assert (
            (
                (abs_data[half_count:-1] - abs_data[half_count + 1 :] > 0)
                | (abs_data[half_count + 1 :] < zero_threshold)
            )
            .mean(axis=1)
            .all()
        )
    except Exception:
        breakpoint()

    return result


def integral(time, length, pi_mass, num_wrappings_squared):
    result = quad_vec(
        lambda momentum: integrand(
            momentum, time, length, pi_mass, num_wrappings_squared**0.5
        ),
        a=0,
        b=np.inf,
    )
    return result[0]


def delta_fv_g_pi_pi(time, length, pi_mass, projected_decay_time):
    """
    The full Eq. (49)
    """
    result = (
        pi_mass**3
        / 3
        * half_infinite_sum(
            [
                wrapping_multiplicities()[num_wrappings_squared]
                * integral(time, length, pi_mass, num_wrappings_squared)
                for num_wrappings_squared in range(1, projected_decay_time)
            ]
        )
    )
    print(result)
    return result


def delta_fv_S(length, pi_mass, inner_projected_decay_time=42):
    """
    Eq. (52), without the $C$ factor.
    """
    projected_decay_time = 100
    data = np.array(
        [
            time**2
            * delta_fv_g_pi_pi(time, length, pi_mass, inner_projected_decay_time)
            for time in range(1, projected_decay_time)
        ]
    )
    return -2 * np.pi * half_infinite_sum(data)


def main():
    args = get_args()

    s_parameter_vp = read_numpy(args.input_file)
    spectrum = pd.read_csv(args.previous_data)
    ensemble = EnsembleManager(s_parameter_vp, spectrum)
    result = delta_fv_S(ensemble.L, ensemble.get_samples("mpi"))
    dump_numpy(
        {
            "delta_fv_S_samples": result,
            "delta_fv_S": jackknife_mean_variance(result),
            "S_infinite_t_samples": s_parameter_vp["S_infinite_t_samples"],
            "S_infinite_t": s_parameter_vp["S_infinite_t"],
            "pi_mass_samples": ensemble.get_samples("mpi"),
            **{key: s_parameter_vp[key] for key in METADATA_KEYS},
        },
        args.output_file,
    )


if __name__ == "__main__":
    main()
