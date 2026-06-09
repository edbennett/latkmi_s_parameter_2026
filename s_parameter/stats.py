"""
Tools for binning and jackknife analysis.
"""

import hashlib

from math import prod
import numpy as np


def jackknife_mean_error(samples, systematic_method=None):
    """
    Given a set of jackknife samples
    (sampled along axis 0, with other axes free),
    compute the mean and standard deviation of the underlying data.

    By default,
    assume that the samples contain only jackknife data;
    if a `systematic_method` is specified,
    then take the appropriate mean of the different systematic models
    before treating the jackknife sampling.
    """
    if systematic_method:
        if systematic_method != "max_deviation":
            raise NotImplementedError(f"{systematic_method} not currently implemented.")
        insufficient_samples = len(samples) <= 2
        samples = samples.mean(axis=0)

    mean = np.mean(samples, axis=0)
    variance = (len(samples) - 1) / len(samples) * ((samples - mean) ** 2).sum(axis=0)
    if systematic_method and insufficient_samples:
        variance = np.nan

    return mean, variance**0.5


def filter_fit_results(samples):
    """
    Apply progressively less strict selection criteria to samples
    until we find a set that we are happy with.
    """
    matched_samples = []
    matched_sample_data = []
    for chisquare_ubound, chisquare_lbound, dof_lbound, quantity_required in [
        (1.5, None, 2, 4),
        (2, None, 2, 2),
        (1.5, 0.5, None, 2),
        (2, 0.5, None, 2),
        (2.5, 0.5, None, 2),
    ]:
        for sample in samples:
            values, errors = sample["fit_result"]
            if (np.abs(errors / values) > 0.99).any():
                continue
            reduced_chisquare = sample["chisquare"][0] / sample["dof"]
            if reduced_chisquare >= chisquare_ubound:
                continue
            if chisquare_lbound is not None and reduced_chisquare <= chisquare_lbound:
                continue
            if dof_lbound is not None and sample["dof"] < dof_lbound:
                continue
            sample_descriptor = (sample["min_timeslice"], sample["max_timeslice"])
            if sample_descriptor in matched_samples:
                continue
            matched_samples.append(sample_descriptor)
            matched_sample_data.append(sample)
        if len(matched_samples) >= quantity_required:
            return matched_sample_data

    return [sorted(samples, key=lambda sample: sample["chisquare"][0] / sample["dof"])][
        0
    ]


def sample_systematics(func, min_timeslice, max_timeslice, min_end_timeslice=0):
    """
    Run `func` for a variety of fit intervals,
    and filter to those results anticipated to be the most reliable.
    `func` must accept two arguments,
    the start and end of the fit range,
    and return a dict containing minimally:
    - `"chisquare"` - the sum of residuals
    - `"dof"` - the number of degrees of freedom in the fit
    """
    fit_samples = []
    for start_timeslice in range(min_timeslice, max_timeslice - 3):
        for end_timeslice in range(
            max(min_end_timeslice, start_timeslice + 4), max_timeslice + 1
        ):
            try:
                fit_samples.append(func(start_timeslice, end_timeslice))
            except RuntimeError:
                # Not all fits work, and that's OK
                continue

    return filter_fit_results(fit_samples)


def jackknife_statistical_intermediary(samples, weighting="flat"):
    """
    Given a set of samples of a quantity estimated with different systematics
    (sampled along axis 0, with other axes free),
    and the final central value estimate,
    return the contribution to the statistics from each jackknife sample,
    weighting for the appropriate systematic sampling `weighting`.

    (Currently only `"flat"` is supported.)
    """
    if weighting != "flat":
        raise NotImplementedError(f"{weighting} not currently implemented.")
    return samples.mean(axis=0)


def jackknife_systematic_intermediary(samples, method="max_deviation"):
    r"""
    Given a set of samples of a quantity estimated with different systematics
    (sampled along axis 0, with other axes free),
    and the final central value estimate,
    return the contribution to the systematic error from each jackknife sample.

    (Currently only `"max_deviation"` is supported.
    This gives $\sqrt(1/N \max (sample - mean(samples)))$
    """
    if method != "max_deviation":
        raise NotImplementedError(f"{method} not currently implemented")

    if samples.shape[0] < 2:
        return np.array([[np.nan]])

    return np.array(
        [
            np.max(np.abs(sample - np.mean(sample, axis=0)), axis=0)
            for sample in samples.swapaxes(0, 1)
        ]
    )


def sample_jackknife(data, free_axes=-1):
    """
    Compute jackknife samples of data,
    for which specified axes are left free (e.g. a correlation function),
    and the first axis is the statistical samples over which to jackknife
    (e.g. Monte Carlo samples).
    Intermediary axes are quasi-independent observations of the same quantity.

    See e.g. https://en.wikipedia.org/wiki/Jackknife_resampling
    """

    num_axes = len(data.shape)
    if isinstance(free_axes, int):
        free_axes = [free_axes]
    if any(axis >= num_axes or axis <= -num_axes for axis in free_axes):
        raise ValueError(f"Invalid axes {free_axes}")

    normalised_free_axes = [axis % num_axes for axis in free_axes]

    sample_axis = 0
    full_axes = tuple(
        axis for axis in range(num_axes) if axis not in normalised_free_axes
    )
    observation_axes = tuple(axis for axis in full_axes if axis != sample_axis)
    observation_count = prod(
        [
            data.shape[axis] - 1 if axis == sample_axis else data.shape[axis]
            for axis in full_axes
        ]
    )
    return (
        data.sum(axis=full_axes) - data.sum(axis=observation_axes)
    ) / observation_count


def sample_jackknife_ratio(numerator, denominator):
    """
    Given two broadcastable Numpy arrays,
    compute the jackknife samples for their ratio.
    """
    assert numerator.shape == denominator.shape
    return sample_jackknife(numerator) / sample_jackknife(denominator)


def bin_data(data, bin_size):
    """
    Given a Numpy array representing observations of a time series,
    bin in all directions except the time direction.
    """
    data_count, *remaining_shape = data.shape
    bin_count = data_count // bin_size

    return (
        data[: bin_size * bin_count]
        .reshape(
            [bin_count, bin_size, *remaining_shape],
        )
        .mean(axis=1)
    )


def product_error_contribution(data, *keys):
    values = data[[f"value_{key}" for key in keys]]
    errors = data[[f"error_{key}" for key in keys]]
    contributions = (errors.to_numpy() / values.to_numpy()) ** 2
    return contributions.sum(axis=1) ** 0.5


def get_rng(data):
    """
    Get an RNG with a consistent seed for a given ensemble,
    so that the data are (more) reproducible,
    rather than introducing large fluctuations each time the samples are regenerated.
    """
    seed_string = "Nf8_mf{mass}_{Nx}x{Ny}x{Nz}x{Nt}_bin{bin_size}".format(**data)
    hash_string = hashlib.md5(seed_string.encode("utf8")).digest()
    seed = abs(int.from_bytes(hash_string, "big"))
    return np.random.default_rng(seed)


def generate_jackknife(mass, mass_error, data, shape):
    """
    The original correlator logs and jackknife samples
    for most of the ensembles we are considering for this work
    are lost to time.
    As such,
    we need to regenerate sample sets
    to propagate the error in the spectrum into the fit,
    since the latter must be done via a bootstrap.
    For consistency,
    we do this for all ensembles,
    even the mf=0.009 ensemble where the bootstrap samples are available
    in the data release to arXiv:2505.08658
    """
    rng = get_rng(data)
    if isinstance(shape, float):
        shape = [shape]
    distribution_std = mass_error / (np.prod(shape) - 1) ** 0.5
    return rng.normal(mass, distribution_std, shape)


def add_quadrature(*values):
    total = 0
    for value in values:
        if isinstance(value, tuple):
            if len(value) == 2:
                # Use unusual ordering so that (value, error) pairs work
                denominator, numerator = value
                total += (numerator / denominator) ** 2
            elif len(value) == 1:
                total += value[0] ** 2
            else:
                raise NotImplementedError("Unsupported error format")
        else:
            total += value**2
    return total**0.5
