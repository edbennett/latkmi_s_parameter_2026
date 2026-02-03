"""
Tools for binning and jackknife analysis.
"""

from math import prod


def jackknife_mean_variance(samples):
    """
    Given a set of jackknife samples
    (sampled along axis 0, with other axes free),
    compute the mean and standard deviation of the underlying data.
    """
    mean = samples.mean(axis=0)
    variance = (len(samples) - 1) / len(samples) * ((samples - mean) ** 2).sum(axis=0)
    return mean, variance**0.5


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
    observation_count = prod([data.shape[axis] for axis in full_axes])

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
