"""
Tools for binning and jackknife analysis.
"""


def jackknife_mean_variance(samples):
    """
    Given a set of jackknife samples,
    compute the mean and standard deviation of the underlying data.
    """
    mean = samples.mean(axis=0)
    variance = (len(samples) - 1) / len(samples) * ((samples - mean) ** 2).sum(axis=0)
    return mean, variance**0.5


def sample_jackknife_ratio(numerator, denominator):
    """
    Given two broadcastable Numpy arrays,
    compute the jackknife samples for their ratio.
    """
    assert numerator.shape == denominator.shape
    full_axes = tuple(range(len(numerator.shape) - 1))
    jackknife_axes = tuple(range(1, len(numerator.shape) - 1))
    samples = (numerator.sum(axis=full_axes) - numerator.sum(axis=jackknife_axes)) / (
        denominator.sum(axis=full_axes) - denominator.sum(axis=jackknife_axes)
    )
    return samples


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
