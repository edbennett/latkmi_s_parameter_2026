#!/usr/bin/env python3

"""
Pade Ansatz f(q^2) = (b0 + b1 q^2)/(1 + c1 q^2 + c2 q^4);
see Eq. (9) of the paper.
"""

import numpy as np


def pade(momentum_squared, b0, b1, c1, c2):
    return (b0 + b1 * momentum_squared) / (
        1 + c1 * momentum_squared + c2 * momentum_squared**2
    )


def get_momentum_filter(
    momentum,
    lattice_volume,
    momentum_squared_upper_bound=1,
    max_momentum_units_per_direction=2,
):
    """
    Return a boolean array of which rows of `momentum` meet the constraints:
    - Each component ≤ max_momentum_units_per_direction
                       times the unit momentum in any direction
    - Magnitude of momentum squared ≤ momentum_squared_upper_bound
    """
    momentum_squared = (momentum**2).sum(axis=1)

    brillouin_zone_size = 2 * np.pi

    momentum_integer_units = np.rint(
        momentum / (brillouin_zone_size / np.array(lattice_volume))
    )
    return (momentum_integer_units <= max_momentum_units_per_direction).all(axis=1) & (
        momentum_squared < momentum_squared_upper_bound
    )
