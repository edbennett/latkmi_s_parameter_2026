#!/usr/bin/env python3

"""
Pade Ansatz f(q^2) = (b0 + b1 q^2)/(1 + c1 q^2 + c2 q^4);
see Eq. (9) of the paper.
"""


def pade(momentum_squared, b0, b1, c1, c2):
    return (b0 + b1 * momentum_squared) / (
        1 + c1 * momentum_squared + c2 * momentum_squared**2
    )
