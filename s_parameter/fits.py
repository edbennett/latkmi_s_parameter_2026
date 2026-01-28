#!/usr/bin/env python3


def fit_form_linear(mf, A, B):
    return A + B * mf


def fit_form_quadratic(mf, A, C):
    return A + C * mf**2
