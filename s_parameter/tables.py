#!/usr/bin/env python3

import numpy as np

from format_multiple_errors import format_multiple_errors


def formatter(value_and_errors, **options):
    options = {
        "abbreviate": True,
        "length_control": "central",
        "significant_figures": 3,
        **options,
    }
    if None in value_and_errors:
        return "---"
    if np.isnan(value_and_errors).any():
        return "nan"
    return format_multiple_errors(*value_and_errors, **options)
