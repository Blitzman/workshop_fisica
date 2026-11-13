"""Shared validation of numerical inputs; no I/O and no plotting."""

import numpy as np


def checked(x, name, positive=False):
    """Return x as a float array, requiring finite and non-negative (or positive) values."""
    x = np.asarray(x, dtype=float)
    if not np.isfinite(x).all() or np.any(x <= 0 if positive else x < 0):
        raise ValueError(
            f'{name}: finite and {"positive" if positive else "non-negative"} values required'
        )
    return x
