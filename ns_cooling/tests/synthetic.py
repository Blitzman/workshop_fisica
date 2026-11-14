"""Synthetic ``Points`` for tests that need a simple star without an equation of state."""

import numpy as np

from background import Points, derived_fermi_quantities, fermi_momentum

PRIMARY = ("rho", "ye", "yn", "yp", "xh", "aa", "zz", "ephi", "tc_crust", "tc_n", "tc_p", "ymu")


def point_fields(table):
    """All the fields of ``Points`` from the primary columns (ymu defaults to 0)."""
    f = {c: np.asarray(table[c], dtype=float) for c in PRIMARY if c in table}
    f.setdefault("ymu", np.zeros(len(f["rho"])))
    values = derived_fermi_quantities(f["rho"], f["ye"], f["yn"], f["yp"])
    values["kf_mu"] = fermi_momentum(f["rho"], f["ymu"])
    return {**f, **values}


def points_from_table(table):
    """``Points`` from the primary columns."""
    return Points(**point_fields(table))
