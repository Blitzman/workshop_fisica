"""D1: radial linearization of the field equations from D0 (moved from axion_ns/model.py::linearize, unchanged).

Input : D0's field_equations.eq.json (loaded by run.py).  Output: the linear parts in the perturbations, with the
e^{-i w t} dependence factored out. This file does not know the equations of the papers.
"""
import sympy as sp

from common.notation import t, Phi, lam, a, n, xi, Phi0, lam0, a0, n0, dPhi, dlam, da1, dn1, xi1, w, eta

LINEAR_KEYS = ('Gtt', 'Gtr', 'Grr', 'KG', 'FLUID_r', 'BAR')


def linearize(expr):
    """Linear part in the perturbations, with the e^{-i w t} dependence factored out."""
    E = sp.exp(-sp.I * w * t)
    sub = {Phi: Phi0 + eta * dPhi * E, lam: lam0 + eta * dlam * E, a: a0 + eta * da1 * E,
           n: n0 + eta * dn1 * E, xi: eta * xi1 * E}
    lin = sp.diff(expr.subs(sub).doit(), eta).subs(eta, 0)
    return sp.expand(sp.powsimp(sp.expand(lin / E)))


def derive_linear(field_equations):
    return {k: linearize(field_equations[k]) for k in LINEAR_KEYS}
