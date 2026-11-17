"""Equations of arXiv:2610.37308 (Paper II) transcribed as typed in the paper, form 'lhs - rhs'.

Moved from axion_ns/nonradial.py (the PAPER dict, unchanged). The shorthand (eps, p, dm, dmm, thp, meff2, dp,
de, ...) comes from the notation contract. These are REFERENCES: no derivation module may import this file.
"""
import sympy as sp

from common.notation import (r, w, l, fa, Phi0, lam0, a0, n0, dn, dth, W, Vt, eps, p, cs2, dm, dmm, thp, meff2, dp,
                             de, d)

PAPER = {
    '(8a) momento radial': d(dp) + (d(Phi0) + l / r) * dp - (
        w**2 * sp.exp(2 * (lam0 - Phi0)) * sp.exp(-lam0) / r * (eps + p) * W - de * d(Phi0)
        - n0 * dm * (d(dth) + l / r * dth) - thp * (dn * dm + n0 * dmm * dth)),
    '(8b) bariones': d((eps + p) * W) + (eps + p) * W * (d(Phi0) + (l + 1) / r)
    + r * sp.exp(lam0) * (dp / cs2 + (eps + p) * Vt * l * (l + 1) / r**2),
    '(8c) momento angular': (eps + p) * Vt + sp.exp(2 * Phi0) / w**2 * (dp + n0 * dm * dth),
    '(8d) KG interior': sp.exp(-2 * lam0) * d(dth, 2)
    + sp.exp(-2 * lam0) * (d(Phi0) - d(lam0) + 2 * (l + 1) / r) * d(dth)
    + (w**2 * sp.exp(-2 * Phi0) - meff2 + l / r * sp.exp(-2 * lam0) * (d(Phi0) - d(lam0))
       + l * (l + 1) / r**2 * (sp.exp(-2 * lam0) - 1)) * dth
    - dn / fa**2 * dm,
}
