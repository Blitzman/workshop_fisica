"""D2: the extended system (Lindblom-Detweiler + axion) written compactly, to be VERIFIED against the derived system.

Moved from axion_ns/fullgr.py (block "sistema extendido con axion"). It is a CANDIDATE (kind = 'candidate'): LD
(corrected, with b and Q read from the metric: R's METFIX set) plus the NEW terms below, which are the new results of
the work. All the NEW terms vanish without axion (theta_0' = 0, dm*/dtheta = 0, V = 0, delta theta = 0). The system is
only published as an artifact after checks_fullgr.extended_checks() shows that every equation follows from the action.

EXT[k] = (LD equation k  [for 'KG': Eq. (8d) of Paper II])  -  NEW[k]
"""
import sympy as sp

from common.notation import (r, w, l, G, V, a0, n0, fa, Phi0, lam0, H0, H1, K, W, Vt, dth, d, ep, eL, dm, dmm, thp, Dtheta,
                             Pp)

V_theta = fa * sp.diff(V(a0), a0)
q1 = n0 * dm
POTENTIAL_TERMS = {
    'H0': 8 * sp.pi * G * r**2 * V_theta * dth,
    "X'": sp.exp(Phi0 + lam0) * q1 * V_theta * W / fa**2,
    'KG': -V_theta * H0 / fa**2,
}

NEW = {   # términos nuevos, en la forma  ecuación_LD = NEW
    "H1'": 8 * sp.pi * G * r**2 * eL * V(a0) * H1,
    "K'": 8 * sp.pi * G * r * fa**2 * thp * dth,
    "W'": sp.S(0),
    "X'": sp.exp(Phi0) * (n0 * dm * r * thp * (H0 / 2 + K) - l * (l + 1) * n0 * dm * thp * Vt / r
                          + n0 * dm * r * d(dth)
                          + (n0 * dm * (l + r * Pp) + n0 * dmm * r * thp - 4 * sp.pi * G * fa**2 * r * thp * ep) * dth
                          + n0 * sp.exp(-lam0) * (dm**2 * eL * n0 / fa**2 - 4 * dm * thp / r + dmm * thp**2) * W),
    'H0': -4 * sp.pi * G * r**2 * sp.exp(-2 * lam0) * d(a0)**2 * H0
    - 8 * sp.pi * G * r * sp.exp(-2 * lam0) * fa * d(a0) * (r * d(dth) + (l + 2 + r * Pp) * dth)
    - 8 * sp.pi * G * r * sp.exp(-lam0) * n0 * dm * thp * W,
    'V': -n0 * dm * sp.exp(2 * Phi0) * Dtheta / (w**2 * ep),
    'KG': -(n0 * dm / fa**2 - 2 * sp.exp(-2 * lam0) * Pp * thp) * H0
    + 16 * sp.pi * G * sp.exp(-2 * lam0) * fa**2 * thp**2 * dth,
}
NEW = {k: v + POTENTIAL_TERMS.get(k, 0) for k, v in NEW.items()}


def build_ext(LD_metfix, P8d):
    """LD_metfix: dict k -> LD equation (R: METFIX/*);  P8d: Eq. (8d) of Paper II (R: paper2_equations)."""
    return {k: (P8d if k == 'KG' else LD_metfix[k]) - v for k, v in NEW.items()}
