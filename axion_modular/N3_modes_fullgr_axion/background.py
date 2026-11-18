"""N3: the 18 background quantities of the numerical matrices, from the state (vectorised, exact arithmetic).

Same quantities and order as axion_ns/axgr.bgvals: (Phi + shift, Phi', Phi'', lambda, lambda', x, x', theta, theta', e0, p, c_s0^2, f_a, Mst, M1, M2, V, V2),
with x = n_b / n_c and e0, p, c_s0^2 zero / 1 outside the fluid. The only difference with the legacy: Phi'' is the analytic derivative of
Phi' = (m + 4 pi r^3 (p - V + kin)) / (r^2 F) (from dm/dr, dp/dr, theta'' of the same state) instead of a central finite difference with h = 1e-5 km, which
is meaningless on tables whose nodes are 1.7 cm apart near the surface. tests/test_n3.py checks it against the finite difference.
"""
import numpy as np


def bg18(st, r, m, p, Phi, t, tp):
    r, m, p, Phi, t, tp = (np.atleast_1d(np.asarray(a, float)) for a in (r, m, p, Phi, t, tp))
    KG, eos = st.KGEO, st.eosv
    pthr = eos.p_min * 1.0001
    ax = bool(st.nu_a)
    pm = p / KG
    inside = pm > pthr
    pms = np.where(inside, pm, pthr)
    e0, nb, cs0 = eos.eps(pms) * KG, eos.nb(pms), eos.cs2(pms)
    x = nb / st.nc
    F = 1 - 2 * m / r
    zero = np.zeros_like(r)
    if ax:
        b, E0, fa2 = st.beta, st.E0, st.fa2
        g = np.sqrt(1 - b * np.sin(t / 2)**2)
        g1 = -b * np.sin(t) / (4 * g)
        g2 = -b * np.cos(t) / (4 * g) - b * b * np.sin(t)**2 / (16 * g**3)
        e, q1 = e0 - E0 * x * (1 - g), E0 * x * g1
        V, kin = E0 * (1 - g), fa2 * tp**2 / 2 * F
    else:
        g = g1 = g2 = zero
        e, q1, V, kin = e0, zero, zero, zero
    pr = np.where(inside, p, 0.0)
    eo = np.where(inside, e, 0.0)
    xi = np.where(inside, x, 0.0)
    dPhi = (m + 4 * np.pi * r**3 * (pr - V + kin)) / (r**2 * F)
    dm = 4 * np.pi * r**2 * (eo + V + kin)
    dp = np.where(inside, -(e + pr) * dPhi - q1 * tp, 0.0)
    if ax:
        tpp = (4 * st.ma**2 / b * (xi - 1) * g1 - 2 * tp / r * (1 - m / r - 2 * np.pi * r**2 * (eo - pr + 2 * V))) / F
        Fp = -(2 / r) * (dm - m / r)
        Pp = dp + E0 * g1 * tp + fa2 * (tp * tpp * F + tp**2 * Fp / 2)             # d/dr (p - V + kin), V' = -E0 g1 theta'
    else:
        Pp = dp
    P = pr - V + kin
    N = m + 4 * np.pi * r**3 * P
    Np = dm + 12 * np.pi * r**2 * P + 4 * np.pi * r**3 * Pp
    D = r**2 - 2 * m * r
    Dp = 2 * r - 2 * m - 2 * r * dm
    P2 = (Np * D - N * Dp) / D**2
    lam, dlam = -0.5 * np.log(F), (dm * r - m) / (r * (r - 2 * m))
    x1 = np.where(inside, dp / np.where(inside, (e0 + p) * cs0 / np.where(x == 0, 1, x), 1.0), 0.0)
    if ax:
        fav = np.sqrt(fa2) * np.ones_like(r)
        axt = (fav, -E0 * (1 - g), E0 * g1, E0 * g2, E0 * (1 - g), -E0 * g2)
    else:
        axt = (np.ones_like(r), zero, zero, zero, zero, zero)
        t = tp = zero
    return (Phi + st.Phi_shift, dPhi, P2, lam, dlam, xi, x1, t, tp, np.where(inside, e0, 0.0), pr, np.where(inside, cs0, 1.0), *axt)


def bgvals(st, rr):
    """The same 18 quantities at one radius, from the table state (plain floats)."""
    m, p, Phi, t, tp = st.y_at(rr)
    return tuple(float(a[0]) for a in bg18(st, rr, m, p, Phi, t, tp))
