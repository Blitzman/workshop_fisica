"""N3: background quantities for the numerical matrices, evaluated on arrays of stellar states.

The order follows systems.BG. Here x = n_b / n_c and V1 and V2 are derivatives
of the potential with respect to theta. Phi'' is obtained by analytically
differentiating Phi' = (m + 4 pi r^3 (p - V + kin)) / (r^2 F).
"""
import numpy as np


def background_arrays(st, r, m, p, Phi, t, tp):
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
        V1 = -E0 * g1
    else:
        g = g1 = g2 = zero
        e, q1, V, V1, kin = e0, zero, zero, zero, zero
    pr = np.where(inside, p, 0.0)
    eo = np.where(inside, e, 0.0)
    xi = np.where(inside, x, 0.0)
    dPhi = (m + 4 * np.pi * r**3 * (pr - V + kin)) / (r**2 * F)
    dm = 4 * np.pi * r**2 * (eo + V + kin)
    dp = np.where(inside, -(e + pr) * dPhi - q1 * tp, 0.0)
    if ax:
        tpp = (4 * st.ma**2 / b * (xi - 1) * g1 - 2 * tp / r * (1 - m / r - 2 * np.pi * r**2 * (eo - pr + 2 * V))) / F
        Fp = -(2 / r) * (dm - m / r)
        Pp = dp - V1 * tp + fa2 * (tp * tpp * F + tp**2 * Fp / 2)  # d/dr (p - V + kin)
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
        axt = (fav, -V, E0 * g1, E0 * g2, V, V1, -E0 * g2)
    else:
        axt = (np.ones_like(r), zero, zero, zero, zero, zero, zero)
        t = tp = zero
    return (Phi + st.Phi_shift, dPhi, P2, lam, dlam, xi, x1, t, tp, np.where(inside, e0, 0.0), pr, np.where(inside, cs0, 1.0), *axt)


def bgvals(st, rr):
    """The matrix background at one radius, from the table state (plain floats)."""
    m, p, Phi, t, tp = st.y_at(rr)
    return tuple(float(a[0]) for a in background_arrays(st, rr, m, p, Phi, t, tp))
