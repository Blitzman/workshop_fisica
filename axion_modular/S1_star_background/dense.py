"""S1: dense STATE table of a star, the hand-over format read by N1.

Only the smooth state y = (ln p, m, Phi, theta, theta') is tabulated against r; everything derived from it (energy density,
sound speed, couplings, Phi', lambda', m') is computed by the consumer from the state and the EOS table, exactly as the
solver does. Tabulating the derived quantities instead was tried and rejected: c_s^2 and friends inherit the EOS structure of
the crust and miss 1e-5 even with the best grid.

Nodes: `n_in` points on [R_MIN, r_k] clustered toward the surface (the pressure varies exponentially in the crust and the
right-hand side error is dominated by the last ~0.4 km), then `n_out` uniform points on [r_k, r_end] (the axion tail).
r_k is the radius where the EOS is switched off (p = 1.0001 p_min), just outside R; it is a node (the last interior one and the
first of the tail). The consumer interpolates every column with TWO cubic splines (not-a-knot), one on each side of r_k
(see split_splines for why): `validate_dense` measures the error of exactly that interpolation against the integrator's dense
solution and the achieved error is written to the metadata.
"""
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

from . import star as S

R_MIN = 1e-3            # km: where the mode integrations start (N1: interior r0 = 1e-3)
CLUSTER = 5.0           # exponent of the surface clustering

COLUMNS = ['r_km', 'ln_p', 'm_km', 'Phi', 'theta', 'thetap']
UNITS = {'r_km': 'km', 'ln_p': '1', 'm_km': 'km', 'Phi': '1', 'theta': 'rad', 'thetap': 'km^-1'}
DESC = {'r_km': 'areal radius (node)', 'ln_p': 'natural log of the pressure in geometric units (km^-2)',
        'm_km': 'gravitational mass inside r (G=c=1)', 'Phi': 'lapse potential, Phi(r_end) unshifted (add Phi_shift from the metadata)',
        'theta': 'axion field a/f_a', 'thetap': 'dtheta/dr'}


def surface_clustered(a, b, n, c=CLUSTER):
    """n nodes on [a, b], strongly clustered toward b."""
    u = np.linspace(0.0, 1.0, n)
    return a + (b - a) * (1 - (np.exp(c * (1 - u)) - 1) / (np.exp(c) - 1))


def kink_radius(st):
    """Radius where the EOS is switched off (p = 1.0001 p_min): the right-hand side has a kink there."""
    thr = S._eos.p_min * 1.0001
    f = lambda r: st.sol.sol(r)[1] / S.KGEO - thr
    rr = np.linspace(st.R, st.sol.t[-1], 4000)
    idx = np.nonzero(st.sol.sol(rr)[1] / S.KGEO - thr < 0)[0]
    if len(idx) == 0:
        return float(st.R + 1e-3)
    return float(brentq(f, rr[idx[0] - 1], rr[idx[0]], xtol=1e-12))


def dense_state(st, n_in=2000, n_out=500, cluster=CLUSTER):
    """Returns (columns, r_k). The node r_k is the LAST interior node and also the first node of the tail."""
    rk = kink_radius(st)
    r_in = surface_clustered(R_MIN, rk, n_in, cluster)
    r_in[-1] = rk
    r = np.concatenate([r_in, np.linspace(rk, st.sol.t[-1], n_out)[1:]])
    m, p, Phi, t, tp = st.sol.sol(r)
    return {'r_km': r, 'ln_p': np.log(p), 'm_km': m, 'Phi': Phi, 'theta': t, 'thetap': tp}, rk


def split_splines(cols, rk, spline=CubicSpline):
    """The interpolation the consumer uses: one spline on [R_MIN, r_k] and one on [r_k, r_end]. They must be separate: ln p has a
    slope of ~ -270 /km just inside r_k and is constant outside, and one spline across that jump rings (p wrong by up to 20x
    in the first tail intervals, where p ~ 1e-15)."""
    r = cols['r_km']
    Y = np.column_stack([cols[c] for c in COLUMNS[1:]])
    i = int(np.searchsorted(r, rk))
    assert r[i] == rk, 'r_k must be a node of the table'
    return spline(r[:i + 1], Y[:i + 1], axis=0), spline(r[i:], Y[i:], axis=0)


def validate_dense(st, cols, rk, seed=11):
    """Relative interpolation error of the stored state (the two cubic splines vs the integrator's dense solution) at ~28,000
    radii (uniform, near the surface, near the center, around the axion wall, in the tail). Floors: 1e-6 of the column maximum."""
    rng = np.random.default_rng(seed)
    s_in, s_out = split_splines(cols, rk)
    wall = np.linspace(R_MIN, st.R, 4000)
    rw = wall[np.argmax(np.abs(st.sol.sol(wall)[4]))]
    top = rk * (1 - 1e-9)
    rend = st.sol.t[-1]
    rt = np.concatenate([rng.uniform(R_MIN, top, 12000), rng.uniform(st.R - 0.6, top, 4000), rng.uniform(R_MIN, 0.5, 2000),
                         np.clip(rw + rng.uniform(-0.5, 0.5, 4000), R_MIN, top), rng.uniform(rk, rend, 6000),
                         rk + 10 ** rng.uniform(-4, np.log10(rend - rk), 4000)])
    m, p, Phi, t, tp = st.sol.sol(rt)
    ref = {'p': p, 'm': m, 'Phi': Phi, 'theta': t, 'theta_p': tp}
    tab = np.where((rt <= rk)[:, None], s_in(np.minimum(rt, rk)), s_out(np.maximum(rt, rk)))
    got = {'p': np.exp(tab[:, 0]), 'm': tab[:, 1], 'Phi': tab[:, 2], 'theta': tab[:, 3], 'theta_p': tab[:, 4]}
    err = {}
    for k in ref:
        scale = np.abs(ref[k]).max()
        err[k] = 0.0 if scale == 0 else float(np.max(np.abs(got[k] - ref[k]) / np.maximum(np.abs(ref[k]), 1e-6 * scale)))
    return err
