"""S1: dense STATE table of a star, the hand-over format read by N1 and N3.

Only the smooth state y = (ln p, m, Phi, theta, theta') is tabulated against r; everything derived from it is computed by the consumer from
the state and the EOS table. Tabulating derived quantities (c_s^2, ...) against r was tried and rejected (they inherit the EOS structure of the crust).

The star is integrated against s = ln p (star_logp.py); the nodes follow the solution, so the table is much finer in r near the surface:
  interior [R_MIN, R_stop]: N_s nodes uniform in s = ln p (dense in r near the surface) + N_r nodes uniform in r + N_geom nodes geometric
      from R_MIN to 1 km (the centre is a regular singular point: N3 tabulates r^2 A and needs them), centre nodes (r < 10 m) from the series;
  tail (R_stop, r_end]: N_tail nodes clustered toward R_stop (the axion tail, no fluid: ln p = -inf, stored as '-inf').
R_stop is where p has dropped to 1e-10 p_c (about 9 m before the EOS floor); it is the last interior node and the first of the tail. Consumers use
TWO cubic splines, one on each side (ln p is not defined in the tail) and `validate_dense` measures exactly that interpolation.
"""
import numpy as np
from scipy.interpolate import CubicSpline

R_MIN = 1e-3            # km: where the mode integrations start

COLUMNS = ['r_km', 'ln_p', 'm_km', 'Phi', 'theta', 'thetap']
UNITS = {'r_km': 'km', 'ln_p': '1', 'm_km': 'km', 'Phi': '1', 'theta': 'rad', 'thetap': 'km^-1'}
DESC = {'r_km': 'areal radius (node)', 'ln_p': 'natural log of the pressure in geometric units (km^-2); -inf in the fluid-free tail',
        'm_km': 'gravitational mass inside r (G=c=1)', 'Phi': 'lapse potential, Phi(r0) = 0 unshifted (add Phi_shift from the metadata)',
        'theta': 'axion field a/f_a', 'thetap': 'dtheta/dr'}


def table_nodes(st, n_s=4000, n_r=2000, n_geom=200, n_tail=1000, tail_c=4.0):
    """Returns (nodes, index of the node R_stop). All nodes strictly increasing."""
    r_s = st.sol_in.sol(np.linspace(st.sol_in.t[0], st.sol_in.t[-1], n_s))[0]
    r_in = np.unique(np.concatenate([r_s, np.linspace(R_MIN, st.R_stop, n_r), np.geomspace(R_MIN, 1.0, n_geom), [R_MIN]]))
    r_in = r_in[(r_in >= R_MIN) & (r_in < st.R_stop - 1e-10)]
    r_in = np.append(r_in, st.R_stop)
    r_end = st.r_end if st.nu_a else st.R_stop + 5.0
    u = np.linspace(0, 1, n_tail)[1:]
    r_tail = st.R_stop + (r_end - st.R_stop) * (np.exp(tail_c * u) - 1) / (np.exp(tail_c) - 1)
    return np.concatenate([r_in, r_tail]), len(r_in) - 1


def dense_state(st, **kw):
    """Returns (columns, index of R_stop)."""
    r, k = table_nodes(st, **kw)
    m, p, Phi, t, tp = st.state_at_r(r)
    with np.errstate(divide='ignore'):
        lnp = np.where(np.arange(len(r)) <= k, np.log(np.where(p > 0, p, 1.0)), -np.inf)
    return {'r_km': r, 'ln_p': lnp, 'm_km': m, 'Phi': Phi, 'theta': t, 'thetap': tp}, k


def split_splines(cols, k, spline=CubicSpline):
    """The interpolation the consumers use: interior [R_MIN, R_stop] (5 columns) and tail [R_stop, r_end] (4 columns: no ln p)."""
    r = cols['r_km']
    Y = np.column_stack([cols[c] for c in COLUMNS[1:]])
    return spline(r[:k + 1], Y[:k + 1], axis=0), spline(r[k:], Y[k:, 1:], axis=0)


def validate_dense(st, cols, k, seed=11):
    """Relative interpolation error of the stored state (the two splines vs the exact state_at_r) at ~28,000 radii (uniform, near the surface,
    near the centre, around the axion wall, in the tail). Floors: 1e-6 of the column maximum."""
    rng = np.random.default_rng(seed)
    s_in, s_out = split_splines(cols, k)
    R = st.R_stop
    wall = np.linspace(R_MIN, R, 4000)
    rw = wall[np.argmax(np.abs(st.state_at_r(wall)[4]))]
    rend = cols['r_km'][-1]
    rt = np.concatenate([rng.uniform(R_MIN, R, 12000), rng.uniform(R - 0.6, R, 4000), rng.uniform(R - 0.05, R, 4000), rng.uniform(R_MIN, 0.5, 2000),
                         np.clip(rw + rng.uniform(-0.5, 0.5, 4000), R_MIN, R), rng.uniform(R, rend, 6000),
                         R + 10 ** rng.uniform(-3, np.log10(rend - R), 4000)])
    m, p, Phi, t, tp = st.state_at_r(rt)
    inside = rt <= R
    tab_in = s_in(np.minimum(rt, R))
    tab_out = s_out(np.maximum(rt, R))
    got = {'p': np.where(inside, np.exp(tab_in[:, 0]), 0.0), 'm': np.where(inside, tab_in[:, 1], tab_out[:, 0]),
           'Phi': np.where(inside, tab_in[:, 2], tab_out[:, 1]), 'theta': np.where(inside, tab_in[:, 3], tab_out[:, 2]),
           'theta_p': np.where(inside, tab_in[:, 4], tab_out[:, 3])}
    ref = {'p': p, 'm': m, 'Phi': Phi, 'theta': t, 'theta_p': tp}
    # floors (fraction of the column maximum below which the error is judged in absolute terms). m, Phi - Phi(r0) and theta' vanish at the centre
    # like r^3, r^2 and r, and the centre series meets the ODE solution at r0 = 10 m with a kink of ~1e-7 relative to those tiny values (about
    # 3e-8 of the column maximum in absolute terms): for them the error is measured against max(|value|, 1% of the maximum). p and theta do not
    # vanish there and use 1e-6 of the maximum.
    floors = {'p': 1e-6, 'm': 1e-2, 'Phi': 1e-2, 'theta': 1e-6, 'theta_p': 1e-2}
    err = {}
    for key in ref:
        sel = inside if key == 'p' else np.ones_like(inside)
        scale = np.abs(ref[key][sel]).max()
        err[key] = 0.0 if scale == 0 else float(np.max(np.abs(got[key][sel] - ref[key][sel]) / np.maximum(np.abs(ref[key][sel]), floors[key] * scale)))
    return err
