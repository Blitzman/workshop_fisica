"""S1: write the EOS and the star profiles as solver-grade tables (CSV + sidecar metadata).

A profile table stores, at a grid of radii, the five state variables y = (m, p, Phi, theta, theta') of the background
ODE AND their derivatives y' = rhs(r, y). A consumer reconstructs y(r) with cubic Hermite interpolation (N1 does this).
The grid is the integrator's own step nodes subdivided K times; K is increased until the Hermite interpolant
reproduces the integrator's dense solution within `tol` at every midpoint (the achieved error is written to the
metadata and checked by run.py). Units: geometric (G = c = 1), lengths in km: m [km], p [km^-2], theta' [km^-1].
"""
import numpy as np
from scipy.interpolate import CubicHermiteSpline

from common import artifacts as A

from . import star as S

COLUMNS = ['r_km', 'm_km', 'p_km2', 'Phi', 'theta', 'thetap', 'dm', 'dp', 'dPhi', 'dtheta', 'dthetap']
UNITS = {'r_km': 'km', 'm_km': 'km', 'p_km2': 'km^-2', 'Phi': '1', 'theta': 'rad', 'thetap': 'km^-1', 'dm': '1',
         'dp': 'km^-3', 'dPhi': 'km^-1', 'dtheta': 'km^-1', 'dthetap': 'km^-2'}
DESC = {'r_km': 'areal radius', 'm_km': 'gravitational mass inside r (G=c=1)', 'p_km2': 'pressure (G=c=1)',
        'Phi': 'lapse potential, Phi(r_end) unshifted (add Phi_shift from the metadata)', 'theta': 'axion field a/f_a',
        'thetap': 'dtheta/dr', 'dm': 'dm/dr', 'dp': 'dp/dr', 'dPhi': 'dPhi/dr', 'dtheta': 'dtheta/dr (= thetap)',
        'dthetap': 'd^2 theta/dr^2'}
FLOORS = np.array([0.0, 0.0, 0.0, 1e-6, 1e-6])      # relative floor (x max|y_j|) for components that can cross zero


def _nodes(t, K):
    t = np.asarray(t)
    a, b = t[:-1], t[1:]
    sub = np.linspace(0.0, 1.0, K, endpoint=False)
    r = (a[:, None] + (b - a)[:, None] * sub[None, :]).ravel()
    return np.concatenate([r, t[-1:]])


def _rhs_at(st, r, Y):
    return np.array([st.rhs(float(ri), Y[:, i]) for i, ri in enumerate(r)]).T


def sample_star(st, tol=1e-6, K0=1, Kmax=128):
    """Returns (columns dict, info dict) with the Hermite-validated grid.

    K doubles from K0 until the interpolation error is <= tol or K reaches Kmax (then the achieved error is reported and
    run.py flags it). With the integrator at rtol = 1e-6 the steps are long, so up to ~64-128 subdivisions are needed
    (at Kmax = 64 the 7 kHz star reached only 1.6e-6)."""
    s = st.sol
    K = K0
    while True:
        r = _nodes(s.t, K)
        Y = s.sol(r)
        dY = _rhs_at(st, r, Y)
        mid = 0.5 * (r[1:] + r[:-1])
        ref = s.sol(mid)
        worst, worst_at = 0.0, None
        scale = np.abs(Y).max(axis=1)
        for j in range(5):
            h = CubicHermiteSpline(r, Y[j], dY[j])(mid)
            den = np.maximum(np.abs(ref[j]), FLOORS[j] * scale[j]) + 1e-300
            e = np.abs(h - ref[j]) / den
            k = int(np.argmax(e))
            if e[k] > worst:
                worst, worst_at = float(e[k]), (j, float(mid[k]))
        if worst <= tol or K >= Kmax:
            break
        K *= 2
    cols = {'r_km': r}
    for j, nm in enumerate(['m_km', 'p_km2', 'Phi', 'theta', 'thetap']):
        cols[nm] = Y[j]
    for j, nm in enumerate(['dm', 'dp', 'dPhi', 'dtheta', 'dthetap']):
        cols[nm] = dY[j]
    info = {'K': K, 'n_rows': int(len(r)), 'hermite_max_rel_err': worst, 'hermite_tol': tol,
            'worst_component_and_r': worst_at, 'n_solver_steps': int(len(s.t) - 1)}
    return cols, info


def star_parameters(st):
    return {'nu_a_kHz': st.nu_a, 'p0_MeV_fm3': st.p0, 'E0_geo_km2': st.E0, 'E0_MeV_fm3': S.E0, 'nc_fm3': st.nc,
            'beta': S.BETA, 'sigma_N_MeV': S.SIGMA_N, 'eps': S.EPS, 'm_pi_MeV': S.MPI, 'f_pi_MeV': S.FPI,
            'ma_per_km': st.ma, 'fa2_geo': (None if not np.isfinite(st.fa2) else st.fa2), 'theta_c': st.thc,
            'R_km': float(st.R), 'M_km': float(st.M), 'M_msun': float(st.M / S.MSUN), 'Phi_shift': float(st.Phi_shift),
            'r_end_km': float(st.sol.t[-1]), 'M_total_km': float(st.sol.y[0, -1]),
            'KGEO': S.KGEO, 'MSUN_km': S.MSUN, 'CKMS': S.CKMS, 'HBARC_MeV_fm': S.HBARC,
            'p_min_MeV_fm3': float(S._eos.p_min),
            'integrator': {'method': 'DOP853', 'rtol': st.rtol, 'atol': st.atol, 'max_step_km': st.max_step, 'start_r_km': 1e-5,
                           'shoot_tol': st.shoot_tol}}


def write_star(out_dir, st, name, module, inputs=(), tol=1e-6, k0=1):
    cols, info = sample_star(st, tol=tol, K0=k0)
    par = star_parameters(st)
    par['hermite'] = info
    meta = A.save_table(out_dir / f'{name}.csv', cols, name=name, units=UNITS, descriptions=DESC, parameters=par,
                        module=module, script='run.py', inputs=inputs,
                        title='Background profile (TOV + Klein-Gordon), state and derivatives')
    return meta, info


def write_eos(out_dir, module):
    e = S._eos
    cols = {'n_fm3': e.n, 'eps_MeV_fm3': e.e, 'p_MeV_fm3': e.p}
    par = {'name': 'BSk26 analytic fits (Potekhin et al., bskfit18.f)', 'n_min': 1e-8, 'n_max': 1.6, 'N_initial': 4000,
           'p_min_MeV_fm3': float(e.p_min), 'note': 'monotonic-in-p subset; consumers rebuild CubicSpline(ln p, ln eps), '
           '(ln p, ln n) exactly as bsk26.EOS does'}
    return A.save_table(out_dir / 'eos_bsk26.csv', cols, name='eos_bsk26',
                        units={'n_fm3': 'fm^-3', 'eps_MeV_fm3': 'MeV fm^-3', 'p_MeV_fm3': 'MeV fm^-3'},
                        descriptions={'n_fm3': 'baryon number density', 'eps_MeV_fm3': 'total energy density',
                                      'p_MeV_fm3': 'pressure'}, parameters=par, module=module, script='run.py',
                        title='BSk26 equation of state')
