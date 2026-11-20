"""S1: write the EOS and the star profiles as tables (CSV + sidecar metadata).

A profile table is the dense STATE table of dense.py (ln p, m, Phi, theta, theta' at ~7,000 nodes built from the ln p solution, finer in r near the
surface); the metadata carry every constant the consumers need and the interpolation error measured at write time. Units: geometric (G = c = 1),
lengths in km: m [km], p [km^-2], theta' [km^-1]. R_km / M_km are the radius and mass at the hand-off point p = 1e-10 p_c (`r_k_km` is an alias
kept for the consumers that split their splines there).
"""
import numpy as np

from common import artifacts as A

from . import star as S
from .dense import COLUMNS, DESC, UNITS, dense_state, validate_dense


def star_parameters(st, rend):
    return {'nu_a_kHz': st.nu_a, 'p0_MeV_fm3': st.p0, 'E0_geo_km2': st.E0, 'E0_MeV_fm3': S.E0, 'nc_fm3': st.nc,
            'beta': S.BETA, 'sigma_N_MeV': S.SIGMA_N, 'eps': S.EPS, 'm_pi_MeV': S.MPI, 'f_pi_MeV': S.FPI,
            'ma_per_km': st.ma, 'fa2_geo': (None if not np.isfinite(st.fa2) else st.fa2), 'theta_c': st.thc,
            'R_km': float(st.R_stop), 'M_km': float(st.M_stop), 'M_msun': float(st.M_stop / S.MSUN), 'Phi_shift': float(st.Phi_shift),
            'r_end_km': float(rend), 'M_total_km': float(st.M_total), 'r_k_km': float(st.R_stop),
            'KGEO': S.KGEO, 'MSUN_km': S.MSUN, 'CKMS': S.CKMS, 'HBARC_MeV_fm': S.HBARC, 'p_min_MeV_fm3': float(S._eos.p_min),
            'integrator': {'variable': 'ln p (structure), r (vacuum tail)', 'method': 'DOP853', 'rtol': st.rtol, 'atol': st.atol,
                           'max_step_km': st.max_step, 'max_step_s': st.max_step_s, 'start_r_km': st.r0, 'p_frac': st.p_frac,
                           'shoot_tol': st.shoot_tol}}


def write_star(out_dir, st, name, module, inputs=(), tol=1e-5, **node_kw):
    cols, k = dense_state(st, **node_kw)
    err = validate_dense(st, cols, k)
    par = star_parameters(st, cols['r_km'][-1])
    par['dense'] = {'n_rows': int(len(cols['r_km'])), 'n_interior': int(k + 1), 'nodes': 'uniform in ln p + uniform in r + geometric near the centre',
                    'interpolation': 'two cubic splines (not-a-knot), split at R_stop; the tail has no ln p',
                    'interp_err': err, 'interp_max_err': max(err.values()), 'interp_tol': tol}
    meta = A.save_table(out_dir / f'{name}.csv', cols, name=name, units=UNITS, descriptions=DESC, parameters=par, module=module, script='run.py',
                        inputs=inputs, title='Background profile (TOV + Klein-Gordon): dense state table from the ln p integration')
    return meta, par['dense']


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
