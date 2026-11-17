"""S1: write the EOS and the star profiles as tables (CSV + sidecar metadata).

A profile table is the dense STATE table of dense.py (ln p, m, Phi, theta, theta' on ~2,500 nodes, about 150 kB). The metadata
carry every constant the consumer needs and the interpolation error measured at write time. Units: geometric (G = c = 1),
lengths in km: m [km], p [km^-2], theta' [km^-1].
"""
import numpy as np

from common import artifacts as A

from . import star as S
from .dense import COLUMNS, DESC, UNITS, dense_state, validate_dense


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


def write_star(out_dir, st, name, module, inputs=(), tol=1e-5, n_in=2000, n_out=500):
    cols, rk = dense_state(st, n_in, n_out)
    err = validate_dense(st, cols, rk)
    par = star_parameters(st)
    par['r_k_km'] = rk
    par['dense'] = {'n_in': n_in, 'n_out': n_out, 'n_rows': int(len(cols['r_km'])), 'interior_grid': 'surface-clustered',
                    'cluster_exponent': 5.0, 'tail_grid': 'uniform', 'interpolation': 'cubic spline, not-a-knot, per column',
                    'interp_err': err, 'interp_max_err': max(err.values()), 'interp_tol': tol}
    meta = A.save_table(out_dir / f'{name}.csv', cols, name=name, units=UNITS, descriptions=DESC, parameters=par, module=module,
                        script='run.py', inputs=inputs, title='Background profile (TOV + Klein-Gordon): dense state table')
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
