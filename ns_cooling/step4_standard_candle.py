"""Step 4 (H4): standard slow cooling against Yakovlev et al. 2011 Eq. (14) and the Ofengeim et al. 2017 MU fit."""
import json
from pathlib import Path

import numpy as np

from ns_cooling import neutrino as NU
from ns_cooling import sly4 as S
from ns_cooling import step3_integrals as I3
from ns_cooling.structure import volume_integral, with_mass
from ns_cooling.tov import C, G, MSUN, Table

RES = Path(__file__).parent / 'results'
YEAR = 3.15576e7
T_C, A14, B14 = 330.0, 3.45e8, 0.12            # Yakovlev et al. 2011, Eq. (14)
TOL_T, TOL_MU = 0.07, 0.66                     # their stated maximum errors
MS = I3.MSTAR


def t_sc(M, R, t_yr):
    x = 2 * G * M * MSUN / (R * C**2)
    return A14 * (1 - x) * (1 + B14 * (R / 1e6)**2) * (T_C / t_yr)**(1 / 6)


def budget(pr, n_cc, fractions, redshift=True, murca=True):
    """C and L_nu^infinity at T~ = 1e9 K; both scale exactly (C ~ T, L ~ T^8)."""
    n = pr['n']
    core = n >= n_cc
    xe, xmu = fractions(np.clip(n, n_cc, None))
    nn, npr, ne, nmu = n * (1 - xe - xmu), n * (xe + xmu), n * xe, n * xmu
    phi = pr['phi'] if redshift else np.zeros_like(n)
    T = 1e9 * np.exp(-phi)
    c = (I3.c_fermi(nn, T, MS * I3.M_N) + I3.c_fermi(npr, T, MS * I3.M_N)
         + I3.c_fermi(ne, T, I3.lepton_mstar(ne, I3.M_E)) + I3.c_fermi(nmu, T, I3.lepton_mstar(nmu, I3.M_MU)))
    parts = {**NU.brems(nn, npr, T / 1e9, MS, MS)}
    if murca:
        parts.update(NU.murca(nn, npr, ne, nmu, T / 1e9, MS, MS))
    L = {k: volume_integral(pr, q * np.exp(2 * phi), core) for k, q in parts.items()}
    return volume_integral(pr, c, core), L


def t_interior(Cap, L, t_yr):
    q = sum(L.values()) / (Cap * (1e9)**7)       # K^-6 s^-1, their Eq. (9) with n = 8
    return (6 * q * t_yr * YEAR)**(-1 / 6)


def run():
    masses = [1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
    tab, fr, ncc = Table(*S.eos()), I3.sly_fractions(), S.composition()[0][0]
    rows = []
    for M in masses:
        pr = with_mass(tab, M, hi=15.44)
        Cap, L = budget(pr, ncc, fr)
        mu = sum(v for k, v in L.items() if k[:2] in ('n_', 'p_'))
        row = {'M': M, 'R_km': pr['R'] / 1e5, 'C': Cap, 'L': L, 'L_mu': mu, 'L_tot': sum(L.values()),
               'mu_ratio': mu / I3.fit('mu', M, pr['R']),
               'T': t_interior(Cap, L, T_C), 'T_ref': t_sc(M, pr['R'], T_C)}
        row['T_brems_only'] = t_interior(*budget(pr, ncc, fr, murca=False), T_C)
        row['T_no_redshift'] = t_interior(*budget(pr, ncc, fr, redshift=False), T_C)
        rows.append(row)
    dT = np.array([r['T'] / r['T_ref'] - 1 for r in rows])
    dmu = np.array([r['mu_ratio'] - 1 for r in rows])
    d1 = np.array([r['T_brems_only'] / r['T_ref'] - 1 for r in rows])
    d2 = np.array([r['T_no_redshift'] / r['T_ref'] - 1 for r in rows])
    pre = 7.5e19 * MS**4 * NU.A_NN * NU.B_NN * NU.N_NU
    checks = [
        {'name': f"Transcription: Eq. (165) with m* = 0.7 m gives {pre:.3e}, Ofengeim et al. quote {I3.FIT['nn'][0]:.2e} (within 1 %)", 'ok': bool(abs(pre / I3.FIT['nn'][0] - 1) <= 0.01), 'note': f"{pre / I3.FIT['nn'][0] - 1:+.2%}"},
        {'name': f"SLy4 interior temperature at {T_C:.0f} yr within {TOL_T:.0%} of Yakovlev et al. 2011 Eq. (14) at every mass", 'ok': bool(np.max(np.abs(dT)) <= TOL_T), 'note': ', '.join(f"{x:+.1%}" for x in dT)},
        {'name': f"Full modified-Urca luminosity within {TOL_MU:.0%} of the Ofengeim et al. fit at every mass", 'ok': bool(np.max(np.abs(dmu)) <= TOL_MU), 'note': ', '.join(f"{r['mu_ratio']:.2f}" for r in rows)},
        {'name': 'NEG1: bremsstrahlung only leaves the 7 % band', 'ok': bool(np.min(np.abs(d1)) > TOL_T), 'note': f"{np.min(d1):+.0%} to {np.max(d1):+.0%}"},
        {'name': 'NEG2: no redshift leaves the 7 % band', 'ok': bool(np.max(np.abs(d2)) > TOL_T), 'note': f"{np.min(d2):+.0%} to {np.max(d2):+.0%}"},
    ]
    out = {'source': 'Yakovlev et al. 2001 Eqs. (140), (142), (165)-(167); Yakovlev et al. 2011 Eq. (14); Ofengeim et al. 2017 Table 2',
           't_yr': T_C, 'rows': rows, 'dT': dT.tolist(), 'checks': checks}
    RES.mkdir(exist_ok=True)
    (RES / 'step4_standard_candle.json').write_text(json.dumps(out, indent=1))
    for c in checks:
        print('PASS' if c['ok'] else 'FAIL', c['name'], f"({c['note']})")
    for r in rows:
        sh = {k: v / r['L_tot'] for k, v in r['L'].items()}
        print(f"  {r['M']:.1f}  T {r['T']:.3e} ref {r['T_ref']:.3e}  shares " + ' '.join(f"{k}:{v:.2f}" for k, v in sh.items()))
    return out


if __name__ == '__main__':
    run()
