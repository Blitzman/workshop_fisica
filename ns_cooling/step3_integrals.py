"""Step 3 (H3): L_nu and C_core of the isothermal core against the fits of Ofengeim et al. 2017 (arXiv:1612.04672)."""
import json
from pathlib import Path

import numpy as np
from scipy import constants as K
from scipy.interpolate import PchipInterpolator

from ns_cooling import bskfit as B
from ns_cooling import sly4 as S
from ns_cooling.step1_structure import table as bsk_table
from ns_cooling.structure import volume_integral, with_mass
from ns_cooling.tov import C, G, MSUN, Table

RES = Path(__file__).parent / 'results'
KB, HBAR = K.k * 1e7, K.hbar * 1e7                       # erg/K, erg s
M_N, M_E, M_MU = K.m_n * 1e3, K.m_e * 1e3, K.physical_constants['muon mass'][0] * 1e3   # g
N0 = 0.16                                               # fm^-3, Ofengeim et al. Eq. (2)
RHO0 = 2.8e14                                           # g/cm^3, their Sect. 4.2
# Ofengeim et al. 2017, Eqs. (2), (4), (7) and Table 2: prefactor, n, k, p, a1, a2, a3, gamma
FIT = {'nn': (1.77e19, 8, 1, 6, 2.03, 1.14, 0.0047, 2.51),
       'mu': (1.75e21, 8, 2, 6, 1.12, 1.14, 0.0060, 2.45),
       'c': (1.12e20, 1, 1, 1, 2.68, 1.14, 0.0174, 2.11)}
MAXERR = {'nn': 0.54, 'c': 0.12}                        # their Table 2, last column
MSTAR = 0.7                                             # m*/m used in their Eqs. (2)-(7)
T9 = 1.0


def fit(case, M, R):
    q0, n, k, p, a1, a2, a3, g = FIT[case]
    xr, xg = M * MSUN / (RHO0 * R**3), 2 * G * M * MSUN / (R * C**2)
    J = a1 * xr**(k / 3) * (1 - xg)**(-p / 2) * (1 + a3 * xr**(g - 1))**((p * g - k / 3) / (g - 1)) / np.sqrt(1 - a2 * xg)
    return q0 * R**3 * T9**n * J


def c_fermi(n_fm3, T, mstar_g):
    """Their Eq. (6): k_B^2 T m* p_F / (3 hbar^3), with p_F from the number density."""
    pf = HBAR * (3 * np.pi**2 * np.maximum(n_fm3, 0.0) * 1e39)**(1 / 3)
    return KB**2 * T * mstar_g * pf / (3 * HBAR**3)


def lepton_mstar(n_fm3, m):
    pf = HBAR * (3 * np.pi**2 * np.maximum(n_fm3, 0.0) * 1e39)**(1 / 3)
    return np.sqrt(m**2 + (pf / C)**2)


def integrals(pr, n_cc, fractions, redshift=True, uniform=False):
    n = pr['n']
    core = n >= n_cc
    if uniform:                                         # negative control: the whole star at the central density
        n, core = np.full_like(n, n[0]), np.ones_like(core)
    xe, xmu = fractions(np.clip(n, n_cc, None))
    nn, npr, ne, nmu = n * (1 - xe - xmu), n * (xe + xmu), n * xe, n * xmu
    phi = pr['phi'] if redshift else np.zeros_like(n)
    T = T9 * 1e9 * np.exp(-phi)
    t9 = T / 1e9
    red = np.exp(2 * phi)
    c = (c_fermi(nn, T, MSTAR * M_N) + c_fermi(npr, T, MSTAR * M_N)
         + c_fermi(ne, T, lepton_mstar(ne, M_E)) + c_fermi(nmu, T, lepton_mstar(nmu, M_MU)))
    return {'c': volume_integral(pr, c, core),
            'nn': volume_integral(pr, FIT['nn'][0] * (nn / N0)**(1 / 3) * t9**8 * red, core),
            'mu': volume_integral(pr, FIT['mu'][0] * (npr / N0)**(1 / 3) * t9**8 * red, core)}


def sly_fractions():
    n, xp, xe, xmu = S.composition()
    fe, fm = PchipInterpolator(n, xe), PchipInterpolator(n, xmu)
    return lambda x: (fe(np.clip(x, n[0], n[-1])), fm(np.clip(x, n[0], n[-1])))


def sweep(tab, n_cc, fractions, masses, hi, **kw):
    rows = []
    for M in masses:
        pr = with_mass(tab, M, hi=hi)                  # hi: log10 of a central density below the maximum-mass one
        I = integrals(pr, n_cc, fractions, **kw)
        rows.append({'M': M, 'R_km': pr['R'] / 1e5, **{k: I[k] for k in I},
                     **{f'{k}_ratio': I[k] / fit(k, M, pr['R']) for k in I}})
    return rows


def run():
    masses = [1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
    sly_tab, sly_fr, sly_ncc = Table(*S.eos()), sly_fractions(), S.composition()[0][0]
    sly = sweep(sly_tab, sly_ncc, sly_fr, masses, 15.44)
    sly_noz = sweep(sly_tab, sly_ncc, sly_fr, masses, 15.44, redshift=False)
    sly_uni = sweep(sly_tab, sly_ncc, sly_fr, masses, 15.44, uniform=True)
    bsk = sweep(bsk_table(24), B.n_cc(24), lambda x: B.fractions(24, x), masses + [2.2], 15.34)
    c0 = float(c_fermi(N0, 1e9, MSTAR * M_N))

    def dev(rows, k):
        return np.array([r[f'{k}_ratio'] - 1 for r in rows])
    stats = {name: {k: {'max': float(np.max(np.abs(dev(rows, k)))), 'rms': float(np.sqrt(np.mean(dev(rows, k)**2)))}
                    for k in ('c', 'nn', 'mu')} for name, rows in (('sly', sly), ('bsk24', bsk))}
    s = stats['sly']
    checks = [
        {'name': f"SLy4 heat capacity within {MAXERR['c']:.0%} of the fit at every mass", 'ok': bool(s['c']['max'] <= MAXERR['c']), 'note': f"max {s['c']['max']:.1%}, rms {s['c']['rms']:.1%}"},
        {'name': f"SLy4 nn-bremsstrahlung luminosity within {MAXERR['nn']:.0%} at every mass, rms at most 20 %", 'ok': bool(s['nn']['max'] <= MAXERR['nn'] and s['nn']['rms'] <= 0.20), 'note': f"max {s['nn']['max']:.1%}, rms {s['nn']['rms']:.1%}"},
        {'name': f"c0 from Eq. (6) with m* = 0.7 m_n at n0: {c0:.3e} vs {FIT['c'][0]:.2e} (within 1 %)", 'ok': bool(abs(c0 / FIT['c'][0] - 1) <= 0.01), 'note': f"{c0 / FIT['c'][0] - 1:+.2%}"},
        {'name': 'RESULT: SLy4 modified-Urca luminosity over the fit, with Omega = 1', 'ok': None, 'note': ', '.join(f"{r['mu_ratio']:.2f}" for r in sly)},
        {'name': 'RESULT: BSk24 (outside their calibration set), maximum deviations', 'ok': None, 'note': f"C {stats['bsk24']['c']['max']:.1%}, nn {stats['bsk24']['nn']['max']:.1%}"},
        {'name': 'NEG1: without redshift, L_nn leaves the 54 % band for at least one mass', 'ok': bool(np.max(np.abs(dev(sly_noz, 'nn'))) > MAXERR['nn']), 'note': f"max {np.max(np.abs(dev(sly_noz, 'nn'))):.0%}"},
        {'name': 'NEG2: a uniform-density star leaves the 12 % band in heat capacity', 'ok': bool(np.max(np.abs(dev(sly_uni, 'c'))) > MAXERR['c']), 'note': f"max {np.max(np.abs(dev(sly_uni, 'c'))):.0%}"},
    ]
    out = {'source': 'Ofengeim et al. 2017 (arXiv:1612.04672), Eqs. (2)-(11), (19)-(20), Table 2', 'T9': T9, 'c0': c0,
           'sly': sly, 'bsk24': bsk, 'stats': stats, 'checks': checks}
    RES.mkdir(exist_ok=True)
    (RES / 'step3_integrals.json').write_text(json.dumps(out, indent=1))
    for c in checks:
        print({True: 'PASS', False: 'FAIL', None: 'RESULT'}[c['ok']], c['name'], f"({c['note']})")
    for r in sly:
        print(f"  SLy {r['M']:.1f} R={r['R_km']:.2f}  C {r['c']:.3e} ({r['c_ratio']:.3f})  Lnn {r['nn']:.3e} ({r['nn_ratio']:.3f})  Lmu ({r['mu_ratio']:.3f})")
    return out


if __name__ == '__main__':
    run()
