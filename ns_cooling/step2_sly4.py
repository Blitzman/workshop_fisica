"""Step 2 (H2): SLy4 structure and direct-Urca thresholds against Douchin & Haensel 2001."""
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq, minimize_scalar

from ns_cooling import sly4 as S
from ns_cooling import tov
from ns_cooling.tov import Table, star

RES = Path(__file__).parent / 'results'


def threshold(lepton, muons=True):
    n, xp, xe, xmu = S.composition()
    if not muons:
        xe = xp
    xl = xe if lepton == 'e' else xmu
    f = PchipInterpolator(n, xp**(1 / 3) + xl**(1 / 3) - (1 - xp)**(1 / 3))
    return brentq(f, 0.2, 1.5)


def run():
    ref, (q_e, q_mu) = S.table6(), S.thresholds_quoted()
    tab = Table(*S.eos())

    def best(**kw):
        opt = minimize_scalar(lambda lr: -star(tab, 10**lr, **kw)[0], bounds=(15.0, 15.75), method='bounded', options={'xatol': 1e-5})
        return 10**opt.x, -opt.fun
    rho_c, m_max = best()
    _, r_max = star(tab, rho_c)
    P_c = float(tab.P_of_rho(rho_c))
    n_c = float(tab.n_of_P(P_c))
    _, m_newt = best(newtonian=True)
    lr2 = brentq(lambda lr: star(tab, 10**lr)[0] - 2.0, 15.0, np.log10(rho_c))
    n_c2 = float(tab.n_of_P(tab.P_of_rho(10**lr2)))
    lr14 = brentq(lambda lr: star(tab, 10**lr)[0] - 1.4, 14.5, 15.2)
    n_e, n_mu, n_e_nomu = threshold('e'), threshold('mu'), threshold('e', muons=False)
    ours = {'m_max': m_max, 'r_max': r_max, 'n_c': n_c, 'rho_c': rho_c, 'p_c': P_c}
    dev = {k: ours[k] / ref[k] - 1 for k in ref}
    checks = [
        {'name': f"Maximum mass {m_max:.3f} Msun vs {ref['m_max']} (within 0.5 %)", 'ok': bool(abs(dev['m_max']) <= 0.005), 'note': f"{100 * dev['m_max']:+.2f} %"},
        {'name': f"Radius at maximum mass {r_max:.2f} km vs {ref['r_max']} (within 0.5 %)", 'ok': bool(abs(dev['r_max']) <= 0.005), 'note': f"{100 * dev['r_max']:+.2f} %"},
        {'name': f"Central baryon density {n_c:.3f} fm-3 vs {ref['n_c']} (within 3 %)", 'ok': bool(abs(dev['n_c']) <= 0.03), 'note': f"{100 * dev['n_c']:+.2f} %"},
        {'name': f"Central pressure {P_c:.3e} vs {ref['p_c']:.2e} dyn/cm2 (within 3 %)", 'ok': bool(abs(dev['p_c']) <= 0.03), 'note': f"{100 * dev['p_c']:+.2f} %"},
        {'name': f"Electron direct-Urca threshold {n_e:.3f} fm-3 vs {q_e} (within 0.02)", 'ok': bool(abs(n_e - q_e) <= 0.02), 'note': ''},
        {'name': f"Muon direct-Urca threshold {n_mu:.3f} fm-3 vs {q_mu} (within 0.02)", 'ok': bool(abs(n_mu - q_mu) <= 0.02), 'note': ''},
        {'name': f"No direct Urca in any stable star: threshold {n_e:.3f} above the maximum central density {n_c:.3f} fm-3", 'ok': bool(n_e > n_c), 'note': f"2 Msun star: {n_c2:.3f} fm-3"},
        {'name': 'NEG1: Newtonian structure changes the maximum mass by more than 10 %', 'ok': bool(abs(m_newt / m_max - 1) > 0.10), 'note': f"{m_newt:.2f} Msun"},
        {'name': 'NEG2: ignoring muons moves the electron threshold by more than 0.05 fm-3', 'ok': bool(abs(n_e_nomu - n_e) > 0.05), 'note': f"{n_e_nomu:.3f} fm-3"},
    ]
    out = {'source': 'sly4.dat; Douchin & Haensel 2001 Tables 4, 6 and Sect. 7.6', 'ref': ref, 'ours': ours, 'dev': dev,
           'quoted_thresholds': {'e': q_e, 'mu': q_mu}, 'thresholds': {'e': n_e, 'mu': n_mu, 'e_no_muons': n_e_nomu},
           'n_c_2msun': n_c2, 'r_14': star(tab, 10**lr14)[1], 'n_c_14': float(tab.n_of_P(tab.P_of_rho(10**lr14))),
           'm_max_newtonian': m_newt, 'checks': checks}
    RES.mkdir(exist_ok=True)
    (RES / 'step2_sly4.json').write_text(json.dumps(out, indent=1))
    for c in checks:
        print('PASS' if c['ok'] else 'FAIL', c['name'], f"({c['note']})" if c['note'] else '')
    return out


if __name__ == '__main__':
    run()
