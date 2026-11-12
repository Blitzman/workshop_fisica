"""Step 1 (H1): BSk24 structure and direct-Urca threshold against the numbers quoted in Marino et al. 2024."""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from ns_cooling import bskfit as B
from ns_cooling.tov import Table, star

RES = Path(__file__).parent / 'results'
# Marino et al. 2024, Methods (paper:1422-1425) and main text (paper:237)
PAPER = {'m_max': 2.279, 'rho_c_max': 2.26e15, 'rho_durca': 8.25e14, 'm_durca': 1.6}


def table(keos):
    n = np.logspace(-9, np.log10(2.0), 6000)
    rho = B.rho_of_n(keos, n)
    return Table(n, rho, 10**B.lgP(keos, np.log10(rho)))


def durca_margin(keos, n, muons=True):
    """> 0 where the Fermi momenta close the triangle p_Fn <= p_Fp + p_Fe."""
    ye, ymu = B.fractions(keos, n)
    yp = ye + ymu if muons else ye
    return yp**(1 / 3) + ye**(1 / 3) - (1 - yp)**(1 / 3)


def analyse(keos, muons=True):
    tab = table(keos)
    opt = minimize_scalar(lambda lr: -star(tab, 10**lr)[0], bounds=(15.0, 15.6), method='bounded', options={'xatol': 1e-5})
    rho_max, m_max = 10**opt.x, -opt.fun
    n_d = brentq(lambda n: durca_margin(keos, n, muons), B.n_cc(keos) * 1.5, 1.5)
    rho_d = float(B.rho_of_n(keos, n_d))
    m_d, r_d = star(tab, rho_d)
    rho14 = brentq(lambda lr: star(tab, 10**lr)[0] - 1.4, 14.5, 15.2)
    return {'m_max': m_max, 'rho_c_max': rho_max, 'n_durca': n_d, 'rho_durca': rho_d, 'm_durca': m_d, 'r_durca': r_d,
            'r_14': star(tab, 10**rho14)[1], 'yp_durca': float(sum(B.fractions(keos, n_d)))}


def run():
    r24, r24_nomu, r26 = analyse(24), analyse(24, muons=False), analyse(26)
    dev = {k: r24[k] / PAPER[k] - 1 for k in ('m_max', 'rho_c_max', 'rho_durca')}
    checks = [
        {'name': f"BSk24 maximum mass {r24['m_max']:.3f} Msun vs {PAPER['m_max']} (within 0.5 %)", 'ok': bool(abs(dev['m_max']) <= 0.005), 'note': f"{100 * dev['m_max']:+.2f} %"},
        {'name': f"Central density at maximum mass {r24['rho_c_max']:.3e} vs {PAPER['rho_c_max']:.2e} g/cm3 (within 3 %)", 'ok': bool(abs(dev['rho_c_max']) <= 0.03), 'note': f"{100 * dev['rho_c_max']:+.2f} %"},
        {'name': f"Direct-Urca threshold density {r24['rho_durca']:.3e} vs {PAPER['rho_durca']:.2e} g/cm3 (within 3 %)", 'ok': bool(abs(dev['rho_durca']) <= 0.03), 'note': f"{100 * dev['rho_durca']:+.2f} %"},
        {'name': f"Direct-Urca threshold mass {r24['m_durca']:.3f} Msun, paper says about {PAPER['m_durca']} (1.55-1.65)", 'ok': bool(1.55 <= r24['m_durca'] <= 1.65), 'note': ''},
        {'name': 'NEG1: without muons in the triangle the threshold density moves more than 3 %', 'ok': bool(abs(r24_nomu['rho_durca'] / r24['rho_durca'] - 1) > 0.03), 'note': f"{100 * (r24_nomu['rho_durca'] / r24['rho_durca'] - 1):+.1f} %"},
        {'name': 'NEG2: BSk26 coefficients miss the BSk24 maximum mass by more than 1 %', 'ok': bool(abs(r26['m_max'] / PAPER['m_max'] - 1) > 0.01), 'note': f"{r26['m_max']:.3f} Msun"},
    ]
    out = {'source': 'bskfit18.f (Pearson et al. 2018); paper numbers: Marino et al. 2024 Methods', 'paper': PAPER,
           'bsk24': r24, 'bsk24_no_muons': r24_nomu, 'bsk26': r26, 'dev': dev, 'checks': checks}
    RES.mkdir(exist_ok=True)
    (RES / 'step1_structure.json').write_text(json.dumps(out, indent=1))
    for c in checks:
        print('PASS' if c['ok'] else 'FAIL', c['name'], f"({c['note']})" if c['note'] else '')
    return out


if __name__ == '__main__':
    run()
