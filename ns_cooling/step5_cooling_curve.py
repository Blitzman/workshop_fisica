"""Step 5 (H5): cooling curve of the isothermal standard candle with the PCY97 envelope."""
import json
from pathlib import Path

import numpy as np
from scipy import constants as K
from scipy.integrate import solve_ivp

from ns_cooling import bskfit as B
from ns_cooling import envelope as E
from ns_cooling import sly4 as S
from ns_cooling import step3_integrals as I3
from ns_cooling.step1_structure import table as bsk_table
from ns_cooling.step4_standard_candle import YEAR, budget, t_interior
from ns_cooling.structure import with_mass
from ns_cooling.tov import C, G, MSUN, Table

HERE = Path(__file__).parent
RES = HERE / 'results'
SIGMA = K.Stefan_Boltzmann * 1e3              # erg cm^-2 s^-1 K^-4
TS_REF, TOL = 1.4e6, 0.07                     # Yakovlev et al. 2011, Sect. 4.3
T_C = 330.0


class Candle:
    def __init__(self, pr, n_cc, fractions):
        self.M, self.R = pr['M'], pr['R']
        self.x = 2 * G * self.M * MSUN / (self.R * C**2)
        self.g14 = G * self.M * MSUN / (self.R**2 * np.sqrt(1 - self.x)) / 1e14
        self.C9, L = budget(pr, n_cc, fractions)
        self.L9, self.Lparts = sum(L.values()), L

    def surface(self, T, env=E.ts_iron, redshift=True):
        """T_s [K] and L_gamma^infinity [erg/s] for interior redshifted temperature T."""
        Tb = T / np.sqrt(1 - self.x) if redshift else T
        ts = env(Tb, self.g14)
        return ts, 4 * np.pi * self.R**2 * SIGMA * ts**4 * (1 - self.x)

    def curve(self, t_yr, env=E.ts_iron, T0=1e10):
        def rhs(t, y):
            T = np.exp(y[0])
            return [-(self.L9 * (T / 1e9)**8 + self.surface(T, env)[1]) / (self.C9 * T / 1e9) / T]
        t = np.asarray(t_yr) * YEAR
        sol = solve_ivp(rhs, (0.0, t[-1]), [np.log(T0)], t_eval=t, method='LSODA', rtol=1e-8, atol=1e-10)
        assert sol.success, sol.message
        T = np.exp(sol.y[0])
        ts, lg = self.surface(T, env)
        return T, ts * np.sqrt(1 - self.x), lg


def run():
    src = json.loads((HERE / 'data' / 'marino2024_sources.json').read_text())
    sly_tab, sly_fr, sly_ncc = Table(*S.eos()), I3.sly_fractions(), S.composition()[0][0]
    rows = []
    for M in (1.0, 1.2, 1.4, 1.6, 1.8, 2.0):
        cd = Candle(with_mass(sly_tab, M, hi=15.44), sly_ncc, sly_fr)
        T, ts_inf, _ = cd.curve([T_C])
        T_an = t_interior(cd.C9, cd.Lparts, T_C)
        rows.append({'M': M, 'x': cd.x, 'g14': cd.g14, 'T': float(T[0]), 'T_analytic': T_an, 'Ts_inf': float(ts_inf[0]),
                     'Ts_inf_no_redshift': float(cd.surface(T[0], redshift=False)[0] * np.sqrt(1 - cd.x)),
                     'Ts_inf_accreted': float(cd.surface(T[0], env=E.ts_accreted)[0] * np.sqrt(1 - cd.x))})
    lo = [r for r in rows if r['x'] < 0.5]
    d = np.array([r['Ts_inf'] / TS_REF - 1 for r in lo])
    dan = np.array([r['T'] / r['T_analytic'] - 1 for r in rows])
    d1 = np.array([r['Ts_inf_no_redshift'] / TS_REF - 1 for r in lo])
    d2 = np.array([r['Ts_inf_accreted'] / TS_REF - 1 for r in lo])
    # curves for the comparison with the sources
    t = np.logspace(0, 6.5, 300)
    ages = sorted({a for o in src['objects'] for a in o['age_yr']})
    stars = {'SLy4 1.4': Candle(with_mass(sly_tab, 1.4, hi=15.44), sly_ncc, sly_fr),
             'BSk24 1.4': Candle(with_mass(bsk_table(24), 1.4, hi=15.34), B.n_cc(24), lambda n: B.fractions(24, n))}
    curves, at_ages = {}, {}
    for name, cd in stars.items():
        for label, env in (('iron', E.ts_iron), ('accreted', E.ts_accreted)):
            T, ts_inf, lg = cd.curve(t, env)
            curves[f'{name} {label}'] = {'t_yr': t.tolist(), 'T': T.tolist(), 'Ts_inf': ts_inf.tolist(), 'L': lg.tolist()}
            at_ages[f'{name} {label}'] = dict(zip(map(str, ages), cd.curve(ages, env)[2].tolist()))
    gaps = []
    for o in src['objects']:
        for name in stars:
            Ls = [at_ages[f'{name} iron'][str(a)] for a in o['age_yr']]
            gaps.append({'object': o['name'], 'star': name, 'L_model': Ls, 'L_obs': o['L32'] * 1e32,
                         'ratio': [L / (o['L32'] * 1e32) for L in Ls],
                         'accreted_over_iron': [at_ages[f'{name} accreted'][str(a)] / at_ages[f'{name} iron'][str(a)] for a in o['age_yr']]})
    checks = [
        {'name': f"SLy4 candles with x < 0.5 at {T_C:.0f} yr: redshifted surface temperature within {TOL:.0%} of 1.4 MK", 'ok': bool(np.max(np.abs(d)) <= TOL), 'note': ', '.join(f"{r['Ts_inf'] / 1e6:.2f}" for r in lo) + ' MK'},
        {'name': 'Numerical solution equals the analytic neutrino-stage one at 330 yr (within 1 %)', 'ok': bool(np.max(np.abs(dan)) <= 0.01), 'note': f"max {np.max(np.abs(dan)):.2%}"},
        {'name': 'NEG1: no redshift at the envelope leaves the 7 % band', 'ok': bool(np.max(np.abs(d1)) > TOL), 'note': f"{np.min(d1):+.0%} to {np.max(d1):+.0%}"},
        {'name': 'NEG2: accreted relation in place of iron leaves the 7 % band', 'ok': bool(np.min(np.abs(d2)) > TOL), 'note': f"{np.min(d2):+.0%} to {np.max(d2):+.0%}"},
    ] + [{'name': f"RESULT: {g['object']} vs {g['star']} iron, model / observed luminosity", 'ok': None,
          'note': ' to '.join(f"{x:.0f}" for x in sorted(set(g['ratio']))) + '; accreted / iron ' + ' to '.join(f"{x:.1f}" for x in sorted(set(g['accreted_over_iron'])))} for g in gaps]
    out = {'source': 'PCY97 Eqs. (A5)-(A10); Yakovlev et al. 2011 Sect. 4.3; Marino et al. 2024 Suppl. Table 6',
           'rows': rows, 'gaps': gaps, 'at_ages': at_ages, 'curves': curves, 'checks': checks}
    RES.mkdir(exist_ok=True)
    (RES / 'step5_cooling_curve.json').write_text(json.dumps(out))
    for c in checks:
        print({True: 'PASS', False: 'FAIL', None: 'RESULT'}[c['ok']], c['name'], f"({c['note']})")
    return out


if __name__ == '__main__':
    run()
