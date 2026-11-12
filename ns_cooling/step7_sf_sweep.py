"""Step 7 (H7): superfluid cooling of SLy4 stars (isothermal interior) for every gap combination, against the three sources."""
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

from ns_cooling import envelope as E
from ns_cooling import neutrino as NU
from ns_cooling import sly4 as S
from ns_cooling import step3_integrals as I3
from ns_cooling import superfluid as SF
from ns_cooling.step4_standard_candle import YEAR
from ns_cooling.step5_cooling_curve import SIGMA
from ns_cooling.structure import volume_integral, with_mass
from ns_cooling.tov import C, G, MSUN, Table

HERE = Path(__file__).parent
RES = HERE / 'results'
MS = I3.MSTAR
VMAX = 200.0                                   # factors are below 1e-80 beyond this; avoids underflow
TGRID = np.logspace(6.0, 10.0, 240)            # redshifted interior temperature [K]


class Star:
    def __init__(self, pr, n_cc, fractions):
        self.pr, self.M, self.R = pr, pr['M'], pr['R']
        self.x = 2 * G * self.M * MSUN / (self.R * C**2)
        self.g14 = G * self.M * MSUN / (self.R**2 * np.sqrt(1 - self.x)) / 1e14
        n = pr['n']
        self.core = n >= n_cc
        xe, xmu = fractions(np.clip(n, n_cc, None))
        self.nn, self.np, self.ne, self.nmu = n * (1 - xe - xmu), n * (xe + xmu), n * xe, n * xmu
        self.kn, self.kp = (3 * np.pi**2 * self.nn)**(1 / 3), (3 * np.pi**2 * self.np)**(1 / 3)

    def tables(self, gap_p=None, gap_n=None, pbf='Y01', pbf_on=True):
        """L_nu^infinity and C on TGRID."""
        s, phi = self, self.pr['phi']
        tcp = SF.tc(s.kp, gap_p, 'A') if gap_p else np.zeros_like(s.kp)
        tcn = SF.tc(s.kn, gap_n, 'B') if gap_n else np.zeros_like(s.kn)
        L, Cap = np.empty_like(TGRID), np.empty_like(TGRID)
        for i, Tt in enumerate(TGRID):
            T = Tt * np.exp(-phi)
            t9 = T / 1e9
            with np.errstate(divide='ignore'):
                vp = np.minimum(SF.v_of_tau(np.where(tcp > 0, T / np.where(tcp > 0, tcp, 1.0), np.inf), 'A'), VMAX)
                vn = np.minimum(SF.v_of_tau(np.where(tcn > 0, T / np.where(tcn > 0, tcn, 1.0), np.inf), 'B'), VMAX)
            mu = NU.murca(s.nn, s.np, s.ne, s.nmu, t9, MS, MS)
            br = NU.brems(s.nn, s.np, t9, MS, MS)
            q = ((mu['n_e'] + mu['n_mu']) * SF.murca_n(vn, vp) + (mu['p_e'] + mu['p_mu']) * SF.murca_p(vn, vp)
                 + br['nn'] * SF.brems_nn(vn) + br['np'] * SF.brems_np(vn, vp) + br['pp'] * SF.brems_pp(vp))
            if pbf_on:
                q = q + np.where(vn > 0, SF.pbf(s.nn, t9, vn, MS, 'n', pbf), 0.0) + np.where(vp > 0, SF.pbf(s.np, t9, vp, MS, 'p', pbf), 0.0)
            c = (I3.c_fermi(s.nn, T, MS * I3.M_N) * SF.heat_factor(vn, 'B') + I3.c_fermi(s.np, T, MS * I3.M_N) * SF.heat_factor(vp, 'A')
                 + I3.c_fermi(s.ne, T, I3.lepton_mstar(s.ne, I3.M_E)) + I3.c_fermi(s.nmu, T, I3.lepton_mstar(s.nmu, I3.M_MU)))
            L[i] = volume_integral(s.pr, q * np.exp(2 * phi), s.core)
            Cap[i] = volume_integral(s.pr, c, s.core)
        return L, Cap

    def photon(self, T):
        ts = E.ts_iron(T / np.sqrt(1 - self.x), self.g14)
        return 4 * np.pi * self.R**2 * SIGMA * ts**4 * (1 - self.x)

    def luminosity(self, ages_yr, L, Cap, T0=1e10):
        lL, lC = CubicSpline(np.log(TGRID), np.log(L)), CubicSpline(np.log(TGRID), np.log(Cap))

        def rhs(t, y):
            T = np.exp(y[0])
            return [-(np.exp(lL(y[0])) + self.photon(T)) / np.exp(lC(y[0])) / T]
        t = np.asarray(ages_yr, float) * YEAR
        sol = solve_ivp(rhs, (0.0, t[-1]), [np.log(T0 * 0.999)], t_eval=t, method='LSODA', rtol=1e-7, atol=1e-9)
        assert sol.success, sol.message
        return self.photon(np.exp(sol.y[0]))


def run(masses=(1.0, 1.2, 1.4, 1.6, 1.8, 2.0)):
    src = json.loads((HERE / 'data' / 'marino2024_sources.json').read_text())
    ages = sorted({a for o in src['objects'] for a in o['age_yr']})
    tab, fr, ncc = Table(*S.eos()), I3.sly_fractions(), S.composition()[0][0]
    ref5 = json.loads((RES / 'step5_cooling_curve.json').read_text())['at_ages']['SLy4 1.4 iron']
    rows, base = [], {}
    for M in masses:
        st = Star(with_mass(tab, M, hi=15.44, N=600), ncc, fr)
        base[M] = st.luminosity(ages, *st.tables())
        for gp, gn in itertools.product([None] + list(SF.GAPS_P), [None] + list(SF.GAPS_NT)):
            if gp is None and gn is None:
                continue
            for pbf in ('Y01', 'updated'):
                Lg = st.luminosity(ages, *st.tables(SF.GAPS_P.get(gp), SF.GAPS_NT.get(gn), pbf))
                rows.append({'M': M, 'p': gp, 'n': gn, 'pbf': pbf, 'L': Lg.tolist()})
        if M == 1.4:
            zero = st.luminosity(ages, *st.tables((0.0,) + SF.GAPS_P['CCDK'][1:], (0.0,) + SF.GAPS_NT['TToa'][1:]))
            nopbf = {gn: st.luminosity(ages, *st.tables(None, SF.GAPS_NT[gn], pbf_on=False)) for gn in SF.GAPS_NT}
            withpbf = {gn: next(np.array(r['L']) for r in rows if r['M'] == M and r['p'] is None and r['n'] == gn and r['pbf'] == 'Y01') for gn in SF.GAPS_NT}
        print(f"  M = {M}: done", flush=True)
    ia = {a: i for i, a in enumerate(ages)}
    summary = []
    for o in src['objects']:
        a = o['age_yr'][1]                                    # the oldest allowed age is the most favourable for the model
        up = (o['L32'] + o['L32_err'][1]) * 1e32
        for pbf in ('Y01', 'updated'):
            rs = [r for r in rows if r['pbf'] == pbf]
            best = min(rs, key=lambda r: r['L'][ia[a]])
            summary.append({'object': o['name'], 'age': a, 'pbf': pbf, 'L_obs': o['L32'] * 1e32, 'L_upper': up,
                            'best': {k: best[k] for k in ('M', 'p', 'n')}, 'best_L': best['L'][ia[a]],
                            'best_ratio': best['L'][ia[a]] / (o['L32'] * 1e32),
                            'n_within': sum(r['L'][ia[a]] <= up for r in rs), 'n_total': len(rs),
                            'frac_colder': float(np.mean([r['L'][ia[a]] < base[r['M']][ia[a]] for r in rs])),
                            'nonsf_ratio': {str(M): float(base[M][ia[a]] / (o['L32'] * 1e32)) for M in masses}})
    d_base = max(abs(base[1.4][ia[a]] / ref5[str(a)] - 1) for a in ages)
    d_zero = float(np.max(np.abs(zero / base[1.4] - 1)))
    d_pbf = max(float(np.max(np.abs(nopbf[g] / withpbf[g] - 1))) for g in SF.GAPS_NT)
    checks = [
        {'name': 'With no gaps the sweep reproduces the step-5 luminosities (within 1 %)', 'ok': bool(d_base <= 0.01), 'note': f"max {d_base:.2%}"},
        {'name': 'Control: gaps scaled to zero give the non-superfluid curve', 'ok': bool(d_zero < 1e-6), 'note': f"max {d_zero:.1e}"},
        {'name': 'Control: switching pair-breaking emission off changes the luminosity by more than 20 % for some triplet model', 'ok': bool(d_pbf > 0.20), 'note': f"max {d_pbf:.0%}"},
    ] + [{'name': f"RESULT: {s['object']} at {s['age']} yr, pair-breaking {s['pbf']}", 'ok': None,
          'note': f"best model/observed {s['best_ratio']:.1f} ({s['best']['M']} Msun, p {s['best']['p']}, n {s['best']['n']}); within the upper error bar {s['n_within']}/{s['n_total']}; colder than non-superfluid {s['frac_colder']:.0%}"} for s in summary]
    out = {'source': 'SLy4, iron envelope, isothermal interior; factors of step 6', 'ages': ages, 'masses': list(masses),
           'base': {str(M): base[M].tolist() for M in masses}, 'rows': rows, 'summary': summary, 'checks': checks}
    (RES / 'step7_sf_sweep.json').write_text(json.dumps(out))
    for c in checks:
        print({True: 'PASS', False: 'FAIL', None: 'RESULT'}[c['ok']], c['name'], f"({c['note']})")
    return out


if __name__ == '__main__':
    run()
