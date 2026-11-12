"""Step 8 (H8): in-medium modified Urca (Shternin, Baldo & Haensel 2018) and its effect on the sweep of step 7."""
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from ns_cooling import neutrino as NU
from ns_cooling import sly4 as S
from ns_cooling import step3_integrals as I3
from ns_cooling import superfluid as SF
from ns_cooling.step7_sf_sweep import MS, Star
from ns_cooling.structure import with_mass
from ns_cooling.tov import Table

HERE = Path(__file__).parent
RES = HERE / 'results'


def quadrature_dev():
    n, xp, xe, _ = S.composition()
    kn, kp, ke = [(3 * np.pi**2 * n * x)**(1 / 3) for x in (1 - xp, xp, xe)]
    dev = []
    for i in np.where(kn > kp + ke)[0][::4]:
        pn, pp, pl = (k[i] * NU.HBARC for k in (kn, kp, ke))
        mu, m = np.hypot(NU.ME_MEV, pl), MS * NU.MN_MEV
        qa = 2 * m**2 * mu**2 / pl * quad(lambda k: 1 / (pp**2 - k**2)**2, pn - pl, pn + pl, epsrel=1e-12)[0]
        qb = 2 * m**2 * mu**2 / (pl * pp) * quad(lambda k: k / (pn**2 - k**2)**2, pp - pl, pp + pl, epsrel=1e-12)[0]
        RA, RB = NU.r_dir(kn[i], kp[i], ke[i], MS, MS, NU.ME_MEV)
        dev += [abs(RA / qa - 1), abs(RB / qb - 1)]
    return float(max(dev))


def exponents(lower=False):
    kp, ke = 0.9, 0.7
    d = np.logspace(-4, -3, 20)
    kn = kp + ke + d
    RA, RB = NU.r_dir(kn, kp, ke, MS, MS, NU.ME_MEV, lower_A=kn if lower else None)
    return float(np.polyfit(np.log(d), np.log(RA), 1)[0]), float(np.polyfit(np.log(d), np.log(RB), 1)[0])


def run(masses=(1.0, 1.2, 1.4, 1.6, 1.8, 2.0)):
    src = json.loads((HERE / 'data' / 'marino2024_sources.json').read_text())
    h7 = json.loads((RES / 'step7_sf_sweep.json').read_text())
    ages = h7['ages']
    ia = {a: i for i, a in enumerate(ages)}
    tab, fr, ncc = Table(*S.eos()), I3.sly_fractions(), S.composition()[0][0]
    qd, (ea, eb), (ea_c, _) = quadrature_dev(), exponents(), exponents(lower=True)
    rows, base, alpha = [], {}, {}
    for M in masses:
        st = Star(with_mass(tab, M, hi=15.44, N=600), ncc, fr)
        ke = (3 * np.pi**2 * st.ne)**(1 / 3)
        au = NU.alpha_u(st.kn[st.core], st.kp[st.core], ke[st.core], MS, MS, NU.ME_MEV)
        alpha[M] = [float(au.min()), float(au.max())]
        base[M] = st.luminosity(ages, *st.tables(inmedium=True))
        for gp, gn in itertools.product([None] + list(SF.GAPS_P), [None] + list(SF.GAPS_NT)):
            if gp is None and gn is None:
                continue
            for pbf in ('Y01', 'updated'):
                Lg = st.luminosity(ages, *st.tables(SF.GAPS_P.get(gp), SF.GAPS_NT.get(gn), pbf, inmedium=True))
                rows.append({'M': M, 'p': gp, 'n': gn, 'pbf': pbf, 'L': Lg.tolist()})
        if M == 1.4:                                         # control: alpha_U forced to the standard value
            keep = NU.alpha_u
            NU.alpha_u = lambda *a: NU.ALPHA_N * NU.BETA_N
            forced = st.luminosity(ages, *st.tables(SF.GAPS_P['CCDK'], SF.GAPS_NT['TToa'], 'Y01', inmedium=True))
            NU.alpha_u = keep
            ref = next(np.array(r['L']) for r in h7['rows'] if r['M'] == M and r['p'] == 'CCDK' and r['n'] == 'TToa' and r['pbf'] == 'Y01')
            d_forced = float(np.max(np.abs(forced / ref - 1)))
        print(f"  M = {M}: done", flush=True)
    summary = []
    for o in src['objects']:
        a = o['age_yr'][1]
        up = (o['L32'] + o['L32_err'][1]) * 1e32
        for pbf in ('Y01', 'updated'):
            rs = [r for r in rows if r['pbf'] == pbf]
            best = min(rs, key=lambda r: r['L'][ia[a]])
            old = next(s for s in h7['summary'] if s['object'] == o['name'] and s['pbf'] == pbf)
            summary.append({'object': o['name'], 'age': a, 'pbf': pbf, 'best': {k: best[k] for k in ('M', 'p', 'n')},
                            'best_ratio': best['L'][ia[a]] / (o['L32'] * 1e32), 'best_ratio_h7': old['best_ratio'],
                            'n_within': sum(r['L'][ia[a]] <= up for r in rs), 'n_total': len(rs),
                            'frac_colder': float(np.mean([r['L'][ia[a]] < base[r['M']][ia[a]] for r in rs])),
                            'nonsf_ratio': {str(M): float(base[M][ia[a]] / (o['L32'] * 1e32)) for M in masses},
                            'nonsf_dimming': {str(M): float(h7['base'][str(M)][ia[a]] / base[M][ia[a]]) for M in masses}})
    dim = [x for s in summary for x in s['nonsf_dimming'].values()]
    checks = [
        {'name': 'Closed forms of Eqs. (8)-(9) agree with quadrature (1e-8)', 'ok': bool(qd < 1e-8), 'note': f"max {qd:.1e}"},
        {'name': 'Near the threshold both factors scale as the inverse distance: exponents within 0.05 of -1', 'ok': bool(abs(ea + 1) <= 0.05 and abs(eb + 1) <= 0.05), 'note': f"A {ea:.3f}, B {eb:.3f}"},
        {'name': 'NEG1: starting integral (8) at p_Fn removes the scaling', 'ok': bool(abs(ea_c + 1) > 0.05), 'note': f"exponent {ea_c:.3f}"},
        {'name': 'NEG2: alpha_U forced to alpha_n beta_n reproduces step 7', 'ok': bool(d_forced < 1e-6), 'note': f"max {d_forced:.1e}"},
        {'name': 'RESULT: alpha_U (electron branch) in the cores of the SLy4 stars', 'ok': None, 'note': ', '.join(f"{M}: {a[0]:.1f}-{a[1]:.1f}" for M, a in alpha.items())},
        {'name': 'RESULT: dimming of the non-superfluid star at the source ages', 'ok': None, 'note': f"factor {min(dim):.2f} to {max(dim):.2f}"},
    ] + [{'name': f"RESULT: {s['object']} at {s['age']} yr, pair-breaking {s['pbf']}", 'ok': None,
          'note': f"best model/observed {s['best_ratio']:.1f} (step 7: {s['best_ratio_h7']:.1f}; {s['best']['M']} Msun, p {s['best']['p']}, n {s['best']['n']}); within the upper error bar {s['n_within']}/{s['n_total']}; colder than non-superfluid {s['frac_colder']:.0%}"} for s in summary]
    out = {'source': 'Shternin, Baldo & Haensel 2018 Eqs. (8), (9), (12), (14); sweep as in step 7', 'ages': ages, 'masses': list(masses),
           'alpha': {str(M): a for M, a in alpha.items()}, 'exponents': [ea, eb], 'base': {str(M): base[M].tolist() for M in masses},
           'rows': rows, 'summary': summary, 'checks': checks}
    (RES / 'step8_inmedium.json').write_text(json.dumps(out))
    for c in checks:
        print({True: 'PASS', False: 'FAIL', None: 'RESULT'}[c['ok']], c['name'], f"({c['note']})")
    return out


if __name__ == '__main__':
    run()
