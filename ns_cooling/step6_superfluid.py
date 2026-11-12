"""Step 6 (H6): superfluid factors against their printed asymptotes and Fig. 19 of Yakovlev et al. 2001."""
import json
from pathlib import Path

import numpy as np

from ns_cooling import sly4 as S
from ns_cooling import superfluid as SF
from ns_cooling.tov import Table

RES = Path(__file__).parent / 'results'
TAU = 0.1
# Y01 low-temperature asymptotes as printed: (name, function of tau, fitted function of v, pairing type)
ASYM = [
    ('(209) Mn branch, protons A', lambda t: 1.166e-4 / t**7.5 * np.exp(-1.764 / t), SF.mn_pA, 'A'),
    ('(210) Mp branch, protons A', lambda t: 0.00764 / t**7 * np.exp(-3.528 / t), SF.mp_pA, 'A'),
    ('(214) Mp branch, neutrons B', lambda t: 3.99e-6 / t**7 * np.exp(-1.188 / t), SF.mp_nB, 'B'),
    ('(218) np bremsstrahlung, protons A', lambda t: 0.910 / t * np.exp(-1.764 / t), SF.np_pA, 'A'),
    ('(219) pp bremsstrahlung, protons A', lambda t: 0.671 / t**2 * np.exp(-3.528 / t), SF.pp_pA, 'A'),
    ('(240) pair-breaking F, singlet', lambda t: 35.5 / t**6.5 * np.exp(-3.528 / t), SF.F_A, 'A'),
    ('(240) pair-breaking F, triplet', lambda t: 1.27 / t**6 * np.exp(-2.376 / t), SF.F_B, 'B'),
]
PI = np.pi
# the same asymptotes in terms of v, as printed (first equality of each equation); single or doubled gap in the exponent
ASYM_V = [
    ('(209)', lambda v: 72 * np.sqrt(2 * PI) / (11513 * PI**8) * v**7.5 * np.exp(-v), SF.mn_pA),
    ('(210)', lambda v: 120960 / (11513 * PI**8) * 0.130 * v**7 * np.exp(-2 * v), SF.mp_pA),
    ('(214)', lambda v: 72 / (11513 * PI**7 * np.sqrt(3)) * v**7 * np.exp(-v), SF.mp_nB),
    ('(218)', lambda v: 945 / (164 * PI**8) * 849 * v * np.exp(-v), SF.np_pA),
    ('(219)', lambda v: 8505 / (41 * PI**6) * v**2 * np.exp(-2 * v), SF.pp_pA),
    ('(240) A', lambda v: np.sqrt(PI) / 2 * v**6.5 * np.exp(-2 * v), SF.F_A),
    ('(240) B', lambda v: PI / (4 * np.sqrt(3)) * v**6 * np.exp(-2 * v), SF.F_B),
]


def limit(f, g):
    r1, r2 = g(100.0) / f(100.0), g(200.0) / f(200.0)
    return float(2 * r2 - r1)


def mn_pA_mutated(v):
    a, b = 0.1477 + SF._s(0.8523, 0.1293, v), 0.1477 + SF._s(0.8523, 0.1297, v)
    return (a**7.5 + b**5.5) / 2 * np.exp(3.4370 - SF._s(3.4370, 1, v))


UNITY = [SF.c_A, SF.c_B, SF.d_A, SF.d_B, SF.mn_pA, SF.mp_pA, SF.mp_nB, SF.np_pA, SF.pp_pA]
FIG19 = {'rho': 2e14, 'Tc': 1e9, 'ms': 0.7, 'peak': 1e21, 'tol_dex': 0.15, 'tau': (0.7, 0.9)}


def asym_dev(swap=False):
    out = []
    for name, f, g, kind in ASYM:
        k = ({'A': 'B', 'B': 'A'}[kind]) if swap else kind
        out.append(float(g(SF.v_of_tau(TAU, k)) / f(TAU) - 1))
    return out


def fig19(n_nu=3):
    tab = Table(*S.eos())
    n = float(np.exp(np.interp(np.log(FIG19['rho']), np.log(tab.rho), np.log(tab.n))))
    xp = float(np.interp(n, *S.composition()[:2]))
    tau = np.linspace(0.3, 0.999, 2000)
    q = SF.pbf(n * (1 - xp), tau * FIG19['Tc'] / 1e9, SF.v_of_tau(tau, 'B'), FIG19['ms'], 'n', n_nu=n_nu)
    qa = SF.pbf(n * (1 - xp), tau * FIG19['Tc'] / 1e9, SF.v_of_tau(tau, 'A'), FIG19['ms'], 'n', n_nu=n_nu) / 4.17   # a_nA = 1
    i = int(np.argmax(q))
    return {'n_n': n * (1 - xp), 'peak': float(q[i]), 'tau_peak': float(tau[i]), 'peak_A': float(qa.max())}


def run():
    dev, dev_swap = asym_dev(), asym_dev(swap=True)
    unity = [float(abs(f(0.0) - 1)) for f in UNITY if f not in (SF.c_A, SF.c_B)]
    v = np.linspace(0, 3, 31)
    red = max(float(np.max(np.abs(SF.d_BA(v, 0 * v) / SF.d_B(v) - 1))), float(np.max(np.abs(SF.d_BA(0 * v, v) / SF.d_A(v) - 1))))
    f19, f19_1 = fig19(), fig19(n_nu=1)
    dex = abs(np.log10(f19['peak'] / FIG19['peak']))
    lim = {name: limit(f, g) for name, f, g in ASYM_V}
    lim_mut = limit(ASYM_V[0][1], mn_pA_mutated)
    tcs = {'p': {k: float(SF.tc(np.linspace(0, 1.5, 600), p, 'A').max()) for k, p in SF.GAPS_P.items()},
           'nt': {k: float(SF.tc(np.linspace(1, 3.5, 1000), p, 'B').max()) for k, p in SF.GAPS_NT.items()}}
    checks = [
        {'name': 'Reduction factors equal 1 at v = 0', 'ok': bool(max(unity) < 2e-3), 'note': f"max deviation {max(unity):.1e}"},
        {'name': f"Fits reproduce their printed asymptotes at T = {TAU} T_c within 5 %", 'ok': bool(max(map(abs, dev)) <= 0.05), 'note': ', '.join(f"{x:+.1%}" for x in dev)},
        {'name': 'H6b: fit over printed asymptote, extrapolated to v = infinity, within 3 % of 1 for all seven', 'ok': bool(max(abs(x - 1) for x in lim.values()) <= 0.03), 'note': ', '.join(f"{x:.3f}" for x in lim.values())},
        {'name': 'H6b NEG: a 10 % change in one coefficient of Eq. (211) takes the limit outside 3 %', 'ok': bool(abs(lim_mut - 1) > 0.03), 'note': f"{lim_mut:.3f}"},
        {'name': 'Fit (205) reduces to the single-superfluid factors (199) for v <= 3 within 5 %', 'ok': bool(red <= 0.05), 'note': f"max {red:.1%}"},
        {'name': f"Fig. 19: triplet pair-breaking peak {f19['peak']:.2e} erg cm-3 s-1 vs 1e21 (within 0.15 dex)", 'ok': bool(dex <= FIG19['tol_dex']), 'note': f"{np.log10(f19['peak'] / FIG19['peak']):+.2f} dex"},
        {'name': f"Fig. 19: peak at T = {f19['tau_peak']:.2f} T_c (0.7 to 0.9)", 'ok': bool(FIG19['tau'][0] <= f19['tau_peak'] <= FIG19['tau'][1]), 'note': ''},
        {'name': 'NEG1: swapping the singlet and triplet gap functions breaks an asymptote', 'ok': bool(max(map(abs, dev_swap)) > 0.05), 'note': f"max {max(map(abs, dev_swap)):.0%}"},
        {'name': 'NEG2: one neutrino flavour takes the Fig. 19 peak outside 0.15 dex', 'ok': bool(abs(np.log10(f19_1['peak'] / FIG19['peak'])) > FIG19['tol_dex']), 'note': f"{np.log10(f19_1['peak'] / FIG19['peak']):+.2f} dex"},
    ]
    out = {'source': 'Yakovlev et al. 2001; Yakovlev, Levenfish & Shibanov 1999; Ho et al. 2015', 'tau': TAU,
           'asymptotes': dict(zip([a[0] for a in ASYM], dev)), 'reduction_205': red, 'limits': lim, 'limit_mutated': lim_mut, 'fig19': f19, 'tc_max': tcs, 'checks': checks}
    RES.mkdir(exist_ok=True)
    (RES / 'step6_superfluid.json').write_text(json.dumps(out, indent=1))
    for c in checks:
        print('PASS' if c['ok'] else 'FAIL', c['name'], f"({c['note']})" if c['note'] else '')
    print('  max Tc [1e9 K]  p:', {k: round(x / 1e9, 2) for k, x in tcs['p'].items()})
    print('  max Tc [1e9 K] nt:', {k: round(x / 1e9, 2) for k, x in tcs['nt'].items()})
    return out


if __name__ == '__main__':
    run()
