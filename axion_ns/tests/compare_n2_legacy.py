"""SLOW, manual:  python -m tests.compare_n2_legacy      (about 3-4 minutes on one core)

Compares the modular N2 (equations from the artifacts, lambdify with cse) with the ORIGINAL axion_ns/fmodes.py on the LMI star (l = 2):
  1. the extended system from D2's artifact vs fullgr.EXT; the star; the compiled matrices at 40 real states;
  2. ONE Cowling evaluation and ONE amplitude evaluation, legacy vs new (value and time);
  3. an rtol scan of the new functions (values relative to rtol = 1e-11, time, right-hand-side calls).
Imports the legacy package on purpose (cross-module test).
"""
import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parent.parent
sys.path[:0] = [str(ROOT), str(ROOT.parent / 'archive')]

from axion_ns import fmodes as FM                      # noqa: E402
from axion_ns import fullgr as FG                      # noqa: E402
from N2_modes_fullgr_noaxion import modes as M2        # noqa: E402
from N2_modes_fullgr_noaxion import systems as S2      # noqa: E402
from N2_modes_fullgr_noaxion.polytrope import Star     # noqa: E402

G, C, MSUN = 6.6743e-8, 2.99792458e10, 1.98841e33
MASS = 1.4 * G * MSUN / C**2 / np.sqrt(1.3346e5 / G)
T = time.perf_counter


def main():
    print('--- setup (symbolic solve + lambdify), one core')
    ext, p2, _ = S2.load_inputs(ROOT / 'D2_derive_nonradial' / 'outputs', ROOT / 'R_paper_refs' / 'outputs')
    t0 = T(); full = S2.full_system(ext, 'energy'); t_full = T() - t0
    t0 = T(); cow = S2.cowling_system(p2, 'energy'); t_cow = T() - t0
    print(f'new   : full_system {t_full:.1f} s, cowling_system {t_cow:.1f} s  (cse=True)', flush=True)
    t0 = T(); Lcow = FM.cowling_system('energy'); Lfull = FM.full_system('energy'); t_leg = T() - t0
    print(f'legacy: both systems {t_leg:.1f} s', flush=True)

    print('\n--- 1. equations, star, matrices')
    same = {k: (ext[k] == FG.EXT[k]) or sp.expand(ext[k] - FG.EXT[k]) == 0 for k in S2.EXT_KEYS}
    print('extended system from D2 artifact == legacy fullgr.EXT:', same, flush=True)
    st, sl = Star('energy', mass=MASS), FM.Star('energy', mass=MASS)
    print(f'star: new M={st.M:.12f} R={st.R:.12f} Phic={st.Phic:.12f} | legacy M={sl.M:.12f} R={sl.R:.12f} Phic={sl.Phic:.12f}', flush=True)
    u = np.sqrt(st.M / st.R**3)
    r0 = 1e-4 * st.R; b0 = st._start(st.nc, r0); b0[2] = st.Phic
    tov = solve_ivp(lambda rr, y: st.bg(rr, *y)[1], (r0, st.R * (1 - 1e-7)), b0, dense_output=True, rtol=1e-12, atol=1e-30)
    rs = np.linspace(0.03, 0.97, 40) * st.R
    worst = {'full mat': 0, 'full ay': 0, 'cow mat': 0, 'cow ay': 0}
    y4, y2 = np.array([1.0, 0.5, 0.3, 0.1]), np.array([1.0, 0.4])
    for x in rs:
        B = st.bg(x, *tov.sol(x))[0]
        for tag, new, leg, y in (('full', full, Lfull, y4), ('cow', cow, Lcow, y2)):
            a = np.array(leg(x, 1.2016 * u, 2, *B), float)
            m = np.asarray(new.mat(x, 1.2016 * u, 2, *B), float)
            ay = np.array(new.ay(x, 1.2016 * u, 2, *B, *y), float)
            worst[tag + ' mat'] = max(worst[tag + ' mat'], np.abs(m - a).max() / np.abs(a).max())
            worst[tag + ' ay'] = max(worst[tag + ' ay'], np.abs(ay - a @ y).max() / np.abs(a @ y).max())
    print('max relative difference of the compiled matrices vs the legacy ones (40 states):', {k: f'{v:.1e}' for k, v in worst.items()}, flush=True)

    print('\n--- 2. ONE evaluation each, legacy vs new (LMI star, l = 2)')
    rows = []
    for label, fl, fn, om in (('Cowling X(R)  omega = 1.3 u', lambda o: FM.cowling_X(sl, Lcow, o, 2), lambda o: M2.cowling_X(st, cow, o, 2), 1.3 * u),
                              ('amplitude     omega = 1.2016 u', lambda o: FM.outgoing_amplitude(sl, Lfull, o, 2), lambda o: M2.outgoing_amplitude(st, full, o, 2), 1.2016 * u)):
        t0 = T(); a = fl(om); ta = T() - t0
        st.nfev = 0
        t0 = T(); b = fn(om); tb = T() - t0
        print(f'{label}: legacy {a:.12e} ({ta:.2f} s)   new {b:.12e} ({tb:.2f} s, {st.nfev} rhs calls)   rel diff {abs(a - b) / abs(a):.1e}   speedup {ta / tb:.1f}x', flush=True)

    print('\n--- 3. rtol scan of the NEW functions (value relative to rtol = 1e-11; time; right-hand-side calls)')
    print(f'{"function @ omega":34s} {"rtol":>7s} {"value":>17s} {"rel diff":>9s} {"time (s)":>9s} {"rhs calls":>10s}')
    for label, fn, om in (('amplitude @ 1.2016 u (resonance)', lambda o, t: M2.outgoing_amplitude(st, full, o, 2, rtol=t), 1.2016 * u),
                          ('amplitude @ 1.45 u (off resonance)', lambda o, t: M2.outgoing_amplitude(st, full, o, 2, rtol=t), 1.45 * u),
                          ('Cowling X(R) @ 1.3 u', lambda o, t: M2.cowling_X(st, cow, o, 2, rtol=t), 1.3 * u)):
        base = None
        for rt in (1e-11, 1e-10, 1e-9, 1e-8, 1e-7):
            st.nfev = 0
            t0 = T(); v = fn(om, rt); tt = T() - t0
            base = v if base is None else base
            print(f'{label:34s} {rt:7.0e} {v:17.9e} {abs(v / base - 1):9.1e} {tt:9.2f} {st.nfev:10d}', flush=True)


if __name__ == '__main__':
    main()
