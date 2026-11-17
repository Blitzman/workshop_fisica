"""N2: validation of the extended Lindblom-Detweiler system (axion switched off) against the literature.

Port of axion_ns/step8_numerics.py: compute() and checks() with the same algorithms, jobs and check texts, minus the axion background
(solve_bvp) and the figures, which are not part of this module (see MIGRATION.md). The numerical systems come from the equation
artifacts (systems.py), the published numbers from R's literature artifact (refs.py).

    * Cowling (8a)-(8c) of Paper II against the 24 Cowling modes of Yoshida and Kojima (1997);
    * the extended system without axion, with the LD centre conditions and Zerilli outside, against the "exact" frequencies of
      Lindblom, Mendell and Ipser (1997) (5 digits) and Yoshida and Kojima (3 digits);
    * Zhao and Lattimer's X' as printed (sign erratum): the f-mode must come out more than 10% off.

Conditioning (measured, MIGRATION.md): |A_out|^2 carries 1e-3..1e-2 relative noise because both interior solutions grow to W ~ 1e6 at the
surface and the X(R) = 0 combination cancels that by a factor ~ 3e5. The defaults rtol = 1e-11, atol = 1e-14 must not be loosened for the
amplitude (rtol 1e-10 already moves the LMI f-mode by 7.7e-5 > 5e-5); the Cowling rtol can be relaxed (option) after checking the roots.
"""
import multiprocessing as mp

import numpy as np
from scipy.optimize import minimize_scalar

from . import modes as M
from .polytrope import Star
from . import systems as S2

G, C, MSUN = 6.6743e-8, 2.99792458e10, 1.98841e33
LUNIT = np.sqrt(1.3346e5 / G)                        # K^{1/2} in cm
MASS_LMI = 1.4 * G * MSUN / C**2 / LUNIT
MODES = ('f', 'p1', 'p2')
LS = (2, 3, 4, 5)
ZL_GUESS = 0.833                                      # the minimum with the Zhao-Lattimer erratum (grid 0.3-2)
_near = lambda w0: w0 * (1 + np.concatenate([-np.geomspace(1e-4, 3e-2, 12), [0], np.geomspace(1e-4, 3e-2, 12)]))
SCAN = np.unique(np.concatenate([np.linspace(0.6, 2.0, 57), _near(1.2016), _near(ZL_GUESS)]))
SCAN_QUICK = np.array([1.1, 1.2016, 1.3])

_W = {}                                               # systems, stars, literature, settings: built BEFORE the fork


def setup(ext, p2, lit, cow_rtol=1e-11, quick=False):
    if not _W:
        _W['cow'] = S2.cowling_system(p2, 'energy')
        _W['full'], _W['ZL'] = S2.full_systems(ext, 'energy')
        _W['stars'] = {'LMI': Star('energy', mass=MASS_LMI)}
        if not quick:
            _W['stars'].update({0.1: Star('energy', 0.1), 0.2: Star('energy', 0.2)})
        _W['lit'], _W['cow_rtol'], _W['quick'] = lit, cow_rtol, quick
    return _W


def _unit(st):
    return np.sqrt(st.M / st.R**3)


def _cowling(key):
    st = _W['stars'][key]
    u = _unit(st)
    out = {}
    for ll in ((2,) if _W['quick'] else LS):
        grid = (np.linspace(0.8, 3.0, 12) if _W['quick'] else np.linspace(0.8, 7.0, 64)) * u
        roots = M.cowling_modes(st, _W['cow'], ll, grid, rtol=_W['cow_rtol'])
        out[ll] = [x / u for x in roots[:3]]
    return key, out


def _full(job):
    key, ll, guess, which = job
    s0, g = M.full_mode(_W['stars'][key], _W[which], ll, guess, width=0.06)
    return job, (s0, g)


def _converged(_):
    """The LMI f-mode l = 2 with another starting point, another surface cut and a farther extraction (as the legacy)."""
    st = _W['stars']['LMI']
    orig = st.integrate
    st.integrate = lambda system, y0, o, ll, r0=None, r1=None, rtol=1e-11, atol=1e-14: orig(
        system, y0, o, ll, r0=1e-5 * st.R, r1=st.R * (1 - 1e-9), rtol=rtol, atol=atol)
    u = _unit(st)
    m = minimize_scalar(lambda s: M.outgoing_amplitude(st, _W['full'], s * u, 2, n_far=120), bracket=(1.18, 1.2016, 1.22), tol=1e-10)
    st.integrate = orig
    return m.x


def _scan(job):
    which, s = job
    st = _W['stars']['LMI']
    return job, M.outgoing_amplitude(st, _W[which], s * _unit(st), 2)


def _map(pool, fn, items):
    return pool.map(fn, items) if pool is not None else [fn(x) for x in items]


def compute(procs=1, quick=False):
    """Returns (cow, full, conv, scan, scan_grid). The pool uses the 'fork' start method (the compiled systems are not picklable)."""
    lit = _W['lit']
    LMI97, YK97 = lit[0], lit[1]
    if quick:
        jobs, scan_grid, cow_keys = [('LMI', 2, LMI97[2][0], 'full')], SCAN_QUICK, ['LMI']
    else:
        jobs = [('LMI', ll, LMI97[ll][0], 'full') for ll in LS]
        jobs += [(MR, ll, YK97[MR][ll][m][0], 'full') for MR in (0.1, 0.2) for ll in LS for m in MODES]
        jobs += [('LMI', 2, ZL_GUESS, 'ZL')]
        scan_grid, cow_keys = SCAN, [0.1, 0.2, 'LMI']
    pool = None
    if procs > 1 and 'fork' in mp.get_all_start_methods():
        pool = mp.get_context('fork').Pool(procs)
    try:
        cow = dict(_map(pool, _cowling, cow_keys))
        full = dict(_map(pool, _full, jobs))
        conv = None if quick else _map(pool, _converged, [0])[0]
        which = ('full',) if quick else ('full', 'ZL')
        sc = dict(_map(pool, _scan, [(w, s) for w in which for s in scan_grid]))
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    scan = {w: np.array([sc[(w, s)] for s in scan_grid]) for w in which}
    return cow, full, conv, scan, scan_grid


def checks(cow, full, conv, lit):
    """The 7 checks of axion_ns/step8_numerics.checks (same texts), for the full run."""
    LMI97, YK97, LMI_R_KM = lit
    res = []
    st = _W['stars']['LMI']
    Rkm = st.R * LUNIT / 1e5
    res.append((f'Fondo: el modelo de 1.4 M_sun con kappa = 1.3346e5 cgs tiene R = {Rkm:.3f} km (LMI: 14.45 km)',
                abs(Rkm - LMI_R_KM) < 0.005, ''))
    pairs = [(cow[MR][ll][i], YK97[MR][ll][m][1], f'{m} l={ll} M/R={MR}') for MR in (0.1, 0.2) for ll in LS for i, m in enumerate(MODES)]
    out = [f'{n}: {a:.4f} vs {b}' for a, b, n in pairs if abs(a - b) > 0.005]
    dmax = max(abs(a - b) for a, b, _ in pairs)
    res.append((f'Cowling, Ecs. (8a)-(8c) del paso 5: los 24 modos f, p1, p2 (l = 2-5, M/R = 0.1 y 0.2) de Yoshida '
                f'y Kojima (1997): {24 - len(out)} dentro de media unidad de su última cifra, máx. dif. {dmax:.4f}',
                dmax < 0.006, '; '.join(out)))
    dre = [abs(full[('LMI', ll, LMI97[ll][0], 'full')][0] - LMI97[ll][0]) for ll in LS]
    res.append((f'Métrica dinámica, sistema extendido sin axión: Re w de los modos f l = 2-5 igual a las 5 cifras de '
                f'Lindblom, Mendell e Ipser (1997) (máx. dif. {max(dre):.5f})', max(dre) <= 5e-5, ''))
    dim = [abs(full[('LMI', ll, LMI97[ll][0], 'full')][1] / LMI97[ll][1] - 1) for ll in LS]
    res.append((f'Im w (emisión de ondas gravitacionales) dentro del 2 % de LMI (máx. {100 * max(dim):.1f} %)', max(dim) < 0.02, ''))
    dyk = [abs(full[(MR, ll, YK97[MR][ll][m][0], 'full')][0] / YK97[MR][ll][m][0] - 1) for MR in (0.1, 0.2) for ll in LS for m in MODES]
    res.append((f'Métrica dinámica: los 24 modos "exactos" de Yoshida y Kojima dentro del {100 * max(dyk):.1f} % '
                '(a 3 cifras; LMI, a 5, coincide entero)', max(dyk) < 0.007, ''))
    s0 = full[('LMI', 2, LMI97[2][0], 'full')][0]
    res.append((f'Convergencia: otro punto de partida, corte en la superficie y extracción mueven Re w en '
                f'{abs(conv - s0):.0e}', abs(conv - s0) < 1e-5, ''))
    zl = full[('LMI', 2, ZL_GUESS, 'ZL')][0]
    res.append((f'Errata de signo de Zhao y Lattimer: con su X\' tal como está impreso, el modo f l = 2 sale en '
                f'{zl:.4f} en vez de 1.2016 ({100 * (zl / LMI97[2][0] - 1):+.0f} %)', abs(zl / LMI97[2][0] - 1) > 0.1, ''))
    return res


def quick_checks(cow, full, lit):
    """Smoke-test checks: the same quantities restricted to what the quick run computes (LMI star, l = 2)."""
    LMI97, _, LMI_R_KM = lit
    st = _W['stars']['LMI']
    Rkm = st.R * LUNIT / 1e5
    s0, g0 = full[('LMI', 2, LMI97[2][0], 'full')]
    return [(f'quick: R = {Rkm:.3f} km (LMI: 14.45 km)', abs(Rkm - LMI_R_KM) < 0.005, ''),
            (f'quick: Re w of the LMI f-mode l = 2 = {s0:.6f} (LMI 1.2016, tolerance 5e-5)', abs(s0 - LMI97[2][0]) <= 5e-5, ''),
            (f'quick: Im w = {g0:.4e} within 2 % of LMI ({LMI97[2][1]})', abs(g0 / LMI97[2][1] - 1) < 0.02, ''),
            (f'quick: Cowling roots of the LMI star for l = 2 in 0.8-3.0 u: {[round(x, 4) for x in cow["LMI"][2]]}', None, 'descriptive')]


def mode_table(cow, full, lit):
    """Long-format table of every mode computed."""
    LMI97, YK97, _ = lit
    rows = []
    for key, byl in cow.items():
        for ll, roots in byl.items():
            for i, x in enumerate(roots):
                ref = YK97[key][ll][MODES[i]][1] if key in YK97 else float('nan')
                rows.append((str(key), ll, MODES[i], 'Cowling (8a-8c)', x, float('nan'), ref, float('nan')))
    for (key, ll, guess, which), (s0, g) in full.items():
        if key == 'LMI':
            mode, ref, refi = 'f', LMI97[ll][0] if which == 'full' else float('nan'), LMI97[ll][1] if which == 'full' else float('nan')
        else:
            mode = next(m for m, v in YK97[key][ll].items() if v[0] == guess)
            ref, refi = guess, float('nan')
        rows.append((str(key), ll, mode, 'extended LD' if which == 'full' else 'ZL as printed', s0, g, ref, refi))
    cols = list(zip(*rows))
    names = ['star', 'l', 'mode', 'system', 're_w', 'im_w', 'ref_re_w', 'ref_im_w']
    return {n: (np.array(c, dtype=str) if n in ('star', 'mode', 'system') else np.array(c, dtype=float)) for n, c in zip(names, cols)}
