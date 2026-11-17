"""N1, part `tables` (old step 9): Cowling modes f, p1, p2 of the BSk26 stars with condensate, Tables I, II and IV of Paper II.

Logic of axion_ns/step9_axion_cowling.py. What changed: the star is a StarTable rebuilt from S1's CSV tables (not an
in-memory Star), the Paper II tables come from R's data file, only the stars present in S1's outputs are computed
(the check texts report how many), and the result is also written as a full-precision table (modes_cowling.csv) instead
of only percent shifts rounded to 4 digits.
"""
import time

import numpy as np
from scipy.optimize import brentq

from . import modes as M
from .star_table import StarTable
from .tables import step9_shapes

REF = {'M': 1.53, 'R': 11.7, 'modes': (2.43, 7.15, 10.18)}
BRACKETS = [(2.2, 3.1), (6.9, 8.6), (9.9, 12.6)]
L = 2
NAMES = ('f', 'p1', 'p2')


def _real_root(fn, lo, hi, n=24):
    xs = np.linspace(lo, hi, n)
    v = [fn(x) for x in xs]
    for i in range(n - 1):
        if np.sign(v[i]) != np.sign(v[i + 1]):
            return brentq(fn, xs[i], xs[i + 1], xtol=1e-7)
    return np.nan


def background_modes(st):
    return [_real_root(lambda x: M.surface_residual(st, L, M.omega_of(x)), *b) for b in BRACKETS]


def coupled(st, guess):
    """Coupled mode near `guess` (kHz). Below the threshold nu_a the root is real (bound state)."""
    det = lambda z: M.coupled_det(st, L, M.omega_of(z))
    if guess < st.nu_a:
        z = _real_root(lambda x: det(x).real, guess - 0.35, min(guess + 0.15, st.nu_a - 1e-4), n=16)
        if np.isfinite(z):
            return complex(z)
    # H29: tol = 1e-9 |z| is < 1e-5 Hz (1e-11 exhausted the 40 steps of Muller)
    return M.muller(det, guess * 0.997, guess, guess * 1.003 - 2e-4j, tol=1e-9)


def row(st, nu, t1_keys, t2_keys):
    t0 = time.time()
    out = {'nu': nu, 'M': st.M / st.MSUN, 'R': st.R, 'thc': st.thc}
    if nu in t1_keys or nu in t2_keys:
        out['T1'] = background_modes(st)
    if nu in t2_keys:
        out['T2'] = []
        for g0 in out['T1']:
            try:
                out['T2'].append(coupled(st, g0))
            except Exception:  # noqa: BLE001
                out['T2'].append(complex(np.nan, np.nan))
    out['secs'] = time.time() - t0
    return out


def compute(star_for, nus, t1_keys, t2_keys, ref_star):
    ref = {'M': ref_star.M / ref_star.MSUN, 'R': ref_star.R, 'modes': background_modes(ref_star)}
    rows = {}
    for nu in nus:
        rows[nu] = row(star_for(nu), nu, t1_keys, t2_keys)
        print(f'  nu_a={nu:g}: M={rows[nu]["M"]:.4f} modes_bg={np.round(rows[nu].get("T1", []), 4)} '
              f'coupled={np.round(rows[nu].get("T2", []), 4)} ({rows[nu]["secs"]:.0f} s)', flush=True)
    return ref, rows


def _shift(v, v0):
    return 100 * (np.asarray(v) / np.asarray(v0) - 1)


def checks(ref, rows, none25, T1p, T2p, T4M):
    """Same checks as old step 9, restricted to the stars computed in this run."""
    res = []
    r0, p0 = np.array(ref['modes']), np.array(REF['modes'])
    dev = 100 * np.abs(r0 / p0 - 1)
    res.append((f"Cowling sin axión: f, p1, p2 = {r0[0]:.3f}, {r0[1]:.3f}, {r0[2]:.3f} kHz (Paper II: 2.43, 7.15, 10.18); "
                f"máx. dif. {dev.max():.2f} %", dev.max() <= 0.3,
                'ajuste analítico de BSk26 frente a la tabla del grupo; el paper da un 0.3 % de acuerdo con otros cálculos'))
    dm = {nu: 100 * (rows[nu]['M'] / m - 1) for nu, m in T4M.items() if nu in rows}
    if dm:
        res.append((f"Tabla IV: masas de las estrellas con axión ({len(dm)} valores) dentro del 0.5 % "
                    f"(máx. {max(abs(v) for v in dm.values()):.2f} %)", max(abs(v) for v in dm.values()) < 0.5, ''))
    else:
        res.append(('Tabla IV: masas de las estrellas con axión', None, 'ninguna estrella de la Tabla IV en esta ejecución'))
    d1 = [_shift(rows[nu]['T1'], r0) - _shift(pv, p0) for nu, pv in T1p.items() if nu in rows and 'T1' in rows[nu]]
    if d1:
        m1 = np.abs(d1).max()
        res.append((f"Tabla I (solo fondo): los {3 * len(d1)} desplazamientos nu/nu_0 - 1 coinciden con el paper dentro de 0.6 "
                    f"puntos (máx. {m1:.2f})", m1 <= 0.6, 'nu_0: la nota sin axión de cada cálculo'))
    t2rows = {nu: pv for nu, pv in T2p.items() if nu in rows and 'T2' in rows[nu]}
    if t2rows:
        d2 = [_shift([z.real for z in rows[nu]['T2']], r0) - _shift([v[0] for v in pv], p0) for nu, pv in t2rows.items()]
        m2 = np.abs(d2).max()
        df = max(abs(d[0]) for d in d2)
        res.append((f"Tabla II (acoplado), parte real: los {3 * len(d2)} desplazamientos dentro de 0.6 puntos "
                    f"(máx. {m2:.2f}; modo f, {df:.2f})", m2 <= 0.6, ''))
        mism, near = [], []
        for nu, pv in t2rows.items():
            for i, z in enumerate(rows[nu]['T2']):
                ours, theirs = abs(z.imag) * 1e3 > 1e-3, pv[i][1] != 0
                if ours != theirs:
                    (near if abs(z.real / nu - 1) < 0.003 else mism).append((nu, NAMES[i]))
        res.append(("Tabla II: los modos amortiguados son los mismos (Im nu = 0 por debajo de nu_a)"
                    + (f" salvo {len(near)} justo en el umbral" if near else ", sin excepciones"), not mism,
                    ', '.join(f'{m} a {nu} kHz: Re nu a menos del 0.3 % de nu_a' for nu, m in near)))
        ratios = [abs(z.imag) * 1e3 / abs(pv[i][1]) for nu, pv in t2rows.items() for i, z in enumerate(rows[nu]['T2'])
                  if abs(pv[i][1]) >= 1]
        im_max = max(abs(z.imag) for nu in t2rows for z in rows[nu]['T2'])
        tau = 1 / (2 * np.pi * im_max * 1e3) * 1e3
        if ratios:
            res.append((f"Tabla II, amortiguamiento: los {len(ratios)} |Im nu| >= 1 Hz del paper, dentro de un factor 2 "
                        f"(cociente {min(ratios):.2f}-{max(ratios):.2f}); tau_a mínimo = {tau:.1f} ms (paper: 12 ms)",
                        0.5 <= min(ratios) and max(ratios) <= 2 and (len(t2rows) < 12 or abs(tau - 12.2) < 1.5),
                        'el amortiguamiento depende de delta theta(R), muy sensible al fondo (Fig. 6 del paper); '
                        'tau_a mínimo solo es comparable con las 12 estrellas de la Tabla II'))
    if none25 is None:
        res.append(('Umbral del condensado (de S1)', None, 'sin dato'))
    return res


def table_columns(ref, rows):
    nus = sorted(rows)
    cols = {k: [] for k in ('nu_a_kHz', 'mode', 'nu_noaxion_kHz', 'nu_bg_kHz', 'nu_cp_re_kHz', 'nu_cp_im_kHz',
                            'M_msun', 'R_km')}
    for nu in nus:
        for i, nm in enumerate(NAMES):
            r = rows[nu]
            cols['nu_a_kHz'].append(nu)
            cols['mode'].append(nm)
            cols['nu_noaxion_kHz'].append(ref['modes'][i])
            cols['nu_bg_kHz'].append(r['T1'][i] if 'T1' in r else np.nan)
            z = r['T2'][i] if 'T2' in r else complex(np.nan, np.nan)
            cols['nu_cp_re_kHz'].append(z.real)
            cols['nu_cp_im_kHz'].append(z.imag)
            cols['M_msun'].append(r['M'])
            cols['R_km'].append(r['R'])
    return {k: np.array(v) for k, v in cols.items()}
