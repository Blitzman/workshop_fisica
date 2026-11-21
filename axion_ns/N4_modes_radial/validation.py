"""N4: validation of the radial modes (old steps 15, 16, 19, 20, 22 and 26; step 17, the Font et al. comparison, is not part of this module by decision).

Port of axion_ns/step15_radial_star.py, step16_radial_stability.py, step19_avoided_crossing.py, step20_damped_count.py, step22_exterior_contour.py and
step26_paper1_tables.py with the same algorithms, starting points and (adapted where the star source changed) check texts. What changed:
  * the stars come from S1's tables (the 12 of Paper II and, for the stars the legacy built on the fly, S1's extra set: nu_a = 3.4, 3.5 and the step 19 values, the
    p_c = 1000 star with condensate, and no-axion stars for p_c = 700..1100 in steps of 10), and the right-hand sides from the coefficient tables (tables.py);
  * step 15, turning point: the legacy located it with a bounded minimiser and a root finder on freshly integrated stars; here M(p_c) and omega_F^2(p_c) are taken on the
    fixed grid of the extra stars and interpolated with cubic splines (same criterion: both at the same p_c within 1%);
  * step 19: the legacy refined the nu_a grid around the minimum of the gap at run time; here the 15 values it ended up with are fixed (all in the extra set);
  * step 15, retired artefact: the window of the 'axion resonance' is 'above the threshold and damped' (the legacy window was tied to its matching radius, and the resonance moves with it);
  * step 15, robustness of Im nu: the start radius variant r0 = 1e-4 km is r0 = 3e-3 km (the tables start at 1e-3 km), and the matching-radius variants are 0.8 and 0.6 of
    the matching radius (the legacy used the first |theta_0| < 1e-5 and 1e-6, which with a tail valid to r_end may not exist);
  * the matching radius is `R + 0.7 (r_end - R)` instead of the first |theta| < 1e-7;
  * the pools are optional (`procs`).
Sampling of the contour (steps 16, 19, 20) is unchanged: every point is one det_entire evaluation (about 0.11 s instead of about 2.7 s).
"""
import copy
import json
import multiprocessing as mp
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

from common.checks import ok_field

from . import modes as MD
from common.star_table import StarTable
from .tables import RadialTables

omega_of = MD.omega_of
MSUN_KM = 1.476625
_W = {}                                  # flat functions, paths: built BEFORE any fork
_CACHE = {}                              # (name) -> (star, tables), per process


def setup(flat_i, flat_e, eos_csv, s1_dir, extra_dir, procs=1):
    _W.update(fi=flat_i, fe=flat_e, eos=eos_csv, s1=Path(s1_dir), extra=Path(extra_dir), procs=procs)
    _CACHE.clear()


def _pool():
    return mp.get_context('fork').Pool(_W['procs']) if (_W['procs'] > 1 and 'fork' in mp.get_all_start_methods()) else None


def star_file(nu=None, p0=None):
    """Path of the table of the star with axion mass nu (kHz; None: no axion) and central pressure p0 (None: 100 MeV/fm^3)."""
    name = ('star_noaxion' if p0 is None else f'star_noaxion_p0_{p0:g}') if nu is None else (f'star_nua{nu:g}' if p0 is None else f'star_nua{nu:g}_p0_{p0:g}')
    for d in (_W['s1'], _W['extra']):
        if (d / f'{name}.csv').exists():
            return d / f'{name}.csv'
    raise FileNotFoundError(f'{name}.csv not found in {_W["s1"]} or {_W["extra"]} (extra stars: python -m S1_star_background.run --extra)')


def get_star(nu=None, p0=None):
    key = (nu, p0)
    if key not in _CACHE:
        st = StarTable(_W['eos'], star_file(nu, p0))
        _CACHE[key] = (st, RadialTables(st, _W['fi'], _W['fe']))
    return _CACHE[key]


def old_lapse(st, tabs=None):
    """Same star with the OLD lapse normalisation (Phi matched at R with m(R), without the tail of the condensate): the positive control of steps 15 and 20."""
    m_R, _, Phi_R = st.y_at(st.R)[:3]
    st2 = copy.copy(st)
    st2.Phi_shift = 0.5 * np.log(1 - 2 * m_R / st.R) - Phi_R
    return st2, RadialTables(st2, _W['fi'], _W['fe'])


def r_theta_below(st, rm, tol):
    """First radius in [R, r_m] with |theta_0| < tol (r_m if none)."""
    rr = np.linspace(st.R, rm, 4000)
    th = np.abs(np.array([st.y_at(x)[3] for x in rr]))
    return float(rr[np.argmax(th < tol)]) if (th < tol).any() else float(rm)


def dump(path, obj):
    Path(path).write_text(json.dumps(obj, indent=1, ensure_ascii=False, default=lambda o: o.item() if hasattr(o, 'item') else str(o)))


def _corners(b):
    x0, x1, y0, y1 = b
    return [complex(x0, y0), complex(x1, y0), complex(x1, y1), complex(x0, y1)]


def _map(pool, fn, items):
    return pool.map(fn, items) if pool is not None else [fn(x) for x in items]


# =============================================================================================== step 22
POINTS22 = {'H1': 8.1192 - 0.0036j, 'artefacto antiguo': 6.585 - 0.563j, 'cero profundo': 6.293 - 4.480j, 'nivel ligado n = 2': 5.944070 + 0j}


def step22():
    res, out = [], {'points': {}}
    st, tabs = get_star(6.0)
    rm, M = float(MD.matching_radius(st)), float(st.Mt)

    def L(nu, **kw):
        kw = {'sheet': 1, 'n_decay': 40.0, **kw}
        ps, dps = MD.exterior_psi(tabs.flat_e, omega_of(nu), M, rm, st.ma, np.sqrt(st.fa2), st.E0, st.beta, rtol=1e-12, **kw)
        return complex(dps / ps)

    for name, nu in POINTS22.items():
        L0 = L(nu)
        rel = lambda z: float(abs(z / L0 - 1))
        try:
            ing = rel(L(nu, sheet=-1))
        except RuntimeError:
            ing = None            # below the threshold the 'ingoing wave' is the growing solution e^{kappa r}: the integration does not finish
        out['points'][name] = {'nu': [nu.real, nu.imag], 'L': [L0.real, L0.imag], 'dphi+0.4': rel(L(nu, dphi=0.4)), 'dphi-0.4': rel(L(nu, dphi=-0.4)),
                               'n_decay60': rel(L(nu, n_decay=60.0)), 'ingoing': ing, 'n_decay1': rel(L(nu, n_decay=1.0))}
    P = out['points']
    ang = max(max(p['dphi+0.4'], p['dphi-0.4']) for p in P.values())
    nd = max(p['n_decay60'] for p in P.values())
    res.append((f'H24(a): girar el rayo complejo ±0.4 rad cambia L = psi\'/psi en r_m como mucho {ang:.1e} (relativo) en ' + ', '.join(P) + ' '
                + '; '.join(f"{k}: {max(p['dphi+0.4'], p['dphi-0.4']):.0e}" for k, p in P.items()), ang < 1e-8, 'la condición saliente no depende del camino en el plano r complejo'))
    res.append((f'H24(b): alejar el punto de partida (n_decay 40 → 60) cambia L como mucho {nd:.1e}: ' + '; '.join(f"{k}: {p['n_decay60']:.0e}" for k, p in P.items()),
                nd < 1e-8, 'el truncado del desarrollo asintótico (con la fase coulombiana) no influye'))
    ing = min(p['ingoing'] for k, p in P.items() if k != 'nivel ligado n = 2')
    bound = P['nivel ligado n = 2']['ingoing']
    res.append((f'Control negativo 1: con la onda entrante, L cambia al menos un {ing:.2g} (relativo) en los puntos por encima del umbral', ing > 1e-2,
                'en el nivel ligado la mutación es la solución creciente: ' + ('la integración no termina y se detecta' if bound is None else f'L cambia {bound:.1e}')))
    n1 = min(p['n_decay1'] for p in P.values())
    res.append((f'Control negativo 2: con n_decay = 1 el cambio mínimo es {n1:.1e}, frente a {max(ang, nd):.1e} en (a)-(b): la prueba detecta un punto de partida demasiado cercano',
                n1 > 100 * max(ang, nd), ''))
    out['checks'] = [{'name': n_, 'ok': ok_field(o_), 'note': t_} for n_, o_, t_ in res]
    return res, out


# =============================================================================================== step 26
NUA26 = (2.8, 3.0, 4.0, 5.0, 6.0, 8.0, 12.0)


def decoupled_dp(nu_a, nu, r0=1e-3, r1=0.5, rtol=1e-10):
    """Delta p(R) with delta theta switched off: the (xi, Delta p) submatrix of the interior matrix (Table 1)."""
    st, tabs = get_star(nu_a)
    om = omega_of(nu)
    M = tabs.mat_int(r0, om, keep=(0, 1)).real
    y0 = np.array([r0, (1.0 - M[0, 0] * r0) / M[0, 1]])
    f = tabs.rhs_int(om, keep=(0, 1))
    s = MD._ok(solve_ivp(f, (r0, r1), y0, method='DOP853', rtol=rtol, atol=1e-40))
    s = MD._ok(solve_ivp(f, (r1, st.R), s.y[:, -1], method='DOP853', rtol=rtol, atol=1e-40))
    return s.y[1, -1].real


def bracket_root(g, seed, rel=0.04, n=9):
    xs = np.linspace(seed * (1 - rel), seed * (1 + rel), n)
    v = [g(x) for x in xs]
    cands = [(xs[i], xs[i + 1]) for i in range(n - 1) if v[i] * v[i + 1] < 0]
    if not cands:
        return None
    a, b = min(cands, key=lambda ab: abs(0.5 * (ab[0] + ab[1]) - seed))
    return brentq(g, a, b, xtol=1e-9)


def _table1(nu_a, n, T):
    return bracket_root(lambda nu: decoupled_dp(nu_a, nu), T['T1'][nu_a][n])


def _table2(nu_a, n, T):
    st, tabs = get_star(nu_a)
    rm = MD.matching_radius(st)
    f = lambda nu: MD.det(st, tabs, omega_of(nu), rm)
    re_, im_ = T['T2'][nu_a][n]
    if im_ == 0.0:
        hi = min(re_ * 1.04, nu_a * (1 - 1e-6))
        xs = np.linspace(re_ * 0.96, hi, 9)
        v = [f(x).real for x in xs]
        c = [i for i in range(8) if v[i] * v[i + 1] < 0]
        if not c:
            return {'paper': [re_, im_]}
        i = min(c, key=lambda i: abs(0.5 * (xs[i] + xs[i + 1]) - re_))
        z = brentq(lambda nu: f(nu).real, xs[i], xs[i + 1], xtol=1e-9)
        return {'paper': [re_, im_], 'nu': [z, 0.0], 'bound': True, 'im_det_rel': float(abs(f(z + 0.01).imag) / abs(f(z + 0.01)))}
    try:
        z = MD.muller(f, complex(re_, im_ * 1e-3))
    except RuntimeError as e:
        return {'paper': [re_, im_], 'err': str(e)[:120]}
    scale = max(abs(f(z + 0.05)), 1e-300)
    return {'paper': [re_, im_], 'nu': [z.real, z.imag], 'bound': bool(z.imag == 0.0),
            'incoming': float(abs(MD.det(st, tabs, omega_of(z), rm, sheet=-1)) / scale)}


def _table3(nu_a, k, T):
    st, tabs = get_star(nu_a)
    rm = MD.matching_radius(st)
    f = lambda nu: MD.det(st, tabs, omega_of(nu), rm)
    re_, im_ = T['T3'][nu_a][k]
    try:
        z = MD.muller(f, complex(re_, im_))
        return {'paper': [re_, im_], 'nu': [z.real, z.imag]}
    except RuntimeError as e:
        return {'paper': [re_, im_], 'nu': None, 'err': str(e)[:120]}


def _task26(args):
    kind, nu_a, n = args
    t = time.time()
    r = {'t1': _table1, 't2': _table2, 't3': _table3}[kind](nu_a, n, _W['T26'])
    return kind, nu_a, n, r, time.time() - t


def load_paper1_tables(path):
    d = json.loads(Path(path).read_text(encoding='utf-8'))
    key = lambda k: None if k == 'none' else float(k)
    return {'NUA': tuple(d['NUA']), 'T1': {key(k): tuple(v) for k, v in d['T1'].items()},
            'T2': {key(k): tuple(tuple(x) for x in v) for k, v in d['T2'].items()}, 'T3': {key(k): tuple(tuple(x) for x in v) for k, v in d['T3'].items()}}


def step26(T):
    _W['T26'] = T
    NUA = T['NUA']
    tasks = ([('t2', x, n) for x in NUA for n in range(4)] + [('t3', x, k) for x in NUA for k in range(2)] + [('t1', x, n) for x in (None,) + NUA for n in range(4)])
    pool = _pool()
    try:
        got = _map(pool, _task26, tasks)
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    D = {'t1': {}, 't2': {}, 't3': {}, 'seconds': {}}
    for kind, nu_a, n, r, sec in got:
        D[kind].setdefault(str(nu_a), {})[n] = r
        D['seconds'][f'{kind}:{nu_a}:{n}'] = sec
    for kind in ('t1', 't2', 't3'):
        D[kind] = {k: [v[i] for i in sorted(v)] for k, v in D[kind].items()}
    D['source'] = 'Paper I (arXiv:2510.16245v2), Tablas 1-3'
    D['t1_paper'] = {str(k): v for k, v in T['T1'].items()}
    T1, T2, T3 = T['T1'], T['T2'], T['T3']
    res = []
    dev1 = [(nu_a, n, (z - T1[nu_a][n]) / T1[nu_a][n]) for nu_a in NUA for n, z in enumerate(D['t1'][str(nu_a)]) if z is not None]
    miss1 = sum(z is None for nu_a in NUA for z in D['t1'][str(nu_a)])
    m1 = max(abs(d) for *_, d in dev1)
    worst1 = max(dev1, key=lambda t: abs(t[2]))
    D['t1_maxdev'] = m1
    res.append((f'Tabla 1 (fondo con axión, δθ apagado): {len(dev1)}/28 frecuencias halladas; desviación máxima {m1 * 100:.2f} % (f{worst1[1]} a {worst1[0]} kHz)', miss1 == 0 and m1 <= 0.01, 'H34(a)'))
    d0 = [(z - T1[None][n]) / T1[None][n] for n, z in enumerate(D['t1']['None'])]
    res.append(('Sin axión frente a la fila sin axión de la Tabla 1: ' + ', '.join(f'{d * 100:+.2f} %' for d in d0), None, 'descriptivo: la prueba previa del cuaderno, ahora en el paso'))
    rows = [(nu_a, n, m) for nu_a in NUA for n, m in enumerate(D['t2'][str(nu_a)])]
    found = [(nu_a, n, m) for nu_a, n, m in rows if 'nu' in m]
    dre = [abs(m['nu'][0] - m['paper'][0]) / m['paper'][0] for *_, m in found]
    cls = all((m['paper'][1] == 0.0) == m['bound'] for *_, m in found)
    bound_real = all(m['im_det_rel'] < 1e-8 for *_, m in found if m['bound'])
    D['t2_maxdev_re'] = max(dre)
    res.append((f'Tabla 2, partes reales: {len(found)}/28 modos hallados; desviación máxima {max(dre) * 100:.2f} %; los {sum(m["bound"] for *_, m in found)} modos con Im = 0 del paper son '
                f'ligados también aquí (frecuencia real)', len(found) == 28 and max(dre) <= 0.01 and cls and bound_real, 'H34(b)'))
    damped = [(nu_a, n, m) for nu_a, n, m in found if not m['bound']]
    ratios = [(nu_a, n, (m['nu'][1] * 1e3) / m['paper'][1]) for nu_a, n, m in damped]
    D['t2_im_ratio'] = [[nu_a, n, q] for nu_a, n, q in ratios]
    out25 = [(nu_a, n, q) for nu_a, n, q in ratios if abs(q - 1) > 0.25]
    res.append((f'Tabla 2, amortiguamientos: {len(ratios) - len(out25)}/{len(ratios)} dentro del 25 %; cociente nuestro/paper entre {min(q for *_, q in ratios):.2f} y {max(q for *_, q in ratios):.2f}; fuera: '
                + (', '.join(f'f{n} a {nu_a} kHz ({q:.2f})' for nu_a, n, q in out25) or 'ninguno'), len(out25) <= 1, 'H34(c): criterio fijado antes de calcular, como mucho un modo fuera'))
    t3 = [(nu_a, k, m) for nu_a in NUA for k, m in enumerate(D['t3'][str(nu_a)])]
    conv = [(nu_a, k, m) for nu_a, k, m in t3 if m['nu']]
    near = [(nu_a, k, m) for nu_a, k, m in conv if abs(complex(*m['nu']) - complex(*m['paper'])) < 0.1 * abs(m['paper'][1])]
    dre3 = [100 * (m['nu'][0] / m['paper'][0] - 1) for *_, m in conv]
    rim3 = [m['nu'][1] / m['paper'][1] for *_, m in conv]
    D['t3_stats'] = {'n_conv': len(conv), 'n': len(t3), 're_pct': [min(dre3), max(dre3)], 'im_ratio': [min(rim3), max(rim3)], 'max_im_conv': max(abs(m['paper'][1]) for *_, m in conv),
                     'min_im_fail': min(abs(m['paper'][1]) for *_, m in t3 if not m['nu'])}
    res.append((f'Tabla 3 (modos del axión): Muller converge en {len(conv)}/{len(t3)} semillas y en {len(near)} queda cerca del valor del paper (a menos de un 10 % de su |Im|). En las que convergen, '
                f'Re ν se desvía entre {min(dre3):+.2f} % y {max(dre3):+.2f} % e Im ν / paper va de {min(rim3):.3f} a {max(rim3):.3f}. Las que no convergen tienen |Im| ≥ '
                f'{D["t3_stats"]["min_im_fail"]:.2f} kHz en el paper; las que sí, ≤ {D["t3_stats"]["max_im_conv"]:.2f}', None, 'H34(d), descriptivo'))
    neg1 = min(abs(D['t1']['None'][n] - T1[x][n]) / T1[x][n] for x in NUA if x >= 4 for n in range(4))
    res.append((f'Control negativo 1: la estrella sin axión se desvía de la Tabla 1 al menos un {neg1 * 100:.1f} % para ν_a ≥ 4 kHz', neg1 > 0.05, ''))
    inc = [m['incoming'] for *_, m in damped]
    res.append((f'Control negativo 2 (onda entrante): en las {len(inc)} raíces amortiguadas el determinante deja de anularse (|det|/escala ≥ {min(inc):.1e})', min(inc) > 1e-3, ''))
    D['checks'] = [{'name': t, 'ok': ok_field(o), 'note': n_} for t, o, n_ in res]
    return res, D


# =============================================================================================== step 15
P0_SCAN = tuple(float(x) for x in range(700, 1101, 10))


def exterior_vs_37(st6, flat_e, drop_mass=False):
    """Derived exterior equation with theta_0 = 0 against Eq. (3.7) of Paper I (verified in the derivation modules)."""
    ma = omega_of(6.0)
    E0 = st6.E0
    fa = np.sqrt(st6.beta * E0 / (4 * ma**2))
    M, r = 2.17, 20.0
    f = 1 - 2 * M / r
    w = omega_of(4.2)
    P1 = M / (r * (r - 2 * M))
    P2 = -2 * (r - M) * M / (r * (r - 2 * M))**2
    G2 = 0.0 if drop_mass else -st6.beta / 4
    bg = (0.5 * np.log(f), P1, P2, -0.5 * np.log(f), -P1, 0, 0, 0.0, 0.0, 0, 0, 1.0, fa, 0.0, 0.0, E0 * G2, 0.0, -E0 * G2)
    A = MD._mat(flat_e, 2, r, w, bg)
    want = np.array([(ma**2 - w**2 / f) / f, -(2 / r) * (1 - M / r) / f])
    return float(max(abs(A[1, 0] - want[0]) / abs(want[0]), abs(A[1, 1] - want[1]) / abs(want[1])))


def wF2(st, tabs, lo=-3e-3, hi=3e-3):
    g = lambda w2: MD.interior(st, tabs, np.sqrt(complex(w2)), 'fluid').y[1, -1].real
    while g(lo) * g(hi) > 0:
        lo, hi = lo * 2, hi * 2
    return brentq(g, lo, hi, xtol=1e-12)


def turning_point():
    """Maximum mass and omega_F^2 = 0 of the star without axion, from the fixed grid of extra stars (cubic splines)."""
    p0s = np.array(P0_SCAN)
    Ms, w2s = [], []
    for p0 in p0s:
        st, tabs = get_star(None, p0)
        Ms.append(st.M)
        w2s.append(wF2(st, tabs))
    cs = CubicSpline(p0s, Ms)
    i = int(np.argmax(Ms))
    pm = brentq(cs.derivative(), p0s[max(i - 1, 0)], p0s[min(i + 1, len(p0s) - 1)])
    cw = CubicSpline(p0s, w2s)
    j = next(k for k in range(len(w2s) - 1) if w2s[k] * w2s[k + 1] < 0)
    pz = brentq(cw, p0s[j], p0s[j + 1])
    return float(pm), float(pz), float(cs(pm) / MSUN_KM), {'p0': p0s.tolist(), 'M_km': Ms, 'wF2': w2s}


def real_root(st, tabs, lo, hi, rm, n=17):
    """Real root of det in [lo, hi] (bound states below the threshold). det is not entire: it has poles where the exterior wave vanishes at r_m (that is why det_entire exists), and
    with a shorter matching radius than the legacy one a pole can fall inside the bracket (the end values then have the same sign). So: scan, take every sign change, refine it
    and keep only those where |det| goes to zero (a pole does not); return the one closest to the middle of the bracket."""
    f = lambda nu: MD.det(st, tabs, omega_of(nu), rm).real
    xs = np.linspace(lo, hi, n)
    v = [f(x) for x in xs]
    scale = float(np.median(np.abs(v)))
    roots = []
    for i in range(n - 1):
        if v[i] * v[i + 1] < 0:
            r = brentq(f, xs[i], xs[i + 1], xtol=1e-9)
            if abs(f(r)) < 1e-4 * scale:
                roots.append(r)
    if not roots:
        raise ValueError(f'real_root: no zero of det in [{lo:.6f}, {hi:.6f}]')
    return min(roots, key=lambda r: abs(r - 0.5 * (lo + hi)))


def complex_root(st, tabs, seed, rm):
    f = lambda nu: MD.det(st, tabs, omega_of(nu), rm)
    z = MD.muller(f, seed)
    return z, abs(f(z))


def step15():
    res, out = [], {}
    st6, tabs6 = get_star(6.0)
    rm6 = MD.matching_radius(st6)
    Mt, fa_v = st6.Mt, np.sqrt(st6.fa2)
    e37, e37_bad = exterior_vs_37(st6, tabs6.flat_e), exterior_vs_37(st6, tabs6.flat_e, drop_mass=True)
    res.append((f'Ecuación exterior derivada (theta_0 = 0) = Ec. (3.7) del Paper I: error relativo {e37:.1e}', e37 < 1e-12, ''))
    res.append((f'Control negativo 1 (sin término de masa, el error real del 2026-11-07): error {e37_bad:.2f}', e37_bad > 0.1, ''))

    pm, pz, Mmax, scan = turning_point()
    out['turning'] = {'p0_Mmax': pm, 'p0_wF2_zero': pz, 'Mmax_msun': Mmax}
    res.append((f'Sin axión (BSk26): masa máxima {Mmax:.4f} M☉ en p_c = {pm:.1f} MeV/fm³; ω_F² = 0 en p_c = {pz:.1f} (diferencia {abs(pm - pz) / pm * 100:.2f} %)',
                abs(pm - pz) / pm < 0.01, 'criterio de Chandrasekhar/punto de retorno; malla fija de estrellas y splines cúbicos'))

    z0, h = omega_of(6.585 - 0.563j), omega_of(1e-3)
    Lc = lambda om: (lambda p: p[1] / p[0])(MD.exterior_psi(tabs6.flat_e, om, Mt, rm6, st6.ma, fa_v, st6.E0, st6.beta))
    Lr = lambda om: MD.exterior_L(st6, 0, om, r_stop=rm6)
    cr = lambda L: abs((L(z0 + h) - L(z0)) / h - (L(z0 + 1j * h) - L(z0)) / (1j * h)) / abs((L(z0 + h) - L(z0)) / h)
    crc = cr(Lc)
    rel = lambda z: abs(Lr(omega_of(z)) - Lc(omega_of(z))) / abs(Lc(omega_of(z)))
    d_strong, d_weak, d_below = rel(6.585 - 0.563j), rel(8.132 - 0.0036j), rel(5.5)
    crr = d_strong
    res.append((f'Exterior en el rayo complejo: cumple Cauchy-Riemann (desajuste {crc:.1e}). El de la recta real (axstar.exterior_L) se desvía {d_strong * 100:.0f} % en la resonancia muy '
                f'amortiguada, {d_weak:.0e} en un modo poco amortiguado y {d_below:.0e} bajo el umbral', crc < 1e-2 and d_strong > 0.5 and d_weak < 1e-3 and d_below < 1e-8,
                'la onda saliente crece con r y la integración hacia dentro es inestable; en los modos poco amortiguados de los pasos 9-10 este error es pequeño (no comprueba otras fuentes de error)'))

    F6 = real_root(st6, tabs6, 3.40, 3.50, rm6)
    fF6 = MD.fluid_fraction(st6, tabs6, omega_of(F6), rm6)
    d45 = MD.det(st6, tabs6, omega_of(4.5), rm6)
    imrel = abs(d45.imag) / abs(d45)
    res.append((f'ν_a = 6 kHz: modo F en {F6:.4f} kHz < ν_a, frecuencia real (|Im det|/|det| = {imrel:.0e} en el eje), fracción fluida {fF6:.3f}', imrel < 1e-10 and fF6 > 0.9, 'H11(b)'))

    alpha_g = Mt * st6.ma
    ladder = []
    for n in range(2, 7):
        pred = 6.0 * (1 - alpha_g**2 / (2 * n * n))
        dn = 6.0 * alpha_g**2 / n**3 * 0.3
        z = real_root(st6, tabs6, pred - dn, min(pred + dn, 6.0 * (1 - 1e-6)), rm6)
        ladder.append((n, z, alpha_g / np.sqrt(2 * (1 - z / 6.0)), MD.fluid_fraction(st6, tabs6, omega_of(z), rm6) if n <= 3 else None))
    ok_l = all(abs(ne - n) < 0.05 for n, _, ne, _ in ladder) and all(fr < 0.05 for *_, fr in ladder if fr is not None)
    res.append((f'Átomo gravitatorio (α_g = {alpha_g:.4f}): n_eff = ' + ', '.join(f'{ne:.3f}' for _, _, ne, _ in ladder) + ' para n = 2…6; carácter axiónico', ok_l, 'H11(c)'))

    H1, rH1 = complex_root(st6, tabs6, 8.13 - 0.004j, rm6)
    fH1 = MD.fluid_fraction(st6, tabs6, omega_of(H1), rm6)
    res.append((f'ν_a = 6 kHz: H1 = {H1.real:.4f} {H1.imag:+.4f}i kHz (fracción fluida {fH1:.3f})', rH1 < 1e-11 and H1.imag < 0 and fH1 > 0.9, 'H11(b)'))

    # retired artefact: the 'axion branch' of very damped resonances was a cavity between R and r_m created by the jump of omega^2 e^{-2 Phi} at the matching with the old lapse
    # the legacy window (6.2 < Re < 7.0, -0.8 < Im < -0.3) was tied to its r_m of 117 km: a cavity resonance moves with r_m (here 6.290 - 0.802i at r_m = 80 km, 8.261 - 1.742i at 64 km). The claim
    # of the control is that the old normalisation creates a damped resonance above the axion threshold that is not a root with the corrected one, so the window is that region.
    win = lambda z: 6.0 < z.real < 12.0 and -1.5 < z.imag < -0.1
    st_old, tabs_old = old_lapse(st6)
    r84 = 0.8 * rm6
    old_a, _ = complex_root(st_old, tabs_old, 6.585 - 0.563j, rm6)
    old_b, _ = complex_root(st_old, tabs_old, 6.585 - 0.563j, r84)
    try:                                    # H25: Muller no longer returns an unconverged approximation
        new_z, _ = complex_root(st6, tabs6, 6.585 - 0.563j, rm6)
    except RuntimeError:
        new_z = complex(np.nan, np.nan)
    moved = abs(old_a - old_b)
    res.append((f'Artefacto retirado: con la normalización antigua del lapso, la «resonancia axiónica» se mueve de {old_a.real:.3f} {old_a.imag:+.3f}i a {old_b.real:.3f} {old_b.imag:+.3f}i al '
                f'cambiar el empalme ({rm6:.0f} → {r84:.0f} km); con la corrección, Muller desde esa semilla sale de la región (hacia {new_z.real:.2f} {new_z.imag:+.2f}i); el recuento de la '
                f'región lo hace el paso 20', moved > 0.1 and win(old_a) and not win(new_z), 'era una cavidad entre R y el radio de empalme (ventana: por encima del umbral y amortiguada, la del original dependía de su r_m; el segundo empalme es 0.8 r_m)'))
    d_ok = abs(MD.det(st6, tabs6, omega_of(H1), rm6, sheet=-1)) / max(abs(MD.det(st6, tabs6, omega_of(H1 + 0.05), rm6)), 1e-300)
    res.append((f'Control negativo 2 (onda entrante): en la raíz H1 el determinante deja de anularse (|det|/escala = {d_ok:.1e})', d_ok > 1e-3, ''))

    st35, tabs35 = get_star(3.5)
    st34, tabs34 = get_star(3.4)
    rm35, rm34 = MD.matching_radius(st35), MD.matching_radius(st34)
    b35 = [real_root(st35, tabs35, 3.38, 3.43, rm35), real_root(st35, tabs35, 3.44, 3.48, rm35)]
    f35 = [MD.fluid_fraction(st35, tabs35, omega_of(z), rm35) for z in b35]
    F34, rF34 = complex_root(st34, tabs34, 3.42 - 0.001j, rm34)
    fF34 = MD.fluid_fraction(st34, tabs34, omega_of(F34), rm34)
    res.append((f'ν_a = 3,5 kHz: dos ligados {b35[0]:.4f} y {b35[1]:.4f} kHz con carácter mezclado ({f35[0]:.2f} y {f35[1]:.2f} de fluido): cruce evitado F – átomo n = 1',
                all(0.1 < x < 0.9 for x in f35), ''))
    res.append((f'ν_a = 3,4 kHz: F ya por encima del umbral, Re ν = {F34.real:.4f} kHz e Im ν = {F34.imag * 1e3:.2f} Hz (fracción fluida {fF34:.3f})',
                rF34 < 1e-11 and F34.real > 3.4 and F34.imag < 0 and fF34 > 0.9, 'el modo fluido pasa a cuasi-normal al cruzar ν_a'))

    stab = {}
    for lab, st_, tabs_, rm_ in (('6', st6, tabs6, rm6), ('3.4', st34, tabs34, rm34)):
        gs = np.geomspace(0.01, 10, 40)
        v = [MD.det(st_, tabs_, omega_of(1j * g), rm_) for g in gs]
        stab[lab] = {'sign_changes': int(sum(v[i].real * v[i + 1].real < 0 for i in range(len(gs) - 1))), 'im_rel': float(max(abs(x.imag) / abs(x) for x in v))}
    res.append(('Estabilidad radial: sin raíces en ω = iγ, γ ∈ [0,01; 10] kHz, para ν_a = 6 y 3,4 kHz (' + ', '.join(f'{k}: {v_["sign_changes"]} cambios de signo' for k, v_ in stab.items()) + ')',
                all(v_['sign_changes'] == 0 and v_['im_rel'] < 1e-10 for v_ in stab.values()), 'H12; un barrido de signo puede no ver pares de raíces muy próximas'))
    robust = {}
    for lab, st_, tabs_, rm_, z0_ in (('H1 (6 kHz)', st6, tabs6, rm6, H1), ('F (3,4 kHz)', st34, tabs34, rm34, F34)):
        rows = {'base': z0_}
        for name, kw in (('rtol 1e-12', {'rtol': 1e-12}), ('r0 3e-3 km', {'r0': 3e-3}), ('exterior ×1,5', {'n_decay': 60.0})):
            rows[name] = MD.muller(lambda nu: MD.det(st_, tabs_, omega_of(nu), rm_, **kw), z0_)
        for fac in (0.8, 0.6):
            rows[f'empalme {fac:g} r_m ({fac * rm_:.0f} km)'] = MD.muller(lambda nu: MD.det(st_, tabs_, omega_of(nu), fac * rm_), z0_)
        robust[lab] = {k: [v_.real, v_.imag] for k, v_ in rows.items()}
    num = max(abs(v_[1] - r_['base'][1]) / abs(r_['base'][1]) for r_ in robust.values() for k, v_ in r_.items() if k in ('rtol 1e-12', 'r0 3e-3 km', 'exterior ×1,5'))
    emp = {lab: abs(list(r_.values())[-1][1] - r_['base'][1]) / abs(r_['base'][1]) for lab, r_ in robust.items()}
    res.append((f'Robustez de Im ν: rtol, r0 y longitud exterior la cambian como mucho {num:.0e} relativo; el radio de empalme (último paso hacia fuera) '
                + ', '.join(f'{v_:.0e} en {k}' for k, v_ in emp.items()), num < 5e-4 and all(v_ < 0.01 for v_ in emp.values()), 'H13; el radio de empalme es la incertidumbre dominante'))

    out.update({'stability': stab, 'robustness': robust,
                'nu6': {'M_msun': st6.M / MSUN_KM, 'R_km': st6.R, 'alpha_g': alpha_g, 'F': F6, 'F_fluid': fF6,
                        'ladder': [{'n': n, 'nu': z, 'n_eff': ne} for n, z, ne, _ in ladder], 'H1': [H1.real, H1.imag], 'H1_fluid': fH1},
                'artifact': {'old_rm': [old_a.real, old_a.imag], 'old_r84': [old_b.real, old_b.imag], 'rm': rm6, 'r84': r84, 'new_escape': [new_z.real, new_z.imag]},
                'nu35': {'bound': b35, 'fluid': f35}, 'nu34': {'F': [F34.real, F34.imag], 'fluid': fF34},
                'exterior': {'cr_complex_ray': crc, 'cr_real_line': crr, 'old_vs_new_strong': d_strong, 'old_vs_new_weak': d_weak},
                'turning_scan': scan, 'checks': [{'name': n_, 'ok': ok_field(o_), 'note': t_} for n_, o_, t_ in res]})
    return res, out


# =============================================================================================== step 19
NU_A19 = [3.46, 3.48, 3.4857142857142858, 3.4914285714285715, 3.4971428571428573, 3.5, 3.5028571428571427, 3.5085714285714285, 3.5142857142857142, 3.52, 3.54, 3.57, 3.6, 3.65, 3.7]


def two_lowest(nu_a, lo=3.25, step=0.002):
    st, tabs = get_star(nu_a)
    f = lambda nu: MD.det_entire(st, tabs, omega_of(nu)).real
    grid = np.arange(lo, nu_a - 0.004, step)
    v = [f(x) for x in grid]
    roots = []
    for i in range(len(grid) - 1):
        if v[i] * v[i + 1] < 0:
            roots.append(brentq(f, grid[i], grid[i + 1], xtol=1e-9))
            if len(roots) == 2:
                break
    fr = [MD.fluid_fraction(st, tabs, omega_of(r)) for r in roots]
    return {'nu_a': nu_a, 'roots': roots, 'fluid': fr, 'M': st.M}


def step19(nus=None):
    pool = _pool()
    try:
        rows = sorted(_map(pool, two_lowest, list(nus or NU_A19)), key=lambda r: r['nu_a'])
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    ok_rows = [r for r in rows if len(r['roots']) == 2]
    gap = [(r['nu_a'], r['roots'][1] - r['roots'][0]) for r in ok_rows]
    i = int(np.argmin([g for _, g in gap]))
    nu_min, g_min = gap[i]
    interior_ = 0 < i < len(gap) - 1
    lo_side, hi_side = ok_rows[0], ok_rows[-1]
    swap = (lo_side['fluid'][0] < 0.5 < lo_side['fluid'][1]) and (hi_side['fluid'][0] > 0.5 > hi_side['fluid'][1])
    res = [(f'H21(a): separación de los dos ligados más bajos con mínimo {1e3 * g_min:.1f} Hz en nu_a = {nu_min:.3f} kHz (interior a la malla, refinada; {len(gap)} valores de nu_a)', interior_ and g_min > 0,
            'la malla fina es la que el original eligió alrededor del mínimo'),
           (f'H21(b): intercambio de carácter: a nu_a = {lo_side["nu_a"]:.2f} kHz fracción fluida {lo_side["fluid"][0]:.2f} (inferior) y {lo_side["fluid"][1]:.2f} (superior); '
            f'a {hi_side["nu_a"]:.2f} kHz, {hi_side["fluid"][0]:.2f} y {hi_side["fluid"][1]:.2f}', swap, '')]
    out = {'rows': rows, 'gap_min_khz': g_min, 'nu_a_min': nu_min, 'checks': [{'name': n_, 'ok': ok_field(o_), 'note': t_} for n_, o_, t_ in res]}
    return res, out


# =============================================================================================== step 20
TILES20 = [(6.05, 9.0, -0.3, -0.002), (9.0, 12.0, -0.3, -0.002), (6.05, 9.0, -1.5, -0.3), (9.0, 12.0, -1.5, -0.3)]
TILES_THIN = [(6.05, 9.0, -0.002, -1e-5), (9.0, 12.0, -0.002, -1e-5)]
TILES_DEEP = [(6.05, 9.0, -5.0, -1.5), (9.0, 12.0, -5.0, -1.5)]
PAPER2_A0 = 9.0 - 3.8j
_S = {}


def _F20(nu):
    return complex(MD.det_entire(_S['st'], _S['tabs'], omega_of(nu), rm=_S.get('rm'), inward=False, sheet=1))


def count20(st, tabs, tiles=None, n0=48, n_max=384):
    _S.update(st=st, tabs=tabs)
    pool = _pool()
    try:
        w = [MD.winding_converged(_F20, _corners(b), n0=n0, n_max=n_max, pool=pool) for b in (tiles or TILES20)]
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    return {'per_tile': [x[0] for x in w], 'tried': [x[4] for x in w], 'points': sum(x[2] for x in w), 'max_dphase': float(max(x[3] for x in w)),
            'total': sum(x[0] for x in w) if all(x[0] is not None for x in w) else None}


def step20(n0=48, n_max=384):
    res, out = [], {}
    st, tabs = get_star(6.0)
    rm = MD.matching_radius(st)
    c = count20(st, tabs, None, n0, n_max)
    found = []
    for seed in (8.119 - 0.0036j, 10.6 - 0.01j, 11.7 - 0.01j):
        try:
            z = MD.muller(_F20, seed, tol=1e-9)
            if 6.05 < z.real < 12 and -1.5 < z.imag < -0.002 and all(abs(z - y) > 1e-3 for y in found):
                found.append(complex(z))
        except RuntimeError:
            pass
    fr = [MD.fluid_fraction(st, tabs, omega_of(z)) for z in found]
    out['nu6'] = {**c, 'fluid_modes': [[z.real, z.imag] for z in found], 'fluid_fraction': fr}
    res.append((f'H22(a): nu_a = 6 kHz, índice {c["total"]} en Re nu in [6.05, 12], Im nu in [-1.5, -0.002] kHz (piezas {c["per_tile"]}, convergidas {c["tried"]}); Muller encuentra '
                f'{len(found)} modos fluidos ahí: ' + ', '.join(f'{z.real:.4f}{z.imag * 1e3:+.2f}i (fluido {f:.2f})' for z, f in zip(found, fr)), c['total'] == len(found),
                'sin ceros adicionales: ninguna familia axiónica en esa región'))
    ct = count20(st, tabs, TILES_THIN, n0, n_max)
    out['thin'] = ct
    res.append((f'H26: franja casi real Re nu in [6.05, 12], Im nu in [-0.002, -1e-5] kHz: índice {ct["total"]} (piezas {ct["per_tile"]}, convergidas {ct["tried"]}; salto de fase máximo '
                f'{ct["max_dphase"]:.2f} rad)', ct['total'] == 0, 'fuera quedan |Im nu| < 0.01 Hz, Re nu < 6.05 kHz y Re nu > 12 kHz'))
    cd = count20(st, tabs, TILES_DEEP, n0, n_max)
    deep = {**cd, 'roots': []}
    if cd['total']:
        r_alt = 0.8 * rm
        seeds = [PAPER2_A0] + [complex(x, y) for x in (7.0, 8.5, 10.0, 11.5) for y in (-2.0, -3.0, -4.0)]
        for sd in seeds:
            try:
                z = complex(MD.muller(_F20, sd, tol=1e-9))
            except RuntimeError:
                continue
            if 6.05 < z.real < 12 and -5 < z.imag < -1.5 and all(abs(z - y['nu']) > 1e-3 for y in deep['roots']):
                _S['rm'] = r_alt
                try:
                    z2 = complex(MD.muller(_F20, z, tol=1e-9))
                except RuntimeError:
                    z2 = complex(np.nan, np.nan)
                _S.pop('rm')
                deep['roots'].append({'nu': z, 'nu_alt_rm': z2, 'r_alt': r_alt, 'shift_pct': float(100 * abs(z2 / z - 1)), 'fluid': MD.fluid_fraction(st, tabs, omega_of(z))})
        deep['roots'] = [{**r_, 'nu': [r_['nu'].real, r_['nu'].imag], 'nu_alt_rm': [r_['nu_alt_rm'].real, r_['nu_alt_rm'].imag]} for r_ in deep['roots']]
    out['deep'] = deep
    robust = [r_ for r_ in deep['roots'] if r_['shift_pct'] < 0.1]
    if not cd['total']:
        txt = 'ningún cero: no hay análogo radial de a0 en esa región'
    elif not deep['roots']:
        txt = f'ningún cero localizado: Muller no converge desde ninguna de las {len(seeds)} semillas'
    else:
        txt = 'ceros: ' + '; '.join(f"{r_['nu'][0]:.3f}{r_['nu'][1]:+.3f}i kHz (fluido {r_['fluid']:.2f}; al cambiar el empalme se mueve un {r_['shift_pct']:.2g} %)" for r_ in deep['roots'])
    res.append((f'H23: región profunda Re nu in [6.05, 12], Im nu in [-5, -1.5] kHz: índice {cd["total"]} (piezas {cd["per_tile"]}, convergidas {cd["tried"]}); {txt}. Secuela, Tabla III '
                f'(l = 2, Cowling): a0 = {PAPER2_A0.real} {PAPER2_A0.imag:+.1f}i kHz', cd['total'] is not None and (cd['total'] == 0 or len(robust) == cd['total']),
                'el resultado decide; no había predicción fuerte'))
    # positive control: the old lapse normalisation (Phi matched at R with m(R), no tail): the count must see the cavity resonances
    st_o, tabs_o = old_lapse(st)
    m_R = st.y_at(st.R)[0]
    out['lapse'] = {'phi_inf_old': float(st_o.Phi_shift - st.Phi_shift), 'tail_mass_frac': float((st.Mt - m_R) / st.Mt)}
    co = count20(st_o, tabs_o, None, n0, n_max)
    out['old_lapse'] = co
    out['old_lapse_deep'] = count20(st_o, tabs_o, TILES_DEEP, n0, n_max)
    res.append((f'H22(b), control positivo: con la normalización antigua del lapso el índice sube a {co["total"]} (piezas {co["per_tile"]}): el recuento ve las resonancias de cavidad',
                co['total'] is not None and c['total'] is not None and co['total'] >= c['total'] + 3, ''))
    out['checks'] = [{'name': n_, 'ok': ok_field(o_), 'note': t_} for n_, o_, t_ in res]
    return res, out


# =============================================================================================== step 16
EPS, YMAX, XC, XMAX, GMAX = 0.05, 12.0, 2.0, 12.0, 20.0
TILES_C = [(-XC, XC, EPS, 4.0), (-XC, XC, 4.0, YMAX)]
TILES_R = [(XC, XMAX, EPS, 4.0), (XC, XMAX, 4.0, YMAX)]


def _F16(nu):
    return complex(MD.det_entire(_S['st'], _S['tabs'], omega_of(nu), **_S['kw']))


def _g16(nu):                              # without axion: Delta p(R) of the fluid seed, entire in omega
    return complex(MD.interior(_S['st'], _S['tabs'], omega_of(nu), 'fluid').y[1, -1])


def count16(st, tabs, f=None, n0=48, n_max=384, eps=EPS, kw=None):
    """Total index in |Re nu| <= XMAX, eps <= Im nu <= YMAX (central strip + 2 x right) and sign changes of F on the imaginary axis."""
    f = f or _F16
    _S.update(st=st, tabs=tabs, kw=kw or {})
    tc = [(a, b, max(c, eps) if c == EPS else c, d) for a, b, c, d in TILES_C]
    tr = [(a, b, max(c, eps) if c == EPS else c, d) for a, b, c, d in TILES_R]
    pool = _pool()
    try:
        wc = [MD.winding_converged(f, _corners(b), n0=n0, n_max=n_max, pool=pool) for b in tc]
        wr = [MD.winding_converged(f, _corners(b), n0=n0, n_max=n_max, pool=pool) for b in tr]
        gs = np.concatenate([np.geomspace(0.01, 1, 60), np.linspace(1, GMAX, 400)[1:]])
        v = np.array(_map(pool, f, list(1j * gs)))
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    conv = all(w[0] is not None for w in wc + wr)
    total = sum(w[0] for w in wc) + 2 * sum(w[0] for w in wr) if conv else None
    sc = [float(gs[i]) for i in range(len(gs) - 1) if v[i].real * v[i + 1].real < 0]
    return {'total': total, 'converged': conv, 'central': [w[0] for w in wc], 'right': [w[0] for w in wr], 'tried': [w[4] for w in wc + wr], 'n0_final': [w[1] for w in wc + wr],
            'points': sum(w[2] for w in wc + wr), 'max_dphase': max(w[3] for w in wc + wr), 'axis_sign_changes': sc, 'axis_imag_rel': float(np.max(np.abs(v.imag) / np.abs(v)))}


def symmetry(st, tabs, sheet=0):
    _S.update(st=st, tabs=tabs, kw={'sheet': sheet})
    pts = [2 + 1j, 5 + 0.3j, 9 + 6j, 11 + 11j]
    return max(abs(_F16(-z.conjugate()) - np.conj(_F16(z))) / abs(_F16(z)) for z in pts)


def step16(n0=48, n_max=384):
    res, out = [], {}
    ctrl = {}
    for p0 in (1000.0, 800.0):
        st_c, tabs_c = get_star(None, p0)
        _S.update(st=st_c, tabs=tabs_c)
        pool = _pool()
        try:
            ctrl[p0] = MD.winding_converged(_g16, _corners((-XMAX, XMAX, EPS, YMAX)), n0=n0, n_max=n_max, pool=pool)[0]
        finally:
            if pool is not None:
                pool.close()
                pool.join()
    res.append((f'Control positivo sin axión: índice {ctrl[1000.0]} más allá de la masa máxima (p_c = 1000) y {ctrl[800.0]} antes (p_c = 800)', ctrl[1000.0] == 1 and ctrl[800.0] == 0, 'H18(b)'))
    out['control_no_axion'] = {str(k): v for k, v in ctrl.items()}
    stc, tabsc = get_star(6.0, 1000.0)
    cc = count16(stc, tabsc, n0=n0, n_max=n_max)
    res.append((f'Control positivo con axión (nu_a = 6 kHz, p_c = 1000): índice {cc["total"]}; cambio de signo en el eje imaginario en gamma = {", ".join(f"{g:.2f}" for g in cc["axis_sign_changes"])} kHz',
                cc['total'] == 1 and len(cc['axis_sign_changes']) == 1, 'H18(c)'))
    out['control_axion'] = cc
    for nu_a in (6.0, 3.4):
        st, tabs = get_star(nu_a)
        sym = symmetry(st, tabs)
        c = count16(st, tabs, n0=n0, n_max=n_max)
        out[f'nu{nu_a}'] = {**c, 'symmetry': sym}
        ok = c['total'] == 0 and not c['axis_sign_changes'] and c['axis_imag_rel'] < 1e-10 and sym < 1e-10
        res.append((f'nu_a = {nu_a} kHz: ningún cero de F con |Re nu| <= {XMAX:.0f}, {EPS} <= Im nu <= {YMAX:.0f} kHz (índice {c["total"]}: central {c["central"]}, derecha x2 {c["right"]}, cada pieza '
                    f'estable al doblar el muestreo, {c["tried"]}; {c["points"]} puntos) y ningún cambio de signo en el eje hasta {GMAX:.0f} kHz; simetría F(-w*) = F(w)* a {sym:.0e}', ok, 'H18(a)'))
        if nu_a == 6.0:
            _S.update(st=st, tabs=tabs, kw={})
            pool = _pool()
            try:
                w2 = MD.winding(_F16, _corners(TILES_C[0]), n0=4 * n0, pool=pool)[0]
                w3 = MD.winding_converged(_F16, _corners((-XC, XC, EPS / 2, 4.0)), n0=n0, n_max=n_max, pool=pool)[0]
                _S.update(kw={'raw': 'inf'})
                wbad = MD.winding(_F16, _corners((0.0, XMAX, EPS, 4.0)), n0=12, pool=pool)[0]
            finally:
                if pool is not None:
                    pool.close()
                    pool.join()
            out['robust6'] = {'n0_96': w2, 'eps_half': w3, 'neg_raw_inf': wbad}
            res.append((f'Robustez (6 kHz, franja central baja): índice {round(w2)} con {4 * n0} puntos por lado y {w3} con eps/2', round(w2) == c['central'][0] and w3 == c['central'][0], 'H18(d)'))
            res.append((f'Control negativo 1: sin quitar el factor e^(i k r_m) y con muestreo grueso, el contorno pierde vueltas (índice {round(wbad)} en vez de 0)', round(wbad) != 0, ''))
            symbad = symmetry(st, tabs, sheet=1)
            out['neg_sheet_principal'] = symbad
            res.append((f'Control negativo 2: con la rama principal de k también para Re omega < 0, la simetría F(-w*) = F(w)* se rompe ({symbad:.1e} frente a {sym:.0e})', symbad > 1e-3, ''))
    out['region'] = {'eps': EPS, 'ymax': YMAX, 'xmax': XMAX, 'gmax': GMAX}
    out['checks'] = [{'name': n_, 'ok': ok_field(o_), 'note': t_} for n_, o_, t_ in res]
    return res, out
