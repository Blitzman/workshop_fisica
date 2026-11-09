"""Paso 26: los modos radiales frente a las Tablas 1-3 del Paper I (hipótesis H34, séptima revisión del informe radial).

Referencia: Gómez-Bañón, Pnigouras y Pons, PRD 113, 023053 (arXiv:2510.16245v2), Sec. V. BSk26, eps = 0,1,
p0 = 100 MeV/fm^3, relatividad general completa, dependencia e^{-i omega t}. Re en kHz; Im en Hz (Tabla 2) o en kHz
(Tabla 3). Los valores del paper solo se usan como semillas y como referencia: cada raíz se acepta si el
determinante se anula (Muller con residuo, o cambio de signo en el eje real).

- Tabla 1: fondo con axión y delta theta apagado. Es la submatriz (xi, Delta p) de las mismas ecuaciones derivadas
  (radial.interior_matrix), igual que la Tabla I del paso 9 para la secuela.
- Tabla 2: problema completo (radial.det).
- Tabla 3: modos del axión, muy amortiguados. Solo descriptivo (en esa región F está en su suelo de ruido, paso 21).
"""
import json
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from . import axgr as AG
from . import axstar as AX
from . import radial as RD
from .pipeline import ok_field, report

OUT = Path(__file__).resolve().parent.parent / 'results'
NUA = (2.8, 3.0, 4.0, 5.0, 6.0, 8.0, 12.0)
# Tabla 1 (kHz): f0, f1, f2, f3; None = sin axión
T1 = {None: (3.10, 7.27, 10.22, 11.63), 2.8: (3.14, 7.32, 10.26, 11.67), 3.0: (3.32, 7.54, 10.46, 11.86),
      4.0: (3.63, 8.03, 11.25, 12.84), 5.0: (3.67, 8.14, 11.56, 13.45), 6.0: (3.68, 8.18, 11.71, 13.87),
      8.0: (3.69, 8.22, 11.88, 14.45), 12.0: (3.72, 8.27, 12.07, 15.16)}
# Tabla 2: (Re kHz, Im Hz) de f0, f1, f2, f3
T2 = {2.8: ((3.13, -8.80), (7.31, -2.60), (10.25, -1.04), (11.66, -0.394)),
      3.0: ((3.26, -43.7), (7.49, -6.10), (10.46, -2.15), (11.84, -0.854)),
      4.0: ((3.49, 0.0), (7.96, -38.9), (11.21, -4.55), (12.82, -1.06)),
      5.0: ((3.50, 0.0), (8.05, -8.38), (11.52, -7.74), (13.44, -4.02)),
      6.0: ((3.47, 0.0), (8.11, -2.05), (11.66, -22.5), (13.85, -2.01)),
      8.0: ((3.42, 0.0), (8.12, -36.2), (11.80, -6.71), (14.42, -5.22)),
      12.0: ((3.36, 0.0), (8.09, 0.0), (11.98, -12.0), (15.08, -11.4))}
# Tabla 3: (Re kHz, Im kHz) de f1a y f2a (f0a: los autores lo dan con una cifra y lo declaran poco fiable)
T3 = {2.8: ((10.50, -4.30), (19.29, -5.82)), 3.0: ((11.23, -4.31), (20.06, -5.94)),
      4.0: ((16.36, -3.45), (22.99, -5.37)), 5.0: ((20.30, -2.22), (25.54, -4.70)),
      6.0: ((23.51, -1.50), (28.22, -3.81)), 8.0: ((29.86, -0.98), (33.55, -2.50)),
      12.0: ((42.97, -1.08), (45.19, -1.50))}


def decoupled_dp(st, nu, r0=1e-3, r1=0.5, rtol=1e-10):
    """Delta p(R) con delta theta apagado: submatriz (xi, Delta p) de radial.interior_matrix (Tabla 1)."""
    A = RD.interior_matrix()
    om = AX.omega_of(nu)
    sub = lambda rr: np.asarray(A(rr, om, *RD._bgvals(st, rr)), dtype=complex)[:2, :2].real
    M = sub(r0)
    y0 = np.array([r0, (1.0 - M[0, 0] * r0) / M[0, 1]])
    f = lambda rr, y: sub(rr) @ y
    s = RD._ok(solve_ivp(f, (r0, r1), y0, method='DOP853', rtol=rtol, atol=1e-40))
    s = RD._ok(solve_ivp(f, (r1, st.R), s.y[:, -1], method='DOP853', rtol=rtol, atol=1e-40))
    return s.y[1, -1]


def bracket_root(g, seed, rel=0.04, n=9):
    """Raíz de g (real) cerca de seed: barrido en seed (1 ± rel) y brentq en el cambio de signo más cercano."""
    xs = np.linspace(seed * (1 - rel), seed * (1 + rel), n)
    v = [g(x) for x in xs]
    cands = [(xs[i], xs[i + 1]) for i in range(n - 1) if v[i] * v[i + 1] < 0]
    if not cands:
        return None
    a, b = min(cands, key=lambda ab: abs(0.5 * (ab[0] + ab[1]) - seed))
    return brentq(g, a, b, xtol=1e-9)


def _table1(nu_a, n):
    return bracket_root(lambda nu: decoupled_dp(_star(nu_a), nu), T1[nu_a][n])


_STARS = {}


def _star(nu_a):
    """Una estrella por proceso y masa: construirla cuesta ~2 min."""
    if nu_a not in _STARS:
        _STARS[nu_a] = AX.Star(nu_a) if nu_a else AX.Star(None)
    return _STARS[nu_a]


def _table2(nu_a, n):
    st = _star(nu_a)
    rm = AG.r_match(st)
    f = lambda nu: RD.det(st, AX.omega_of(nu), rm)
    re_, im_ = T2[nu_a][n]
    if im_ == 0.0:
        hi = min(re_ * 1.04, nu_a * (1 - 1e-6))
        xs = np.linspace(re_ * 0.96, hi, 9)
        v = [f(x).real for x in xs]
        c = [i for i in range(8) if v[i] * v[i + 1] < 0]
        if not c:
            return {'paper': [re_, im_]}
        i = min(c, key=lambda i: abs(0.5 * (xs[i] + xs[i + 1]) - re_))
        z = brentq(lambda nu: f(nu).real, xs[i], xs[i + 1], xtol=1e-9)
        return {'paper': [re_, im_], 'nu': [z, 0.0], 'bound': True,
                'im_det_rel': float(abs(f(z + 0.01).imag) / abs(f(z + 0.01)))}
    try:
        z = RD.muller(f, complex(re_, im_ * 1e-3))
    except RuntimeError as e:
        return {'paper': [re_, im_], 'err': str(e)[:120]}
    scale = max(abs(f(z + 0.05)), 1e-300)
    return {'paper': [re_, im_], 'nu': [z.real, z.imag], 'bound': bool(z.imag == 0.0),
            'incoming': float(abs(RD.det(st, AX.omega_of(z), rm, sheet=-1)) / scale)}


def _table3(nu_a, k):
    st = _star(nu_a)
    rm = AG.r_match(st)
    f = lambda nu: RD.det(st, AX.omega_of(nu), rm)
    re_, im_ = T3[nu_a][k]
    try:
        z = RD.muller(f, complex(re_, im_))
        return {'paper': [re_, im_], 'nu': [z.real, z.imag]}
    except RuntimeError as e:
        return {'paper': [re_, im_], 'nu': None, 'err': str(e)[:120]}


def _task(args):
    kind, nu_a, n = args
    t = time.time()
    r = {'t1': _table1, 't2': _table2, 't3': _table3}[kind](nu_a, n)
    return kind, nu_a, n, r, time.time() - t


def run():
    # Ordenadas por masa para que cada proceso reutilice su estrella; las caras (t2, t3) primero
    tasks = ([('t2', x, n) for x in NUA for n in range(4)] + [('t3', x, k) for x in NUA for k in range(2)]
             + [('t1', x, n) for x in (None,) + NUA for n in range(4)])
    with Pool(6) as p:
        got = p.map(_task, tasks, chunksize=1)
    D = {'t1': {}, 't2': {}, 't3': {}, 'seconds': {}}
    for kind, nu_a, n, r, sec in got:
        D[kind].setdefault(str(nu_a), {})[n] = r
        D['seconds'][f'{kind}:{nu_a}:{n}'] = sec
    for kind in ('t1', 't2', 't3'):
        D[kind] = {k: [v[i] for i in sorted(v)] for k, v in D[kind].items()}
    D['source'] = 'Paper I (arXiv:2510.16245v2), Tablas 1-3'
    D['t1_paper'] = {str(k): v for k, v in T1.items()}
    res = []

    # (a) Tabla 1
    dev1 = [(nu_a, n, (z - T1[nu_a][n]) / T1[nu_a][n]) for nu_a in NUA for n, z in enumerate(D['t1'][str(nu_a)])
            if z is not None]
    miss1 = sum(z is None for nu_a in NUA for z in D['t1'][str(nu_a)])
    m1 = max(abs(d) for *_, d in dev1)
    worst1 = max(dev1, key=lambda t: abs(t[2]))
    D['t1_maxdev'] = m1
    res.append((f'Tabla 1 (fondo con axión, δθ apagado): {len(dev1)}/28 frecuencias halladas; desviación máxima '
                f'{m1 * 100:.2f} % (f{worst1[1]} a {worst1[0]} kHz)', miss1 == 0 and m1 <= 0.01, 'H34(a)'))
    d0 = [(z - T1[None][n]) / T1[None][n] for n, z in enumerate(D['t1']['None'])]
    res.append((f'Sin axión frente a la fila sin axión de la Tabla 1: ' + ', '.join(f'{d * 100:+.2f} %' for d in d0),
                None, 'descriptivo: la prueba previa del cuaderno, ahora en el paso'))

    # (b) y (c) Tabla 2
    rows = [(nu_a, n, m) for nu_a in NUA for n, m in enumerate(D['t2'][str(nu_a)])]
    found = [(nu_a, n, m) for nu_a, n, m in rows if 'nu' in m]
    dre = [abs(m['nu'][0] - m['paper'][0]) / m['paper'][0] for *_, m in found]
    cls = all((m['paper'][1] == 0.0) == m['bound'] for *_, m in found)
    bound_real = all(m['im_det_rel'] < 1e-8 for *_, m in found if m['bound'])
    D['t2_maxdev_re'] = max(dre)
    res.append((f'Tabla 2, partes reales: {len(found)}/28 modos hallados; desviación máxima {max(dre) * 100:.2f} %; '
                f'los {sum(m["bound"] for *_, m in found)} modos con Im = 0 del paper son ligados también aquí '
                f'(frecuencia real)', len(found) == 28 and max(dre) <= 0.01 and cls and bound_real, 'H34(b)'))
    damped = [(nu_a, n, m) for nu_a, n, m in found if not m['bound']]
    ratios = [(nu_a, n, (m['nu'][1] * 1e3) / m['paper'][1]) for nu_a, n, m in damped]
    D['t2_im_ratio'] = [[nu_a, n, q] for nu_a, n, q in ratios]
    out25 = [(nu_a, n, q) for nu_a, n, q in ratios if abs(q - 1) > 0.25]
    res.append((f'Tabla 2, amortiguamientos: {len(ratios) - len(out25)}/{len(ratios)} dentro del 25 %; cociente '
                f'nuestro/paper entre {min(q for *_, q in ratios):.2f} y {max(q for *_, q in ratios):.2f}; fuera: '
                + (', '.join(f'f{n} a {nu_a} kHz ({q:.2f})' for nu_a, n, q in out25) or 'ninguno'),
                len(out25) <= 1, 'H34(c): criterio fijado antes de calcular, como mucho un modo fuera'))

    # (d) Tabla 3
    t3 = [(nu_a, k, m) for nu_a in NUA for k, m in enumerate(D['t3'][str(nu_a)])]
    conv = [(nu_a, k, m) for nu_a, k, m in t3 if m['nu']]
    near = [(nu_a, k, m) for nu_a, k, m in conv if abs(complex(*m['nu']) - complex(*m['paper'])) < 0.1 * abs(m['paper'][1])]
    dre3 = [100 * (m['nu'][0] / m['paper'][0] - 1) for *_, m in conv]
    rim3 = [m['nu'][1] / m['paper'][1] for *_, m in conv]
    D['t3_stats'] = {'n_conv': len(conv), 'n': len(t3), 're_pct': [min(dre3), max(dre3)], 'im_ratio': [min(rim3), max(rim3)],
                     'max_im_conv': max(abs(m['paper'][1]) for *_, m in conv),
                     'min_im_fail': min(abs(m['paper'][1]) for *_, m in t3 if not m['nu'])}
    res.append((f'Tabla 3 (modos del axión): Muller converge en {len(conv)}/{len(t3)} semillas y en {len(near)} queda '
                f'cerca del valor del paper (a menos de un 10 % de su |Im|). En las que convergen, Re ν se desvía entre '
                f'{min(dre3):+.2f} % y {max(dre3):+.2f} % e Im ν / paper va de {min(rim3):.3f} a {max(rim3):.3f}. Las que no convergen '
                f'tienen |Im| ≥ {D["t3_stats"]["min_im_fail"]:.2f} kHz en el paper; las que sí, ≤ {D["t3_stats"]["max_im_conv"]:.2f}', None, 'H34(d), descriptivo'))

    # Controles negativos
    neg1 = min(abs(D['t1']['None'][n] - T1[x][n]) / T1[x][n] for x in NUA if x >= 4 for n in range(4))
    res.append((f'Control negativo 1: la estrella sin axión se desvía de la Tabla 1 al menos un {neg1 * 100:.1f} % '
                f'para ν_a ≥ 4 kHz', neg1 > 0.05, ''))
    inc = [m['incoming'] for *_, m in damped]
    res.append((f'Control negativo 2 (onda entrante): en las {len(inc)} raíces amortiguadas el determinante deja de '
                f'anularse (|det|/escala ≥ {min(inc):.1e})', min(inc) > 1e-3, ''))

    D['checks'] = [{'name': t, 'ok': ok_field(o), 'note': n_} for t, o, n_ in res]
    (OUT / 'radial_paper1.json').write_text(json.dumps(D, indent=1, ensure_ascii=False, default=lambda o: o.item()))
    return report('Paso 26: modos radiales frente a las Tablas 1-3 del Paper I', res)


if __name__ == '__main__':
    run()
