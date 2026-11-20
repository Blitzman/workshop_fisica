"""Paso 25: rama del exterior y convergencia del paso 10 en sus modos más delicados (revisión del informe principal).

(a) H31. `axstar.exterior_L` arrancaba con la rama principal de k = sqrt(w^2 - m_a^2). Por debajo del umbral, con el modo
    amortiguado, esa rama tiene Im k < 0: la solución que crece hacia fuera. Ahora se toma la que decae. Se comprueba
    que L = psi'/psi no cambiaba: la integración hacia dentro recorre ~40 longitudes de decaimiento y filtra la creciente.
    Controles negativos: (1) arrancando a solo 2 longitudes de decaimiento la rama sí cambia L; (2) por encima del
    umbral, la rama opuesta (onda entrante) cambia L.
(b) H33. La prueba de convergencia del paso 10 solo miraba el modo f a 6 kHz. Aquí se repite (tolerancia, radio de
    empalme, distancia de extracción) en los modos más exigentes: los más cercanos al umbral por debajo y por encima,
    el de menor amortiguamiento que radia axiones y el modo f con el axión más pesado. Controles negativos: empalmar
    en la superficie (sin la cola del campo) y extraer la onda a una sola longitud de onda.
"""
import json
import multiprocessing as mp
from pathlib import Path

import numpy as np

from . import axgr as AG
from . import axstar as AX
from . import step10_axion_gr as S10
from .pipeline import ok_field, report

OUT = Path(__file__).resolve().parent.parent / 'results'
L2 = 2
# (a) frecuencias de prueba (kHz) para nu_a = 6 kHz: por debajo del umbral con distinto amortiguamiento, y por encima
BELOW = (2.102 - 0.001j, 2.4 - 0.01j, 5.0 - 0.1j, 5.9 - 0.5j)
ABOVE = (7.4 - 0.001j, 11.2 - 0.01j)
# (b) modos del paso 10: (nombre, índice del modo, nu_a)
DELICATE = (('p1, justo por debajo del umbral (menor amortiguamiento)', 1, 8.0),
            ('p2, justo por debajo del umbral', 2, 12.0),
            ('p1, justo por encima del umbral', 1, 7.0),
            ('p2, justo por encima del umbral', 2, 11.0),
            ('p2, el que radia axiones con menor amortiguamiento', 2, 8.0),
            ('f, axión más pesado', 0, 15.0))
VARIANTS = {'base': {}, 'rtol': {'rtol': 1e-11}, 'rm': 'rm08', 'n_far': {'n_far': 80}}
NEG_MODE = 2                                          # índice en DELICATE para los controles negativos (p1, 7 kHz)
NEGS = {'rm=R': 'rmR', 'n_far=1': {'n_far': 1}}


def branch_test():
    st = AX.Star(6.0)
    rm = AG.r_match(st)
    L = lambda nu, **kw: complex(AX.exterior_L(st, L2, AX.omega_of(nu), r_stop=rm, **kw))
    rel = lambda a, b: abs(a / b - 1)
    out = {'below': [], 'above': []}
    for nu in BELOW:
        k = np.sqrt(complex(AX.omega_of(nu)**2 - st.ma**2))
        out['below'].append({'nu': [nu.real, nu.imag], 'k_principal': [k.real, k.imag],
                             'rel': rel(L(nu, branch='principal'), L(nu)),
                             'rel_short': rel(L(nu, branch='principal', n_decay=2), L(nu))})
    for nu in ABOVE:
        out['above'].append({'nu': [nu.real, nu.imag], 'rel': rel(L(nu, branch='principal'), L(nu)),
                             'rel_flipped': rel(L(nu, branch='flipped'), L(nu))})
    return out


def _guess(i, nu_a):
    d = json.loads((OUT / 'fig_axgr.json').read_text())
    j = d['nus'].index(nu_a)
    m = d['modes'][i]
    return complex(m['re'][j], -m['im'][j] * 1e-3)


def _task(args):
    i, nu_a, kw = args
    st = AX.Star(nu_a)
    z0 = _guess(i, nu_a)
    if kw == 'rm08':
        kw = {'rm': 0.8 * AG.r_match(st)}
    elif kw == 'rmR':
        kw = {'rm': float(st.R)}
    try:
        z = S10._mode(st, z0.real, dyn=True, **kw)    # el mismo arranque y la misma tolerancia que el paso 10
    except RuntimeError:
        z = complex(np.nan, np.nan)
    return z


def run():
    res, out = [], {}
    # (a) rama
    B = branch_test()
    out['branch'] = B
    worst = max(b['rel'] for b in B['below'])
    res.append((f'H31: por debajo del umbral (nu_a = 6 kHz) la rama principal de k tiene Im k < 0 en '
                f'{sum(b["k_principal"][1] < 0 for b in B["below"])}/{len(B["below"])} frecuencias, pero L = psi\'/psi coincide '
                f'con la de la rama que decae a {worst:.0e} relativo como mucho: '
                + ', '.join(f'{b["nu"][0]:.3g}{b["nu"][1]:+.3g}i kHz: {b["rel"]:.0e}' for b in B['below']),
                worst < 1e-10, 'la integración hacia dentro (40 longitudes de decaimiento) filtra la solución creciente'))
    above = max(b['rel'] for b in B['above'])
    res.append((f'H31: por encima del umbral la rama no cambia (diferencia {above:.0e})', above == 0, ''))
    short = min(b['rel_short'] for b in B['below'])
    res.append((f'Control negativo 1: arrancando a 2 longitudes de decaimiento, la rama principal cambia L al menos un '
                f'{short:.0e} relativo (frente a {worst:.0e} con 40)', short > 1e-4, ''))
    flip = min(b['rel_flipped'] for b in B['above'])
    res.append((f'Control negativo 2: por encima del umbral, la onda entrante cambia L al menos un {flip:.1g} relativo',
                flip > 1e-2, ''))
    # (b) convergencia de los modos delicados
    tasks = [(i, nu_a, kw) for _, i, nu_a in DELICATE for kw in VARIANTS.values()]
    _, ni, nn = DELICATE[NEG_MODE]
    tasks += [(ni, nn, kw) for kw in NEGS.values()]
    with mp.get_context('fork').Pool(6) as pool:
        Z = pool.map(_task, tasks)
    nv = len(VARIANTS)
    rows = []
    for k, (name, i, nu_a) in enumerate(DELICATE):
        zs = dict(zip(VARIANTS, Z[k * nv:(k + 1) * nv]))
        z0 = zs['base']
        rows.append({'name': name, 'mode': ['f', 'p1', 'p2'][i], 'nu_a': nu_a, 'nu': [z0.real, z0.imag],
                     'dre_rel': {v: abs(zs[v].real / z0.real - 1) for v in VARIANTS if v != 'base'},
                     'dim_Hz': {v: abs(zs[v].imag - z0.imag) * 1e3 for v in VARIANTS if v != 'base'}})
    neg = dict(zip(NEGS, Z[len(DELICATE) * nv:]))
    out['delicate'] = rows
    # tolerancia de la prueba: 1e-7 relativo en Re; en Im, 1 % o la resolución de Muller (1e-9 |z|), lo que sea mayor
    res_im = lambda r: max(0.01 * abs(r['nu'][1]) * 1e3, 1e-9 * abs(complex(*r['nu'])) * 1e3)
    finite = all(np.isfinite(z) for z in Z[:len(DELICATE) * nv])
    ok_re = finite and all(max(r['dre_rel'].values()) < 1e-7 for r in rows)
    ok_im = finite and all(max(r['dim_Hz'].values()) < res_im(r) for r in rows)
    res.append((f'H33: en los {len(rows)} modos delicados del paso 10, tolerancia, radio de empalme y distancia de extracción '
                f'mueven Re nu como mucho {max(max(r["dre_rel"].values()) for r in rows):.0e} relativo: '
                + '; '.join(f'{r["mode"]} a {r["nu_a"]:g} kHz {max(r["dre_rel"].values()):.0e}' for r in rows), ok_re, ''))
    res.append(('H33: e Im nu como mucho: ' + '; '.join(
        f'{r["mode"]} a {r["nu_a"]:g} kHz (Im {-r["nu"][1] * 1e3:.3g} Hz) {max(r["dim_Hz"].values()):.1e} Hz' for r in rows),
        ok_im, 'criterio: 1 % de Im nu o la resolución de Muller (1e-9 |nu|), lo que sea mayor'))
    r0 = rows[NEG_MODE]
    z0 = complex(*r0['nu'])
    dneg = {k: (abs(v.imag - z0.imag) * 1e3 if np.isfinite(v) else np.inf) for k, v in neg.items()}
    out['neg'] = {'mode': r0['name'], 'nu_a': r0['nu_a'], 'dim_Hz': dneg}
    base = max(r0['dim_Hz'].values())
    res.append((f'Control negativo 3: empalmando en la superficie (sin la cola del campo), Im nu de {r0["mode"]} a '
                f'{r0["nu_a"]:g} kHz se mueve {dneg["rm=R"]:.2g} Hz, frente a {base:.1e} Hz con las variantes de arriba',
                dneg['rm=R'] > 10 * max(base, res_im(r0)), ''))
    res.append((f'Control negativo 4: extrayendo la onda a una sola longitud de onda se mueve {dneg["n_far=1"]:.2g} Hz',
                dneg['n_far=1'] > 10 * max(base, res_im(r0)), ''))
    out['checks'] = [{'name': n_, 'ok': ok_field(o_), 'note': t_} for n_, o_, t_ in res]
    OUT.mkdir(exist_ok=True)
    (OUT / 'exterior_convergence.json').write_text(json.dumps(out, indent=1, default=str))
    return report('Paso 25: rama del exterior y convergencia de los modos delicados', res)


if __name__ == '__main__':
    run()
