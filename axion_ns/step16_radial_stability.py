"""Paso 16: estabilidad radial de la estrella con condensado por el principio del argumento (hipótesis H18).

Responde a la crítica A de las revisiones externas: el barrido de signo en el eje imaginario (paso 15) no excluye
raíces fuera del eje ni pares próximos. Aquí se cuentan los ceros de F(omega) = radial.det_entire en el semiplano
inestable Im omega > 0, donde la solución exterior decae (rama Im k > 0) y F es analítica y sin polos.

- Cada índice se acepta solo si no cambia al doblar el muestreo (radial.winding_converged): cerca del umbral la
  escalera de estados ligados del átomo gravitatorio, que se acumula en nu_a, hace girar la fase más deprisa de lo
  que resuelve un muestreo fijo (a 3,4 kHz una pieza daba 1 con 48 puntos por lado y 0 con 96).
- Contorno: piezas |Re nu| <= 12 kHz, eps <= Im nu <= 12 kHz. Una franja central (|Re nu| <= 2) que contiene el eje
  y dos piezas a la derecha; las de la izquierda valen lo mismo por la simetría F(-omega*) = F(omega)*, que se
  comprueba. Un único rectángulo grande pierde vueltas: su borde superior es largo y la fase gira deprisa
  (lo vimos: 2, 4 u 8 «ceros» según el muestreo, que desaparecen al trocear).
- En el eje imaginario F es real: se cuentan además sus cambios de signo hasta 20 kHz.
- Controles positivos: estrellas más allá de la masa máxima (sin axión y con axión) deben dar índice 1.
- Controles negativos: sin quitar el factor e^{i k r_m} el contorno pierde vueltas, y con la rama principal de k
  para Re omega < 0 se rompe la simetría; las comprobaciones lo detectan.
"""
import json
import multiprocessing as mp
from pathlib import Path

import numpy as np

from . import axstar as AX
from . import radial as RD
from .pipeline import report

OUT = Path(__file__).resolve().parent.parent / 'results'
EPS, YMAX, XC, XMAX, GMAX = 0.05, 12.0, 2.0, 12.0, 20.0
TILES_C = [(-XC, XC, EPS, 4.0), (-XC, XC, 4.0, YMAX)]
TILES_R = [(XC, XMAX, EPS, 4.0), (XC, XMAX, 4.0, YMAX)]
_S = {}


def _F(nu):
    st, kw = _S['st'], _S['kw']
    return complex(RD.det_entire(st, AX.omega_of(nu), **kw))


def _g(nu):                              # sin axión: Delta p(R) de la semilla fluida, entera en omega
    return complex(RD.interior(_S['st'], AX.omega_of(nu), 'fluid').y[1, -1])


def _corners(b):
    x0, x1, y0, y1 = b
    return [complex(x0, y0), complex(x1, y0), complex(x1, y1), complex(x0, y1)]


def count(st, f=_F, n0=48, eps=EPS, kw=None):
    """Índice total en |Re nu| <= XMAX, eps <= Im nu <= YMAX (franja central + 2 x derecha) y cambios de signo de F
    en el eje imaginario. Devuelve un dict con todo."""
    _S.update(st=st, kw=kw or {})
    tc = [(a, b, max(c, eps) if c == EPS else c, d) for a, b, c, d in TILES_C]
    tr = [(a, b, max(c, eps) if c == EPS else c, d) for a, b, c, d in TILES_R]
    with mp.get_context('fork').Pool(6) as pool:
        wc = [RD.winding_converged(f, _corners(b), n0=n0, pool=pool) for b in tc]
        wr = [RD.winding_converged(f, _corners(b), n0=n0, pool=pool) for b in tr]
        gs = np.concatenate([np.geomspace(0.01, 1, 60), np.linspace(1, GMAX, 400)[1:]])
        v = np.array(pool.map(f, list(1j * gs)))
    conv = all(w[0] is not None for w in wc + wr)
    total = sum(w[0] for w in wc) + 2 * sum(w[0] for w in wr) if conv else None
    sc = [float(gs[i]) for i in range(len(gs) - 1) if v[i].real * v[i + 1].real < 0]
    return {'total': total, 'converged': conv, 'central': [w[0] for w in wc], 'right': [w[0] for w in wr],
            'tried': [w[4] for w in wc + wr], 'n0_final': [w[1] for w in wc + wr],
            'points': sum(w[2] for w in wc + wr), 'max_dphase': max(w[3] for w in wc + wr),
            'axis_sign_changes': sc, 'axis_imag_rel': float(np.max(np.abs(v.imag) / np.abs(v)))}


def symmetry(st, sheet=0):
    _S.update(st=st, kw={'sheet': sheet})
    pts = [2 + 1j, 5 + 0.3j, 9 + 6j, 11 + 11j]
    return max(abs(_F(-z.conjugate()) - np.conj(_F(z))) / abs(_F(z)) for z in pts)


def run():
    res, out = [], {}
    # 1. Controles positivos sin axión: F = Delta p(R), el punto de retorno está en p_c ~ 905 MeV/fm^3 (paso 15)
    ctrl = {}
    for p0 in (1000.0, 800.0):
        _S.update(st=AX.Star(None, p0=p0))
        with mp.get_context('fork').Pool(6) as pool:
            ctrl[p0] = RD.winding_converged(_g, _corners((-XMAX, XMAX, EPS, YMAX)), pool=pool)[0]
    res.append((f'Control positivo sin axión: índice {ctrl[1000.0]} más allá de la masa máxima (p_c = 1000) y '
                f'{ctrl[800.0]} antes (p_c = 800)', ctrl[1000.0] == 1 and ctrl[800.0] == 0, 'H18(b)'))
    out['control_no_axion'] = {str(k): v for k, v in ctrl.items()}

    # 2. Control positivo con axión: estrella con condensado más allá de su masa máxima
    stc = AX.Star(6.0, p0=1000.0)
    cc = count(stc)
    res.append((f'Control positivo con axión (nu_a = 6 kHz, p_c = 1000): índice {cc["total"]}; cambio de signo en el '
                f'eje imaginario en gamma = {", ".join(f"{g:.2f}" for g in cc["axis_sign_changes"])} kHz',
                cc['total'] == 1 and len(cc['axis_sign_changes']) == 1, 'H18(c)'))
    out['control_axion'] = cc

    # 3. Las estrellas del paso 15
    for nu_a in (6.0, 3.4):
        st = AX.Star(nu_a)
        sym = symmetry(st)
        c = count(st)
        out[f'nu{nu_a}'] = {**c, 'symmetry': sym}
        ok = c['total'] == 0 and not c['axis_sign_changes'] and c['axis_imag_rel'] < 1e-10 and sym < 1e-10
        res.append((f'nu_a = {nu_a} kHz: ningún cero de F con |Re nu| <= {XMAX:.0f}, {EPS} <= Im nu <= {YMAX:.0f} kHz '
                    f'(índice {c["total"]}: central {c["central"]}, derecha x2 {c["right"]}, cada pieza estable al doblar '
                    f'el muestreo, {c["tried"]}; {c["points"]} puntos) y '
                    f'ningún cambio de signo en el eje hasta {GMAX:.0f} kHz; simetría F(-w*) = F(w)* a {sym:.0e}',
                    ok, 'H18(a)'))
        if nu_a == 6.0:
            # robustez: doble muestreo y eps/2 en la pieza central baja (la más delicada)
            _S.update(st=st, kw={})
            with mp.get_context('fork').Pool(6) as pool:
                w2 = RD.winding(_F, _corners(TILES_C[0]), n0=192, pool=pool)[0]
                w3 = RD.winding_converged(_F, _corners((-XC, XC, EPS / 2, 4.0)), pool=pool)[0]
                # control negativo 1: sin quitar e^{i k r_m}, muestreo grueso en una pieza ancha
                _S.update(kw={'raw': 'inf'})
                wbad = RD.winding(_F, _corners((0.0, XMAX, EPS, 4.0)), n0=12, pool=pool)[0]
            out['robust6'] = {'n0_96': w2, 'eps_half': w3, 'neg_raw_inf': wbad}
            res.append((f'Robustez (6 kHz, franja central baja): índice {round(w2)} con 192 puntos por lado y {w3} con eps/2',
                        round(w2) == c['central'][0] and w3 == c['central'][0], 'H18(d)'))
            res.append((f'Control negativo 1: sin quitar el factor e^(i k r_m) y con muestreo grueso, el contorno pierde '
                        f'vueltas (índice {round(wbad)} en vez de 0)', round(wbad) != 0, ''))
            symbad = symmetry(st, sheet=1)
            out['neg_sheet_principal'] = symbad
            res.append((f'Control negativo 2: con la rama principal de k también para Re omega < 0, la simetría '
                        f'F(-w*) = F(w)* se rompe ({symbad:.1e} frente a {sym:.0e})', symbad > 1e-3, ''))
    OUT.mkdir(exist_ok=True)
    out['region'] = {'eps': EPS, 'ymax': YMAX, 'xmax': XMAX, 'gmax': GMAX}
    out['checks'] = [{'name': n_, 'ok': bool(o_), 'note': t_} for n_, o_, t_ in res]
    (OUT / 'radial_stability.json').write_text(json.dumps(out, indent=1))
    return report('Paso 16: estabilidad radial por el principio del argumento', res)


if __name__ == '__main__':
    run()
