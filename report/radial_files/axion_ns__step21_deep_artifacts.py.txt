"""Paso 21: ¿son modos los ceros profundos del paso 20? (hipótesis H23b).

El paso 20 encuentra dos ceros de F en Re nu in [6.05, 9], Im nu in [-5, -1.5] kHz (nu_a = 6 kHz). No tienen
fluido y se mueven mucho al cambiar el radio de empalme. Aquí se comprueba si son artefactos del empalme:
  (i) control: los modos fluidos H1 y H2 deben quedarse quietos (< 0.1 %) con el mismo cambio de radio;
  (ii) los ceros profundos deben moverse > 5 % también con un tercer radio (donde theta_0 = 1e-6).
Controles negativos:
  1. con onda entrante (sheet = -1), H1 deja de ser un cero;
  2. empalmando en la superficie (r_m = R, sin la cola del campo), la prueba debe notar que se ha quitado algo:
     H1 se mueve más que con el cambio de radio de (i).
"""
import json
from pathlib import Path

import numpy as np

from . import axgr as AG
from . import axstar as AX
from . import radial as RD
from .pipeline import report

OUT = Path(__file__).resolve().parent.parent / 'results'
_S = {}


def _F(nu):
    return complex(RD.det_entire(_S['st'], AX.omega_of(nu), rm=_S.get('rm'), inward=False, sheet=_S.get('sheet', 1)))


def _root(z, rm=None):
    _S['rm'] = rm
    try:
        return complex(RD.muller(_F, z, tol=1e-9))
    except RuntimeError:
        return complex(np.nan, np.nan)
    finally:
        _S.pop('rm')


def _radius(st, th):
    rr = np.linspace(st.R, AG.r_match(st), 4000)
    return float(rr[np.argmax(np.abs(st.sol.sol(rr)[3]) < th)])


def run():
    res, out = [], {}
    D = json.loads((OUT / 'radial_damped.json').read_text())
    st = AX.Star(6.0)
    _S['st'] = st
    r0, r5, r6 = float(AG.r_match(st)), _radius(st, 1e-5), _radius(st, 1e-6)
    out['radii'] = {'default': r0, 'theta_1e-5': r5, 'theta_1e-6': r6, 'R': float(st.R)}
    # (i) control: modos fluidos
    fl = []
    for z in (complex(*v) for v in D['nu6']['fluid_modes']):
        z0 = _root(z)
        z5, z6, zR = _root(z0, r5), _root(z0, r6), _root(z0, float(st.R))
        fl.append({'nu': [z0.real, z0.imag], 'shift_pct': {k: float(100 * abs(w / z0 - 1)) for k, w in
                   (('theta_1e-5', z5), ('theta_1e-6', z6), ('R', zR))},
                   'dIm_Hz': {k: float(1e3 * (w - z0).imag) for k, w in (('theta_1e-5', z5), ('theta_1e-6', z6), ('R', zR))}})
    # (ii) ceros profundos, sin duplicados (Muller converge a puntos a 0.002-0.004 kHz desde semillas distintas)
    deep = []
    for r_ in D['deep']['roots']:
        z = complex(*r_['nu'])
        if all(abs(z - complex(*d['nu'])) > 0.05 for d in deep):
            z0 = _root(z)
            z5, z6 = _root(z0, r5), _root(z0, r6)
            deep.append({'nu': [z0.real, z0.imag], 'nu_1e-5': [z5.real, z5.imag], 'nu_1e-6': [z6.real, z6.imag],
                         'shift_pct': {'theta_1e-5': float(100 * abs(z5 / z0 - 1)), 'theta_1e-6': float(100 * abs(z6 / z0 - 1))},
                         'fluid': r_['fluid']})
    out['fluid'], out['deep'] = fl, deep
    ok_i = all(max(f['shift_pct']['theta_1e-5'], f['shift_pct']['theta_1e-6']) < 0.1 for f in fl)
    res.append(('H23b(i), control: los modos fluidos casi no se mueven al cambiar el radio de empalme: ' + '; '.join(
        f"(Re {f['nu'][0]:.4f} kHz, Im {f['nu'][1] * 1e3:.2f} Hz): {f['shift_pct']['theta_1e-5']:.2g} % y {f['shift_pct']['theta_1e-6']:.2g} % "
        f"(Im: {f['dIm_Hz']['theta_1e-5']:+.2g} y {f['dIm_Hz']['theta_1e-6']:+.2g} Hz)" for f in fl), ok_i, ''))
    ok_ii = bool(deep) and all(min(d['shift_pct'].values()) > 5 for d in deep)
    res.append(('H23b(ii): ' + ('ningún cero profundo que mover: Muller no converge en la región (paso 20); la prueba no '
                                'decide' if not deep else 'los ceros profundos se mueven con los dos radios alternativos: ') + '; '.join(
        f"{d['nu'][0]:.3f}{d['nu'][1]:+.3f}i kHz -> {d['nu_1e-5'][0]:.3f}{d['nu_1e-5'][1]:+.3f}i ({d['shift_pct']['theta_1e-5']:.2g} %) "
        f"y {d['nu_1e-6'][0]:.3f}{d['nu_1e-6'][1]:+.3f}i ({d['shift_pct']['theta_1e-6']:.2g} %)" for d in deep),
        ok_ii, 'si (i) y (ii) pasan: artefactos del empalme, no modos'))
    # control negativo 1: onda entrante
    h1 = complex(*fl[0]['nu'])
    _S['sheet'] = -1
    ratio_in = abs(_F(h1)) / abs(_F(h1 + 0.05))
    _S['sheet'] = 1
    ratio_out = abs(_F(h1)) / abs(_F(h1 + 0.05))
    out['neg1'] = {'ratio_ingoing': float(ratio_in), 'ratio_outgoing': float(ratio_out)}
    res.append((f'Control negativo 1: con onda entrante, H1 deja de ser cero (|F(H1)|/|F(H1 + 50 Hz)| = {ratio_in:.2g}, '
                f'frente a {ratio_out:.1e} con la saliente)', ratio_in > 1e-2 and ratio_out < 1e-4, ''))
    # control negativo 2: empalme en la superficie
    sR = fl[0]['shift_pct']['R']
    s_alt = max(fl[0]['shift_pct']['theta_1e-5'], fl[0]['shift_pct']['theta_1e-6'])
    out['neg2'] = {'shift_R_pct': sR, 'shift_alt_pct': s_alt}
    res.append((f'Control negativo 2: empalmando en la superficie (sin la cola del campo), H1 se mueve un {sR:.2g} %, '
                f'más que con los radios alternativos ({s_alt:.2g} %)', sR > 3 * s_alt, ''))
    out['checks'] = [{'name': n_, 'ok': bool(o_), 'note': t_} for n_, o_, t_ in res]
    OUT.mkdir(exist_ok=True)
    (OUT / 'radial_deep.json').write_text(json.dumps(out, indent=1))
    return report('Paso 21: ceros profundos frente al radio de empalme', res)


if __name__ == '__main__':
    run()
