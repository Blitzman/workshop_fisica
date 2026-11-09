"""Paso 22: la condición exterior no depende del rayo complejo ni del punto de partida (hipótesis H24).

radial.exterior_psi impone la onda saliente e^{i(kr + eta ln r)}/r (con la fase coulombiana eta) en
r_m + S e^{i phi}, |k| S = n_decay, y la integra hacia dentro por el rayo phi = pi/2 - arg k, donde la onda
entrante se atenúa hacia dentro. Se mide L = psi'/psi en r_m:
  (a) girando el rayo dphi = +-0.4 rad;
  (b) alejando el punto de partida (n_decay 40 -> 60).
Controles negativos: (1) onda entrante (sheet = -1); (2) n_decay = 1, donde el truncado del desarrollo asintótico
ya no se suprime.
"""
import json
from pathlib import Path

import numpy as np

from . import axgr as AG
from . import axstar as AX
from . import radial as RD
from .pipeline import ok_field, report

OUT = Path(__file__).resolve().parent.parent / 'results'
# H1 (paso 20), el antiguo artefacto (paso 15), un cero profundo (paso 21) y el nivel ligado n = 2 (paso 18), en kHz
POINTS = {'H1': 8.1192 - 0.0036j, 'artefacto antiguo': 6.585 - 0.563j, 'cero profundo': 6.293 - 4.480j,
          'nivel ligado n = 2': 5.944070 + 0j}


def run():
    res, out = [], {'points': {}}
    st = AX.Star(6.0)
    rm, M = float(AG.r_match(st)), float(st.sol.y[0, -1])

    def L(nu, **kw):
        kw = {'sheet': 1, 'n_decay': 40.0, **kw}
        ps, dps = RD.exterior_psi(AX.omega_of(nu), M, rm, st.ma, np.sqrt(st.fa2), st.E0, rtol=1e-12, **kw)
        return complex(dps / ps)

    for name, nu in POINTS.items():
        L0 = L(nu)
        rel = lambda z: float(abs(z / L0 - 1))
        try:
            ing = rel(L(nu, sheet=-1))
        except RuntimeError:
            # bajo el umbral la «onda entrante» es la solución creciente e^{kappa r}: la integración no termina
            # (antes se usaba el resultado a medias sin mirar status). Este punto no entra en el control negativo 1.
            ing = None
        out['points'][name] = {'nu': [nu.real, nu.imag], 'L': [L0.real, L0.imag],
                               'dphi+0.4': rel(L(nu, dphi=0.4)), 'dphi-0.4': rel(L(nu, dphi=-0.4)),
                               'n_decay60': rel(L(nu, n_decay=60.0)), 'ingoing': ing,
                               'n_decay1': rel(L(nu, n_decay=1.0))}
    P = out['points']
    ang = max(max(p['dphi+0.4'], p['dphi-0.4']) for p in P.values())
    nd = max(p['n_decay60'] for p in P.values())
    res.append((f'H24(a): girar el rayo complejo ±0.4 rad cambia L = psi\'/psi en r_m como mucho {ang:.1e} (relativo) en '
                + ', '.join(P) + ' ' + '; '.join(f"{k}: {max(p['dphi+0.4'], p['dphi-0.4']):.0e}" for k, p in P.items()),
                ang < 1e-8, 'la condición saliente no depende del camino en el plano r complejo'))
    res.append((f'H24(b): alejar el punto de partida (n_decay 40 → 60) cambia L como mucho {nd:.1e}: '
                + '; '.join(f"{k}: {p['n_decay60']:.0e}" for k, p in P.items()), nd < 1e-8,
                'el truncado del desarrollo asintótico (con la fase coulombiana) no influye'))
    ing = min(p['ingoing'] for k, p in P.items() if k != 'nivel ligado n = 2')
    bound = P['nivel ligado n = 2']['ingoing']
    res.append((f'Control negativo 1: con la onda entrante, L cambia al menos un {ing:.2g} (relativo) en los puntos '
                f'por encima del umbral', ing > 1e-2,
                'en el nivel ligado la mutación es la solución creciente: ' +
                ('la integración no termina y se detecta' if bound is None else f'L cambia {bound:.1e}')))
    n1 = min(p['n_decay1'] for p in P.values())
    res.append((f'Control negativo 2: con n_decay = 1 el cambio mínimo es {n1:.1e}, frente a {max(ang, nd):.1e} en (a)-(b): '
                'la prueba detecta un punto de partida demasiado cercano', n1 > 100 * max(ang, nd), ''))
    out['checks'] = [{'name': n_, 'ok': ok_field(o_), 'note': t_} for n_, o_, t_ in res]
    OUT.mkdir(exist_ok=True)
    (OUT / 'radial_contour.json').write_text(json.dumps(out, indent=1))
    return report('Paso 22: condición exterior frente al rayo complejo', res)


if __name__ == '__main__':
    run()
