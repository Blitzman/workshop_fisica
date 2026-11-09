"""Paso 20: recuento de modos en el semiplano amortiguado de la estrella realista (hipótesis H22).

¿Existe una familia física de resonancias axiónicas muy amortiguadas? Las que mostraba una versión anterior eran
artefactos del lapso (paso 15, H17), pero eso no prueba que no exista ninguna. Aquí se cuentan los ceros de
F(omega) = radial.det_entire(inward=False, sheet=1) en Re nu in [6.05, 12], Im nu in [-1.5, -0.002] kHz
(nu_a = 6 kHz): onda saliente (rama principal, Re k > 0, Im k < 0), semillas integradas hacia fuera (lo estable
cuando la solución exterior crece) y Psi normalizada en el infinito (analítica, sin polos). Cada índice se acepta
solo si no cambia al doblar el muestreo.
Control positivo: con la normalización antigua del lapso, las resonancias de cavidad (artefactos) deben aparecer
en el recuento.
"""
import json
import multiprocessing as mp
from pathlib import Path

import numpy as np

from . import axstar as AX
from . import radial as RD
from .pipeline import report

OUT = Path(__file__).resolve().parent.parent / 'results'
TILES = [(6.05, 9.0, -0.3, -0.002), (9.0, 12.0, -0.3, -0.002), (6.05, 9.0, -1.5, -0.3), (9.0, 12.0, -1.5, -0.3)]
# región profunda (H23): ahí está a0 = 9.0 - 3.8i kHz de la Tabla III de la secuela (l = 2, Cowling, nu_a = 6 kHz)
# franja casi real (H26): el texto decía Im nu < 0 y las piezas de arriba paraban en -0.002 kHz
TILES_THIN = [(6.05, 9.0, -0.002, -1e-5), (9.0, 12.0, -0.002, -1e-5)]
TILES_DEEP = [(6.05, 9.0, -5.0, -1.5), (9.0, 12.0, -5.0, -1.5)]
PAPER2_A0 = 9.0 - 3.8j
_S = {}


def _F(nu):
    return complex(RD.det_entire(_S['st'], AX.omega_of(nu), rm=_S.get('rm'), inward=False, sheet=1))


def _corners(b):
    x0, x1, y0, y1 = b
    return [complex(x0, y0), complex(x1, y0), complex(x1, y1), complex(x0, y1)]


def count(st, tiles=TILES):
    _S['st'] = st
    with mp.get_context('fork').Pool(6) as pool:
        w = [RD.winding_converged(_F, _corners(b), pool=pool) for b in tiles]
    return {'per_tile': [x[0] for x in w], 'tried': [x[4] for x in w], 'points': sum(x[2] for x in w),
            'max_dphase': float(max(x[3] for x in w)),
            'total': sum(x[0] for x in w) if all(x[0] is not None for x in w) else None}


def run():
    res, out = [], {}
    st = AX.Star(6.0)
    c = count(st)
    # los modos fluidos cuasi-normales de la región (semillas: H1 del paso 15 y sus sobretonos sin axión)
    found = []
    for seed in (8.119 - 0.0036j, 10.6 - 0.01j, 11.7 - 0.01j):
        try:
            z = RD.muller(_F, seed, tol=1e-9)
            if 6.05 < z.real < 12 and -1.5 < z.imag < -0.002 and all(abs(z - y) > 1e-3 for y in found):
                found.append(complex(z))
        except RuntimeError:
            pass
    fr = [RD.fluid_fraction(st, AX.omega_of(z)) for z in found]
    out['nu6'] = {**c, 'fluid_modes': [[z.real, z.imag] for z in found], 'fluid_fraction': fr}
    res.append((f'H22(a): nu_a = 6 kHz, índice {c["total"]} en Re nu in [6.05, 12], Im nu in [-1.5, -0.002] kHz '
                f'(piezas {c["per_tile"]}, convergidas {c["tried"]}); Muller encuentra {len(found)} modos fluidos ahí: '
                + ', '.join(f'{z.real:.4f}{z.imag * 1e3:+.2f}i (fluido {f:.2f})' for z, f in zip(found, fr)),
                c['total'] == len(found), 'sin ceros adicionales: ninguna familia axiónica en esa región'))
    # H26: franja casi real
    ct = count(st, TILES_THIN)
    out['thin'] = ct
    res.append((f'H26: franja casi real Re nu in [6.05, 12], Im nu in [-0.002, -1e-5] kHz: índice {ct["total"]} (piezas '
                f'{ct["per_tile"]}, convergidas {ct["tried"]}; salto de fase máximo {ct["max_dphase"]:.2f} rad)',
                ct['total'] == 0, 'fuera quedan |Im nu| < 0.01 Hz, Re nu < 6.05 kHz y Re nu > 12 kHz'))
    # H23: región profunda, donde la secuela tiene su modo axiónico a0 (l = 2, Cowling)
    cd = count(st, TILES_DEEP)
    deep = {**cd, 'roots': []}
    if cd['total']:
        rr = np.linspace(st.R, RD.AG.r_match(st), 4000)
        r_alt = float(rr[np.argmax(np.abs(st.sol.sol(rr)[3]) < 1e-5)])
        seeds = [PAPER2_A0] + [complex(x, y) for x in (7.0, 8.5, 10.0, 11.5) for y in (-2.0, -3.0, -4.0)]
        for sd in seeds:
            try:
                z = complex(RD.muller(_F, sd, tol=1e-9))
            except RuntimeError:
                continue
            if 6.05 < z.real < 12 and -5 < z.imag < -1.5 and all(abs(z - y['nu']) > 1e-3 for y in deep['roots']):
                _S['rm'] = r_alt
                try:
                    z2 = complex(RD.muller(_F, z, tol=1e-9))
                except RuntimeError:
                    z2 = complex(np.nan, np.nan)
                _S.pop('rm')
                deep['roots'].append({'nu': z, 'nu_alt_rm': z2, 'r_alt': r_alt, 'shift_pct': float(100 * abs(z2 / z - 1)),
                                      'fluid': RD.fluid_fraction(st, AX.omega_of(z))})
        deep['roots'] = [{**r_, 'nu': [r_['nu'].real, r_['nu'].imag], 'nu_alt_rm': [r_['nu_alt_rm'].real, r_['nu_alt_rm'].imag]}
                         for r_ in deep['roots']]
    out['deep'] = deep
    robust = [r_ for r_ in deep['roots'] if r_['shift_pct'] < 0.1]
    if not cd['total']:
        txt = 'ningún cero: no hay análogo radial de a0 en esa región'
    elif not deep['roots']:
        txt = f'ningún cero localizado: Muller no converge desde ninguna de las {len(seeds)} semillas'
    else:
        txt = ('ceros: ' + '; '.join(f"{r_['nu'][0]:.3f}{r_['nu'][1]:+.3f}i kHz (fluido {r_['fluid']:.2f}; al cambiar el "
                                     f"empalme se mueve un {r_['shift_pct']:.2g} %)" for r_ in deep['roots']))
    res.append((f'H23: región profunda Re nu in [6.05, 12], Im nu in [-5, -1.5] kHz: índice {cd["total"]} (piezas '
                f'{cd["per_tile"]}, convergidas {cd["tried"]}); {txt}. Secuela, Tabla III (l = 2, Cowling): a0 = '
                f'{PAPER2_A0.real} {PAPER2_A0.imag:+.1f}i kHz', cd['total'] is not None and
                (cd['total'] == 0 or len(robust) == cd['total']), 'el resultado decide; no había predicción fuerte'))
    # control positivo: normalización antigua del lapso (Phi ajustado en R con m(R), sin la cola)
    m_R, Phi_R = st.sol.sol(st.R)[0], st.sol.sol(st.R)[2]
    Ps_new = st.Phi_shift
    st._Ps = 0.5 * np.log(1 - 2 * m_R / st.R) - Phi_R
    M_tot = st.sol.y[0, -1]
    out['lapse'] = {'phi_inf_old': float(st._Ps - Ps_new), 'tail_mass_frac': float((M_tot - m_R) / M_tot)}
    co = count(st)
    out['old_lapse'] = co
    out['old_lapse_deep'] = count(st, TILES_DEEP)
    res.append((f'H22(b), control positivo: con la normalización antigua del lapso el índice sube a {co["total"]} '
                f'(piezas {co["per_tile"]}): el recuento ve las resonancias de cavidad', co['total'] is not None
                and c['total'] is not None and co['total'] >= c['total'] + 3, ''))
    out['checks'] = [{'name': n_, 'ok': bool(o_), 'note': t_} for n_, o_, t_ in res]
    OUT.mkdir(exist_ok=True)
    (OUT / 'radial_damped.json').write_text(json.dumps(out, indent=1))
    return report('Paso 20: recuento en el semiplano amortiguado', res)


if __name__ == '__main__':
    run()
