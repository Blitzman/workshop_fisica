"""Paso 23: sensibilidad a la inconsistencia termodinámica de los ajustes BSk26 (hipótesis H27).

Los ajustes publicados eps(n) (BSkEofN) y P(rho) (BSKfit) son independientes y no cumplen exactamente
d eps/dn = (eps + p)/n. El código toma c_s^2 = dp/deps sobre el ajuste y de ahí delta n/n = delta p/(c_s^2 (eps + p)),
que supone la identidad. Aquí se repiten F y H1 (nu_a = 6 kHz) y F (3,4 kHz) con la otra definición,
c_s^2 = n (dp/dn)/(eps + p), y se mide cuánto se mueven.
Controles negativos: (1) la identidad evaluada con el ajuste de P desplazado un 1 % en lg rho debe empeorar mucho;
(2) aplicar el cambio de c_s^2 dos veces (la media geométrica de las dos definiciones) debe dar un desplazamiento
intermedio: si no, el desplazamiento medido no viene de c_s^2.
"""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from . import axgr as AG
from . import axstar as AX
from . import bsk26 as B
from . import radial as RD
from .pipeline import report

OUT = Path(__file__).resolve().parent.parent / 'results'
E = AX._eos
_cs2_eps = E.cs2


def _cs2_n(p):          # n dp/dn/(eps + p) = p/(eps + p) / (d ln n / d ln p)
    return p / (E.eps(p) + p) / E._ln(np.log(p), 1)


def _cs2_mid(p):
    return np.sqrt(_cs2_eps(p) * _cs2_n(p))


def mismatch(p, eos=E):
    lp = np.log(p)
    e = eos.eps(p)
    return float(eos._le(lp, 1) / eos._ln(lp, 1) * e / (e + p) - 1)


def modes(R0):
    st6, st34 = AX.Star(6.0), AX.Star(3.4)
    rm6, rm34 = AG.r_match(st6), AG.r_match(st34)
    F6 = brentq(lambda nu: RD.det(st6, AX.omega_of(nu), rm6).real, R0['F6'] - 0.01, R0['F6'] + 0.01, xtol=1e-10)
    H1 = complex(RD.muller(lambda nu: RD.det(st6, AX.omega_of(nu), rm6), R0['H1']))
    F34 = complex(RD.muller(lambda nu: RD.det(st34, AX.omega_of(nu), rm34), R0['F34']))
    return {'F6': F6, 'H1': H1, 'F34': F34, 'pc6': float(st6.sol.y[1, 0] / AX.KGEO)}


def run():
    res, out = [], {}
    S = json.loads((OUT / 'radial_star.json').read_text())
    R0 = {'F6': S['nu6']['F'], 'H1': complex(*S['nu6']['H1']), 'F34': complex(*S['nu34']['F'])}
    out['mismatch'] = {str(p): mismatch(p) for p in (1, 10, 100, 300, 500, 900)}
    base = modes(R0)
    out['pc6_MeVfm3'] = base['pc6']
    out['mismatch_center'] = mismatch(base['pc6'])
    E.cs2 = _cs2_n
    alt = modes(R0)
    E.cs2 = _cs2_mid
    mid = modes(R0)
    E.cs2 = _cs2_eps
    cmp_ = lambda a, b: {'dRe_rel': float(abs((a.real - b.real) / b.real)),
                         'dIm_rel': float(abs((a.imag - b.imag) / b.imag)) if b.imag else 0.0,
                         'dRe_Hz': float(1e3 * (a.real - b.real)), 'dIm_Hz': float(1e3 * (a.imag - b.imag))}
    out['base'] = {k: [complex(v).real, complex(v).imag] for k, v in base.items() if k != 'pc6'}
    out['alt'] = {k: [complex(v).real, complex(v).imag] for k, v in alt.items() if k != 'pc6'}
    out['shift'] = {k: cmp_(complex(alt[k]), complex(base[k])) for k in ('F6', 'H1', 'F34')}
    out['shift_mid'] = {k: cmp_(complex(mid[k]), complex(base[k])) for k in ('F6', 'H1', 'F34')}
    sh = out['shift']
    res.append((f'H27: inconsistencia de los ajustes BSk26, n deps/dn/(eps + p) - 1 = '
                + ', '.join(f'{v:+.1e} (p = {k})' for k, v in out['mismatch'].items())
                + f' MeV/fm3; en el centro de la estrella de 6 kHz (p_c = {base["pc6"]:.0f}): {out["mismatch_center"]:+.1e}', True,
                'medido; no es un fallo del código sino de los ajustes publicados'))
    big = max(sh[k]['dRe_rel'] for k in sh)
    res.append((f'H27: con c_s^2 = n dp/dn/(eps + p) en vez de dp/deps, Re nu se mueve '
                + ', '.join(f"{k} {sh[k]['dRe_rel']:.1e} ({sh[k]['dRe_Hz']:+.2f} Hz)" for k in sh)
                + f"; Im nu de H1 {sh['H1']['dIm_rel']:.1e} y de F (3,4 kHz) {sh['F34']['dIm_rel']:.1e} relativo",
                1e-5 < big < 1e-2, 'mucho más que el 1e-6 de los controles del integrador: es la incertidumbre dominante del sector fluido'))
    # controles negativos
    Bsh = B.EOS.__new__(B.EOS)
    n = E.n
    e = B.eps_of_n(n)
    p = 10**B.lgP_of_lgrho(np.log10(e * B.MEV_FM3_TO_GCC) * 1.01) / B.MEV_FM3_TO_DYN
    ok = np.concatenate([[True], np.diff(p) > 0])
    from scipy.interpolate import CubicSpline
    Bsh.n, Bsh.e, Bsh.p = n[ok], e[ok], p[ok]
    Bsh._le, Bsh._ln = CubicSpline(np.log(Bsh.p), np.log(Bsh.e)), CubicSpline(np.log(Bsh.p), np.log(Bsh.n))
    m_bad = abs(mismatch(100.0, Bsh))
    res.append((f'Control negativo 1: con P(rho) desplazado un 1 % en lg rho, el desajuste a 100 MeV/fm3 sube a {m_bad:.1e} '
                f'(frente a {abs(out["mismatch"]["100"]):.1e})', m_bad > 10 * abs(out['mismatch']['100']), ''))
    sm = out['shift_mid']
    inter = all(sm[k]['dRe_rel'] < sh[k]['dRe_rel'] for k in sh)
    res.append(('Control negativo 2: con la media geométrica de las dos definiciones, el desplazamiento es menor: '
                + ', '.join(f"{k} {sm[k]['dRe_rel']:.1e}" for k in sm), inter, 'el efecto viene de c_s^2'))
    out['checks'] = [{'name': n_, 'ok': bool(o_), 'note': t_} for n_, o_, t_ in res]
    OUT.mkdir(exist_ok=True)
    (OUT / 'radial_eos.json').write_text(json.dumps(out, indent=1))
    return report('Paso 23: sensibilidad a la consistencia termodinámica de la EOS', res)


if __name__ == '__main__':
    run()
