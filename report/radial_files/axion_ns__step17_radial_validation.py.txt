"""Paso 17: validación externa de los modos radiales en GR sin axión (hipótesis H19, regla del §4 de AGENTS.md).

Referencia: Font et al. 2002, PRD 65, 084024 (arXiv:gr-qc/0110047), Tabla 2: politropo N = 1, K = 100
(G = c = M_sun = 1), rho_c = 1.28e-3 (densidad de masa en reposo), M = 1.4 M_sun, R = 14.15 km (Sec. IV.1);
frecuencias de la teoría de perturbaciones lineales en evolución acoplada: F, H1, H2, H3 = 1.442, 3.955, 5.916,
7.776 kHz. Segunda fuente: tabla del grupo de Garching (rel_hydro/rotating_neutron_stars/tables.shtml), modelo B0:
F = 1.458, H1 = 3.971 kHz, sin indicar el método.

Se usa radial.interior sin tocarlo (con nu_a = None no hay axión): los modos son los ceros de Delta p(R)(omega^2).
"""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from . import axstar as AX
from . import radial as RD
from .pipeline import ok_field, report

OUT = Path(__file__).resolve().parent.parent / 'results'
FONT = {'F': 1.442, 'H1': 3.955, 'H2': 5.916, 'H3': 7.776}          # kHz, Tabla 2
GARCHING = {'F': 1.458, 'H1': 3.971}                                 # kHz, modelo B0
K_CODE, RHO_C, GAMMA = 100.0, 1.28e-3, 2.0
M_REF, R_REF = 1.4, 14.15                                            # M_sun, km (Sec. IV.1)


class Poly:
    """p = K rho^2, eps = rho + p/(Gamma - 1), en MeV/fm^3. gamma1: índice adiabático de las perturbaciones
    (Gamma salvo en el control negativo). internal=False: eps = rho (control negativo)."""
    p_min = 1e-12

    def __init__(self, gamma1=GAMMA, internal=True):
        self.Km = K_CODE * AX.MSUN**2 * AX.KGEO                      # (MeV/fm^3)^-1
        self.gamma1, self.internal = gamma1, internal

    def rho(self, p):
        return np.sqrt(p / self.Km)

    def eps(self, p):
        return self.rho(p) + (p / (GAMMA - 1) if self.internal else 0.0)

    def nb(self, p):
        return self.rho(p)

    def cs2(self, p):
        return self.gamma1 * p / (self.eps(p) + p)


def p_central():
    rho_geo = RHO_C / AX.MSUN**2                                     # km^-2
    return K_CODE * AX.MSUN**2 * rho_geo**2 / AX.KGEO                # MeV/fm^3


def modes(eos, numax=8.5, n=4):
    keep, AX._eos = AX._eos, eos
    try:
        st = AX.Star(None, p0=p_central())
        g = lambda nu: RD.interior(st, AX.omega_of(nu), 'fluid').y[1, -1].real
        nus = np.linspace(0.3, numax, 140)
        v = [g(x) for x in nus]
        roots = [brentq(g, nus[i], nus[i + 1], xtol=1e-10) for i in range(len(nus) - 1) if v[i] * v[i + 1] < 0]
        return st.M / AX.MSUN, st.R, roots[:n]
    finally:
        AX._eos = keep


def run():
    res = []
    M, R, nu = modes(Poly())
    res.append((f'Fondo: M = {M:.4f} M_sun, R = {R:.3f} km (Font et al.: 1.4, 14.15)',
                abs(M - M_REF) < 0.005 and abs(R - R_REF) < 0.01, 'H19(a)'))
    dev = {k: 100 * (v / FONT[k] - 1) for k, v in zip(FONT, nu)}
    res.append(('Modos radiales sin axión frente a Font et al. 2002, Tabla 2: ' +
                ', '.join(f'{k} = {v:.4f} kHz (ref. {FONT[k]}, {dev[k]:+.3f} %)' for k, v in zip(FONT, nu)),
                len(nu) == 4 and all(abs(d) < 0.1 for d in dev.values()), 'H19(b)'))
    dg = {k: 100 * (v / GARCHING[k] - 1) for k, v in zip(GARCHING, nu)}
    res.append(('Frente a la tabla de Garching (modelo B0): ' +
                ', '.join(f'{k} {dg[k]:+.2f} %' for k in GARCHING) + ' (método no indicado en la tabla)', None, 'H19(c)'))
    _, _, nu_g = modes(Poly(gamma1=2.1), n=1)
    _, _, nu_e = modes(Poly(internal=False), n=1)
    d1, d2 = 100 * (nu_g[0] / FONT['F'] - 1), 100 * (nu_e[0] / FONT['F'] - 1)
    res.append((f'Control negativo 1: Gamma_1 = 2.1 en las perturbaciones mueve F un {d1:+.2f} %', abs(d1) > 0.1, ''))
    res.append((f'Control negativo 2: eps = rho (sin energía interna) mueve F un {d2:+.2f} %', abs(d2) > 0.1, ''))
    out = {'M': M, 'R': R, 'modes': dict(zip(FONT, nu)), 'font': FONT, 'garching': GARCHING, 'dev_font_pct': dev,
           'dev_garching_pct': dg, 'neg_gamma1': d1, 'neg_no_internal': d2,
           'checks': [{'name': n_, 'ok': ok_field(o_), 'note': t_} for n_, o_, t_ in res]}
    OUT.mkdir(exist_ok=True)
    (OUT / 'radial_validation.json').write_text(json.dumps(out, indent=1))
    return report('Paso 17: validación externa de los modos radiales (Font et al. 2002)', res)


if __name__ == '__main__':
    run()
