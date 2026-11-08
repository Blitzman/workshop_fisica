"""Paso 18: la escalera del átomo gravitatorio sigue una ley de Rydberg con defecto cuántico constante (hipótesis H20).

Lee los niveles del paso 15 (results/radial_star.json, nu_a = 6 kHz) y no ajusta nada para la prueba principal.
El parámetro de Coulomb relativista es el de la ecuación de Whittaker exterior, derivado y verificado en el paso 5
(Ec. A14 del Paper II): kappa = alpha (1 - 2 q^2)/q, con q^2 = 1 - (nu/nu_a)^2 y alpha = G M m_a. Si la región
interior (estrella + campo fuerte) solo añade un desfase independiente de la energía, kappa_n = n - delta con delta
constante (y correcciones de Rydberg-Ritz delta_2/(n - delta_0)^2).
"""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from .pipeline import report

RES = Path(__file__).resolve().parent.parent / 'results'


def kappa_whittaker(alpha, x):
    q = np.sqrt(1 - x * x)
    return alpha * (1 - 2 * q * q) / q


def n_hydrogen(alpha, x):                      # la fórmula no relativista usada hasta ahora
    return alpha / np.sqrt(2 * (1 - x))


def kappa_newton(alpha, x):                    # control negativo: Coulomb newtoniano con la q exacta, sin el 2 omega^2
    return alpha / np.sqrt(1 - x * x)


def run():
    S = json.loads((RES / 'radial_star.json').read_text())
    a, nua = S['nu6']['alpha_g'], 6.0
    ns = np.array([l['n'] for l in S['nu6']['ladder']], dtype=float)
    xs = np.array([l['nu'] for l in S['nu6']['ladder']]) / nua
    d_w = ns - kappa_whittaker(a, xs)
    d_h = ns - n_hydrogen(a, xs)
    d_n = ns - kappa_newton(a, xs)
    rng = lambda d: float(d.max() - d.min())
    fit = least_squares(lambda p: p[0] + p[1] / (ns - p[0])**2 - d_w, [0.05, 0.0])
    res = [
        (f'H20(a): con el kappa de Whittaker (Ec. A14), n - kappa_n = {", ".join(f"{v:.4f}" for v in d_w)} para '
         f'n = {int(ns[0])}…{int(ns[-1])}: constante dentro de {rng(d_w):.4f}', rng(d_w) < 0.005, 'sin ajuste'),
        (f'H20(b): Rydberg-Ritz delta_n = delta_0 + delta_2/(n - delta_0)^2 con delta_0 = {fit.x[0]:.4f}, '
         f'delta_2 = {fit.x[1]:.4f}; residuo máximo {abs(fit.fun).max():.0e}', abs(fit.fun).max() < 1e-4, '2 parámetros, 5 niveles'),
        (f'Control negativo 1: con la fórmula hidrogenoide no relativista, n - n_eff varía {rng(d_h):.3f} '
         f'(lo que parecía un defecto creciente)', rng(d_h) > 0.02, ''),
        (f'Control negativo 2: con el Coulomb newtoniano alpha/q (sin el factor 2 omega^2 - m_a^2), n - kappa varía '
         f'{rng(d_n):.3f}', rng(d_n) > 0.02, ''),
    ]
    out = {'alpha_g': a, 'n': ns.tolist(), 'x': xs.tolist(), 'n_minus_kappa': d_w.tolist(),
           'n_minus_neff': d_h.tolist(), 'n_minus_kappa_newton': d_n.tolist(),
           'ritz': {'delta0': float(fit.x[0]), 'delta2': float(fit.x[1]), 'resid': float(abs(fit.fun).max())},
           'checks': [{'name': n_, 'ok': bool(o_), 'note': t_} for n_, o_, t_ in res]}
    (RES / 'radial_rydberg.json').write_text(json.dumps(out, indent=1))
    return report('Paso 18: ley de Rydberg del átomo gravitatorio', res)


if __name__ == '__main__':
    run()
