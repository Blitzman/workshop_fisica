"""D2: nonradial perturbations in the Cowling approximation, derived from the action (Paper II setting).

Moved from axion_ns/nonradial.py (everything except the equations typed in the paper, which now live in R).
Same action as D0, with the background metric frozen (Cowling: no metric perturbations) and fields depending on
(t, r, theta). Polar perturbations in the convention of the paper (factor (r/R)^l extracted), axisymmetric m = 0:

    a     = a0 + f_a (r/R)^l dtheta(r) Y e^{-i w t}
    n     = n0 +     (r/R)^l dn(r)     Y e^{-i w t}
    xi^r  =  (r/R)^l e^{-lambda} W(r) / r     Y  e^{-i w t}
    xi^th = -(r/R)^l V(r) / r^2         dY/dth e^{-i w t}

with Y(th) obeying Legendre's equation  Y'' + cot(th) Y' + l(l+1) Y = 0.
This file knows nothing about the equations of the papers.
"""
import sympy as sp

from common.geometry import christoffel, div_mixed
from common.notation import (t, r, th, X, V, mstar, eps_of, p_of, Phi0, lam0, a0, n0, w, l, Rs, fa, eta, A, N, XR, XT,
                             dth, dn, W, Vt, Y)

g = sp.diag(sp.exp(2 * Phi0), -sp.exp(2 * lam0), -r**2, -r**2 * sp.sin(th)**2)
sqrtg = sp.exp(Phi0 + lam0) * r**2 * sp.sin(th)


def derive():
    ginv = g.inv()
    Gam = christoffel(g, ginv, X)
    dA = [sp.diff(A, x) for x in X]

    L_a = sp.Rational(1, 2) * sum(ginv[i, i] * dA[i]**2 for i in range(4)) - V(A) - eps_of(N, A)
    KG = -sp.euler_equations(sqrtg * L_a, A, [t, r, th])[0].lhs / sqrtg

    h = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f'h{i}{j}'))
    L_h = sp.Rational(1, 2) * sum(h[i, j] * dA[i] * dA[j] for i in range(4) for j in range(4)) - V(A)
    L_val = L_h.subs({h[i, j]: ginv[i, j] for i in range(4) for j in range(4)})
    T_ax = sp.Matrix(4, 4, lambda i, j: 2 * sp.diff(L_h, h[i, j]) - g[i, j] * L_val)

    vr, vt = sp.diff(XR, t), sp.diff(XT, t)
    ut = 1 / sp.sqrt(sp.exp(2 * Phi0) - sp.exp(2 * lam0) * vr**2 - r**2 * vt**2)
    u_up = sp.Matrix([ut, ut * vr, ut * vt, 0])
    u_dn = g * u_up
    e, p = eps_of(N, A), p_of(N)
    T_fl = sp.Matrix(4, 4, lambda i, j: (e + p) * u_dn[i] * u_dn[j] - p * g[i, j])

    mixed = lambda Tdn: ginv * Tdn
    out = {'KG': KG}
    for nu, name in [(1, 'FLUID_r'), (2, 'FLUID_th')]:
        divT = div_mixed(mixed(T_fl), Gam, sqrtg, X, nu) + div_mixed(mixed(T_ax), Gam, sqrtg, X, nu)
        out[name] = divT - KG * sp.diff(A, X[nu])
    out['BAR'] = sum(sp.diff(sqrtg * N * u_up[m], X[m]) for m in range(4)) / (r**2 * sp.sin(th))
    # T^r_t del axión (flujo de energía) para el apéndice B del paper
    out['Trt_ax'] = mixed(T_ax)[1, 0]
    return out


def _legendre(expr):
    """Elimina Y'' y Y''' con la ecuación de Legendre."""
    Y2 = -sp.cos(th) / sp.sin(th) * sp.diff(Y, th) - l * (l + 1) * Y
    Y3 = sp.diff(Y2, th).subs(sp.diff(Y, th, 2), Y2)
    return expr.subs(sp.diff(Y, th, 3), Y3).subs(sp.diff(Y, th, 2), Y2)


def _divide(lin, ang):
    """lin / ang, ambos lineales en Y e Y'. Sin `simplify`: su trigsimp puede no terminar según
    PYTHONHASHSEED. Cociente de polinomios en (cos, sin), reducido módulo cos^2 + sin^2 = 1."""
    y, yp, c, s = sp.symbols('y yp c s')
    Yp = sp.diff(Y, th)
    trig = {sp.tan(th): sp.sin(th) / sp.cos(th), sp.cot(th): sp.cos(th) / sp.sin(th)}
    to = lambda e: sp.expand(e.subs(trig).subs(Yp, yp).subs(Y, y).subs({sp.cos(th): c, sp.sin(th): s}))
    num, den = to(lin), to(ang)
    v = y if den.coeff(y) != 0 else yp
    q = sp.cancel(num.coeff(v) / den.coeff(v))
    qn, qd = sp.fraction(q)
    one = [c**2 + s**2 - 1]
    qn, qd = (sp.reduced(sp.expand(e), one, c, s)[1] if e.has(c, s) else e for e in (qn, qd))
    q = sp.cancel(qn / qd)
    if sp.expand(sp.reduced(sp.numer(sp.together(num - q * den)), one, c, s, y, yp)[1]) != 0:
        raise ValueError('la dependencia angular no se factoriza')
    return q.subs({c: sp.cos(th), s: sp.sin(th)})


def linearize(expr, angular):
    """Parte lineal; divide por e^{-i w t} (r/R)^l y por Y (angular='Y') o dY/dth (angular='dY')."""
    E = sp.exp(-sp.I * w * t)
    rl = (r / Rs)**l
    sub = {A: a0 + eta * fa * rl * dth * Y * E,
           N: n0 + eta * rl * dn * Y * E,
           XR: eta * rl * sp.exp(-lam0) * W / r * Y * E,
           XT: -eta * rl * Vt / r**2 * sp.diff(Y, th) * E}
    lin = sp.diff(expr.subs(sub).doit(), eta).subs(eta, 0)
    lin = _legendre(sp.expand(lin / E / rl))
    div = Y if angular == 'Y' else sp.diff(Y, th)
    lin = sp.expand(sp.powdenest(sp.powsimp(sp.expand(_divide(lin, div)), force=True), force=True))
    if lin.has(th):
        raise ValueError('la dependencia angular no se factoriza')
    return lin
