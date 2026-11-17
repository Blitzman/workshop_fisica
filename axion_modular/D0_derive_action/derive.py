"""D0: derivation of the field equations from the action (moved from axion_ns/model.py, logic unchanged).

    S = Int d^4x sqrt(-g) [ -R/(16 pi G) + 1/2 g^{mu nu} d_mu a d_nu a - V(a) - eps(n, a) ],
    eps(n, a) = eps0(n) + n [m*(a) - m_N]   =>   p = n d eps/dn - eps = p0(n).

Everything is derived from this action: Einstein (geometry + T_{mu nu} by metric variation),
Klein-Gordon (Euler-Lagrange) and the fluid (conservation of total T + baryon number).
Fields are general functions of (t, r); the static background and the perturbations are obtained
afterwards by substitution (background() here; linearization lives in D1).

This file knows nothing about the equations of the papers (rule of axion_ns/AGENTS.md section 1).
"""
import sympy as sp

from common.geometry import einstein, div_mixed
from common.notation import (t, r, th, X, G, Phi, lam, a, n, xi, V, mstar, Phi0, lam0, a0, n0, eps_of, p_of)

g = sp.diag(sp.exp(2 * Phi), -sp.exp(2 * lam), -r**2, -r**2 * sp.sin(th)**2)
sqrtg = sp.exp(Phi + lam) * r**2 * sp.sin(th)


def derive_field_equations():
    """Dict with the complete field equations (= 0), depending on (t, r)."""
    Gmn, Gam, ginv = einstein(g, X)
    da = [sp.diff(a, x) for x in X]

    # --- Klein-Gordon by Euler-Lagrange of the density sqrt(-g) L_a (n fixed)
    L_a = sp.Rational(1, 2) * sum(ginv[i, j] * da[i] * da[j] for i in range(4) for j in range(4)) \
        - V(a) - eps_of(n, a)
    KG_el = sp.euler_equations(sqrtg * L_a, a, [t, r])[0].lhs
    # normalization: KG = box a + V'(a) + n m*'(a)
    KG = sp.simplify(-KG_el / sqrtg)

    # --- axion T_{mu nu} by variation with respect to the inverse metric:  T = 2 dL/dg^{mu nu} - g L
    h = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f'h{i}{j}'))
    L_h = sp.Rational(1, 2) * sum(h[i, j] * da[i] * da[j] for i in range(4) for j in range(4)) - V(a)
    L_val = L_h.subs({h[i, j]: ginv[i, j] for i in range(4) for j in range(4)})
    T_ax = sp.Matrix(4, 4, lambda i, j: 2 * sp.diff(L_h, h[i, j]) - g[i, j] * L_val)

    # --- perfect fluid: u^r = u^t d_t xi, normalized
    ut = 1 / sp.sqrt(sp.exp(2 * Phi) - sp.exp(2 * lam) * sp.diff(xi, t)**2)
    u_up = sp.Matrix([ut, ut * sp.diff(xi, t), 0, 0])
    u_dn = g * u_up
    e, p = eps_of(n, a), p_of(n)
    T_fl = sp.Matrix(4, 4, lambda i, j: (e + p) * u_dn[i] * u_dn[j] - p * g[i, j])

    T = T_fl + T_ax
    EIN = {(i, j): Gmn[i, j] - 8 * sp.pi * G * T[i, j] for (i, j) in [(0, 0), (0, 1), (1, 1)]}

    # --- conservation: nabla_mu T^mu_nu(total) = 0, eliminating KG
    mixed = lambda Tdn: ginv * Tdn  # T^mu_nu = g^{mu a} T_{a nu}
    divT_fl = div_mixed(mixed(T_fl), Gam, sqrtg, X, 1)
    divT_ax = div_mixed(mixed(T_ax), Gam, sqrtg, X, 1)
    # identity: nabla T_ax^mu_r = (box a + V') d_r a  =>  with KG: = -n m*' d_r a
    ident_ax = sp.simplify(divT_ax - (KG - n * sp.diff(mstar(a), a)) * sp.diff(a, r))
    FLUID_r = divT_fl + divT_ax - KG * sp.diff(a, r)

    # --- baryon number conservation: d_mu (sqrt(-g) n u^mu) = 0
    BAR = sum(sp.diff(sqrtg * n * u_up[m], X[m]) for m in range(4)) / (r**2 * sp.sin(th))

    return dict(Gtt=EIN[(0, 0)], Gtr=EIN[(0, 1)], Grr=EIN[(1, 1)], KG=KG,
                FLUID_r=FLUID_r, BAR=BAR, ident_ax=ident_ax)


def background(expr):
    """Static configuration (order eta^0)."""
    return sp.simplify(expr.subs({Phi: Phi0, lam: lam0, a: a0, n: n0, xi: 0}).doit())
