"""D2: nonradial polar perturbations WITH the perturbed metric (beyond Cowling), derived from the action.

Moved from axion_ns/fullgr.py (the derivation half: perturbed metric, linearized geometry, field equations,
`linear_system`). The Lindblom-Detweiler transcriptions moved to R; the extended system (LD + axion) is assembled in
ext_system.py and verified in checks_fullgr.py. Convention of Lindblom and Detweiler (LD), in signature (+,-,-,-)
and e^{-i w t}:

    ds^2 = e^{2 Phi}(1 + H0 rY) dt^2 - 2 i w r H1 rY dt dr - e^{2 lambda}(1 - H2 rY) dr^2
           - r^2 (1 - K rY) dOmega^2,          rY = (r/R)^l Y(th) e^{-i w t}

Geometry is linearized by hand (Gamma, Ricci and Einstein to first order). With H0 = H1 = H2 = K = 0 the Cowling
case is recovered. This file knows nothing about the equations of the papers.
"""
import sympy as sp

from common.notation import (t, r, th, X, G, V, eps_of, p_of, Phi0, lam0, a0, n0, w, l, Rs, fa, eta, A, N, XR, XT,
                             dth, dn, W, Vt, Y, H0, H1, H2, K)

from .derive_cowling import _divide, _legendre

E = sp.exp(-sp.I * w * t)
rl = (r / Rs)**l

g0 = sp.diag(sp.exp(2 * Phi0), -sp.exp(2 * lam0), -r**2, -r**2 * sp.sin(th)**2)
g0inv = sp.diag(sp.exp(-2 * Phi0), -sp.exp(-2 * lam0), -1 / r**2, -1 / (r**2 * sp.sin(th)**2))
h = sp.zeros(4, 4)
h[0, 0] = sp.exp(2 * Phi0) * H0
h[0, 1] = h[1, 0] = -sp.I * w * r * H1
h[1, 1] = sp.exp(2 * lam0) * H2
h[2, 2] = r**2 * K
h[3, 3] = r**2 * sp.sin(th)**2 * K
h = h * rl * Y * E
g1inv = -g0inv * h * g0inv                      # g^{mu nu} a primer orden
sqrtg0 = sp.exp(Phi0 + lam0) * r**2 * sp.sin(th)
sqrtg1 = sqrtg0 * sum(g0inv[i, i] * h[i, i] for i in range(4)) / 2


def _gamma(gup, gdn_list):
    """Gamma^a_bc = 1/2 g^{ad}(d_b g_dc + d_c g_db - d_d g_bc), con g^{ad} diagonal o no."""
    Gam = [[[0] * 4 for _ in range(4)] for _ in range(4)]
    for a in range(4):
        for b in range(4):
            for c in range(b, 4):
                s = 0
                for gu, gd in gdn_list(gup):
                    s += sum(gu[a, d_] * (sp.diff(gd[d_, c], X[b]) + sp.diff(gd[d_, b], X[c])
                                          - sp.diff(gd[b, c], X[d_])) for d_ in range(4)) / 2
                Gam[a][b][c] = Gam[a][c][b] = sp.expand(s)
    return Gam


def _ricci(G0, G1=None):
    """R_bd = d_a Gam^a_db - d_d Gam^a_ab + Gam^a_ae Gam^e_db - Gam^a_de Gam^e_ab (parte lineal si G1)."""
    Ric = sp.zeros(4, 4)
    for b in range(4):
        for d_ in range(b, 4):
            if G1 is None:
                s = sum(sp.diff(G0[a][d_][b], X[a]) - sp.diff(G0[a][a][b], X[d_]) for a in range(4))
                s += sum(G0[a][a][e] * G0[e][d_][b] - G0[a][d_][e] * G0[e][a][b]
                         for a in range(4) for e in range(4))
            else:
                s = sum(sp.diff(G1[a][d_][b], X[a]) - sp.diff(G1[a][a][b], X[d_]) for a in range(4))
                s += sum(G0[a][a][e] * G1[e][d_][b] + G1[a][a][e] * G0[e][d_][b]
                         - G0[a][d_][e] * G1[e][a][b] - G1[a][d_][e] * G0[e][a][b]
                         for a in range(4) for e in range(4))
            Ric[b, d_] = Ric[d_, b] = sp.expand(s)
    return Ric


def linear_einstein():
    """G_{mu nu} a primer orden: (G0, G1)."""
    Gam0 = _gamma(None, lambda _: [(g0inv, g0)])
    Gam1 = _gamma(None, lambda _: [(g0inv, h), (g1inv, g0)])
    Ric0, Ric1 = _ricci(Gam0), _ricci(Gam0, Gam1)
    R0 = sum(g0inv[i, i] * Ric0[i, i] for i in range(4))
    R1 = sum(g0inv[i, i] * Ric1[i, i] for i in range(4)) + \
        sum(g1inv[i, j] * Ric0[i, j] for i in range(4) for j in range(4))
    G0 = Ric0 - R0 * g0 / 2
    G1 = Ric1 - (R1 * g0 + R0 * h) / 2
    return G0, G1, Gam0, Gam1


def derive():
    """Ecuaciones de campo con eta marcando el orden en la métrica; los campos materiales son
    generales (t, r, th) y se linealizan después con `linearize`."""
    G0, G1, Gam0, Gam1 = linear_einstein()
    g = g0 + eta * h
    ginv = g0inv + eta * g1inv
    sqrtg = sqrtg0 + eta * sqrtg1
    Gam = [[[Gam0[a][b][c] + eta * Gam1[a][b][c] for c in range(4)] for b in range(4)] for a in range(4)]
    dA = [sp.diff(A, x) for x in X]

    # Klein-Gordon por Euler-Lagrange (n fijo); normalización box a + V' + d eps/da
    L_a = sp.Rational(1, 2) * sum(ginv[i, j] * dA[i] * dA[j] for i in range(4) for j in range(4)) \
        - V(A) - eps_of(N, A)
    KG = -sp.euler_equations(sqrtg * L_a, A, [t, r, th])[0].lhs / sqrtg0

    # T_{mu nu} del axión por variación respecto a g^{mu nu}
    hh = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f'h{i}{j}'))
    L_h = sp.Rational(1, 2) * sum(hh[i, j] * dA[i] * dA[j] for i in range(4) for j in range(4)) - V(A)
    L_val = L_h.subs({hh[i, j]: ginv[i, j] for i in range(4) for j in range(4)})
    T_ax = sp.Matrix(4, 4, lambda i, j: 2 * sp.diff(L_h, hh[i, j]) - g[i, j] * L_val)

    # fluido perfecto: u^mu = u^t (1, d_t xi^r, d_t xi^th, 0), normalizada con la métrica perturbada
    v = sp.Matrix([1, sp.diff(XR, t), sp.diff(XT, t), 0])
    ut = 1 / sp.sqrt((v.T * g * v)[0, 0])
    u_up = ut * v
    u_dn = g * u_up
    e, p = eps_of(N, A), p_of(N)
    T_fl = sp.Matrix(4, 4, lambda i, j: (e + p) * u_dn[i] * u_dn[j] - p * g[i, j])
    T = T_fl + T_ax

    out = {'KG': KG}
    for (i, j), name in [((0, 0), 'Ett'), ((0, 1), 'Etr'), ((0, 2), 'Etth'), ((1, 1), 'Err'),
                         ((1, 2), 'Erth'), ((2, 2), 'Ethth'), ((3, 3), 'Ephph')]:
        out[name] = (G0[i, j] + eta * G1[i, j]) - 8 * sp.pi * G * T[i, j]

    def div(nu):  # nabla_mu T^mu_nu con la conexión a primer orden
        M = ginv * T
        s = sum(sp.diff(sqrtg * M[m, nu], X[m]) for m in range(4)) / sqrtg
        return s - sum(Gam[a][m][nu] * M[m, a] for m in range(4) for a in range(4))

    for nu, name in [(1, 'FLUID_r'), (2, 'FLUID_th')]:
        out[name] = div(nu) - KG * sp.diff(A, X[nu])
    out['BAR'] = sum(sp.diff(sqrtg * N * u_up[m], X[m]) for m in range(4)) / (r**2 * sp.sin(th))
    return out


# factores angulares tras eliminar Y'' con Legendre
Yp = sp.diff(Y, th)
ANG = {'Y': Y, 'dY': Yp, 'traceless': -2 * sp.cos(th) / sp.sin(th) * Yp - l * (l + 1) * Y}


def linearize(expr, angular):
    """Parte lineal (métrica y materia), dividida por e^{-i w t} (r/R)^l y el factor angular."""
    sub = {A: a0 + eta * fa * rl * dth * Y * E,
           N: n0 + eta * rl * dn * Y * E,
           XR: eta * rl * sp.exp(-lam0) * W / r * Y * E,
           XT: -eta * rl * Vt / r**2 * Yp * E}
    lin = sp.diff(expr.subs(sub).doit(), eta).subs(eta, 0)
    lin = _legendre(sp.expand(lin / E / rl))
    lin = sp.expand(sp.powdenest(sp.powsimp(sp.expand(_divide(lin, ANG[angular])), force=True), force=True))
    if lin.has(th):
        raise ValueError('la dependencia angular no se factoriza')
    return lin


def linear_system():
    """Todas las ecuaciones linealizadas, con su factor angular."""
    D = derive()
    out = {}
    for k, ang in [('Ett', 'Y'), ('Etr', 'Y'), ('Err', 'Y'), ('Etth', 'dY'), ('Erth', 'dY'),
                   ('KG', 'Y'), ('FLUID_r', 'Y'), ('FLUID_th', 'dY'), ('BAR', 'Y')]:
        out[k] = linearize(D[k], ang)
    s = sp.sin(th)**2
    out['Eang'] = linearize(D['Ethth'] + D['Ephph'] / s, 'Y')        # traza angular
    out['Etl'] = linearize(D['Ethth'] - D['Ephph'] / s, 'traceless')  # parte sin traza
    return out


