"""Paso 5 (Paper II, arXiv:2610.37308): perturbaciones no radiales en Cowling, desde la acción.

Ecuaciones (8a)-(8d), exterior (13), ajuste (15), regularidad en el origen (10), (16), (17)
y las transformaciones del apéndice A (A2-A18), Ecs. (21)-(23) y los flujos/energías del apéndice B.
"""
from functools import lru_cache

import sympy as sp

from . import nonradial as NR
from .compare import check_modulo, is_zero
from .model import r, t, w, V, mstar, eps0, Phi0, lam0, a0, n0
from .pipeline import background_rules, report

d = lambda f, k=1: sp.diff(f, r, k)
l, fa, Rs = NR.l, NR.fa, NR.Rs

MATCH = {
    '(8a) momento radial': ('FLUID_r', 'Y', NR.W),
    '(8b) bariones': ('BAR', 'Y', NR.dn),
    '(8c) momento angular': ('FLUID_th', 'dY', NR.Vt),
    '(8d) KG interior': ('KG', 'Y', d(NR.dth, 2)),
}


@lru_cache
def derived():
    D = NR.derive()
    return D, {name: NR.linearize(D[der], ang) for name, (der, ang, _) in MATCH.items()}


def equations_checks():
    _, L = derived()
    R = background_rules()
    res = []
    for name, (_, _, piv) in MATCH.items():
        ok, _ = check_modulo(NR.PAPER[name], L[name], piv, R)
        res.append((f'Paper II {name} == linealización de la acción (Cowling)', ok, ''))

    ep = NR.eps + NR.p
    bad = {
        '(8a) sin el término l/r en (dtheta\' + l dtheta/r)': (
            '(8a) momento radial', -n0 * NR.dm * l / r * NR.dth),
        '(8b) sin Phi\' en el lado izquierdo': ('(8b) bariones', -ep * NR.W * d(Phi0)),
        '(8d) con 2l/r en vez de 2(l+1)/r': ('(8d) KG interior', -sp.exp(-2 * lam0) * 2 / r * d(NR.dth)),
    }
    for label, (name, delta) in bad.items():
        ok, _ = check_modulo(NR.PAPER[name] + delta, L[name], MATCH[name][2], R)
        res.append((f'Control negativo: {label} -> se detecta', not ok, ''))

    # exterior (n_b = 0): (8d) con dtheta = (R/r)^l dtheta_ext, por e^{2 lambda} -> Ec. (13)
    u = sp.Function('u')(r)
    ext = NR.PAPER['(8d) KG interior'].subs(n0, 0).subs({NR.dn: 0, NR.dth: (Rs / r)**l * u}).doit()
    m2 = sp.diff(V(a0), a0, 2)
    eq13 = d(u, 2) + (d(Phi0) - d(lam0) + 2 / r) * d(u) \
        - sp.exp(2 * lam0) * (m2 + l * (l + 1) / r**2 - w**2 * sp.exp(-2 * Phi0)) * u
    res.append(('Exterior: (8d) con n_b=0 y dtheta=(R/r)^l dtheta_ext es la Ec. (13)',
                sp.simplify(ext * sp.exp(2 * lam0) * (r / Rs)**l - eq13) == 0, ''))
    # Ec. (15): ajuste de derivadas logarítmicas
    dth_in = (Rs / r)**l * u
    res.append(('Ec. (15): dtheta\'/dtheta = L_ext - l/R',
                sp.simplify(d(dth_in) / dth_in - (d(u) / u - l / r)) == 0, ''))
    return res


def origin_checks():
    """Desarrollo en serie en r=0 de (8a)-(8d) con fondo regular (par en r) genérico."""
    Pc, P2, L2, nc, N2, ac, A2 = sp.symbols('Phi_c Phi_2 lambda_2 n_c n_2 a_c a_2')
    e0, e1, e2, e3, m0, m1, m2, v0, v1, v2 = sp.symbols('e0:4 m0:3 v0:3')
    x = sp.Symbol('x')
    N0, N1, T0, T1, T2, X0, X1 = sp.symbols('N0 N1 T0 T1 T2 X0 X1')
    funcs = {eps0: sp.Lambda(x, e0 + e1 * (x - nc) + e2 * (x - nc)**2 / 2 + e3 * (x - nc)**3 / 6),
             mstar: sp.Lambda(x, m0 + m1 * (x - ac) + m2 * (x - ac)**2 / 2),
             V: sp.Lambda(x, v0 + v1 * (x - ac) + v2 * (x - ac)**2 / 2)}
    bg = {Phi0: Pc + P2 * r**2 / 2, lam0: L2 * r**2 / 2, n0: nc + N2 * r**2 / 2, a0: ac + A2 * r**2 / 2}

    def concrete(e):
        return e.subs(funcs).doit().subs(bg).doit()

    ep = concrete(NR.eps + NR.p)
    pert = {NR.dn: N0 + N1 * r, NR.dth: T0 + T1 * r + T2 * r**2,
            NR.W: (X0 + X1 * r) / ep}
    Vsol = sp.solve(NR.PAPER['(8c) momento angular'], NR.Vt)[0]

    def ser(name, order=1):
        e = NR.PAPER[name].subs(NR.Vt, Vsol)
        e = concrete(e).subs(pert).doit()
        sr = sp.series(e * r, r, 0, order + 1).removeO()   # r * ecuación = c_{-1} + c_0 r + ...
        return [sp.simplify(sr.subs(r, 0)), sp.simplify(sp.diff(sr, r).subs(r, 0))]

    mom, bar, kg = ser('(8a) momento radial'), ser('(8b) bariones'), ser('(8d) KG interior')
    res = []
    sT1 = sp.solve(kg[0], T1)
    res.append(('Ec. (10b): el orden 1/r de (8d) impone dtheta\'(0) = 0', sT1 == [0], ''))
    sX0 = sp.solve(mom[0].subs(T1, 0), X0)[0]
    dp0 = concrete(NR.dpdn).subs(r, 0) * N0
    dmc = fa * m1
    eq10a = sp.exp(-2 * Pc) / l * w**2 * sX0 - nc * dmc * T0 - dp0
    res.append(('Ec. (10a): delta p(0) del orden 1/r de (8a)', sp.simplify(eq10a) == 0, ''))
    res.append(('El orden 1/r de (8b) es compatible con (10a) (misma condición)',
                sp.simplify(bar[0].subs({T1: 0, X0: sX0})) == 0, ''))
    sol = sp.solve([mom[1].subs({T1: 0, X0: sX0}), bar[1].subs({T1: 0, X0: sX0})], [N1, X1], dict=True)
    ok16 = len(sol) == 1 and sp.simplify(sol[0][N1]) == 0 and sp.simplify(sol[0][X1]) == 0
    res.append(('Ecs. (16a,b): delta p\'(0) = 0 y [(eps+p)W]\'(0) = 0', ok16, ''))
    sT2 = sp.solve(kg[1].subs(T1, 0), T2)[0]
    meff0 = v2 + nc * m2
    paper16c = (T0 * (-w**2 * sp.exp(-2 * Pc) - l * P2 + l * (l + 2) * L2 + meff0)
                + N0 / fa**2 * dmc) / (2 * l + 3)
    res.append(('Ec. (16c): dtheta\'\'(0)', sp.simplify(2 * sT2 - paper16c) == 0, ''))
    eq17 = concrete(NR.cs2 * (NR.eps + NR.p) / n0).subs(r, 0) * N0 - dp0
    res.append(('Ec. (17): delta p(0) = c_s^2 (eps+p) delta n_b(0)/n_b es exacta con p = p(n)',
                sp.simplify(eq17) == 0, ''))
    return res


def appendix_checks():
    """Apéndice A y Ecs. (21)-(23), en unidades m_a = 1 (r = r~, k = k~, GM = alpha)."""
    x, al, k = sp.symbols('x alpha k', positive=True)
    kc = sp.Symbol('k')  # k~ complejo
    res = []
    F = 1 - 2 * al / x
    u = sp.Function('u')(x)
    w2 = 1 + kc**2
    # Ec. (13) en Schwarzschild sin axión de fondo (m_eff = m_a = 1)
    ext = sp.diff(u, x, 2) + (sp.diff(F, x) / F + 2 / x) * sp.diff(u, x) - (1 + l * (l + 1) / x**2 - w2 / F) / F * u
    Psi = sp.Function('Psi')(x)
    U = (kc**2 + 2 * al / x - l * (l + 1) / x**2 * F + (al / x**2)**2) / F**2
    lio = ext.subs(u, F**sp.Rational(-1, 2) * Psi / x).doit()
    lio = lio.subs(sp.diff(Psi, x, 2), -U * Psi)
    res.append(('Ecs. (A2)-(A4): Liouville elimina la derivada primera y da el potencial U', is_zero(lio), ''))

    Q, wf = sp.Function('Q')(x), sp.Function('w')(x)
    Psip = Q**sp.Rational(-1, 2) * sp.exp(sp.I * wf)
    milne = 3 / (2 * Q) * sp.diff(Q, x)**2 + 2 * U * Q - 2 * Q**3
    e = (sp.diff(Psip, x, 2) + U * Psip).subs(sp.diff(wf, x, 2), sp.diff(Q, x)).subs(sp.diff(wf, x), Q)
    e = e.subs(sp.diff(Q, x, 2), milne)
    res.append(('Ec. (A6): Psi = Q^{-1/2} e^{iw} resuelve (A3) si Q cumple la ecuación de Milne', is_zero(e), ''))
    Lp = sp.diff(Psip, x).subs(sp.diff(wf, x), Q) / Psip
    Lm = sp.diff(Q**sp.Rational(-1, 2) * sp.exp(-sp.I * wf), x).subs(sp.diff(wf, x), Q) \
        / (Q**sp.Rational(-1, 2) * sp.exp(-sp.I * wf))
    okA7 = is_zero(Lp - (sp.I * Q - sp.diff(Q, x) / (2 * Q))) and is_zero(Lm - (-sp.I * Q - sp.diff(Q, x) / (2 * Q)))
    okA18 = is_zero((Lp - Lm) / (2 * sp.I) - Q) and is_zero(-Q * (Lp + Lm) - sp.diff(Q, x))
    res.append(('Ecs. (A7) y (A18): derivadas logarítmicas L+- y su inversión para (Q, Q\')', okA7 and okA18, ''))
    th_ext = F**sp.Rational(-1, 2) * Psi / x
    LA10 = sp.diff(Psi, x) / Psi - 1 / x - al / (x**2 * F)
    res.append(('Ec. (A10): L_ext a partir de L+', is_zero(sp.diff(th_ext, x) / th_ext - LA10), ''))

    s = sp.Symbol('s', positive=True)
    Us = sp.series(U.subs(x, 1 / s), s, 0, 3).removeO()
    coul = kc**2 + 2 * al * (2 * kc**2 + 1) * s + (4 * al**2 * (3 * kc**2 + 2) - l * (l + 1)) * s**2
    res.append(('Ec. (A13): U desarrollado hasta O(1/r~^2)', sp.simplify(sp.expand(Us - coul)) == 0, ''))

    # Psi(x) = W(z), z = -2ik~x  =>  Psi'' = (-2ik~)^2 W''(z), con W'' = (1/4 - kappa/z - (1/4 - mu^2)/z^2) W
    kap = sp.I * al * (2 * kc**2 + 1) / kc
    mu2 = (l + sp.Rational(1, 2))**2 - 4 * al**2 * (3 * kc**2 + 2)
    zz = -2 * sp.I * kc * x
    whit = sp.Rational(1, 4) - kap / zz - (sp.Rational(1, 4) - mu2) / zz**2
    eqc = sp.simplify((-2 * sp.I * kc)**2 * whit + coul.subs(s, 1 / x))
    res.append(('Ec. (A14): z = -2ik~r~, kappa, mu reducen (A13) a la ecuación de Whittaker', eqc == 0, ''))

    # Ec. (21)-(22): forma de Schrödinger en coordenada tortuga (d/dr* = F d/dr), psi = r dtheta_ext
    psi = x * u
    dV = -2 * al / x + F * (l * (l + 1) / x**2 + 2 * al / x**3)
    schr = F * sp.diff(F * sp.diff(psi, x), x) + (kc**2 - dV) * psi
    res.append(('Ecs. (21)-(22): psi = r dtheta_ext en coordenada tortuga y potencial Delta V',
                is_zero(schr / (x * F**2) - ext), ''))
    # Ec. (23): niveles hidrogénicos n = l+1 (sin nodos) y n = l+2 (un nodo)
    c = sp.Symbol('c')
    hyd = []
    for nn, poly in [(l + 1, 1), (l + 2, 1 - c * x)]:
        ps = x**(l + 1) * poly * sp.exp(-al / nn * x)
        e = sp.diff(ps, x, 2) + (-(al / nn)**2 + 2 * al / x - l * (l + 1) / x**2) * ps
        hyd.append(sp.expand(sp.simplify(e / (x**(l - 1) * sp.exp(-al / nn * x)))))
    sol = sp.solve(sp.Poly(hyd[1], x).coeffs(), c, dict=True)
    ok = hyd[0] == 0 and len(sol) == 1
    res.append(('Ec. (23): k~_n^2 = -alpha^2/n^2 para el problema de Coulomb (n = l+1, l+2)', ok, ''))
    return res


def energy_checks():
    """Apéndice B: flujo (B4), luminosidad media (B6) y energías (B7) desde T_{mu nu} de la acción."""
    D, _ = derived()
    res = []
    th_f = sp.Function('theta')(t, r)
    Trt = D['Trt_ax'].subs(NR.A, fa * th_f).doit()
    sqrtg = sp.exp(Phi0 + lam0) * r**2
    B4 = -fa**2 * r**2 * sp.exp(Phi0 - lam0) * sp.diff(th_f, r) * sp.diff(th_f, t)
    res.append(('Ec. (B4): sqrt(-g) T^r_t del axión', sp.simplify(sqrtg * Trt - B4) == 0, ''))

    om = sp.Symbol('omega', positive=True)
    uR, uI = sp.Function('u_R', real=True)(r), sp.Function('u_I', real=True)(r)
    uc = uR + sp.I * uI
    th = sp.re(sp.expand_complex(uc * sp.exp(-sp.I * om * t)))
    avg = lambda e: sp.integrate(sp.expand(e), (t, 0, 2 * sp.pi / om)) * om / (2 * sp.pi)
    La = -fa**2 * r**2 * sp.exp(Phi0 - lam0) * avg(sp.diff(th, r) * sp.diff(th, t))
    ImL_abs2 = sp.im(sp.expand_complex(d(uc) * sp.conjugate(uc)))  # Im(L_ext)|u|^2 = Im(u' u*)
    ImL_abs2 = ImL_abs2.subs({sp.re(d(uR)): d(uR), sp.im(d(uR)): 0, sp.re(d(uI)): d(uI), sp.im(d(uI)): 0})
    B6 = fa**2 / 2 * r**2 * sp.exp(Phi0 - lam0) * om * ImL_abs2
    res.append(('Ec. (B6): luminosidad promediada <L> = f_a^2/2 r^2 e^{Phi-lambda} Re w Im L_ext |dtheta_ext|^2',
                sp.simplify(La - B6) == 0, ''))

    # (B7b): 2 x energía cinética media del axión = (sqrt(-g) T^t_t)_cin
    Ea = 2 * avg(sqrtg * fa**2 * sp.exp(-2 * Phi0) * sp.diff(th, t)**2 / 2)
    B7b = om**2 / 2 * fa**2 * sp.exp(lam0 - Phi0) * (uR**2 + uI**2) * r**2
    res.append(('Ec. (B7b): energía del modo axiónico', sp.simplify(Ea - B7b) == 0, ''))
    # (B7a): término cuadrático en velocidades de sqrt(-g) T^t_t del fluido, con xi^r, xi^th del paper
    vr, vt, ee, pp = sp.symbols('v_r v_t e p')
    ut2 = 1 / (sp.exp(2 * Phi0) - sp.exp(2 * lam0) * vr**2 - r**2 * vt**2)
    Ttt = (ee + pp) * sp.exp(2 * Phi0) * ut2 - pp
    sc = sp.Symbol('sc')
    kin = sp.series(Ttt.subs({vr: sc * vr, vt: sc * vt}), sc, 0, 3).removeO().subs(sc, 1) - ee
    Wc, Vc = sp.symbols('W V')
    # amplitudes de |xi^r|^2 y |xi^th|^2 ya integradas en ángulos: <Y^2> = 1, <(dY)^2> = l(l+1)
    rl2 = (r / Rs)**(2 * l)
    sub = {vr**2: om**2 / 2 * rl2 * sp.exp(-2 * lam0) * Wc**2 / r**2,
           vt**2: om**2 / 2 * rl2 * Vc**2 / r**4 * l * (l + 1)}
    Ef = sp.expand(sqrtg * kin).subs(sub)
    B7a = om**2 / 2 * (ee + pp) * sp.exp(lam0 - Phi0) * rl2 * (Wc**2 + l * (l + 1) * Vc**2)
    res.append(('Ec. (B7a): energía del modo fluido (término cuadrático de T^t_t)', sp.simplify(Ef - B7a) == 0, ''))
    return res


def run():
    res = equations_checks() + origin_checks() + appendix_checks() + energy_checks()
    return report('Paso 5: Paper II, perturbaciones no radiales (simbólico)', res)


if __name__ == '__main__':
    run()
