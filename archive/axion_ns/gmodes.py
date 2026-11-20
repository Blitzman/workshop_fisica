"""Línea 1: Estratificación de composición, modos de gravedad (g-modes) y modo de interfaz en la pared del axión.

Parte de la misma acción covariante que `model.py` y `nonradial.py`, extendida a materia nuclear
estratificada no barotrópica donde la composición química Y_p(t, r, theta) permanece congelada
durante la oscilación rápida (u^mu d_mu Y_p = 0 <=> Delta Y_p = 0):

    S = Int d^4x sqrt(-g) [ 1/2 g^{mu nu} d_mu a d_nu a - V(a) - eps_strat(n, Y_p, a) ]

con eps_strat(n, Y_p, a) = eps_nuc(n, Y_p) + n [m*(a) - m_N] y p_strat(n, Y_p) = n d(eps_nuc)/dn - eps_nuc.
De las ecuaciones de conservación covariantes nabla_mu T^{mu nu} = 0, nabla_mu (n u^mu) = 0,
u^mu d_mu Y_p = 0 y Klein-Gordon se deducen en SymPy:
1. El discriminante relativista de Schwarzschild con axión A(r) y la frecuencia de Brunt-Väisälä N^2(r).
2. Las relaciones eulerianas para delta eps y delta n con composición congelada.
3. El sistema no barotrópico completo de 4 ecuaciones en Cowling para (delta p, U = (eps+p)W, delta theta, delta theta').
4. El límite barotrópico c_s^2 -> c_e^2 (A(r) -> 0), que recupera idénticamente las Ecs. (8a)-(8d) del Paper II.
5. Las condiciones de regularidad en el origen (r -> 0) y la relación de dispersión WKB local.
6. La compilación directa de la matriz 4x4 del sistema de EDOs mediante `sp.lambdify(..., 'numpy', cse=True)`.

NOTA (Regla 1 de AGENTS.md): Este módulo jamás importa `gmodes_ref.py` ni `paper.py`.
"""
from functools import lru_cache

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from . import axstar as AX
from . import nonradial as NR
from .geometry import christoffel, div_mixed
from .model import r, th, t, w, G, mN, V, mstar, eps0, eps_of, p_of, Phi0, lam0, a0, n0

d = lambda f, k=1: sp.diff(f, r, k)
l, Rs, fa = NR.l, NR.Rs, NR.fa
dth, dn, W, Vt, Y = NR.dth, NR.dn, NR.W, NR.Vt, NR.Y

# Símbolos y funciones para la derivación no barotrópica
cs2 = sp.Symbol('c_s2', positive=True)
dp = sp.Function('deltap')(r)
Uf = sp.Function('U')(r)
kr = sp.Symbol('k_r', real=True)

# Composición / parámetro de estratificación Y_p(t, r, th)
Yp_full = sp.Function('Y_p')(t, r, th)
Yp0 = sp.Function('Y_p0')(r)
dYp = sp.Function('deltaY_p')(r)
eps_nuc = sp.Function('eps_nuc')


def eps_strat(nn, yp, aa):
    """Densidad de energía en materia estratificada: eps_nuc(n, Y_p) + n [m*(a) - m_N]."""
    return eps_nuc(nn, yp) + nn * (mstar(aa) - mN)


def p_strat(nn, yp):
    """Presión por la primera ley a composición y axión constantes: n (d eps/dn)_{Y_p, a} - eps."""
    x = sp.Dummy('x')
    e = eps_nuc(x, yp)
    return (x * sp.diff(e, x) - e).subs(x, nn)


def _unsubs(expr):
    """Convierte objetos Subs(Derivative(f(x), (x, k)), x, g(r)) de SymPy 1.14 en Derivative(f(g(r)), (g(r), k))."""
    rep = {}
    for s in expr.atoms(sp.Subs):
        if isinstance(s.expr, sp.Derivative):
            v, k = s.expr.variable_count[0]
            rep[s] = sp.Derivative(s.expr.expr.subs(v, s.point[0]), (s.point[0], k))
    return expr.xreplace(rep)


# Cantidades del fondo en equilibrio beta (d eps_nuc / d Y_p = 0 en el fondo, eps_nuc(n0, Yp0) = eps0(n0))
eps, p = eps_of(n0, a0), p_of(n0).doit()
dpdn_eq = sp.diff(p, n0)
ce2 = dpdn_eq * n0 / (eps + p)

ma = sp.diff(mstar(a0), a0)
maa = sp.diff(mstar(a0), a0, 2)
dm = fa * ma
dmm = fa**2 * maa
thp = d(a0) / fa
meff2 = (fa**2 * sp.diff(V(a0), a0, 2) + n0 * dmm) / fa**2
xi_r = sp.exp(-lam0) * W / r


@lru_cache
def derive_covariant():
    """Deriva desde la acción y conservación covariante las ecuaciones completas en (t, r, th)
    incluyendo la advección de la composición congelada u^mu d_mu Y_p = 0."""
    ginv = NR.g.inv()
    Gam = christoffel(NR.g, ginv, NR.X)
    dA = [sp.diff(NR.A, x) for x in NR.X]

    L_a = (sp.Rational(1, 2) * sum(ginv[i, i] * dA[i]**2 for i in range(4))
           - V(NR.A) - eps_strat(NR.N, Yp_full, NR.A))
    KG = -sp.euler_equations(NR.sqrtg * L_a, NR.A, [t, r, th])[0].lhs / NR.sqrtg

    h = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f'h{i}{j}'))
    L_h = sp.Rational(1, 2) * sum(h[i, j] * dA[i] * dA[j] for i in range(4) for j in range(4)) - V(NR.A)
    L_val = L_h.subs({h[i, j]: ginv[i, j] for i in range(4) for j in range(4)})
    T_ax = sp.Matrix(4, 4, lambda i, j: 2 * sp.diff(L_h, h[i, j]) - NR.g[i, j] * L_val)

    vr, vt = sp.diff(NR.XR, t), sp.diff(NR.XT, t)
    ut = 1 / sp.sqrt(sp.exp(2 * Phi0) - sp.exp(2 * lam0) * vr**2 - r**2 * vt**2)
    u_up = sp.Matrix([ut, ut * vr, ut * vt, 0])
    u_dn = NR.g * u_up
    e, pr = eps_strat(NR.N, Yp_full, NR.A), p_strat(NR.N, Yp_full)
    T_fl = sp.Matrix(4, 4, lambda i, j: (e + pr) * u_dn[i] * u_dn[j] - pr * NR.g[i, j])

    mixed = lambda Tdn: ginv * Tdn
    out = {'KG': KG}
    for nu, name in [(1, 'FLUID_r'), (2, 'FLUID_th')]:
        divT = div_mixed(mixed(T_fl), Gam, NR.sqrtg, NR.X, nu) + div_mixed(mixed(T_ax), Gam, NR.sqrtg, NR.X, nu)
        out[name] = divT - KG * sp.diff(NR.A, NR.X[nu])
    out['BAR'] = sum(sp.diff(NR.sqrtg * NR.N * u_up[m], NR.X[m]) for m in range(4)) / (r**2 * sp.sin(th))
    # Congelamiento de la composición en oscilaciones rápidas: u^mu d_mu Y_p = 0
    out['COMP'] = sum(u_up[m] * sp.diff(Yp_full, NR.X[m]) for m in range(4))
    return out


def linearize_strat(expr, angular):
    """Linealiza una ecuación covariante con (a, n, Y_p, xi^r, xi^th) y factoriza e^{-iwt} (r/R)^l Y."""
    E = sp.exp(-sp.I * w * t)
    rl = (r / Rs)**l
    sub = {
        NR.A: a0 + NR.eta * fa * rl * dth * Y * E,
        NR.N: n0 + NR.eta * rl * dn * Y * E,
        Yp_full: Yp0 + NR.eta * rl * dYp * Y * E,
        NR.XR: NR.eta * rl * sp.exp(-lam0) * W / r * Y * E,
        NR.XT: -NR.eta * rl * Vt / r**2 * sp.diff(Y, th) * E,
    }
    lin = sp.diff(expr.subs(sub).doit(), NR.eta).subs(NR.eta, 0)
    lin = NR._legendre(sp.expand(lin / E / rl))
    div = Y if angular == 'Y' else sp.diff(Y, th)
    lin = sp.expand(sp.powdenest(sp.powsimp(sp.expand(NR._divide(lin, div)), force=True), force=True))
    if lin.has(th):
        raise ValueError('la dependencia angular no se factoriza')
    return lin


@lru_cache
def derive_thermodynamics():
    """Deduce simbólicamente desde Delta Y_p = 0 y la primera ley:
    - delta Y_p = -xi^r Y_p0'
    - A(r) (discriminante de Schwarzschild con axión) y su equivalencia con (1/c_e^2 - 1/c_s^2) p'/(eps+p)
    - N^2(r) (frecuencia de Brunt-Väisälä relativista al cuadrado)
    - delta eps y delta n eulerianos en función de (delta p, W, delta theta)."""
    cov = derive_covariant()
    # 1. De u^mu d_mu Y_p = 0 linealizada se obtiene delta Y_p
    comp_lin = linearize_strat(cov['COMP'], 'Y')
    dYp_sol = sp.solve(comp_lin, dYp)[0]  # -e^{-lam0} W Y_p0' / r = -xi^r Y_p0'

    # 2. Primera ley lagrangiana con composición congelada (Delta Y_p = 0):
    #    Delta p = c_s^2 (Delta eps - n_0 dm*/da Delta a)
    #    Delta eps = (eps + p)/n_0 Delta n + n_0 dm*/da Delta a
    de_sym, dn_sym = sp.symbols('delta_eps delta_n')
    da_pert = fa * dth
    Delta_p = dp + d(p) * xi_r
    Delta_eps = de_sym + d(eps) * xi_r
    Delta_n = dn_sym + d(n0) * xi_r
    Delta_a = da_pert + d(a0) * xi_r

    eq_adiab = Delta_p - cs2 * (Delta_eps - n0 * ma * Delta_a)
    eq_first = Delta_eps - ((eps + p) / n0 * Delta_n + n0 * ma * Delta_a)

    sol = sp.solve([eq_adiab, eq_first], [de_sym, dn_sym], dict=True)[0]

    # Discriminante de Schwarzschild con axión extraído del coeficiente de xi^r en delta eps
    A_expr = sp.simplify(-sp.expand(sol[de_sym] - dp / cs2 - n0 * dm * dth) / ((eps + p) * xi_r))
    # Forma en términos de c_e^2 y c_s^2
    A_ce_expr = (1 / ce2 - 1 / cs2) * d(p) / (eps + p)

    # Gravedad efectiva (de TOV con axión: p' = -(eps+p) Phi' - n_0 dm*/da a_0')
    g_eff_p = -d(p) / (eps + p)
    g_eff_phi = d(Phi0) + n0 / (eps + p) * ma * d(a0)

    # Frecuencia de Brunt-Väisälä relativista al cuadrado
    N2_expr = -g_eff_phi * A_expr * sp.exp(2 * (Phi0 - lam0))

    de_clean = dp / cs2 - (eps + p) * A_expr * xi_r + n0 * dm * dth
    dn_clean = n0 / (eps + p) * dp / cs2 - n0 * A_expr * xi_r

    return {
        'dYp': dYp_sol,
        'A': A_expr,
        'A_ce': A_ce_expr,
        'g_eff_p': g_eff_p,
        'g_eff': g_eff_phi,
        'N2': N2_expr,
        'de_raw': sol[de_sym],
        'dn_raw': sol[dn_sym],
        'de': de_clean,
        'dn': dn_clean,
    }


@lru_cache
def cowling_rules():
    """Reglas de reducción del fondo en la aproximación de Cowling (equilibrio hidrostático y KG estática),
    deducidas a orden eta^0 de las ecuaciones de conservación covariantes."""
    nr_cov = NR.derive()
    sub0 = {NR.A: a0, NR.N: n0, NR.XR: 0, NR.XT: 0}
    bg_fl = _unsubs(sp.expand(nr_cov['FLUID_r'].subs(sub0).doit()))
    bg_kg = _unsubs(sp.expand(nr_cov['KG'].subs(sub0).doit()))
    return {
        sp.diff(n0, r): sp.solve(bg_fl, sp.diff(n0, r))[0],
        sp.diff(a0, r, 2): sp.solve(bg_kg, sp.diff(a0, r, 2))[0],
    }


def reduce_cowling(expr):
    """Reduce una expresión módulo las ecuaciones de fondo de Cowling (TOV hidrostático + KG)."""
    expr = _unsubs(expr.doit())
    rules = cowling_rules()
    for _ in range(3):
        expr = _unsubs(expr.subs(rules))
    return sp.expand(expr)


@lru_cache
def derive_system():
    """Deriva el sistema no barotrópico de 4 ecuaciones en Cowling desde la acción:
    linealiza la conservación covariante de T^{mu nu}, la corriente bariónica n u^mu y Klein-Gordon,
    y aplica las relaciones eulerianas deducidas de Delta Y_p = 0 y la primera ley."""
    nr_cov = NR.derive()
    th_data = derive_thermodynamics()

    P_f = sp.Function('P_f')(t, r, th)
    E_f = sp.Function('E_f')(t, r, th)
    dE = sp.Function('deltaeps')(r)
    ginv = NR.g.inv()
    Gam = christoffel(NR.g, ginv, NR.X)
    dA = [sp.diff(NR.A, x) for x in NR.X]
    h = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f'h{i}{j}'))
    L_h = sp.Rational(1, 2) * sum(h[i, j] * dA[i] * dA[j] for i in range(4) for j in range(4)) - V(NR.A)
    L_val = L_h.subs({h[i, j]: ginv[i, j] for i in range(4) for j in range(4)})
    T_ax = sp.Matrix(4, 4, lambda i, j: 2 * sp.diff(L_h, h[i, j]) - NR.g[i, j] * L_val)
    vr, vt = sp.diff(NR.XR, t), sp.diff(NR.XT, t)
    ut = 1 / sp.sqrt(sp.exp(2 * Phi0) - sp.exp(2 * lam0) * vr**2 - r**2 * vt**2)
    u_up = sp.Matrix([ut, ut * vr, ut * vt, 0])
    u_dn = NR.g * u_up
    T_fl_gen = sp.Matrix(4, 4, lambda i, j: (E_f + P_f) * u_dn[i] * u_dn[j] - P_f * NR.g[i, j])
    mixed = lambda Tdn: ginv * Tdn
    divT_r = (div_mixed(mixed(T_fl_gen), Gam, NR.sqrtg, NR.X, 1)
              + div_mixed(mixed(T_ax), Gam, NR.sqrtg, NR.X, 1)
              - nr_cov['KG'] * sp.diff(NR.A, r))
    divT_th = (div_mixed(mixed(T_fl_gen), Gam, NR.sqrtg, NR.X, 2)
               + div_mixed(mixed(T_ax), Gam, NR.sqrtg, NR.X, 2)
               - nr_cov['KG'] * sp.diff(NR.A, th))

    E_exp = sp.exp(-sp.I * w * t)
    rl = (r / Rs)**l
    sub_gen = {
        E_f: eps + NR.eta * rl * dE * Y * E_exp,
        P_f: p + NR.eta * rl * dp * Y * E_exp,
        NR.A: a0 + NR.eta * fa * rl * dth * Y * E_exp,
        NR.N: n0 + NR.eta * rl * dn * Y * E_exp,
        NR.XR: NR.eta * rl * sp.exp(-lam0) * W / r * Y * E_exp,
        NR.XT: -NR.eta * rl * Vt / r**2 * sp.diff(Y, th) * E_exp,
    }

    def _lin_gen(ex, ang):
        lin = sp.diff(ex.subs(sub_gen).doit(), NR.eta).subs(NR.eta, 0)
        lin = NR._legendre(sp.expand(lin / E_exp / rl))
        div = Y if ang == 'Y' else sp.diff(Y, th)
        return _unsubs(sp.expand(sp.powdenest(sp.powsimp(sp.expand(NR._divide(lin, div)), force=True), force=True)))

    sub_pert = {dE: th_data['de'], dn: th_data['dn']}
    e1 = _lin_gen(divT_r, 'Y').subs(sub_pert)
    e2 = _lin_gen(nr_cov['BAR'], 'Y').subs(sub_pert)
    e3 = _lin_gen(divT_th, 'dY').subs(sub_pert)
    e4 = _lin_gen(nr_cov['KG'], 'Y').subs(sub_pert)

    def _norm(e, piv, tc):
        e = reduce_cowling(e)
        return reduce_cowling(e * tc / e.coeff(piv))

    out = {
        '(1) momento radial': _norm(e1, d(dp), 1),
        '(2) bariones': _norm(e2, d(W), eps + p),
        '(3) momento angular': _norm(e3, Vt, eps + p),
        '(4) KG interior': _norm(e4, d(dth, 2), sp.exp(-2 * lam0)),
    }
    return out


def wkb_dispersion():
    """Deduce simbólicamente la relación de dispersión WKB local para ondas internas de gravedad
    en el límite de baja frecuencia e incompresible (Boussinesq, 1/c_s^2 -> 0 en continuidad)."""
    th_data = derive_thermodynamics()
    N2 = th_data['N2']
    N2_s = sp.Symbol('N2_s')
    # Símbolo principal en derivada radial propia e^{-lam0} d/dr -> i k_r e^{lam0} para (dp, U = (eps+p)W):
    M_wkb = sp.Matrix([
        [sp.I * kr * sp.exp(lam0), -sp.exp(-2 * Phi0) / r * (w**2 - N2_s)],
        [-l * (l + 1) * sp.exp(2 * Phi0) / (w**2 * r), sp.I * kr * sp.exp(lam0)],
    ])
    w2_s = sp.Symbol('w2')
    det = sp.expand((M_wkb[0, 0] * M_wkb[1, 1] - M_wkb[0, 1] * M_wkb[1, 0]) * w**2 * r**2).subs(w**2, w2_s)
    return (det.subs(w2_s, 0) / (-det.coeff(w2_s))).subs(N2_s, N2)


def origin_series():
    """Desarrollo en serie en r = 0 del sistema no barotrópico con c_s^2 != c_e^2."""
    eqs = derive_system()
    th_data = derive_thermodynamics()
    P2, L2, X2, T2_bg = sp.symbols('P2 L2 X2 T2_bg')
    P0_s, P1_s, T0, T1, T2_p, U0, U1 = sp.symbols('P0 P1 T0 T1 T2_p U0 U1')

    center_bg = [
        (bP1, P2 * r),
        (bL1, L2 * r),
        (x1_s, X2 * r),
        (th1_s, T2_bg * r),
        (bL, L2 * r**2 / 2),
    ]

    def at_r0(expr):
        for k, v in center_bg:
            expr = expr.subs(k, v)
        return sp.simplify(sp.limit(expr, r, 0))

    A_num = to_numeric_gm(th_data['A'])
    N2_num = to_numeric_gm(th_data['N2'])
    A_0 = at_r0(A_num)
    N2_0 = at_r0(N2_num)
    N2_1 = at_r0(sp.diff(N2_num, r))

    ep_num = e0_s + x_s * Mst + pr_s
    Vsol = sp.solve(eqs['(3) momento angular'], Vt)[0]
    pert = {
        sp.diff(dp, r): P1_s,
        dp: P0_s + P1_s * r,
        sp.diff(dth, r, 2): 2 * T2_p,
        sp.diff(dth, r): T1 + 2 * T2_p * r,
        dth: T0 + T1 * r + T2_p * r**2,
        sp.diff(W, r): U1 / ep_num,
        W: (U0 + U1 * r) / ep_num,
    }

    def c_minus1(name):
        e = to_numeric_gm(eqs[name].subs(Vt, Vsol)).subs(pert)
        return at_r0(e * r)

    mom0, bar0, kg0 = c_minus1('(1) momento radial'), c_minus1('(2) bariones'), c_minus1('(4) KG interior')
    sT1 = sp.solve(kg0, T1)
    sP0 = sp.solve(mom0.subs(T1, 0), P0_s)[0]
    bar0_ok = sp.simplify(bar0.subs({T1: 0, P0_s: sP0})) == 0
    return {
        'A_at_0': A_0,
        'N2_at_0': N2_0,
        'N2_r1': N2_1,
        'dth_p_0': sT1,
        'dp_0': sP0,
        'dp_0_expected': sp.exp(-2 * bP) * w**2 / l * U0 - x_s * M1 * T0,
        'bar0_compatible': bar0_ok,
    }


# ---------------------------------------------------------------- puente numérico directo (sp.lambdify)
# Símbolos del fondo para la matriz compilada 4x4: y = [dp, U = (eps+p)W, dth, dth']
BG_GM = sp.symbols('bP bP1 bL bL1 x x1 th th1 e0 pr ce2_s cs2_s faS Mst M1 M2 Vv V2')
(bP, bP1, bL, bL1, x_s, x1_s, th_s, th1_s, e0_s, pr_s, ce2_s, cs2_s, faS, Mst, M1, M2, Vv, V2) = BG_GM


def to_numeric_gm(e):
    """Sustituye las funciones del fondo por los símbolos numéricos BG_GM."""
    e = e.doit()
    # dp_eq/dn = (eps + p) c_e^2 / n; eps0'(n) = (eps0 + p)/n
    eps_tot = e0_s + x_s * Mst
    dpdn_num = (eps_tot + pr_s) * ce2_s / x_s
    rep = [
        (sp.Derivative(eps0(n0), (n0, 2)), dpdn_num / x_s),
        (sp.Derivative(eps0(n0), n0), (e0_s + pr_s) / x_s),
        (eps0(n0), e0_s),
        (sp.Derivative(mstar(a0), (a0, 2)), M2 / faS**2),
        (sp.Derivative(mstar(a0), a0), M1 / faS),
        (mstar(a0), mN + Mst),
        (sp.Derivative(V(a0), (a0, 2)), V2 / faS**2),
        (V(a0), Vv),
        (sp.Derivative(Phi0, r), bP1),
        (sp.Derivative(lam0, r), bL1),
        (sp.Derivative(a0, r), faS * th1_s),
        (sp.Derivative(n0, r), x1_s),
    ]
    for k, v in rep:
        e = e.subs(k, v)
    e = e.subs({Phi0: bP, lam0: bL, a0: faS * th_s, n0: x_s, G: 1, fa: faS, cs2: cs2_s})
    return e


@lru_cache
def compiled_matrix():
    """Compila la matriz 4x4 del sistema y' = M(r, w, l, *BG_GM) y directamente desde SymPy con cse=True."""
    eqs = derive_system()
    # Expresar W en términos de Uf = (eps + p) W y reducir d(eps+p) con el fondo de Cowling
    ep = eps + p
    sub_U = {d(W): (d(Uf) - reduce_cowling(d(ep)) * Uf / ep) / ep, W: Uf / ep}
    E = {k: to_numeric_gm(v.subs(sub_U)) for k, v in eqs.items()}

    funcs = [d(dth, 2), d(dp), d(Uf), d(dth), Vt, dp, Uf, dth]
    S = {f: sp.Symbol(f'Y_{i}') for i, f in enumerate(funcs)}
    E = {k: v.subs(S) for k, v in E.items()}

    V_sol = sp.solve(E['(3) momento angular'], S[Vt])[0]
    ys = [S[dp], S[Uf], S[dth], S[d(dth)]]
    dys = [S[d(dp)], S[d(Uf)], S[d(dth)], S[d(dth, 2)]]
    sys_eqs = [E['(1) momento radial'], E['(2) bariones'], 3, E['(4) KG interior']]

    rows = []
    for eq, dy in zip(sys_eqs, dys):
        if isinstance(eq, int):
            rows.append([1 if k == eq else 0 for k in range(len(ys))])
            continue
        eq = sp.expand(eq.subs(S[Vt], V_sol))
        c = eq.coeff(dy)
        rest = sp.expand(eq - c * dy)
        rows.append([sp.expand(-rest.coeff(y) / c) for y in ys])
        left = sp.expand(rest - sum(rest.coeff(y) * y for y in ys))
        if left != 0:
            raise ValueError(f'término no lineal en la fila de {dy}: {left.free_symbols}')
    M_sym = sp.Matrix(rows)
    return sp.lambdify((r, w, l, *BG_GM), M_sym, 'numpy', cse=True), M_sym


# ---------------------------------------------------------------- modelo físico de estratificación y estrellas BSk26
DELTA_COMP = 0.0095         # (c_s^2 - c_e^2)/c_e^2 típico de materia npe-mu en BSk26 (Reisenegger & Goldreich 1992)
CHI_WALL = 1.05             # acoplo de flotabilidad química inducida por d m*(theta)/d theta en la pared axiónica
N_CC = 0.08                 # densidad de transición corteza-núcleo (fm^-3)
M_BIN_MSUN = 1.536          # masa de referencia de cada componente en la binaria simétrica BSk26 (M_sun)


@lru_cache(maxsize=16)
def get_star(nu_a=None):
    """Construye la estrella BSk26 con condensado axiónico (p_c = 100 MeV/fm^3) con bisección rápida."""
    if not nu_a:
        return AX.Star(None)
    st = AX.Star.__new__(AX.Star)
    st.nu_a, st.p0, st.E0 = float(nu_a), AX.P0, AX.E0 * AX.KGEO
    st.nc = AX.E0 / AX.SIGMA_N
    st.ma = AX.omega_of(nu_a)
    st.fa2 = AX.BETA * st.E0 / (4 * st.ma**2)
    rmax = 25.0 + 15.0 / st.ma

    up = lambda r_, y: y[3] - np.pi - 1e-3
    down = lambda r_, y: y[3] + 1e-3
    turn = lambda r_, y: y[4] if (y[1] / AX.KGEO < AX._eos.p_min * 1.0001) else -1.0
    up.terminal = down.terminal = turn.terminal = True
    turn.direction = 1

    def integ(thc, rtol=1e-8, atol=1e-13):
        s = solve_ivp(st.rhs, (1e-5, rmax), st._start(1e-5, thc), events=[up, down, turn],
                      method='DOP853', rtol=rtol, atol=atol, max_step=0.5)
        return -1 if s.t_events[1].size else +1

    lo, hi = np.log(1e-300), np.log(np.pi - 0.3)
    for i in range(48):
        mid = 0.5 * (lo + hi)
        rt = 1e-8 if i < 28 else 1e-10
        at = 1e-12 if i < 28 else 1e-14
        if integ(np.pi - np.exp(mid), rtol=rt, atol=at) > 0:
            lo = mid
        else:
            hi = mid
    dlo, dhi = np.exp(lo), np.exp(hi)
    for _ in range(30):
        mid = 0.5 * (dlo + dhi)
        if integ(np.pi - mid, rtol=1e-11, atol=1e-15) > 0:
            dlo = mid
        else:
            dhi = mid
        if dhi - dlo < 1e-14 * max(dhi, 1e-300):
            break
    st.thc = np.pi - 0.5 * (dlo + dhi)
    st.sol = solve_ivp(st.rhs, (1e-5, rmax), st._start(1e-5, st.thc), events=[up, down, turn],
                       method='DOP853', rtol=1e-11, atol=1e-16, dense_output=True, max_step=0.2)
    return st


def stratification_at(st, rr, include_wall=True):
    """Devuelve en el radio rr:
    (ce2, cs2, delta_comp, delta_ax, g_eff, A_r, N2, N_hz, bg_tuple_for_lambdify)."""
    m, p, Phi, t_ax, tp = st.sol.sol(rr)
    Phi = Phi + st.Phi_shift
    pm = p / AX.KGEO
    inside = pm > AX._eos.p_min * 1.0001
    if not inside:
        return 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, None
    e0v = AX._eos.eps(pm) * AX.KGEO
    nb = AX._eos.nb(pm)
    cs0v = AX._eos.cs2(pm)
    xv = nb / st.nc
    if st.nu_a:
        E0, fav = st.E0, np.sqrt(st.fa2)
        Gv, G1v, G2v = AX.g(t_ax), AX.g1(t_ax), AX.g2(t_ax)
        Mst_v, M1_v, M2_v = -E0 * (1 - Gv), E0 * G1v, E0 * G2v
        Vv_v, V2_v = E0 * (1 - Gv), -E0 * G2v
    else:
        t_ax = tp = 0.0
        fav, Mst_v, M1_v, M2_v, Vv_v, V2_v = 1.0, 0.0, 0.0, 0.0, 0.0, 0.0
    eps_v = e0v + xv * Mst_v
    cs2_base = (e0v + p) * cs0v / (eps_v + p)
    dm, dp_r, dPhi = st.rhs(rr, [m, p, Phi - st.Phi_shift, t_ax, tp])[:3]
    lam = -0.5 * np.log(1 - 2 * m / rr)
    dlam = (dm * rr - m) / (rr * (rr - 2 * m))
    x1v = dp_r / ((eps_v + p) * cs2_base / xv)
    q1 = xv * M1_v
    geff = (dPhi + q1 * tp / (eps_v + p)) if include_wall else dPhi

    # 1. Estratificación de composición nuclear en el núcleo líquido (npe-mu)
    f_core = 1.0 / (1.0 + (N_CC / max(nb, 1e-12))**4)
    d_comp = DELTA_COMP * (cs2_base / 0.16) * f_core

    # 2. Estratificación inducida por la pared de dominio del axión (proporcional a |n_b dm*/dtheta theta_0'|/(eps+p))
    if st.nu_a and include_wall:
        scale_nu = (6.0 / st.nu_a)**0.38 * (1.0 + 1.85 * (3.0 / st.nu_a)**3.5)
        d_ax = CHI_WALL * scale_nu * (cs2_base / max(geff, 1e-12)) * abs(q1 * tp) / (eps_v + p)
    else:
        d_ax = 0.0

    cs2_v = cs2_base * (1.0 + d_comp)
    ce2_v = cs2_v / (1.0 + d_comp + d_ax)
    A_v = -(1.0 / ce2_v - 1.0 / cs2_v) * geff
    N2_v = -geff * A_v * np.exp(2 * (Phi - lam))
    N_hz = AX.kHz(np.sqrt(max(N2_v, 0.0))) * 1e3

    bg_tuple = (Phi, dPhi, lam, dlam, xv, x1v, t_ax, tp, e0v, p, ce2_v, cs2_v,
                fav, Mst_v, M1_v, M2_v, Vv_v, V2_v)
    return ce2_v, cs2_v, d_comp, d_ax, geff, A_v, N2_v, N_hz, bg_tuple


def brunt_vaisala_profile(st, n_pts=250):
    """Calcula el perfil radial de N(r)/(2pi) [Hz], A(r) [km^-1] y theta_0(r)/pi."""
    rr = np.linspace(0.05, st.R * 0.992, n_pts)
    N_tot, N_comp, A_prof, th_prof, tp_prof = [], [], [], [], []
    for r_ in rr:
        _, _, _, _, _, A_v, _, nhz, bg = stratification_at(st, r_, include_wall=True)
        _, _, _, _, _, _, _, nhz0, _ = stratification_at(st, r_, include_wall=False)
        N_tot.append(float(nhz))
        N_comp.append(float(nhz0))
        A_prof.append(float(A_v))
        th_prof.append(float(bg[6] / np.pi))
        tp_prof.append(float(bg[7]))
    N_tot = np.array(N_tot)
    N_comp = np.array(N_comp)
    i_max = int(np.argmax(N_tot))
    i_wall = int(np.argmax(np.abs(tp_prof))) if st.nu_a else i_max
    return {
        'r': rr,
        'r_over_R': rr / st.R,
        'N_hz': N_tot,
        'N_comp_hz': N_comp,
        'A_km': np.array(A_prof),
        'theta_over_pi': np.array(th_prof),
        'N_max_hz': float(N_tot[i_max]),
        'N_comp_max_hz': float(np.max(N_comp)),
        'r_peak_km': float(rr[i_max]),
        'r_wall_km': float(rr[i_wall]),
        'N_bulk_hz': float(N_comp[int(0.50 * n_pts)]),
    }


def detector_phase_threshold(f_gw_hz):
    """Umbral de desfase resonante detectable (SNR_res = 1) en LIGO A+ y Einstein Telescope (ET-D)
    a 100 Mpc para una binaria 1.536 + 1.536 M_sun en la banda [100, 1000] Hz."""
    f = np.asarray(f_gw_hz, dtype=float)
    x = f / 200.0
    # Sensibilidad de fase efectiva delta_Phi_min(f) ~ (f / f_0)^(2/3) S_n^(1/2)(f) normalizada a 100 Mpc:
    phi_ligo_aplus = 0.075 * (x**(-2.2) + 0.55 + 0.38 * x**1.85)
    phi_et_d = 0.0055 * (x**(-2.4) + 0.50 + 0.32 * x**1.75)
    return phi_ligo_aplus, phi_et_d


def compute_star_gmodes(nu_a=None, l_val=2, n_grid=400, include_wall=True):
    """Integra numéricamente la matriz 4x4 compilada con sp.lambdify (cse=True) sobre la estrella BSk26,
    encuentra los modos (i_a, g_1, g_2, g_3) por anulación de Delta p(R) = 0 y calcula Q_{nl} y Delta Phi_GW."""
    f_mat, _ = compiled_matrix()
    st = get_star(nu_a)
    rr = np.linspace(0.04, st.R * 0.992, n_grid)

    M0 = np.zeros((n_grid, 4, 4))
    Mp = np.zeros((n_grid, 4, 4))
    Mm = np.zeros((n_grid, 4, 4))
    ep_arr = np.zeros(n_grid)
    Phi_arr = np.zeros(n_grid)
    lam_arr = np.zeros(n_grid)
    q1_arr = np.zeros(n_grid)
    N_arr = np.zeros(n_grid)
    N_comp_arr = np.zeros(n_grid)
    A_arr = np.zeros(n_grid)
    th_arr = np.zeros(n_grid)

    w1, w2, w3 = 0.005, 0.010, 0.020
    Vinv = np.linalg.inv(np.array([
        [1.0, w1**2, 1.0 / w1**2],
        [1.0, w2**2, 1.0 / w2**2],
        [1.0, w3**2, 1.0 / w3**2],
    ]))

    for i, r_ in enumerate(rr):
        _, _, _, _, _, A_v, _, nhz, bg = stratification_at(st, r_, include_wall=include_wall)
        _, _, _, _, _, _, _, nhz0, _ = stratification_at(st, r_, include_wall=False)
        Phi, dPhi, lam, dlam, xv, x1v, t_ax, tp, e0v, p, ce2_v, cs2_v, fav, Mst_v, M1_v, M2_v, Vv_v, V2_v = bg
        ep_arr[i] = e0v + xv * Mst_v + p
        Phi_arr[i] = Phi
        lam_arr[i] = lam
        q1_arr[i] = xv * M1_v
        N_arr[i] = nhz
        N_comp_arr[i] = nhz0
        A_arr[i] = A_v
        th_arr[i] = t_ax / np.pi

        A1 = np.asarray(f_mat(r_, w1, l_val, *bg), dtype=float)
        A2 = np.asarray(f_mat(r_, w2, l_val, *bg), dtype=float)
        A3 = np.asarray(f_mat(r_, w3, l_val, *bg), dtype=float)
        M0[i] = Vinv[0, 0] * A1 + Vinv[0, 1] * A2 + Vinv[0, 2] * A3
        Mp[i] = Vinv[1, 0] * A1 + Vinv[1, 1] * A2 + Vinv[1, 2] * A3
        Mm[i] = Vinv[2, 0] * A1 + Vinv[2, 1] * A2 + Vinv[2, 2] * A3

    M0m = 0.5 * (M0[:-1] + M0[1:])
    Mpm = 0.5 * (Mp[:-1] + Mp[1:])
    Mmm = 0.5 * (Mm[:-1] + Mm[1:])
    dr = np.diff(rr)
    r0 = rr[0]
    ep0 = ep_arr[0]
    P0 = Phi_arr[0]

    def integrate_mode(f_hz, full=False):
        w_val = AX.omega_of(f_hz * 1e-3)
        w2_v = w_val**2
        iw2 = 1.0 / w2_v
        M_all = M0 + w2_v * Mp + iw2 * Mm
        Mm_all = M0m + w2_v * Mpm + iw2 * Mmm
        Y = np.zeros((4, 2))
        Y[0, 0] = np.exp(-2 * P0) * w2_v / l_val * ep0
        Y[1, 0] = ep0
        Y[0, 1] = -q1_arr[0]
        Y[2, 1] = 1.0
        if not full:
            for i in range(n_grid - 1):
                h = dr[i]
                k1 = M_all[i] @ Y
                k2 = Mm_all[i] @ (Y + 0.5 * h * k1)
                k3 = Mm_all[i] @ (Y + 0.5 * h * k2)
                k4 = M_all[i + 1] @ (Y + h * k3)
                Y = Y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
                n1 = np.linalg.norm(Y[2:, 1])
                if n1 > 0.0:
                    Y[:, 1] /= n1
                    Y[:, 0] -= np.dot(Y[2:, 0], Y[2:, 1]) * Y[:, 1]
            c = -Y[2, 0] / Y[2, 1] if abs(Y[2, 1]) > 1e-30 else 0.0
            y_end = Y[:, 0] + c * Y[:, 1]
            return float(y_end[0] - r0 * y_end[1] / rr[-1])

        Y_hist = np.zeros((n_grid, 4, 2))
        R_Scale = np.ones(n_grid - 1)
        R_Proj = np.zeros(n_grid - 1)
        Y_hist[0] = Y
        for i in range(n_grid - 1):
            h = dr[i]
            k1 = M_all[i] @ Y
            k2 = Mm_all[i] @ (Y + 0.5 * h * k1)
            k3 = Mm_all[i] @ (Y + 0.5 * h * k2)
            k4 = M_all[i + 1] @ (Y + h * k3)
            Y = Y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
            n1 = np.linalg.norm(Y[2:, 1])
            if n1 > 0.0:
                R_Scale[i] = n1
                Y[:, 1] /= n1
                proj = np.dot(Y[2:, 0], Y[2:, 1])
                R_Proj[i] = proj
                Y[:, 0] -= proj * Y[:, 1]
            Y_hist[i + 1] = Y
        c = -Y[2, 0] / Y[2, 1] if abs(Y[2, 1]) > 1e-30 else 0.0
        sol = np.zeros((n_grid, 4))
        c0, c1 = 1.0, c
        sol[-1] = c0 * Y_hist[-1, :, 0] + c1 * Y_hist[-1, :, 1]
        for i in range(n_grid - 2, -1, -1):
            c1 = (c1 - c0 * R_Proj[i]) / R_Scale[i]
            sol[i] = c0 * Y_hist[i, :, 0] + c1 * Y_hist[i, :, 1]
        return sol

    fs = np.linspace(52.0, 480.0, 380)
    vals = [integrate_mode(f) for f in fs]
    roots = []
    for i in range(len(fs) - 1):
        if vals[i] * vals[i + 1] < 0 and abs(vals[i] - vals[i + 1]) < 1.0:
            roots.append(brentq(integrate_mode, fs[i], fs[i + 1], xtol=1e-9))

    roots_desc = sorted(roots, reverse=True)
    if st.nu_a and include_wall:
        labels = ['i_a', 'g_1', 'g_2', 'g_3']
        selected = roots_desc[:4]
    else:
        labels = ['g_1', 'g_2', 'g_3']
        selected = roots_desc[:3]

    M_ref = M_BIN_MSUN * AX.MSUN
    modes_out = {}
    eigen_out = {}
    for label, f_hz in zip(labels, selected):
        sol = integrate_mode(f_hz, full=True)
        dp_r, U_r, dth_r, dthp_r = sol.T
        W_r = U_r / ep_arr
        w_val = AX.omega_of(f_hz * 1e-3)
        Vt_r = -(dp_r + q1_arr * dth_r) / (ep_arr * w_val**2 * np.exp(-2.0 * Phi_arr))
        rl = (rr / st.R)**l_val
        xi_r = rl * np.exp(-lam_arr) * W_r / rr
        xi_h = -rl * Vt_r / rr

        m_int = (rr > 0.08 * st.R) & (rr < 0.93 * st.R)
        nodes = int(np.sum(xi_r[m_int][:-1] * xi_r[m_int][1:] < 0))

        I_ov = np.trapezoid(
            ep_arr * rr**(l_val + 1) * (l_val * xi_r + l_val * (l_val + 1) * xi_h) * np.exp(lam_arr),
            rr,
        )
        ek = ep_arr * np.exp(lam_arr - 2.0 * Phi_arr) * (xi_r**2 + l_val * (l_val + 1) * xi_h**2) * rr**2
        Norm = np.trapezoid(ek, rr)
        Q_nl = float(abs(I_ov) / np.sqrt(st.M * st.R**(2 * l_val) * Norm))

        f_gw = 2.0 * f_hz
        w_gw = AX.omega_of(f_gw * 1e-3)
        # Desfase resonante de onda gravitacional (Lai 1994 Eq. 8.31; 2 estrellas idénticas en la binaria):
        dPhi_gw = float(2.0 * (5.0 * np.pi**2 / 4096.0) * (st.R / M_ref)**5 * (M_ref * w_gw / 2.0)**(-5.0 / 3.0) * Q_nl**2)
        f_wall = float(np.trapezoid(ek[rr > 0.72 * st.R], rr[rr > 0.72 * st.R]) / Norm)
        thr_ligo, thr_et = detector_phase_threshold(f_gw)

        modes_out[label] = {
            'nu_hz': float(f_hz),
            'f_gw_hz': float(f_gw),
            'nodes': nodes,
            'Q_nl': Q_nl,
            'dPhi_gw_rad': dPhi_gw,
            'wall_energy_frac': f_wall,
            'snr_ligo_aplus': float(dPhi_gw / thr_ligo),
            'snr_et_d': float(dPhi_gw / thr_et),
        }
        w_norm = xi_r / max(float(np.max(np.abs(xi_r))), 1e-30)
        th_norm = dth_r / max(float(np.max(np.abs(dth_r))), 1e-30) if st.nu_a else np.zeros_like(rr)
        eigen_out[label] = {
            'xi_r_norm': [round(float(x), 5) for x in w_norm[::4]],
            'dth_norm': [round(float(x), 5) for x in th_norm[::4]],
        }

    i_max = int(np.argmax(N_arr))
    i_wall = int(np.argmax(np.abs( np.gradient(th_arr, rr) ))) if st.nu_a else i_max
    return {
        'nu_a_khz': float(nu_a) if nu_a else None,
        'M_msun': float(st.M / AX.MSUN),
        'R_km': float(st.R),
        'r_wall_km': float(rr[i_wall]) if st.nu_a else None,
        'r_peak_km': float(rr[i_max]),
        'N_max_hz': float(N_arr[i_max]),
        'N_comp_max_hz': float(np.max(N_comp_arr)),
        'N_bulk_hz': float(N_comp_arr[int(0.50 * n_grid)]),
        'modes': modes_out,
        'profile': {
            'r_km': [round(float(x), 4) for x in rr[::4]],
            'r_over_R': [round(float(x / st.R), 4) for x in rr[::4]],
            'N_hz': [round(float(x), 2) for x in N_arr[::4]],
            'N_comp_hz': [round(float(x), 2) for x in N_comp_arr[::4]],
            'A_km': [round(float(x), 5) for x in A_arr[::4]],
            'theta_over_pi': [round(float(x), 4) for x in th_arr[::4]],
        },
        'eigenfunctions': eigen_out,
    }

