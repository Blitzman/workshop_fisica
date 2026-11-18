"""N4: the radial (l = 0) system with the axion condensate, as numerical matrices, from D0's and D1's equation artifacts.

Port of axion_ns/radial.py (reduced_equations, first_order_system, interior_matrix, exterior_matrix). Nothing is typed by hand: D1's radial_linear.eq.json (the
linearised field, fluid, baryon and Klein-Gordon equations) and D0's background.eq.json are reduced exactly as the legacy did with `pipeline.linear_equations` and
`pipeline.background_equations`:

    y = (xi, DP, dtheta, dtheta')   with DP = Delta p = delta p + xi p_0'(r)  (Lagrangian pressure); delta n = (DP - xi T)/(dp/dn), T = p_0'(r) from the background fluid equation
    algebraic constraints: delta lambda (G_tr) and delta Phi' (G_rr); delta a = f_a dtheta;  interior 4 x 4 (12 non-zero entries), exterior (dtheta, dtheta') 2 x 2 (3 non-zero)

Both are compiled with lambdify(cse=True) as FLAT lists of entries (16 and 4). A depends on omega only through omega^2 and, for this system, A = A_0(r) + omega^2 A_1(r) exactly
(residual 2.5e-16; no 1/omega^2 term): tables.py exploits it. `to_numeric` and the background symbols are the same as N3's (copied: modules never import each other).
"""
import sympy as sp

from common import artifacts as A
from common.modulo import Rules
from common.notation import (r, w, G, V, mstar, mN, eps0, a0, n0, Phi0, lam0, d, fa, dPhi, dlam, da1, dn1, xi1, p_of)


# background symbols, in this order in the numerical functions (as N3)
BG = sp.symbols('bP bP1 bP2 bL bL1 x x1 th th1 e0 pr cs0 fa Mst M1 M2 Vv V2')
(bP, bP1, bP2, bL, bL1, x, x1, th, th1, e0, pr, cs0, faS, Mst, M1, M2, Vv, V2) = BG


def to_numeric(e):
    """Background functions -> symbols. eps0'(n) = (eps0 + p)/n (first law) and eps0''(n) = (dp/dn)/n with dp/dn = (eps0 + p) c_s0^2 / n."""
    e = e.doit()
    dpdn = (e0 + pr) * cs0 / x
    rep = [(sp.Derivative(eps0(n0), (n0, 2)), dpdn / x), (sp.Derivative(eps0(n0), n0), (e0 + pr) / x), (eps0(n0), e0),
           (sp.Derivative(mstar(a0), (a0, 2)), M2 / faS**2), (sp.Derivative(mstar(a0), a0), M1 / faS), (mstar(a0), mN + Mst),
           (sp.Derivative(V(a0), (a0, 2)), V2 / faS**2), (V(a0), Vv),
           (sp.Derivative(Phi0, (r, 2)), bP2), (sp.Derivative(Phi0, r), bP1), (sp.Derivative(lam0, r), bL1),
           (sp.Derivative(a0, r), faS * th1), (sp.Derivative(n0, r), x1)]
    for k, v in rep:
        e = e.subs(k, v)
    return e.subs({Phi0: bP, lam0: bL, a0: faS * th, n0: x, G: 1, fa: faS})


def _states(funcs):
    return {f: sp.Symbol(f'S{i}') for i, f in enumerate(funcs)}


DP = sp.Function('DeltaP')(r)
dth = sp.Function('dtheta')(r)
_nn = sp.Symbol('nn', positive=True)


def load_inputs(d0_dir, d1_dir):
    """(linear equations of D1, background equations of D0, input records)."""
    pb, pl = f'{d0_dir}/background.eq.json', f'{d1_dir}/radial_linear.eq.json'
    b = A.load_equations(pb, expect_name='background')
    l = A.load_equations(pl, expect_name='radial_linear')
    return dict(l.eq), dict(b.eq), [A.input_record(pb), A.input_record(pl)]


def reduced_equations(L, B, rules):
    """Linear equations in (xi, DP, dth) with the background already reduced. Returns (eqs, pn, T)."""
    pn = sp.diff(p_of(_nn), _nn).subs(_nn, n0)                      # dp/dn at n0
    n0p = sp.solve(B['FLUID_r'], d(n0))[0]                           # T = p0'(r) from the background fluid equation
    T = sp.expand(pn * n0p)
    a0pp = sp.solve(B['KG'], d(a0, 2))[0]
    dn_of = (DP - xi1 * T) / pn
    sub = {da1: fa * dth, dn1: dn_of}
    out = {}
    for k in ('Gtr', 'Grr', 'KG', 'FLUID_r', 'BAR'):
        e = L[k].subs(sub).doit()
        e = e.subs(d(a0, 3), d(a0pp)).subs(d(a0, 2), a0pp).doit()
        e = sp.expand(e.subs(d(a0, 2), a0pp))                        # the one that reintroduces d(a0pp)
        # delta Phi (undifferentiated) only appears in baryon conservation and its coefficient is zero modulo the background
        c = e.coeff(dPhi)
        if c != 0:
            if sp.simplify(rules.reduce(c)) != 0:
                raise ValueError(f'{k}: undifferentiated delta Phi with a non-zero coefficient')
            e = sp.expand(e - c * dPhi)
        out[k] = e
    return out, pn, T


def first_order_system(L, B, rules):
    """First-order system (xi', DP', dth'') in (xi, DP, dth, dth'), all symbolic."""
    E, _, _ = reduced_equations(L, B, rules)
    a0pp = sp.solve(B['KG'], d(a0, 2))[0]
    fix = lambda e: sp.expand(e.doit().subs(d(a0, 2), a0pp))
    dlam_s = sp.solve(E['Gtr'], dlam)[0]
    rep1 = lambda e: fix(e.subs(d(dlam), d(dlam_s)).subs(dlam, dlam_s))
    dPhip_s = sp.solve(rep1(E['Grr']), d(dPhi))[0]
    xip_s = sp.solve(rep1(E['BAR']), d(xi1))[0]
    rep2 = lambda e: sp.expand(rep1(e).subs(d(dPhi), dPhip_s).subs(d(xi1), xip_s))
    DPp_s = sp.solve(rep2(E['FLUID_r']), d(DP))[0]
    thpp_s = sp.solve(rep2(E['KG']), d(dth, 2))[0]
    return {'xi': xip_s, 'DP': DPp_s, 'thpp': thpp_s}


def _to_num(e):
    # V'(a0) = -M1/f_a in this model: V = E0 (1 - g) and n (m* - m_N) = -x E0 (1 - g). First V'' (otherwise substituting V' inside V'' = d(V')/da gives 0).
    e = e.doit().subs(sp.Derivative(V(a0), (a0, 2)), V2 / faS**2)
    e = e.subs(sp.Derivative(V(a0), a0), -M1 / faS)
    return sp.expand(to_numeric(e))


def interior_matrix_symbolic(L, B, rules):
    F = first_order_system(L, B, rules)
    S = _states([d(dth), xi1, DP, dth])
    ys = [S[xi1], S[DP], S[dth], S[d(dth)]]
    rows = []
    for key in ('xi', 'DP', None, 'thpp'):
        if key is None:
            rows.append([0, 0, 0, 1])
            continue
        e = _to_num(F[key]).subs(S)
        left = sp.expand(e - sum(e.coeff(y) * y for y in ys))
        if left != 0 or e.atoms(sp.Derivative):
            raise ValueError(f'{key}: leftover {left.free_symbols} / {e.atoms(sp.Derivative)}')
        rows.append([sp.expand(e.coeff(y)) for y in ys])
    return sp.Matrix(rows)


def exterior_matrix_symbolic(L, B, rules):
    """Outside the star (no fluid: xi = DP = 0, x = e0 = p = 0) but with the tail of the condensate: delta lambda and delta Phi' stay coupled to the axion. y = (dth, dth')."""
    E, _, _ = reduced_equations(L, B, rules)
    a0pp = sp.solve(B['KG'], d(a0, 2))[0]
    fix = lambda e: sp.expand(e.doit().subs({xi1: 0, DP: 0}).doit().subs(d(a0, 2), a0pp))
    Ex = {k: fix(v) for k, v in E.items()}
    dlam_s = sp.solve(Ex['Gtr'], dlam)[0]
    rep1 = lambda e: sp.expand(e.subs(d(dlam), d(dlam_s)).subs(dlam, dlam_s).doit().subs(d(a0, 2), a0pp))
    dPhip_s = sp.solve(rep1(Ex['Grr']), d(dPhi))[0]
    thpp = sp.solve(sp.expand(rep1(Ex['KG']).subs(d(dPhi), dPhip_s)), d(dth, 2))[0]
    zero = {x: 0, e0: 0, pr: 0, x1: 0}
    S = _states([d(dth), dth])
    e = sp.expand(_to_num(thpp).subs(zero).subs(S))
    ys = [S[dth], S[d(dth)]]
    left = sp.expand(e - sum(e.coeff(y) * y for y in ys))
    if left != 0 or e.atoms(sp.Derivative):
        raise ValueError(f'exterior: leftover {left.free_symbols} / {e.atoms(sp.Derivative)}')
    return sp.Matrix([[0, 1], [sp.expand(e.coeff(ys[0])), sp.expand(e.coeff(ys[1]))]])


def compile_flat(L, B):
    """(interior_flat, exterior_flat): functions (r, w, *background) -> list of 16 / 4 entries (some scalars)."""
    rules = Rules(B)
    Ai, Ae = interior_matrix_symbolic(L, B, rules), exterior_matrix_symbolic(L, B, rules)
    args = (r, w, *BG)
    return sp.lambdify(args, list(Ai), 'numpy', cse=True), sp.lambdify(args, list(Ae), 'numpy', cse=True), Ai, Ae
