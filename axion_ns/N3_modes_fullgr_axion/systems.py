"""N3: the extended Lindblom-Detweiler system with the axion, as numerical matrices, from D2's equation artifact.

Port of the symbolic half of axion_ns/axgr.py. D2's fullgr_extended.eq.json (EXT: LD + axion, verified against the action) is turned into
y' = A(r, w, l, background) y by replacing the background functions with symbols (to_numeric) and solving the algebraic constraints:

    interior: y = (H1, K, W, X, dtheta, dtheta')   (H0 and V algebraic)           -> 6 x 6, 31 non-zero entries
    exterior: y = (H1, K, dtheta, dtheta')         (no fluid, the axion tail stays coupled) -> 4 x 4

Both are compiled with lambdify(cse=True) as FLAT lists of entries (36 and 16): tables.py evaluates them vectorised at the table nodes. The matrices depend on
omega only through omega^2, A = A_-1/w^2 + A_0 + w^2 A_1 (measured, residual 4e-15; exact at complex w to 1.5e-16), which tables.py exploits.
"""
from pathlib import Path

import sympy as sp

from common import artifacts as A
from common import compiled as C
from common.notation import (r, w, l, G, V, mstar, mN, eps0, dth, a0, n0, Phi0, lam0, d, fa, Xf, X_def, dn, Vt, H0, H1, K, W, eps, p, cs2, dm, dmm,
                             meff2, thp)

# background symbols, in this order in the numerical functions
BG = sp.symbols('bP bP1 bP2 bL bL1 x x1 th th1 e0 pr cs0 fa Mst M1 M2 Vv V2')
(bP, bP1, bP2, bL, bL1, x, x1, th, th1, e0, pr, cs0, faS, Mst, M1, M2, Vv, V2) = BG


def load_ext(d2_dir):
    """The extended system from D2's artifact: dict key -> expression, and the input record."""
    path = f'{d2_dir}/fullgr_extended.eq.json'
    art = A.load_equations(path, expect_name='fullgr_extended')
    return dict(art.eq), A.input_record(path)


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


def _linear(eqs, ys, dys, alg):
    """Matrix A with dys = A ys from linear equations; `alg` = [(eq, var)] algebraic constraints, in order."""
    sol = {}
    for eq, var in alg:
        sol[var] = sp.solve(eq.subs(sol), var)[0]
    rows = []
    for eq, dy in zip(eqs, dys):
        if isinstance(eq, int):                     # identity row: (y_j)' = y_eq
            rows.append([1 if k == eq else 0 for k in range(len(ys))])
            continue
        eq = sp.expand(eq.subs(sol))
        c = eq.coeff(dy)
        rest = sp.expand(eq - c * dy)
        rows.append([sp.expand(-rest.coeff(y) / c) for y in ys])
        left = sp.expand(rest - sum(rest.coeff(y) * y for y in ys))
        if left != 0:
            raise ValueError(f'nonlinear term or unsolved variable: {left.free_symbols}')
    return sp.Matrix(rows)


def _states(funcs):
    return {f: sp.Symbol(f'S{i}') for i, f in enumerate(funcs)}


def interior_matrix_symbolic(ext):
    dn_of = sp.solve(Xf - X_def, dn)[0]
    E = {k: to_numeric(v.subs(dn, dn_of)) for k, v in ext.items()}
    S = _states([d(dth, 2), d(H1), d(K), d(W), d(Xf), d(dth), H0, Vt, H1, K, W, Xf, dth])
    E = {k: v.subs(S) for k, v in E.items()}
    ys = [S[H1], S[K], S[W], S[Xf], S[dth], S[d(dth)]]
    dys = [S[d(H1)], S[d(K)], S[d(W)], S[d(Xf)], S[d(dth)], S[d(dth, 2)]]
    eqs = [E["H1'"], E["K'"], E["W'"], E["X'"], 5, E['KG']]
    return _linear(eqs, ys, dys, [(E['H0'], S[H0]), (E['V'], S[Vt])])


def exterior_matrix_symbolic(ext):
    zero = {x: 0, e0: 0, pr: 0, x1: 0}
    E = {}
    for k in ("H1'", "K'", 'KG', 'H0'):
        e = to_numeric(ext[k].subs({W: 0, Xf: 0, Vt: 0, dn: 0}).doit())
        E[k] = sp.expand(e.subs(zero))
    S = _states([d(dth, 2), d(H1), d(K), d(dth), H0, H1, K, dth])
    E = {k: v.subs(S) for k, v in E.items()}
    ys = [S[H1], S[K], S[dth], S[d(dth)]]
    dys = [S[d(H1)], S[d(K)], S[d(dth)], S[d(dth, 2)]]
    eqs = [E["H1'"], E["K'"], 3, E['KG']]
    return _linear(eqs, ys, dys, [(E['H0'], S[H0])])


def compile_flat(ext):
    """(interior_flat, exterior_flat): functions (r, w, l, *background) -> list of 36 / 16 entries (some scalars). From the compiled cache
    (common/compiled.py, compiled/extended.json; --recompile rebuilds it)."""
    here = Path(__file__).resolve().parent
    args = (r, w, l, *BG)
    key = {'equations': C.expr_hash(ext), 'code': C.code_hash(here / 'systems.py', here.parent / 'common' / 'notation.py')}
    f = C.compiled(here / 'compiled' / 'extended.json', key,
                   lambda: {'interior': (args, list(interior_matrix_symbolic(ext)), 'numpy'),
                            'exterior': (args, list(exterior_matrix_symbolic(ext)), 'numpy')})
    return f['interior'], f['exterior']


def substitution_checks():
    """The axion substitutions: quantities of Paper II converted to the numerical symbols against what is expected."""
    eps_ = e0 + x * Mst
    cases = [('n dm*/dtheta', n0 * dm, x * M1), ('n d2m*/dtheta2', n0 * dmm, x * M2), ('m_eff^2', meff2, (V2 + x * M2) / faS**2),
             ('eps(n, theta)', eps, eps_), ('p(n)', p, pr), ('c_s^2', cs2, (e0 + pr) * cs0 / (eps_ + pr)),
             ("theta_0'", thp, th1), ('V(a_0)', V(a0), Vv)]
    return [(name, sp.simplify(to_numeric(e) - want) == 0) for name, e, want in cases]
