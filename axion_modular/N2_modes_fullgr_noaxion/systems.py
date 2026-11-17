"""N2: numerical perturbation systems compiled from the equation ARTIFACTS (not from sibling-module code).

Reads: D2's fullgr_extended.eq.json (the extended Lindblom-Detweiler system, verified against the action) and R's
paper2_equations.eq.json ((8a)-(8c) of Paper II). With the axion switched off and a concrete polytrope EOS, the algebraic
constraints are solved symbolically (as axion_ns/fmodes.py does) and the matrix A(r, w, l, background) in y' = A y is compiled with
lambdify(cse=True). Two evaluators per system: `mat` (4x4 / 2x2 array) and `ay` (the product A y, as a list of floats, fastest).
A depends on w only through w^2: A = A_-1/w^2 + A_0 + w^2 A_1 (measured, residual 1e-15), not exploited here.
"""
import numpy as np
import sympy as sp

from common import artifacts as A
from common.notation import (r, w, l, G, V, mstar, mN, eps0, dth, a0, n0, Phi0, lam0, d, Xf, X_def, dn, Vt, H0, H1, K, W)

from .polytrope import EOS, _x, eos_num

_B = sp.symbols('bn bn1 bn2 bP bP1 bP2 bL bL1')     # n, n', n'', Phi, Phi', Phi'', lambda, lambda'
EXT_KEYS = ("H1'", "K'", "W'", "X'", 'H0', 'V')
COWLING_KEYS = ('(8a) momento radial', '(8b) bariones', '(8c) momento angular')


class Sys:
    def __init__(self, mat, ay, n):
        self.mat, self.ay, self.n = mat, ay, n


def load_inputs(d2_dir, r_dir):
    """Returns (EXT dict, Paper II dict, input records)."""
    ep, pp = f'{d2_dir}/fullgr_extended.eq.json', f'{r_dir}/paper2_equations.eq.json'
    ext = A.load_equations(ep, expect_name='fullgr_extended')
    p2 = A.load_equations(pp, expect_name='paper2_equations')
    return dict(ext.eq), dict(p2.eq), [A.input_record(ep), A.input_record(pp)]


def _numeric(e, kind):
    """Without axion, with the concrete EOS, and the background functions replaced by symbols."""
    e = e.subs(dth, 0).subs(a0, 0).doit()
    e = e.replace(lambda x: isinstance(x, sp.Subs), lambda x: 0)
    e = e.subs({mstar(0): mN, V(0): 0}).subs(eps0, sp.Lambda(_x, EOS[kind])).doit()
    bn, bn1, bn2, bP, bP1, bP2, bL, bL1 = _B
    for k, v in [(d(Phi0, 2), bP2), (d(n0, 2), bn2), (d(Phi0), bP1), (d(lam0), bL1), (d(n0), bn1)]:
        e = e.subs(k, v)
    return e.subs({Phi0: bP, lam0: bL, n0: bn, G: 1})


def _compile(Amat, n):
    args = (r, w, l, *_B)
    mat = sp.lambdify(args, Amat, 'numpy', cse=True)
    ys = sp.symbols(f'y0:{n}')
    ay = sp.lambdify(list(args) + list(ys), list(Amat * sp.Matrix(ys)), 'math', cse=True)
    return Sys(lambda *a: np.asarray(mat(*a), dtype=float), ay, n)


def _matrix(eqs, ys, alg, kind):
    """y' = A(r, w, l, background) y for the system `eqs` (derivatives of `ys`, algebraic constraints `alg`)."""
    eqs = {k: _numeric(v, kind) for k, v in eqs.items()}
    sol = {}
    for name, var in alg:                       # algebraic constraints, in order
        sol[var] = sp.solve(eqs[name].subs(sol), var)[0]
    ders = [d(y) for y in ys]
    rest = [eqs[k].subs(sol) for k in eqs if k not in dict(alg)]
    rhs = sp.solve(rest, ders, dict=True)[0]
    Amat = sp.Matrix([[sp.factor(sp.cancel(sp.diff(rhs[dy], y))) for y in ys] for dy in ders])
    return Amat


def _zl_from(sysm, kind):
    """Zhao & Lattimer as printed: adds the sign error of 2 w^2 e^{lambda-2Phi} W in X' to the compiled matrix."""
    eos = eos_num(kind)[0]

    def mat_zl(rr, omega, ll, n, n1, n2, P, P1, P2, L, L1):
        M = np.array(sysm.mat(rr, omega, ll, n, n1, n2, P, P1, P2, L, L1), dtype=float)
        e, p = eos(n)[:2]
        M[3, 2] += 2 * omega**2 * (e + p) * np.exp(L - P) / rr
        return M
    return Sys(mat_zl, None, 4)


def full_systems(ext, kind):
    """(EXT, ZL): the extended system without axion and the ZL variant, from ONE symbolic solve (the legacy code solved twice)."""
    Amat = _matrix({k: ext[k] for k in EXT_KEYS}, [H1, K, W, Xf], [('H0', H0), ('V', Vt)], kind)
    sysm = _compile(Amat, 4)
    return sysm, _zl_from(sysm, kind)


def full_system(ext, kind, which='EXT'):
    """LD in variables (H1, K, W, X). 'EXT': the extended system without axion. 'ZL': Zhao & Lattimer as printed (differs only in
    the sign of 2 w^2 e^{lambda-2Phi} W inside X')."""
    sysm, zl = full_systems(ext, kind)
    return sysm if which == 'EXT' else zl


def cowling_system(paper2, kind):
    """(8a)-(8c) of Paper II without axion, in variables (W, X) with X = -e^Phi Delta p."""
    dn_of = sp.solve(Xf - X_def, dn)[0]
    eqs = {k: paper2[k].subs(dn, dn_of).doit() for k in COWLING_KEYS}
    return _compile(_matrix(eqs, [W, Xf], [('(8c) momento angular', Vt)], kind), 2)
