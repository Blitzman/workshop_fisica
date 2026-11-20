"""T: numerical bridge between the derived equation (D1) and the hand-coded Eq. (4.13) used here.

toy.f_toy and branches.f_k are typed by hand. D1 publishes the matching determinant derived from Eqs. (3.3)
(`f_4_13_tan` in toy_model.eq.json, symbols k_p, k_m, k_a, k_out, x_0). This module compiles that expression with
lambdify (mpmath) and compares it with both hand-coded forms at random complex frequencies. A mismatch means the
numerics no longer implement the verified equation.
"""
import mpmath as mp
import numpy as np
import sympy as sp

from common import artifacts as A

from .branches import f_k
from .toy import f_toy


def _kvars(w, lam2, alpha, cs2, mass2_in):
    ka2, ks2 = w**2 - mass2_in, w**2 / cs2
    Aa = alpha**2 * lam2
    sD = mp.sqrt((ka2 - ks2 + Aa)**2 + 4 * Aa * ks2)
    kp, km = mp.sqrt((ka2 + ks2 + Aa + sD) / 2), mp.sqrt((ka2 + ks2 + Aa - sD) / 2)
    return kp, km, mp.sqrt(ka2)


def bridge_checks(eq_path, n_points=24, seed=7):
    art = A.load_equations(eq_path, keys=['f_4_13_tan'], expect_name='toy_model')
    f_tan = art['f_4_13_tan']
    sy = {s.name: s for s in f_tan.free_symbols}
    core = ['k_p', 'k_m', 'k_a', 'k_out', 'x_0']
    extra = [k for k in ('alpha', 'lambda2') if k in sy]       # they cancel analytically in the stored (unsimplified) determinant
    unknown = set(sy) - set(core) - set(extra)
    if not all(k in sy for k in core) or unknown:
        return [('Puente D1 -> T: la expresión derivada no contiene los símbolos esperados', False, str(sorted(sy)))]
    fn = sp.lambdify([sy[k] for k in core + extra], f_tan, modules='mpmath')
    rng = np.random.default_rng(seed)
    e_tan, e_ent = 0.0, 0.0
    for _ in range(n_points):
        lam2, alpha, cs2, x0 = 10**rng.uniform(-4, -2), rng.uniform(2, 12), rng.uniform(0.15, 0.6), rng.uniform(0.3, 1.2)
        mass2 = 1 - alpha
        w = mp.mpc(rng.uniform(0.5, 6.0), -rng.uniform(0.001, 1.0))
        kp, km, kq = _kvars(w, lam2, alpha, cs2, mass2)
        kout = w * mp.sqrt(1 - 1 / w**2)
        mine = f_toy(w, lam2, alpha, cs2, x0, entire=False)
        deriv = fn(kp, km, kq, kout, x0, *[{'alpha': alpha, 'lambda2': lam2}[k] for k in extra])
        e_tan = max(e_tan, float(abs(mine - deriv) / abs(deriv)))
        # entire form in the k variable (branches.f_k): w^2 = 1 + k^2, k_out = k
        k = mp.mpc(rng.uniform(0.1, 3.0), -rng.uniform(0.001, 0.5))
        w2 = 1 + k**2
        kp2, km2, kq2 = _kvars(mp.sqrt(w2), lam2, alpha, cs2, mass2)
        mine_k = f_k(k, lam2, alpha, cs2, x0)
        deriv_k = fn(kp2, km2, kq2, k, x0, *[{'alpha': alpha, 'lambda2': lam2}[k_] for k_ in extra]) * mp.cos(kp2 * x0) * mp.cos(km2 * x0)
        e_ent = max(e_ent, float(abs(mine_k - deriv_k) / abs(deriv_k)))
    return [
        (f'Puente D1 -> T: f_toy (Ec. 4.13, forma con tangentes) coincide con la f(omega) derivada en D1 '
         f'({n_points} puntos aleatorios, diferencia relativa máx. {e_tan:.1e})', e_tan < 1e-15, 'lambdify(mpmath) de f_4_13_tan'),
        (f'Puente D1 -> T: f_k (forma entera en k, branches.py) = f_4_13_tan x cos(k_p x0) cos(k_m x0) '
         f'(diferencia relativa máx. {e_ent:.1e})', e_ent < 1e-15, ''),
    ]
