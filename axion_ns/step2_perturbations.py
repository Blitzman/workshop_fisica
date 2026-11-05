"""Paso 2: perturbaciones radiales lineales y comparación con Ec. (3.3) y (3.7)."""
import sympy as sp

from . import paper
from .compare import check_modulo, is_zero
from .model import r, G, w, V, mstar, Phi0, lam0, a0, n0, dPhi, dlam, da1, dn1, xi1
from .pipeline import linear_equations, background_rules, report

d = lambda f, k=1: sp.diff(f, r, k)

# ecuación derivada que corresponde a cada una del paper, y variable para eliminar
MATCH = {
    '(3.3a) KG': ('KG', d(da1, 2)),
    '(3.3b) rr': ('Grr', d(dPhi)),
    '(3.3c) tr': ('Gtr', dlam),
    '(3.3d) Euler': ('FLUID_r', xi1),
    '(3.3e) continuidad': ('BAR', dn1),
}


def run():
    L, R = linear_equations(), background_rules()
    res = []
    for name, (der, piv) in MATCH.items():
        ok, _ = check_modulo(paper.PERTURBATIONS[name], L[der], piv, R)
        res.append((f'Paper {name} == linealización de la acción', ok, ''))

    # control negativo: errores plantados deben detectarse
    bad = {
        '(3.3b) con 8 pi en vez de 4 pi': ('(3.3b) rr', lambda e: e.subs(sp.pi, 2 * sp.pi)),
        '(3.3a) con signo cambiado en el término 2 delta lambda': ('(3.3a) KG', lambda e: e + 4 * dlam * sp.diff(V(a0), a0)),
        '(3.3e) sin el término delta lambda': ('(3.3e) continuidad', lambda e: e.subs(dlam, 0)),
    }
    for label, (name, mutate) in bad.items():
        der, piv = MATCH[name]
        ok, _ = check_modulo(mutate(paper.PERTURBATIONS[name]), L[der], piv, R)
        res.append((f'Control negativo: {label} -> se detecta', not ok, ''))

    # redundancia: la componente tt linealizada es consecuencia de tr (identidad de Bianchi)
    # (primero delta n de la conservación bariónica, que reintroduce delta lambda; luego delta lambda de G_tr)
    sol_dn = sp.solve(L['BAR'], dn1)[0]
    tt = L['Gtt'].subs({d(dn1): d(sol_dn)}).subs({dn1: sol_dn}).doit()
    sol_dlam = sp.solve(L['Gtr'], dlam)[0]
    tt = R.reduce(tt.subs({d(dlam): d(sol_dlam)}).subs({dlam: sol_dlam}).doit())
    res.append(('G_tt linealizada es redundante (Bianchi + conservación bariónica)', is_zero(tt), ''))

    # exterior: Schwarzschild + axión en vacío -> Ec. (3.7) del paper
    Mc, ma = sp.symbols('M m_a', positive=True)
    x = sp.Symbol('x')
    f = 1 - 2 * G * Mc / r
    ext = paper.PERTURBATIONS['(3.3a) KG'].subs({V: sp.Lambda(x, ma**2 * x**2 / 2), mstar: sp.Lambda(x, 0)}).doit()
    ext = ext.subs({n0: 0, a0: 0, Phi0: sp.log(f) / 2, lam0: -sp.log(f) / 2}).doit()
    eq11 = f * d(da1, 2) + 2 / r * (1 - G * Mc / r) * d(da1) - (ma**2 - w**2 / f) * da1
    res.append(('Exterior: (3.3a) -> KG en Schwarzschild, Ec. (3.7)', sp.simplify(ext - eq11) == 0, ''))
    return report('Paso 2: perturbaciones radiales', res)


if __name__ == '__main__':
    run()
