"""Comparación simbólica 'módulo ecuaciones': ¿es la ecuación del paper consecuencia de las derivadas?"""
import sympy as sp

from .model import r, Phi0, lam0, a0, n0


def symbolize(expr):
    """Sustituye funciones y derivadas por símbolos para que simplify trabaje con polinomios."""
    atoms = sorted(expr.atoms(sp.Derivative), key=lambda e: -sp.count_ops(e))
    atoms += [f for f in expr.atoms(sp.core.function.AppliedUndef)]
    rep = {at: sp.Symbol(f'S{i}') for i, at in enumerate(atoms)}
    return expr.xreplace(rep)


def is_zero(expr):
    s = symbolize(sp.expand(expr))
    return sp.simplify(sp.together(sp.expand(s))) == 0


class Rules:
    """Reglas de sustitución para derivadas del fondo, obtenidas de las ecuaciones de fondo derivadas."""

    def __init__(self, bg):
        self.first = {
            sp.diff(lam0, r): sp.solve(bg['Gtt'], sp.diff(lam0, r))[0],
            sp.diff(Phi0, r): sp.solve(bg['Grr'], sp.diff(Phi0, r))[0],
            sp.diff(n0, r): sp.solve(bg['FLUID_r'], sp.diff(n0, r))[0],
            sp.diff(a0, r, 2): sp.solve(bg['KG'], sp.diff(a0, r, 2))[0],
        }
        self.first = {k: self._apply(v, self.first) for k, v in self.first.items()}
        # derivadas de orden superior por diferenciación + resustitución
        self.higher = {}
        for f, k0 in [(lam0, 1), (Phi0, 1), (n0, 1), (a0, 2)]:
            for k in (k0 + 1, k0 + 2):
                prev = self.higher.get(sp.diff(f, r, k - 1), self.first.get(sp.diff(f, r, k - 1)))
                self.higher[sp.diff(f, r, k)] = self.reduce(sp.diff(prev, r))

    @staticmethod
    def _apply(expr, rules):
        for _ in range(6):
            new = expr.subs(rules)
            if new == expr:
                break
            expr = new
        return expr

    def reduce(self, expr):
        expr = expr.subs({k: v for k, v in sorted(getattr(self, 'higher', {}).items(),
                                                  key=lambda kv: -sp.count_ops(kv[0]))})
        return self._apply(expr, self.first)


def check_modulo(paper_expr, derived_expr, pivot, rules):
    """Resuelve la ecuación derivada para `pivot`, sustituye en la del paper y reduce con el fondo."""
    sol = sp.solve(derived_expr, pivot)
    if len(sol) != 1:
        raise ValueError(f'pivot {pivot}: {len(sol)} soluciones')
    res = rules.reduce(paper_expr.subs(pivot, sol[0]))
    return is_zero(res), res
