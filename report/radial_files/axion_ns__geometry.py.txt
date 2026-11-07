"""Utilidades de geometría diferencial (componentes explícitas, sin sympy.diffgeom).

Convenios: Gamma^a_{bc} = 1/2 g^{ad}(d_b g_{dc} + d_c g_{db} - d_d g_{bc}),
R^a_{bcd} = d_c Gamma^a_{db} - d_d Gamma^a_{cb} + Gamma^a_{ce}Gamma^e_{db} - Gamma^a_{de}Gamma^e_{cb},
R_{bd} = R^a_{bad}.  Con estos convenios G_{mu nu} = 8 pi G T_{mu nu} vale tanto en
signatura (+,-,-,-) como (-,+,+,+), con T_{mu nu} = (e+p) u_mu u_nu - p g_{mu nu} en (+,-,-,-).
"""
import sympy as sp


def christoffel(g, ginv, X):
    n = len(X)
    dg = [[[sp.diff(g[i, j], X[k]) for k in range(n)] for j in range(n)] for i in range(n)]
    Gam = [[[0] * n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(b, n):
                s = sum(ginv[a, d] * (dg[d][c][b] + dg[d][b][c] - dg[b][c][d]) for d in range(n)) / 2
                s = sp.simplify(s)
                Gam[a][b][c] = Gam[a][c][b] = s
    return Gam


def ricci(Gam, X):
    n = len(X)
    Ric = sp.zeros(n, n)
    for b in range(n):
        for d in range(b, n):
            s = 0
            for a in range(n):
                s += sp.diff(Gam[a][d][b], X[a]) - sp.diff(Gam[a][a][b], X[d])
                for e in range(n):
                    s += Gam[a][a][e] * Gam[e][d][b] - Gam[a][d][e] * Gam[e][a][b]
            Ric[b, d] = Ric[d, b] = sp.simplify(s)
    return Ric


def einstein(g, X):
    """Devuelve (G_{mu nu}, Gamma, ginv)."""
    ginv = sp.simplify(g.inv())
    Gam = christoffel(g, ginv, X)
    Ric = ricci(Gam, X)
    R = sp.simplify(sum(ginv[i, j] * Ric[i, j] for i in range(len(X)) for j in range(len(X))))
    return sp.simplify(Ric - R * g / 2), Gam, ginv


def div_mixed(T_up_down, Gam, sqrtg, X, nu):
    """nabla_mu T^mu_nu = (1/sqrt(-g)) d_mu (sqrt(-g) T^mu_nu) - Gamma^l_{mu nu} T^mu_l."""
    n = len(X)
    s = sum(sp.diff(sqrtg * T_up_down[m, nu], X[m]) for m in range(n)) / sqrtg
    s -= sum(Gam[l][m][nu] * T_up_down[m, l] for m in range(n) for l in range(n))
    return s
