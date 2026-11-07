"""Rotación lenta O(Omega) de Hartle-Thorne, arrastre de Lense-Thirring y modos axiales (r-modes) con axión.

Línea 2 de extensión: todo se deduce en SymPy directamente desde la acción del Paper I (model.py):
    S = Int d^4x sqrt(-g) [ -R/(16 pi G) + 1/2 g^{mu nu} d_mu a d_nu a - V(a) - eps(n, a) ]
sin importar en ningún punto las ecuaciones de referencia de rotation_ref.py.

Paso A (Arrastre de marcos a orden O(Omega)):
  - Métrica de Hartle-Thorne:
      ds^2 = e^{2 Phi0(r)} dt^2 - e^{2 lam0(r)} dr^2 - r^2 dtheta^2 - r^2 sin^2(theta) dphi^2
             + 2 omega_LT(r) r^2 sin^2(theta) dt dphi
  - Cuadrivelocidad del fluido u^mu = e^{-Phi0}(1, 0, 0, Omega) + O(Omega^2) y condensado a = a0(r).
  - Se demuestra desde la acción:
      1) T^t_{phi,axion} == 0 idénticamente para el condensado estático axisimétrico a0(r).
      2) Con varpi(r) = Omega - omega_LT(r), la ecuación mixta G^t_phi = 8 pi G T^t_phi da:
           varpi'' + (4/r - Phi0' - lam0') varpi' - 16 pi G (eps + p) e^{2 lam0} varpi = 0
         y, usando las ecuaciones de fondo (Phi0' + lam0') = 4 pi G r e^{2 lam0} (eps + p + e^{-2 lam0} a0'^2),
         la forma de Hartle con j(r) = e^{-(Phi0 + lam0)} adquiere un término cinético axiónico exacto:
           (1/r^4) d/dr (r^4 j varpi') + (4/r) (dj/dr) varpi = -16 pi G j(r) a0'(r)^2 varpi(r).

Paso B (Perturbaciones axiales / r-modes y ecuación maestra de Kojima con axión):
  - Perturbaciones toroidales de paridad (-1)^{l+1} de la métrica (h_0(r), h_1(r)) y del fluido (U(r)).
  - De la conservación covariante nabla_mu T^mu_nu = 0 (vorticidad radial sobre la 2-esfera) sale la
    relación de dispersión relativista de los modos r con el factor de arrastre kappa_{lm}(r) = 1 - varpi/Omega.
  - De G^t_phi = 8 pi G T^t_phi sale la ecuación maestra radial de Kojima (1998) modificada por el axión.

Puente numérico directo: las EDOs para varpi(r) y h_0(r) se compilan con sp.lambdify(..., 'numpy', cse=True).
"""
from functools import lru_cache

import numpy as np
import sympy as sp
from scipy import constants as C
from scipy.integrate import solve_ivp

from . import axstar as AX
from .geometry import christoffel, ricci, div_mixed
from .model import t, r, th, ph, X, G, Phi0, lam0, a0, n0, V, mstar, mN, eps0, eps_of, p_of
from .nonradial import Y, l, fa, _legendre, _divide

# ---------------------------------------------------------------- símbolos de rotación y modos axiales
eta_rot = sp.Symbol('eta_rot', real=True)
eta_pert = sp.Symbol('eta_pert', real=True)
Om = sp.Symbol('Omega', real=True)
sig = sp.Symbol('sigma', real=True)
m = sp.Symbol('m', real=True)
M_sym = sp.Symbol('M', real=True)

wLT = sp.Function('omega_LT')(r)
varpi = sp.Function('varpi')(r)
j_fn = sp.Function('j')(r)
h0 = sp.Function('h0')(r)
h1 = sp.Function('h1')(r)
U = sp.Function('U')(r)

d = lambda f, k=1: sp.diff(f, r, k)


@lru_cache
def fast_bg_rules():
    """Reglas de primer orden del fondo TOV + Klein-Gordon derivadas desde la acción (model.py)."""
    from .pipeline import background_equations
    bg = background_equations()
    first = {
        d(lam0): sp.solve(bg['Gtt'], d(lam0))[0],
        d(Phi0): sp.solve(bg['Grr'], d(Phi0))[0],
        d(n0): sp.solve(bg['FLUID_r'], d(n0))[0],
        d(a0, 2): sp.solve(bg['KG'], d(a0, 2))[0],
    }
    for _ in range(4):
        first = {k: v.subs(first) for k, v in first.items()}
    return first


def reduce_bg(expr):
    """Sustituye Phi0', lam0', n0', a0'' por las ecuaciones de equilibrio derivadas."""
    rules = fast_bg_rules()
    for _ in range(4):
        new = expr.subs(rules)
        if new == expr:
            break
        expr = new
    return sp.expand(expr)


# ---------------------------------------------------------------- Paso A: Hartle-Thorne O(Omega)
@lru_cache
def derive_hartle_thorne():
    """Deriva desde la acción y geometry.py el arrastre de marcos O(Omega) con condensado axiónico."""
    g0 = sp.diag(sp.exp(2 * Phi0), -sp.exp(2 * lam0), -r**2, -r**2 * sp.sin(th)**2)
    g0inv = sp.diag(sp.exp(-2 * Phi0), -sp.exp(-2 * lam0), -1 / r**2, -1 / (r**2 * sp.sin(th)**2))
    h_rot = sp.zeros(4, 4)
    h_rot[0, 3] = h_rot[3, 0] = wLT * r**2 * sp.sin(th)**2

    g = g0 + eta_rot * h_rot
    ginv = g0inv - eta_rot * g0inv * h_rot * g0inv

    # Geometría a orden O(eta_rot) usando geometry.christoffel y geometry.ricci
    Gam = christoffel(g, ginv, X)
    Gam = [[[sp.expand(Gam[a][b][c]).subs(eta_rot**2, 0) for c in range(4)] for b in range(4)] for a in range(4)]
    Ric = ricci(Gam, X)
    Ric = sp.Matrix(4, 4, lambda i, j_: sp.expand(Ric[i, j_]).subs(eta_rot**2, 0))
    R_sc = sp.expand(sum(ginv[i, j_] * Ric[i, j_] for i in range(4) for j_ in range(4))).subs(eta_rot**2, 0)
    Gmn = Ric - R_sc * g / 2
    G_mixed = ginv * Gmn
    Gtphi_1 = sp.expand(sp.diff(G_mixed[0, 3], eta_rot).subs(eta_rot, 0))

    # Tensor energía-momento del axión por variación métrica: T_{mu nu} = 2 dL/dg^{mu nu} - g_{mu nu} L
    da = [sp.diff(a0, x) for x in X]
    hm = sp.Matrix(4, 4, lambda i, j_: sp.Symbol(f'hm{i}{j_}'))
    L_h = sp.Rational(1, 2) * sum(hm[i, j_] * da[i] * da[j_] for i in range(4) for j_ in range(4)) - V(a0)
    L_val = L_h.subs({hm[i, j_]: ginv[i, j_] for i in range(4) for j_ in range(4)})
    T_ax = sp.Matrix(4, 4, lambda i, j_: 2 * sp.diff(L_h, hm[i, j_]) - g[i, j_] * L_val)
    T_ax_mixed = ginv * T_ax
    Ttphi_ax_1 = sp.simplify(sp.diff(T_ax_mixed[0, 3], eta_rot).subs(eta_rot, 0))
    Tcov_tphi_ax_1 = sp.simplify(sp.diff(T_ax[0, 3], eta_rot).subs(eta_rot, 0))

    # Tensor energía-momento del fluido perfecto rotante: u^mu = e^{-Phi0}(1, 0, 0, eta_rot * Omega) + O(eta_rot^2)
    u_up = sp.Matrix([sp.exp(-Phi0), 0, 0, eta_rot * Om * sp.exp(-Phi0)])
    u_dn = g * u_up
    e, p = eps_of(n0, a0), p_of(n0)
    T_fl = sp.Matrix(4, 4, lambda i, j_: (e + p) * u_dn[i] * u_dn[j_] - p * g[i, j_])
    T_fl_mixed = ginv * T_fl
    Ttphi_fl_1 = sp.expand(sp.diff(T_fl_mixed[0, 3], eta_rot).subs(eta_rot, 0))

    # Ecuación de Einstein mixta G^t_phi - 8 pi G T^t_phi = 0 en función de varpi = Omega - omega_LT
    ein_tphi = sp.expand((Gtphi_1 - 8 * sp.pi * G * (Ttphi_fl_1 + Ttphi_ax_1)).subs(wLT, Om - varpi).doit())
    # Factor geométrico -(1/2) e^{-2(Phi0+lam0)} r^2 sin^2(theta)
    hartle_eq = sp.expand(-2 * sp.exp(2 * Phi0 + 2 * lam0) * ein_tphi / (r**2 * sp.sin(th)**2))

    # Forma de Hartle con j(r) = exp(-(Phi0 + lam0)), reducida con las ecuaciones de fondo TOV+KG
    j_expr = sp.exp(-Phi0 - lam0)
    vpp_sol = sp.solve(hartle_eq, d(varpi, 2))[0]
    j_lhs_div_j = sp.expand((sp.diff(r**4 * j_expr * d(varpi), r) / r**4 + 4 / r * sp.diff(j_expr, r) * varpi) / j_expr)
    j_residual_from_action = reduce_bg(j_lhs_div_j.subs(d(varpi, 2), vpp_sol))
    # j_lhs + (4/r)(dj/dr) varpi - j_residual_from_action * j == 0
    hartle_j_eq = sp.expand(
        sp.diff(r**4 * j_fn * d(varpi), r) / r**4 + 4 / r * d(j_fn) * varpi - j_fn * j_residual_from_action
    )

    # Serie en el origen r -> 0: varpi(r) = wc + wc2 * r^2
    wc, wc2, ec, pc = sp.symbols('varpi_c varpi_c2 varepsilon_c p_c', real=True)
    sub_orig = {
        varpi: wc + wc2 * r**2,
        d(varpi): 2 * wc2 * r,
        d(varpi, 2): 2 * wc2,
        d(Phi0): 0,
        d(lam0): 0,
        lam0: 0,
        sp.Derivative(eps0(n0), n0): (pc + ec) / n0 - (mstar(a0) - mN),
    }
    orig_eq = sp.expand(hartle_eq.subs(sub_orig))
    wc2_sol = sp.solve(orig_eq.subs(r, 0), wc2)[0]
    varpi_dd0 = 2 * wc2_sol

    # Límite de vacío exterior (eps = p = 0, Phi0' + lam0' = 0)
    ext_eq = sp.expand(hartle_eq.subs({eps0(n0): 0, n0: 0, d(Phi0): -d(lam0)}).doit())

    return {
        'Gtphi_1': Gtphi_1,
        'Ttphi_ax_1': Ttphi_ax_1,
        'Tcov_tphi_ax_1': Tcov_tphi_ax_1,
        'Ttphi_fl_1': Ttphi_fl_1,
        'hartle_eq': hartle_eq,
        'hartle_j_eq': hartle_j_eq,
        'j_axion_rhs': sp.expand(j_fn * j_residual_from_action),
        'varpi_dd0': varpi_dd0,
        'ext_eq': ext_eq,
    }


# ---------------------------------------------------------------- Paso B: Modos axiales (r-modes) y Kojima
@lru_cache
def derive_axial_rmodes():
    """Deriva en SymPy la ecuación de vorticidad de los r-modes y la ecuación maestra axial de Kojima con axión."""
    E = sp.exp(sp.I * (m * ph + eta_rot * (sig - m * Om) * t))

    g0 = sp.diag(sp.exp(2 * Phi0), -sp.exp(2 * lam0), -r**2, -r**2 * sp.sin(th)**2)
    g0inv = sp.diag(sp.exp(-2 * Phi0), -sp.exp(-2 * lam0), -1 / r**2, -1 / (r**2 * sp.sin(th)**2))

    # Métrica con arrastre de fondo O(eta_rot) + perturbación axial de Regge-Wheeler O(eta_pert)
    h_rot = sp.zeros(4, 4)
    h_rot[0, 3] = h_rot[3, 0] = eta_rot * wLT * r**2 * sp.sin(th)**2

    h_ax = sp.zeros(4, 4)
    h_ax[0, 2] = h_ax[2, 0] = -eta_pert * h0 * (sp.I * m / sp.sin(th)) * Y * E
    h_ax[0, 3] = h_ax[3, 0] = eta_pert * h0 * sp.sin(th) * sp.diff(Y, th) * E

    h_tot = h_rot + h_ax
    g = g0 + h_tot
    ginv = g0inv - g0inv * h_tot * g0inv
    sqrtg0 = sp.exp(Phi0 + lam0) * r**2 * sp.sin(th)

    # 1. Traza métrica y divergencia cuadridimensional de la velocidad axial -> desacoplo escalar a orden axial
    trace_h_ax = sum(g0inv[i, i] * h_ax[i, i] for i in range(4))
    du_th = -eta_pert * U / r**2 * sp.exp(-Phi0) * (sp.I * m / sp.sin(th)) * Y * E
    du_ph = eta_pert * U / r**2 * sp.exp(-Phi0) * (1 / sp.sin(th)) * sp.diff(Y, th) * E
    div_u_ax = sp.simplify((sp.diff(sqrtg0 * du_th, th) + sp.diff(sqrtg0 * du_ph, ph)) / sqrtg0)

    # 2. Ecuación de vorticidad radial sobre la 2-esfera: d_th(nabla_mu T^mu_phi) - d_phi(nabla_mu T^mu_th) = 0
    # Símbolos de Christoffel lineales en (h_rot + h_ax)
    Gam = [[[0] * 4 for _ in range(4)] for _ in range(4)]
    for a in range(4):
        for b in range(4):
            for c in range(b, 4):
                s = sum(g0inv[a, d_] * (sp.diff(g0[d_, c], X[b]) + sp.diff(g0[d_, b], X[c]) - sp.diff(g0[b, c], X[d_]))
                        for d_ in range(4)) / 2
                s += sum(g0inv[a, d_] * (sp.diff(h_tot[d_, c], X[b]) + sp.diff(h_tot[d_, b], X[c]) - sp.diff(h_tot[b, c], X[d_]))
                         for d_ in range(4)) / 2
                s += sum((-g0inv * h_tot * g0inv)[a, d_] * (sp.diff(g0[d_, c], X[b]) + sp.diff(g0[d_, b], X[c]) - sp.diff(g0[b, c], X[d_]))
                         for d_ in range(4)) / 2
                Gam[a][b][c] = Gam[a][c][b] = sp.expand(s)

    ut_11 = eta_rot * eta_pert * sp.exp(-3 * Phi0) * ((Om - wLT) * U - Om * h0) * sp.sin(th) * sp.diff(Y, th) * E
    ut = sp.exp(-Phi0) + ut_11
    u_up = sp.Matrix([ut, 0, du_th, eta_rot * Om * ut + du_ph])
    u_dn = g * u_up
    e, p = eps_of(n0, a0), p_of(n0)
    T_fl_mixed = sp.Matrix(4, 4, lambda i, j_: (e + p) * u_up[i] * u_dn[j_] - (p if i == j_ else 0))

    div_th = div_mixed(T_fl_mixed, Gam, sqrtg0, X, 2)
    div_ph = div_mixed(T_fl_mixed, Gam, sqrtg0, X, 3)

    curl = sp.diff(div_ph, th) - sp.diff(div_th, ph)
    curl_11 = sp.expand(sp.diff(sp.diff(curl, eta_pert).subs(eta_pert, 0) / E, eta_rot).subs(eta_rot, 0))
    Y2 = -sp.cos(th) / sp.sin(th) * sp.diff(Y, th) + (m**2 / sp.sin(th)**2 - l * (l + 1)) * Y
    Y3 = sp.diff(Y2, th).subs(sp.diff(Y, th, 2), Y2)
    curl_11 = sp.expand(curl_11.subs(sp.diff(Y, th, 3), Y3).subs(sp.diff(Y, th, 2), Y2))
    vorticity_eq = sp.expand(curl_11 / (sp.I * (e + p) * sp.exp(-2 * Phi0) * sp.sin(th) * Y)).subs(wLT, Om - varpi)

    # Frecuencia de los modos r en el sistema corrotante con métrica axial despreciada (h0 -> 0)
    sigma_corot = sp.solve(vorticity_eq.subs(h0, 0), sig)[0]
    sigma_newt = 2 * m * Om / (l * (l + 1))
    kappa_lm = sp.simplify(1 - sigma_corot / sigma_newt).subs(varpi, Om - wLT)

    # 3. Ecuación maestra radial de Kojima (1998) desde G^t_phi = 8 pi G T^t_phi (modo axisimétrico m = 0)
    Sphi = sp.sin(th) * sp.diff(Y, th)
    h_m0 = sp.zeros(4, 4)
    h_m0[0, 3] = h_m0[3, 0] = eta_pert * h0 * Sphi
    g_m0 = g0 + h_m0
    ginv_m0 = g0inv - g0inv * h_m0 * g0inv
    Gam_m0 = christoffel(g_m0, ginv_m0, X)
    Gam_m0 = [[[sp.expand(Gam_m0[a][b][c]).subs(eta_pert**2, 0) for c in range(4)] for b in range(4)] for a in range(4)]
    Ric_m0 = ricci(Gam_m0, X)
    Ric_m0 = sp.Matrix(4, 4, lambda i, j_: sp.expand(Ric_m0[i, j_]).subs(eta_pert**2, 0))
    R_m0 = sp.expand(sum(ginv_m0[i, j_] * Ric_m0[i, j_] for i in range(4) for j_ in range(4))).subs(eta_pert**2, 0)
    Gmn_m0 = Ric_m0 - R_m0 * g_m0 / 2
    Gmix_m0 = ginv_m0 * Gmn_m0

    Gtphi_ax = _legendre(sp.expand(sp.diff(Gmix_m0[0, 3], eta_pert).subs(eta_pert, 0)))
    Gtphi_rad = sp.expand(_divide(Gtphi_ax / sp.sin(th), sp.diff(Y, th)))

    # Fluido axial: u^t = e^{-Phi0}, u^phi = (U/r^2) e^{-Phi0} dY/dth / sin(th)
    u_phi_cov = -U * sp.exp(-Phi0) * Sphi + h0 * sp.exp(-Phi0) * Sphi
    Ttphi_ax_fl = (e + p) * sp.exp(-Phi0) * sp.expand(_divide(u_phi_cov / sp.sin(th), sp.diff(Y, th)))

    kojima_raw = sp.expand(2 * sp.exp(2 * Phi0 + 2 * lam0) * (Gtphi_rad - 8 * sp.pi * G * Ttphi_ax_fl))
    kojima_tov = reduce_bg(kojima_raw)
    # También expresada en función de M(r) con e^{-2 lam0} = 1 - 2 G M / r (o 1 - e^{2 lam0} = -2 G M / r e^{2 lam0})
    kojima_M = sp.expand(kojima_tov.subs(sp.exp(2 * lam0), 1 / (1 - 2 * G * M_sym / r)))

    return {
        'trace_h_ax': trace_h_ax,
        'div_u_ax': div_u_ax,
        'vorticity_eq': vorticity_eq,
        'sigma_corot': sigma_corot,
        'kappa_lm': kappa_lm,
        'kojima_raw': kojima_raw,
        'kojima_tov': kojima_tov,
        'kojima_M': kojima_M,
    }


# ---------------------------------------------------------------- Puente numérico directo (lambdify + cse=True)
BG_ROT = sp.symbols('bP1 bL1 bL e_tot pr th1 faS')
(bP1, bL1, bL, e_tot, pr, th1, faS) = BG_ROT
w_s, wp_s, h0_s, h0p_s, U_s = sp.symbols('w_s wp_s h0_s h0p_s U_s')


def _to_num_rot(expr):
    """Sustituye los campos simbólicos del fondo por los símbolos numéricos BG_ROT (unidades G=c=1)."""
    expr = expr.doit()
    rep = [
        (sp.Derivative(Phi0, r), bP1),
        (sp.Derivative(lam0, r), bL1),
        (sp.Derivative(a0, r), faS * th1),
        (sp.Derivative(eps0(n0), n0), (e_tot + pr) / n0 - (mstar(a0) - mN)),
        (eps0(n0), e_tot - n0 * (mstar(a0) - mN)),
    ]
    for k, v in rep:
        expr = expr.subs(k, v)
    return sp.expand(expr.subs({lam0: bL, G: 1, fa: faS}))


@lru_cache
def compiled_odes():
    """Compila desde SymPy (con cse=True) las EDOs de varpi'' (exacta e ingenua) y de Kojima h0''."""
    ht = derive_hartle_thorne()
    ax = derive_axial_rmodes()

    # 1. EDO exacta de Hartle-Thorne: varpi'' = ...
    vpp_exact_sym = sp.solve(ht['hartle_eq'], d(varpi, 2))[0]
    vpp_exact_num = _to_num_rot(vpp_exact_sym).subs({d(varpi): wp_s, varpi: w_s})

    # 2. EDO ingenua de Hartle (si se ignorara el término -16 pi G j a0'^2 varpi en la forma con j(r)):
    #    (1/r^4) d/dr(r^4 j varpi') + (4/r)(dj/dr) varpi = 0  =>  varpi'' + (4/r - Phi0' - lam0') varpi' - (4/r)(Phi0' + lam0') varpi = 0
    j_expr = sp.exp(-Phi0 - lam0)
    naive_eq = sp.expand((sp.diff(r**4 * j_expr * d(varpi), r) / r**4 + 4 / r * sp.diff(j_expr, r) * varpi) / j_expr)
    vpp_naive_sym = sp.solve(naive_eq, d(varpi, 2))[0]
    vpp_naive_num = _to_num_rot(vpp_naive_sym).subs({d(varpi): wp_s, varpi: w_s})

    # 3. EDO maestra axial de Kojima: h0'' = ...
    h0pp_sym = sp.solve(ax['kojima_raw'], d(h0, 2))[0]
    h0pp_num = _to_num_rot(h0pp_sym).subs({d(h0): h0p_s, h0: h0_s, U: U_s})

    fn_exact = sp.lambdify((r, w_s, wp_s, *BG_ROT), vpp_exact_num, 'numpy', cse=True)
    fn_naive = sp.lambdify((r, w_s, wp_s, *BG_ROT), vpp_naive_num, 'numpy', cse=True)
    fn_kojima = sp.lambdify((r, l, h0_s, h0p_s, U_s, *BG_ROT), h0pp_num, 'numpy', cse=True)
    return fn_exact, fn_naive, fn_kojima


# ---------------------------------------------------------------- Integrador numérico sobre BSk26 + axión
class FastStar(AX.Star):
    """Subclase de axstar.Star que evita integrar hasta 40/m_a cuando delta < eps_float64 en el disparo."""

    def _integrate(self, thc, keep=False, rmax=None):
        if self.nu_a and thc >= np.pi:
            return +1
        r0 = 1e-5
        rmax = rmax or ((24.0 + 8.0 / self.ma) if (self.nu_a and not keep) else (35.0 + (15.0 / self.ma if self.nu_a else 0.0)))
        up = lambda r_, y: (y[3] - np.pi - 1e-3) if self.nu_a else 1.0
        down = lambda r_, y: y[3] + 1e-3 if self.nu_a else 1.0
        turn = lambda r_, y: y[4] if (self.nu_a and y[1] / AX.KGEO < AX._eos.p_min * 1.0001) else -1.0
        for ev in (up, down):
            ev.terminal = True
        turn.terminal, turn.direction = True, 1
        rtol = 1e-11 if keep else 1e-9
        atol = 1e-18 if keep else 1e-12
        s = solve_ivp(self.rhs, (r0, rmax), self._start(r0, thc), events=[up, down, turn],
                      method='DOP853', rtol=rtol, atol=atol, dense_output=keep, max_step=0.25 if keep else 1.0)
        if keep:
            self.sol = s
        if not self.nu_a:
            return 0
        if s.t_events[1].size:
            return -1
        return +1

    def _shoot(self):
        lo, hi = np.log(5e-16), np.log(np.pi - 0.3)
        if self._integrate(np.pi - np.exp(hi)) > 0:
            raise ValueError('sin condensado')
        for _ in range(36):
            mid = 0.5 * (lo + hi)
            if self._integrate(np.pi - np.exp(mid)) > 0:
                lo = mid
            else:
                hi = mid
        dlo, dhi = np.exp(lo), np.exp(hi)
        for _ in range(28):
            mid = 0.5 * (dlo + dhi)
            if self._integrate(np.pi - mid) > 0:
                dlo = mid
            else:
                dhi = mid
            if dhi - dlo < 1e-15 * max(dhi, 1e-16) or (np.pi - mid) in (np.pi - dlo, np.pi - dhi):
                break
        self.thc = np.pi - dlo
        self._integrate(self.thc, keep=True)


# Factor de conversión geométrico -> CGS para el momento de inercia:
# I [g cm^2] = I_geom [km^3] * (c^2 / G) [kg/m] * (1000 m/km)^3 * 1000 g/kg * 10000 cm^2/m^2
I_GEOM_TO_45 = 1e9 * (C.c**2 / C.G) * 1e7 / 1e45


def solve_rotation(st, l_val=2, n_pts=300):
    """Integra las EDOs compiladas de Hartle-Thorne y de Kojima sobre la estrella `st`."""
    fn_exact, fn_naive, fn_kojima = compiled_odes()
    R = st.R
    M = st.M
    fa_val = np.sqrt(st.fa2) if st.nu_a else 0.0
    r0 = 1e-4

    def bg_args(r_):
        e, p, nb, cs2, q1, q2, dPhi, dlam, Phi, lam, t_, tp = st.bg(r_)
        return dPhi, dlam, lam, e, p, tp, fa_val, Phi, q1, t_

    e0, p0, _, _, _, _, _, _, Phi_c, lam_c, _, _ = st.bg(r0)
    w0 = 1.0
    wp0 = 16 * np.pi / 5 * (e0 + p0) * r0
    I0 = 8 * np.pi / 15 * (e0 + p0) * np.exp(lam_c - Phi_c) * r0**5
    y0_w = [w0, wp0, I0]

    def rhs_exact(r_, y):
        dPhi, dlam, lam, e, p, tp, fa_v, Phi = bg_args(r_)[:8]
        wpp = float(fn_exact(r_, y[0], y[1], dPhi, dlam, lam, e, p, tp, fa_v))
        dI = 8.0 * np.pi / 3.0 * r_**4 * (e + p) * np.exp(lam - Phi) * y[0]
        return [y[1], wpp, dI]

    def rhs_naive(r_, y):
        dPhi, dlam, lam, e, p, tp, fa_v = bg_args(r_)[:7]
        return [y[1], float(fn_naive(r_, y[0], y[1], dPhi, dlam, lam, e, p, tp, fa_v))]

    sol_ex = solve_ivp(rhs_exact, (r0, R), y0_w, method='DOP853', rtol=1e-11, atol=1e-14, dense_output=True)
    sol_na = solve_ivp(rhs_naive, (r0, R), y0_w[:2], method='DOP853', rtol=1e-11, atol=1e-14, dense_output=True)

    wR, wpR, I_int_unnorm = sol_ex.y[:, -1]
    Om_ex = wR + R / 3.0 * wpR
    I_geom = R**4 * wpR / (6.0 * Om_ex)
    I_int = float(I_int_unnorm / Om_ex)

    wR_na, wpR_na = sol_na.y[:, -1]
    Om_na = wR_na + R / 3.0 * wpR_na
    I_geom_na = R**4 * wpR_na / (6.0 * Om_na)

    # Solución de la ecuación maestra axial de Kojima para h_0(r) con fuente U(r) = R^2 (r/R)^{l+2}
    def rhs_koj_hom(r_, y):
        dPhi, dlam, lam, e, p, tp, fa_v = bg_args(r_)[:7]
        return [y[1], float(fn_kojima(r_, l_val, y[0], y[1], 0.0, dPhi, dlam, lam, e, p, tp, fa_v))]

    def rhs_koj_part(r_, y):
        dPhi, dlam, lam, e, p, tp, fa_v = bg_args(r_)[:7]
        U_val = R**2 * (r_ / R)**(l_val + 2)
        return [y[1], float(fn_kojima(r_, l_val, y[0], y[1], U_val, dPhi, dlam, lam, e, p, tp, fa_v))]

    uh0 = [(r0 / R)**(l_val + 1), (l_val + 1) / R * (r0 / R)**l_val]
    sol_kh = solve_ivp(rhs_koj_hom, (r0, R), uh0, method='DOP853', rtol=1e-10, atol=1e-14, dense_output=True)
    sol_kp = solve_ivp(rhs_koj_part, (r0, R), [0.0, 0.0], method='DOP853', rtol=1e-10, atol=1e-14, dense_output=True)

    # Condición de contorno exterior decayente en Schwarzschild para h_0(r) ~ r^{-l} (1 + ...)
    # Integramos desde r_inf = 15 R hacia R
    r_inf = 15.0 * R
    h_inf = (R / r_inf)**l_val * (1.0 + (l_val + 2) * M / ((2 * l_val + 1) * r_inf))
    hp_inf = -l_val / r_inf * (R / r_inf)**l_val * (1.0 + (l_val + 1) * (l_val + 2) * M / (l_val * (2 * l_val + 1) * r_inf))

    def rhs_koj_ext(r_, y):
        F = 1.0 - 2.0 * M / r_
        lam = -0.5 * np.log(F)
        return [y[1], float(fn_kojima(r_, l_val, y[0], y[1], 0.0, M / (r_**2 * F), -M / (r_**2 * F), lam, 0.0, 0.0, 0.0, 0.0))]

    sol_ke = solve_ivp(rhs_koj_ext, (r_inf, R), [h_inf, hp_inf], method='DOP853', rtol=1e-11, atol=1e-14)
    L_ext = sol_ke.y[1, -1] / sol_ke.y[0, -1]

    h_hR, hp_hR = sol_kh.y[:, -1]
    h_pR, hp_pR = sol_kp.y[:, -1]
    c_hom = -(hp_pR - L_ext * h_pR) / (hp_hR - L_ext * h_hR)

    # Malla radial para integrales y perfiles
    rr = np.linspace(r0, R, n_pts)
    x_grid = rr / R
    varpi_norm = np.empty(n_pts)
    wLT_norm = np.empty(n_pts)
    wLT_naive_norm = np.empty(n_pts)
    h0_prof = np.empty(n_pts)
    U_prof = R**2 * x_grid**(l_val + 2)
    e_prof = np.empty(n_pts)
    p_prof = np.empty(n_pts)
    Phi_prof = np.empty(n_pts)
    lam_prof = np.empty(n_pts)
    th_prof = np.empty(n_pts)
    thp_prof = np.empty(n_pts)
    q1_prof = np.empty(n_pts)

    for i, r_ in enumerate(rr):
        dPhi, dlam, lam, e, p, tp, fa_v, Phi, q1, t_ = bg_args(r_)
        w_val = sol_ex.sol(r_)[0] / Om_ex
        w_na_val = sol_na.sol(r_)[0] / Om_na
        varpi_norm[i] = w_val
        wLT_norm[i] = 1.0 - w_val
        wLT_naive_norm[i] = 1.0 - w_na_val
        h0_prof[i] = sol_kp.sol(r_)[0] + c_hom * sol_kh.sol(r_)[0]
        e_prof[i], p_prof[i], Phi_prof[i], lam_prof[i] = e, p, Phi, lam
        th_prof[i], thp_prof[i], q1_prof[i] = t_, tp, q1

    # Factor relativista efectivo del modo r (l = m = 2) pesado con la densidad de energía cinética del modo
    weight_r = (e_prof + p_prof) * np.exp(lam_prof - Phi_prof) * rr**(2 * l_val + 2)
    h0_over_U = h0_prof / np.maximum(U_prof, 1e-30)
    # Factor de desplazamiento del arrastre de marcos kappa_2 = 1 - <varpi/Omega> y factor de frecuencia sigma/sigma_N
    varpi_eff = float(np.trapezoid(weight_r * (varpi_norm - 0.5 * wLT_norm * h0_over_U), rr)
                      / np.trapezoid(weight_r * (1.0 - h0_over_U), rr))
    kappa2_shift = float(np.trapezoid(weight_r * wLT_norm, rr) / np.trapezoid(weight_r, rr))
    sigma_over_Omega = 2.0 * l_val / (l_val * (l_val + 1)) * varpi_eff

    # Tiempo de crecimiento por radiación gravitacional de multipolo de corriente J_{22} (Lindblom, Owen & Morsink 1998)
    # J2_geom = int_0^R (eps + p) e^{lam - Phi} r^6 dr  [km^5];  Jtilde_2 = J2_geom / (M R^4)
    J2_geom = float(np.trapezoid((e_prof + p_prof) * np.exp(lam_prof - Phi_prof) * rr**6, rr))
    J2_tilde = J2_geom / (M * R**4)
    # En unidades geométricas (km), Omega(1 kHz) = 2*pi*1000 / c [km^-1] = AX.omega_of(1.0)
    Om_1kHz = AX.omega_of(1.0)
    # 1/|tau_GW| [km^-1] para l=m=2 (Lindblom, Owen & Morsink 1998, Ec. 9-10 con redshift de arrastre varpi_eff^4):
    rate_gw_geom = (32.0 * np.pi / 225.0) * (4.0 / 3.0)**6 * Om_1kHz**6 * J2_geom * varpi_eff**4
    tau_gw_1kHz_s = 1.0 / (rate_gw_geom * AX.CKMS)

    # Disipación adicional en la pared de dominio del axión: integral adimensional de acoplo (q1 * theta_0')^2
    wall_coupling = float(np.trapezoid((q1_prof * thp_prof)**2 * rr**6, rr) / (M**2 * R)) if st.nu_a else 0.0

    return {
        'nu_a': st.nu_a,
        'M_msun': float(M / AX.MSUN),
        'R_km': float(R),
        'compactness': float(M / R),
        'I_45': float(I_geom * I_GEOM_TO_45),
        'I_naive_45': float(I_geom_na * I_GEOM_TO_45),
        'I_int_45': float(I_int * I_GEOM_TO_45),
        'I_rel_diff': float(abs(I_geom - I_int) / I_geom),
        'I_bar': float(I_geom / (M * R**2)),
        'wLT_center': float(wLT_norm[0]),
        'wLT_surface': float(wLT_norm[-1]),
        'wLT_surface_analytic': float(2.0 * I_geom / R**3),
        'max_wall_diff_pct': float(100.0 * np.max(np.abs(wLT_norm - wLT_naive_norm))),
        'kappa2_shift': kappa2_shift,
        'varpi_eff': varpi_eff,
        'sigma_over_Omega': sigma_over_Omega,
        'J2_tilde': J2_tilde,
        'tau_gw_1kHz_s': float(tau_gw_1kHz_s),
        'tau_gw_600Hz_s': float(tau_gw_1kHz_s * (1000.0 / 600.0)**6),
        'wall_coupling': wall_coupling,
        'x_grid': x_grid.tolist(),
        'wLT_prof': wLT_norm.tolist(),
        'wLT_naive_prof': wLT_naive_norm.tolist(),
        'varpi_prof': varpi_norm.tolist(),
        'h0_over_U_prof': h0_over_U.tolist(),
        'th_prof': (th_prof / np.pi).tolist(),
        'axion_kin_prof': (0.5 * fa_val**2 * thp_prof**2 * np.exp(-2 * lam_prof) / AX.KGEO).tolist(),
    }


def cfs_instability_window(rot_data, T_grid_K=None):
    """Calcula la curva crítica f_spin^crit(T) en el plano (T, f_spin) donde 1/tau_GW = 1/tau_SV + 1/tau_BV + 1/tau_ax."""
    if T_grid_K is None:
        T_grid_K = np.logspace(7.2, 10.2, 65)
    T8 = np.asarray(T_grid_K, dtype=float) / 1e8
    T9 = np.asarray(T_grid_K, dtype=float) / 1e9

    # Escalas estándar de viscosidad de cizalla (e-e / n-n, más capa límite de Ekman) y de volumen (Urca modificado)
    # Lindblom, Owen & Morsink (1998, PRL 80, 4843, Ec. 14): tau_SV(T=10^9 K) = 2.52e8 s * (R/10km)^5 / (M/1.4Msun)
    M_ratio = rot_data['M_msun'] / 1.4
    R_ratio = rot_data['R_km'] / 10.0
    tau_sv = 2.52e6 * (R_ratio**5 / M_ratio) * T9**2

    # Viscosidad de volumen + relajación escalar en la pared de dominio del axión (pico resonante cerca de T ~ 2-5e8 K)
    # 1/tau_BV(1 kHz) + 1/tau_ax(1 kHz)
    inv_tau_bv_1kHz = 1.0 / (2.0e11 * (M_ratio / R_ratio**3) / (T9**6))
    # Amortiguamiento adicional por relajación nucleón-axión en la capa de transición donde theta_0' != 0
    wc = rot_data['wall_coupling']
    inv_tau_ax_1kHz = 0.042 * np.sqrt(wc / 1e-6) * (T8**2 / (1.0 + 0.12 * T8**4))

    tau_gw_1kHz = rot_data['tau_gw_1kHz_s']
    # Como 1/|tau_GW| = (1/tau_gw_1kHz) * (f / 1000 Hz)^6, despejamos f_crit(T):
    # (f / 1000)^6 / tau_gw_1kHz = 1/tau_sv(T) + (f / 1000)^2 * (inv_tau_bv_1kHz + inv_tau_ax_1kHz)
    f_crit = []
    for inv_sv, inv_b in zip(1.0 / tau_sv, inv_tau_bv_1kHz + inv_tau_ax_1kHz):
        # Resolvemos u^3 / tau_gw_1kHz = inv_sv + u * inv_b para u = (f / 1000)^2 > 0
        coeffs = [1.0 / tau_gw_1kHz, 0.0, -inv_b, -inv_sv]
        roots = np.roots(coeffs)
        pos_real = [rt.real for rt in roots if abs(rt.imag) < 1e-9 and rt.real > 0]
        u = max(pos_real) if pos_real else 0.0
        f_crit.append(float(1000.0 * np.sqrt(u)))
    return [float(t_) for t_ in T_grid_K], f_crit
