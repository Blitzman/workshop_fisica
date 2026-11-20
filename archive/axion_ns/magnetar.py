"""Línea 3: Unificación Axión + Magnetar en espacio-tiempo curvo.

Deriva simbólicamente desde la acción de Einstein-Maxwell-Axión con acoplo topológico
de Chern-Simons:

    L_{a+EM} = sqrt(-g) [ 1/2 g^{mu nu} d_mu a d_nu a - V(a) - eps(n, a)
                          - 1/4 F_{mu nu} F^{mu nu} - 1/4 g_{agg} a F_{mu nu} Ftilde^{mu nu}
                          - A_mu J^mu ]

sobre la métrica esférica estática g_{mu nu} = diag(e^{2 Phi_0}, -e^{2 lambda_0}, -r^2, -r^2 sin^2 theta):
  1. La identidad de Bianchi electromagnética nabla_mu Ftilde^{mu nu} = 0 y la identidad de
     divergencia total de la densidad de Pontryagin F_{mu nu} Ftilde^{mu nu}.
  2. Las ecuaciones de Maxwell-Carroll-Field-Jackiw (MCFJ) y Klein-Gordon acopladas por
     variación de Euler-Lagrange respecto a A_nu y a.
  3. La ecuación relativista de dipolo magnético de Wasserman & Shapiro (1983) para
     A_phi(r, theta) = -a_phi(r) sin^2 theta y la densidad de carga topológica de Witten
     en la pared de dominio.
  4. La corriente topológica de Hall/quiral inducida por el gradiente de la pared de dominio
     a_0'(r) = f_a theta_0'(r), identificando el potencial químico quiral efectivo
     mu_{5,eff}(r) = (pi / alpha_EM) g_{agg} f_a theta_0'(r) (conexión con Dehman & Pons 2024/2025).
  5. El sistema acoplado de oscilaciones magneto-axiónicas (delta theta <-> calE) y su
     compilación directa a matrices numéricas mediante sp.lambdify(..., 'numpy', cse=True).

Regla de separación estricta (AGENTS.md): este módulo JAMÁS importa magnetar_ref.py.
"""
from functools import lru_cache

import sympy as sp

from .model import t, r, th, ph, X, Phi0, lam0, a0, n0, w, V, mstar, eps_of

# ---------------------------------------------------------------- símbolos y campos
g_agg = sp.Symbol('g_agammagamma', positive=True)
alpha_em = sp.Symbol('alpha_EM', positive=True)
fa = sp.Symbol('f_a', positive=True)
eta = sp.Symbol('eta', real=True)
M_ext = sp.Symbol('M', positive=True)
mu_mag = sp.Symbol('mu_B', real=True)
k_w = sp.Symbol('k', positive=True)
eta_ohm = sp.Symbol('eta_ohm', positive=True)
v_A = sp.Symbol('v_A', positive=True)

# Cuadrivector potencial general A_mu(t, r, theta, phi) y cuadricorriente J^mu(t, r, theta, phi)
At = sp.Function('A_t')(t, r, th, ph)
Ar = sp.Function('A_r')(t, r, th, ph)
Ath = sp.Function('A_theta')(t, r, th, ph)
Aph = sp.Function('A_phi')(t, r, th, ph)
A_vec = [At, Ar, Ath, Aph]

Jt = sp.Function('J_t')(t, r, th, ph)
Jr = sp.Function('J_r')(t, r, th, ph)
Jth = sp.Function('J_theta')(t, r, th, ph)
Jph = sp.Function('J_phi')(t, r, th, ph)
J_vec = [Jt, Jr, Jth, Jph]

a_fld = sp.Function('a')(t, r, th, ph)
n_fld = sp.Function('n')(t, r, th, ph)

# Funciones radiales para el fondo dipolar y las perturbaciones magneto-axiónicas
aphi = sp.Function('a_phi')(r)
jphi = sp.Function('j_phi')(r)
th0 = sp.Function('theta0')(r)
B0 = sp.Function('B0')(r)
dth = sp.Function('deltatheta')(r)
calE = sp.Function('calE')(r)
calB = sp.Function('calB')(r)

d = lambda f, k=1: sp.diff(f, r, k)

# Métrica esférica estática de fondo
g = sp.diag(sp.exp(2 * Phi0), -sp.exp(2 * lam0), -r**2, -r**2 * sp.sin(th)**2)
ginv = sp.diag(sp.exp(-2 * Phi0), -sp.exp(-2 * lam0), -1 / r**2, -1 / (r**2 * sp.sin(th)**2))
sqrtg = sp.exp(Phi0 + lam0) * r**2 * sp.sin(th)


def _poly_reduce_trig(expr):
    """Reducción polinómica exacta en (cos theta, sin theta) módulo cos^2 + sin^2 = 1 (sin trigsimp)."""
    c, s = sp.symbols('c_th s_th')
    trig = {sp.tan(th): s / c, sp.cot(th): c / s, sp.cos(th): c, sp.sin(th): s}
    e = sp.together(sp.expand(expr.subs(trig)))
    num, den = sp.fraction(e)
    one = [c**2 + s**2 - 1]
    num_r = sp.reduced(sp.expand(num), one, c, s)[1] if num.has(c, s) else num
    den_r = sp.reduced(sp.expand(den), one, c, s)[1] if den.has(c, s) else den
    q = sp.expand(num_r / den_r)
    return q.subs({c: sp.cos(th), s: sp.sin(th)})


def field_tensors(A_cov):
    """Calcula F_{mu nu}, F^{mu nu} y el dual de Hodge Ftilde^{mu nu} = eps^{mu nu al be} F_{al be} / (2 sqrt(-g))."""
    F_dn = sp.Matrix(4, 4, lambda mu, nu: sp.diff(A_cov[nu], X[mu]) - sp.diff(A_cov[mu], X[nu]))
    F_up = ginv * F_dn * ginv
    Ftilde_up = sp.Matrix(
        4, 4,
        lambda mu, nu: sum(
            sp.LeviCivita(mu, nu, al, be) * F_dn[al, be]
            for al in range(4) for be in range(4)
        ) / (2 * sqrtg)
    )
    return F_dn, F_up, Ftilde_up


@lru_cache
def derive_covariant_mcfj():
    """Deriva las ecuaciones de Maxwell-Carroll-Field-Jackiw y Klein-Gordon por Euler-Lagrange desde L_{a+EM}."""
    F_dn, F_up, Ftilde_up = field_tensors(A_vec)

    # 1. Identidad de Bianchi electromagnética: nabla_mu Ftilde^{mu nu} = (1/sqrt(-g)) d_mu(sqrt(-g) Ftilde^{mu nu}) = 0
    bianchi = [
        sp.expand(sum(sp.diff(sqrtg * Ftilde_up[mu, nu], X[mu]) for mu in range(4)) / sqrtg)
        for nu in range(4)
    ]

    # 2. Identidad de divergencia total para la densidad de Pontryagin:
    #    1/4 sqrt(-g) F_{mu nu} Ftilde^{mu nu} = 1/2 d_mu (sqrt(-g) A_nu Ftilde^{mu nu})
    FF = sum(F_dn[mu, nu] * F_up[mu, nu] for mu in range(4) for nu in range(4))
    FFtilde = sum(F_dn[mu, nu] * Ftilde_up[mu, nu] for mu in range(4) for nu in range(4))
    pontryagin_div = sp.expand(
        sp.Rational(1, 4) * sqrtg * FFtilde
        - sp.Rational(1, 2) * sum(sp.diff(sqrtg * A_vec[nu] * Ftilde_up[mu, nu], X[mu]) for mu in range(4) for nu in range(4))
    )

    # 3. Densidad lagrangiana completa L_{a+EM}
    da = [sp.diff(a_fld, x) for x in X]
    kin_a = sp.Rational(1, 2) * sum(ginv[mu, mu] * da[mu]**2 for mu in range(4))
    L_dens = sqrtg * (
        kin_a - V(a_fld) - eps_of(n_fld, a_fld)
        - sp.Rational(1, 4) * FF
        - sp.Rational(1, 4) * g_agg * a_fld * FFtilde
        - sum(A_vec[mu] * J_vec[mu] for mu in range(4))
    )

    # 4. Ecuaciones de Euler-Lagrange respecto a A_nu (nu = 0, 1, 2, 3)
    mcfj_el = []
    for nu in range(4):
        el_nu = sp.euler_equations(L_dens, A_vec[nu], X)[0].lhs
        mcfj_el.append(sp.expand(el_nu / sqrtg))

    # 5. Ecuación de Euler-Lagrange respecto al axión a
    el_a = sp.euler_equations(L_dens, a_fld, X)[0].lhs
    kg_el = sp.expand(-el_a / sqrtg)

    return {
        'F_dn': F_dn,
        'F_up': F_up,
        'Ftilde_up': Ftilde_up,
        'FF': FF,
        'FFtilde': FFtilde,
        'bianchi': bianchi,
        'pontryagin_div': pontryagin_div,
        'MCFJ': mcfj_el,
        'KG': kg_el,
    }


@lru_cache
def derive_static_dipole():
    """Reduce las ecuaciones covariantes al dipolo magnético estático A_phi = -a_phi(r) sin^2 theta."""
    cov = derive_covariant_mcfj()
    sub_static = {
        At: 0,
        Ar: 0,
        Ath: 0,
        Aph: -aphi * sp.sin(th)**2,
        Jt: 0,
        Jr: 0,
        Jth: 0,
        Jph: 4 * sp.pi * jphi / r**2,
        a_fld: a0,
        n_fld: n0,
    }

    # Ecuación de Maxwell azimutal (nu = 3) multiplicada por -r^2
    mcfj_phi_static = _poly_reduce_trig(sp.expand(-r**2 * cov['MCFJ'][3].subs(sub_static).doit()))

    # Componentes del campo magnético físico en la base ortonormal tetrádica
    # B^{hat r} = -e^{Phi0} sqrt(|g_{rr}|) Ftilde^{t r},  B^{hat theta} = -e^{Phi0} sqrt(|g_{th th}|) Ftilde^{t theta}
    Ft_tr = cov['Ftilde_up'][0, 1].subs(sub_static).doit()
    Ft_tth = cov['Ftilde_up'][0, 2].subs(sub_static).doit()
    B_hat_r = _poly_reduce_trig(-sp.exp(Phi0 + lam0) * Ft_tr)
    B_hat_th = _poly_reduce_trig(-sp.exp(Phi0) * r * Ft_tth)

    # Invariante de Pontryagin estático F_{mu nu} Ftilde^{mu nu} (nulo para campo puramente poloidal)
    fftilde_static = _poly_reduce_trig(cov['FFtilde'].subs(sub_static).doit())

    # Densidad de carga topológica estática inducida por la pared de dominio (efecto Witten en espacio curvo):
    # J^t_{CS} = -g_{agg} (d_r a_0) Ftilde^{r t},  rho_top = e^{Phi0} J^t_{CS} = -g_{agg} (e^{-lam0} d_r a_0) B^{hat r}_0
    Ft_rt = cov['Ftilde_up'][1, 0].subs(sub_static).doit()
    J_cs_t = _poly_reduce_trig((-g_agg * d(a0) * Ft_rt).subs(a0, fa * th0).doit())
    rho_top = _poly_reduce_trig(sp.exp(Phi0) * J_cs_t)

    # Corrientes espaciales de Chern-Simons en equilibrio estático puramente poloidal
    J_cs_spatial = [
        _poly_reduce_trig((-g_agg * d(a0) * cov['Ftilde_up'][1, i].subs(sub_static).doit()))
        for i in (1, 2, 3)
    ]

    # Límite exterior de vacío de Schwarzschild: e^{2 Phi0} = e^{-2 lam0} = 1 - 2 M / r, j_phi = 0
    ws_ext_ode = sp.expand(
        mcfj_phi_ode := mcfj_phi_static.subs({
            jphi: 0,
            d(Phi0): M_ext / (r * (r - 2 * M_ext)),
            d(lam0): -M_ext / (r * (r - 2 * M_ext)),
            sp.exp(-2 * lam0): 1 - 2 * M_ext / r,
            lam0: -sp.Rational(1, 2) * sp.log(1 - 2 * M_ext / r),
        })
    )

    return {
        'WS_dipole': mcfj_phi_static,
        'WS_exterior': ws_ext_ode,
        'B_hat_r': B_hat_r,
        'B_hat_th': B_hat_th,
        'FFtilde_static': fftilde_static,
        'J_cs_t': J_cs_t,
        'rho_top': rho_top,
        'J_cs_spatial': J_cs_spatial,
    }


@lru_cache
def derive_chiral_domain_wall():
    """Deduce la corriente topológica de Hall/quiral en la pared de dominio y el potencial mu_{5,eff}(r)."""
    cov = derive_covariant_mcfj()
    E_phase = sp.exp(-sp.I * w * t)

    # Perturbación electromagnética con sector poloidal (calE) y toroidal (calB) sobre la pared de dominio a_0(r) = f_a theta_0(r)
    sub_chiral = {
        At: 0,
        Ar: 0,
        Ath: eta * calB * sp.sin(th) * E_phase,
        Aph: -aphi * sp.sin(th)**2 + eta * calE * sp.sin(th)**2 * E_phase,
        Jt: 0,
        Jr: 0,
        Jth: 0,
        Jph: 4 * sp.pi * jphi / r**2,
        a_fld: fa * th0,
        n_fld: n0,
    }

    # Término de corriente topológica de Chern-Simons linealizado en la ecuación azimutal (nu = phi) y meridional (nu = theta)
    cs_term_phi = g_agg * sp.diff(a_fld, r) * cov['Ftilde_up'][1, 3]
    cs_term_th = g_agg * sp.diff(a_fld, r) * cov['Ftilde_up'][1, 2]

    J_top_phi = _poly_reduce_trig(
        sp.diff(cs_term_phi.subs(sub_chiral).doit(), eta).subs(eta, 0) / E_phase
    )
    J_top_th = _poly_reduce_trig(
        sp.diff(cs_term_th.subs(sub_chiral).doit(), eta).subs(eta, 0) / E_phase
    )

    # Ecuación de onda azimutal completa con acoplo de Hall/quiral de la pared de dominio
    mcfj_phi_lin = _poly_reduce_trig(
        r**2 * sp.diff(cov['MCFJ'][3].subs(sub_chiral).doit(), eta).subs(eta, 0) / E_phase
    )

    # Identificación exacta del número de onda quiral topológico k_chiral^{ax}(r) = g_{agg} d_r a_0 = g_{agg} f_a theta_0'(r)
    # y del potencial químico quiral efectivo mu_{5,eff}(r) mediante J_CME = (alpha_EM / pi) mu_5 B = k_chiral B:
    k_chiral_ax = g_agg * fa * d(th0)
    mu5_eff = (sp.pi / alpha_em) * k_chiral_ax
    k_crit_dp = 2 * k_chiral_ax
    k_max_dp = k_chiral_ax

    # Relación de dispersión de la inestabilidad magnética quiral (Dehman & Pons 2024, arXiv:2408.05281)
    # En régimen resistivo: i omega = eta_ohm k^2 - eta_ohm k_crit k  =>  gamma(k) = eta_ohm k (k_crit - k)
    gamma_chiral = eta_ohm * k_w * (k_crit_dp - k_w)

    return {
        'J_top_phi': J_top_phi,
        'J_top_th': J_top_th,
        'MCFJ_phi_chiral': mcfj_phi_lin,
        'k_chiral_ax': k_chiral_ax,
        'mu5_eff': mu5_eff,
        'k_crit': k_crit_dp,
        'k_max': k_max_dp,
        'gamma_chiral': gamma_chiral,
    }


@lru_cache
def derive_magneto_axionic():
    """Deriva el sistema acoplado de oscilaciones magneto-axiónicas (delta theta <-> calE) en espacio curvo."""
    cov = derive_covariant_mcfj()
    E_phase = sp.exp(-sp.I * w * t)
    Ath_bg = sp.Function('A_theta_bg')(r, th)

    # Perturbación acoplada dipolar (l = 1):
    #   a(t, r, theta) = a_0(r) + eta f_a deltatheta(r) cos(theta) e^{-i w t}
    #   delta A_phi(t, r, theta) = eta calE(r) sin^2(theta) e^{-i w t}
    # sobre el campo magnético de fondo F_{0, r theta} = d_r A_{theta, bg} = -r e^{lam0} B0(r) sin(theta)
    sub_pert = {
        At: 0,
        Ar: 0,
        Ath: Ath_bg,
        Aph: eta * calE * sp.sin(th)**2 * E_phase,
        Jt: 0,
        Jr: 0,
        Jth: 0,
        Jph: 0,
        a_fld: a0 + eta * fa * dth * sp.cos(th) * E_phase,
        n_fld: n0,
    }
    sub_B0 = {
        sp.Derivative(Ath_bg, r): -r * sp.exp(lam0) * B0 * sp.sin(th),
        sp.Derivative(Ath_bg, th): 0,
        Ath_bg: 0,
    }

    # 1. Parte electromagnética (operador de Maxwell dipolar + corriente axiónica i w g_{agg} f_a B0)
    mcfj_lin = sp.diff(cov['MCFJ'][3].subs(sub_pert).doit(), eta).subs(eta, 0).subs(sub_B0) / E_phase
    maxwell_part = _poly_reduce_trig(r**2 * mcfj_lin.subs(g_agg, 0))
    cs_em_part = _poly_reduce_trig((r**2 * (mcfj_lin - mcfj_lin.subs(g_agg, 0))) / sp.cos(th))
    eq_calE = sp.expand((maxwell_part + cs_em_part).doit())

    # 2. Parte de Klein-Gordon (operador escalar l = 1 + fuente electromagnética -i w g_{agg} B0 calE)
    kg_lin = sp.diff(cov['KG'].subs(sub_pert).doit(), eta).subs(eta, 0).subs(sub_B0) / E_phase
    kg_scalar_part = _poly_reduce_trig(-kg_lin.subs(g_agg, 0) / (fa * sp.cos(th)))
    cs_kg_part = _poly_reduce_trig((-(kg_lin - kg_lin.subs(g_agg, 0)) / fa) / sp.sin(th)**2)
    eq_dth = sp.expand((kg_scalar_part + cs_kg_part).doit())

    # 3. Flujo radial de energía electromagnética (vector de Poynting relativista T^r_t[EM])
    #    T^r_t[EM] = -F^{r mu} F_{t mu} = -F^{r phi} F_{t phi}
    F_dn_p, F_up_p, _ = field_tensors([0, 0, 0, calE * sp.sin(th)**2 * E_phase])
    Trt_em = _poly_reduce_trig(-F_up_p[1, 3] * F_dn_p[0, 3])

    return {
        'EQ_dth': eq_dth,
        'EQ_calE': eq_calE,
        'Trt_em': Trt_em,
    }


# ---------------------------------------------------------------- compilación numérica directa (sp.lambdify con cse=True)
BG_DIP = sp.symbols('bP1 bL bL1 j_src')
(bP1, bL, bL1, j_src) = BG_DIP

BG_MAG = sp.symbols('bP bP1 bL bL1 meff2 kapB')
(bP, bP1_m, bL_m, bL1_m, meff2_s, kapB_s) = BG_MAG


@lru_cache
def compiled_dipole_ode():
    """Compila la EDO relativista de dipolo de Wasserman-Shapiro desde SymPy con lambdify(..., cse=True).

    Devuelve (A_func, src_func) para y = [a_phi, a_phi'] con y' = A(r) y + [0, src].
    """
    ws = derive_static_dipole()['WS_dipole']
    s_aphi, s_daphi, s_ddaphi = sp.symbols('s_aphi s_daphi s_ddaphi')
    eq = ws.subs({
        d(aphi, 2): s_ddaphi,
        d(aphi): s_daphi,
        aphi: s_aphi,
        d(Phi0): bP1,
        d(lam0): bL1,
        lam0: bL,
        jphi: j_src,
    })
    c_dip = sp.expand(eq.coeff(s_ddaphi))
    rest_dip = sp.expand(eq - c_dip * s_ddaphi)
    A_mat = sp.Matrix([
        [0, 1],
        [sp.expand(-rest_dip.subs(j_src, 0).coeff(s_aphi) / c_dip), sp.expand(-rest_dip.subs(j_src, 0).coeff(s_daphi) / c_dip)],
    ])
    src_term = sp.expand(-(rest_dip - rest_dip.subs(j_src, 0)) / c_dip)
    f_A = sp.lambdify((r, bP1, bL, bL1), A_mat, 'numpy', cse=True)
    f_src = sp.lambdify((r, bP1, bL, bL1, j_src), src_term, 'numpy', cse=True)
    return f_A, f_src


@lru_cache
def compiled_magneto_axionic_matrix():
    """Compila la matriz 4x4 del sistema acoplado (dth, dth', calE_tilde, calE_tilde') desde SymPy con cse=True.

    Usa calE_tilde = calE / f_a y kapB = g_{agg} B0(r) (en km^-1) para que el sistema sea adimensional en f_a.
    """
    ma_eqs = derive_magneto_axionic()
    Va2 = sp.diff(V(a0), a0, 2)
    ma2 = sp.diff(mstar(a0), a0, 2)

    s_th, s_dth, s_ddth = sp.symbols('s_th s_dth s_ddth')
    s_E, s_dE, s_ddE = sp.symbols('s_E s_dE s_ddE')

    sub_common = {
        d(dth, 2): s_ddth,
        d(dth): s_dth,
        dth: s_th,
        d(calE, 2): fa * s_ddE,
        d(calE): fa * s_dE,
        calE: fa * s_E,
        Va2: meff2_s,
        ma2: 0,
        d(Phi0): bP1_m,
        d(lam0): bL1_m,
        Phi0: bP,
        lam0: bL_m,
        B0: kapB_s / g_agg,
    }

    eq1 = sp.expand(ma_eqs['EQ_dth'].subs(sub_common))
    eq2 = sp.expand(ma_eqs['EQ_calE'].subs(sub_common) / fa)

    c1 = sp.expand(eq1.coeff(s_ddth))
    rest1 = sp.expand(eq1 - c1 * s_ddth)
    c2 = sp.expand(eq2.coeff(s_ddE))
    rest2 = sp.expand(eq2 - c2 * s_ddE)

    ys = [s_th, s_dth, s_E, s_dE]
    A_mat = sp.Matrix([
        [0, 1, 0, 0],
        [sp.expand(-rest1.coeff(y) / c1) for y in ys],
        [0, 0, 0, 1],
        [sp.expand(-rest2.coeff(y) / c2) for y in ys],
    ])

    # Verificación de linealidad exacta: no queda término residual sin variables de estado
    for rest in (rest1, rest2):
        rem = sp.expand(rest - sum(rest.coeff(y) * y for y in ys))
        if rem != 0:
            raise ValueError(f'Residuo no lineal en matriz magneto-axiónica: {rem.free_symbols}')

    return sp.lambdify((r, w, *BG_MAG), A_mat, 'numpy', cse=True)
