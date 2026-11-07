"""Ecuaciones de referencia de la literatura para la Línea 3 (Axión + Magnetar).

Referencias transcritas de forma independiente (forma 'lhs - rhs' = 0):
  1. S. M. Carroll, G. B. Field & R. Jackiw, Phys. Rev. D 41, 1231 (1990):
     Electrodinámica de Maxwell-Chern-Simons en espacio-tiempo curvo:
       (1 / sqrt(-g)) d_mu (sqrt(-g) F^{mu nu}) + g_{agg} (d_mu a) Ftilde^{mu nu} - J^nu = 0
     y ecuación de Klein-Gordon acoplada a la densidad de Pontryagin:
       Box a + V'(a) + n m*'(a) + 1/4 g_{agg} F_{mu nu} Ftilde^{mu nu} = 0.
  2. I. Wasserman & S. L. Shapiro, Astrophys. J. 265, 1036 (1983):
     Ecuación relativista del dipolo magnético estático A_phi = -a_phi(r) sin^2 theta:
       e^{-2 lambda_0} a_phi'' + e^{-2 lambda_0} (Phi_0' - lambda_0') a_phi' - 2 / r^2 a_phi + 4 pi j_phi(r) = 0
     y solución analítica exterior exacta en la métrica de Schwarzschild (1 - 2M/r):
       a_phi^{ext}(r) = -(3 mu_B r^2) / (8 M^3) [ ln(1 - 2M/r) + 2M/r + 2M^2/r^2 ].
  3. C. Dehman & J. A. Pons, Astrophys. J. 975, 100 (2024, arXiv:2408.05281) y arXiv:2505.06196 (2025):
     Inestabilidad magnética quiral (Chiral Magnetic Effect, J_CME = (alpha_EM / pi) mu_5 B)
     y correspondencia exacta con la corriente topológica en la pared de dominio del axión:
       k_chiral^{ax}(r) = g_{agg} f_a theta_0'(r),
       mu_{5,eff}(r)    = (pi / alpha_EM) g_{agg} f_a theta_0'(r),
       k_crit           = (2 alpha_EM / pi) mu_{5,eff} = 2 g_{agg} f_a theta_0'(r),
       k_max            = k_crit / 2,
       gamma(k)         = eta_ohm k (k_crit - k).
  4. Sistema acoplado de oscilaciones magneto-axiónicas (conversión axión-fotón / Raffelt & Stodolsky
     1988, Phys. Rev. D 37, 1237, en métrica esférica TOV):
       e^{-2 lambda_0} dth'' + e^{-2 lambda_0} (Phi_0' - lambda_0' + 2/r) dth'
         + (omega^2 e^{-2 Phi_0} - m_eff^2 - 2/r^2) dth - i omega (g_{agg} B_0 / (f_a r)) e^{-Phi_0} calE = 0,
       e^{-2 lambda_0} calE'' + e^{-2 lambda_0} (Phi_0' - lambda_0') calE'
         + (omega^2 e^{-2 Phi_0} - 2/r^2) calE + i omega g_{agg} f_a B_0 r e^{-Phi_0} dth = 0.
"""
import sympy as sp

from .model import r, th, X, Phi0, lam0, a0, n0, w, V, mstar
from .magnetar import (
    g_agg, alpha_em, fa, M_ext, mu_mag, k_w, eta_ohm,
    A_vec, J_vec, a_fld, n_fld, aphi, jphi, th0, B0, dth, calE, calB,
    ginv, sqrtg, field_tensors, d
)


def reference_mcfj():
    """Ecuaciones covariantes de Maxwell-Carroll-Field-Jackiw (1990) y Klein-Gordon acoplada."""
    F_dn, F_up, Ftilde_up = field_tensors(A_vec)
    da = [sp.diff(a_fld, x) for x in X]

    # (1/sqrt(-g)) d_mu(sqrt(-g) F^{mu nu}) + g_{agg} (d_mu a) Ftilde^{mu nu} - J^nu = 0
    mcfj_ref = []
    for nu in range(4):
        div_F = sum(sp.diff(sqrtg * F_up[mu, nu], X[mu]) for mu in range(4)) / sqrtg
        cs_cur = g_agg * sum(da[mu] * Ftilde_up[mu, nu] for mu in range(4))
        mcfj_ref.append(sp.expand(div_F + cs_cur - J_vec[nu]))

    # Box a + V'(a) + n dm*/da + 1/4 g_{agg} F_{mu nu} Ftilde^{mu nu} = 0
    box_a = sum(sp.diff(sqrtg * ginv[mu, mu] * da[mu], X[mu]) for mu in range(4)) / sqrtg
    fftilde = sum(F_dn[mu, nu] * Ftilde_up[mu, nu] for mu in range(4) for nu in range(4))
    kg_ref = sp.expand(
        box_a + sp.diff(V(a_fld), a_fld) + n_fld * sp.diff(mstar(a_fld), a_fld)
        + sp.Rational(1, 4) * g_agg * fftilde
    )

    return {'MCFJ': mcfj_ref, 'KG': kg_ref}


# ---- Wasserman & Shapiro (1983), ApJ 265, 1036: dipolo relativista estático
e2l = sp.exp(-2 * lam0)
e2P = sp.exp(-2 * Phi0)

WASSERMAN_SHAPIRO = {
    'dipole_ode': e2l * d(aphi, 2) + e2l * (d(Phi0) - d(lam0)) * d(aphi) - 2 * aphi / r**2 + 4 * sp.pi * jphi,
    'exterior_ode': (1 - 2 * M_ext / r) * d(aphi, 2) + 2 * M_ext / r**2 * d(aphi) - 2 * aphi / r**2,
    'exterior_analytic': -3 * mu_mag * r**2 / (8 * M_ext**3) * (
        sp.log(1 - 2 * M_ext / r) + 2 * M_ext / r + 2 * M_ext**2 / r**2
    ),
    'B_hat_r': 2 * aphi * sp.cos(th) / r**2,
    'B_hat_th': -sp.exp(-lam0) * d(aphi) * sp.sin(th) / r,
    'rho_top_witten': -g_agg * fa * sp.exp(-lam0) * d(th0) * (2 * aphi * sp.cos(th) / r**2),
}

# ---- Dehman & Pons (2024/2025), arXiv:2408.05281: corriente quiral efectiva en la pared de dominio
CHIRAL_DOMAIN_WALL = {
    'J_top_phi': sp.I * w * g_agg * fa * d(th0) * calB / (sp.exp(Phi0 + lam0) * r**2),
    'J_top_th': -sp.I * w * g_agg * fa * d(th0) * calE * sp.sin(th) / (sp.exp(Phi0 + lam0) * r**2),
    'MCFJ_phi_chiral': (
        e2l * d(calE, 2) + e2l * (d(Phi0) - d(lam0)) * d(calE)
        + (w**2 * e2P - 2 / r**2) * calE
        + sp.I * w * sp.exp(-Phi0 - lam0) * g_agg * fa * d(th0) * calB
    ),
    'k_chiral_ax': g_agg * fa * d(th0),
    'mu5_eff': (sp.pi / alpha_em) * g_agg * fa * d(th0),
    'k_crit': 2 * g_agg * fa * d(th0),
    'k_max': g_agg * fa * d(th0),
    'gamma_chiral': eta_ohm * k_w * (2 * g_agg * fa * d(th0) - k_w),
}

# ---- Sistema acoplado de oscilaciones magneto-axiónicas en espacio-tiempo curvo
Va2 = sp.diff(V(a0), a0, 2)
ma2 = sp.diff(mstar(a0), a0, 2)
meff2_ref = Va2 + n0 * ma2

MAGNETO_AXIONIC = {
    'EQ_dth': (
        e2l * d(dth, 2) + e2l * (d(Phi0) - d(lam0) + 2 / r) * d(dth)
        + (w**2 * e2P - meff2_ref - 2 / r**2) * dth
        - sp.I * w * g_agg * B0 * sp.exp(-Phi0) * calE / (fa * r)
    ),
    'EQ_calE': (
        e2l * d(calE, 2) + e2l * (d(Phi0) - d(lam0)) * d(calE)
        + (w**2 * e2P - 2 / r**2) * calE
        + sp.I * w * g_agg * fa * B0 * r * sp.exp(-Phi0) * dth
    ),
}
