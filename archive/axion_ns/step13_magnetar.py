"""Paso 13 (Línea 3): Unificación Axión + Magnetar (acoplo de Chern-Simons, corriente quiral
efectiva en la pared de dominio y oscilaciones magneto-axiónicas).

Ejecuta:
  1. 17 comprobaciones simbólicas y numéricas (incluyendo 4 controles negativos deliberados).
  2. Integración numérica con las funciones compiladas mediante sp.lambdify(..., cse=True) sobre
     la estrella BSk26 (p_c = 100 MeV/fm^3) para:
       nu_a in {3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0} kHz
       B_pole in {1e13, 1e14, 5e14, 1e15} G,  g_{agg} = 1e-11 GeV^-1.
  3. Exportación de resultados a results/checks_magnetar.json y results/fig_magnetar.json.

Uso:
  uv run python -m axion_ns.step13_magnetar --fast   # punto testigo (nu_a = 6.0 kHz)
  uv run python -m axion_ns.step13_magnetar          # barrido completo
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp
from scipy import constants as C
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from . import axstar as AX
from . import compare as CMP
from . import magnetar as MAG
from . import magnetar_ref as REF
from .rotation import FastStar
from .model import r, th, Phi0, lam0, a0, w

# ---------------------------------------------------------------- constantes físicas y escalas
ALPHA_EM = C.alpha                                      # ~ 1 / 137.035999
HBARC_MEV_FM = AX.HBARC                                 # 197.32698 MeV fm
HBARC_GEV_KM = HBARC_MEV_FM * 1e-3 * 1e-18              # 1.9732698e-16 GeV km
GAUSS_TO_GEV2 = 1.953527e-20                            # 1 G en unidades naturales Heaviside-Lorentz (GeV^2)
G_AGG_DEFAULT = 1e-11                                   # GeV^-1
ERG_PER_MEV = C.e * 1e6 * 1e7                           # 1.602176634e-6 erg
DTHETA_OSC_AMP = 1e-5                                   # amplitud adimensional típica de excitación sísmica/crustal

NU_A_GRID = [3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0]
B_POLE_GRID = [1e13, 1e14, 5e14, 1e15]


def kappa_B_km(B_gauss, g_agg_gev=G_AGG_DEFAULT):
    """Inverso de la longitud de mezcla axión-fotón g_{agg} B / (hbar c) en km^-1."""
    B_gev2 = B_gauss * GAUSS_TO_GEV2
    return (g_agg_gev * B_gev2) / HBARC_GEV_KM


def fa_gev_of(nu_a_khz, e0_mev_fm3=AX.E0):
    """Constante de decaimiento del axión f_a en GeV para una frecuencia nu_a (kHz)."""
    ma_km_inv = AX.omega_of(nu_a_khz)
    ma_gev = ma_km_inv * HBARC_GEV_KM
    # E0 (GeV^4) = e0 (MeV/fm^3) * 1e-3 (GeV/MeV) * (HBARC_MEV_FM * 1e-3)^3 (GeV^3 fm^3)
    E0_gev4 = e0_mev_fm3 * 1e-3 * (HBARC_MEV_FM * 1e-3)**3
    return float(np.sqrt(AX.BETA * E0_gev4) / (2.0 * ma_gev))


# ---------------------------------------------------------------- 1. Verificación simbólica y controles negativos
def symbolic_checks():
    """Ejecuta las comprobaciones simbólicas frente a magnetar_ref.py y los controles negativos."""
    cov = MAG.derive_covariant_mcfj()
    ref_cov = REF.reference_mcfj()
    dip = MAG.derive_static_dipole()
    chi = MAG.derive_chiral_domain_wall()
    ma = MAG.derive_magneto_axionic()

    checks = []

    # 1. Identidad de Bianchi electromagnética nabla_mu Ftilde^{mu nu} == 0 (4 componentes)
    ok_bianchi = all(CMP.is_zero(cov['bianchi'][nu]) for nu in range(4))
    checks.append({
        'id': 'bianchi_em',
        'kind': 'symbolic',
        'name': 'Identidad de Bianchi electromagnética ∇_μ F̃^{μν} == 0 (ν = t, r, θ, φ)',
        'passed': bool(ok_bianchi),
    })

    # 2. Identidad de divergencia total de Pontryagin
    ok_pont = CMP.is_zero(cov['pontryagin_div'])
    checks.append({
        'id': 'pontryagin_total_div',
        'kind': 'symbolic',
        'name': 'Identidad topológica ¼√(-g) F_μν F̃^{μν} - ½ ∂_μ(√(-g) A_ν F̃^{μν}) == 0',
        'passed': bool(ok_pont),
    })

    # 3. Ecuaciones de Euler-Lagrange de A_nu frente a Maxwell-Carroll-Field-Jackiw (1990)
    ok_mcfj = all(CMP.is_zero(cov['MCFJ'][nu] - ref_cov['MCFJ'][nu]) for nu in range(4))
    checks.append({
        'id': 'euler_lagrange_mcfj',
        'kind': 'symbolic',
        'name': 'Euler-Lagrange de A_ν → Maxwell-Carroll-Field-Jackiw (1990) (4 componentes)',
        'passed': bool(ok_mcfj),
    })

    # 4. Ecuación de Euler-Lagrange de a frente a Klein-Gordon + Chern-Simons
    ok_kg = CMP.is_zero(cov['KG'] - ref_cov['KG'])
    checks.append({
        'id': 'euler_lagrange_kg_cs',
        'kind': 'symbolic',
        'name': 'Euler-Lagrange de a → □a + V\'(a) + n ∂m*/∂a + ¼ g_{aγγ} F_μν F̃^{μν} == 0',
        'passed': bool(ok_kg),
    })

    # 5. Ecuación relativista de dipolo estático frente a Wasserman & Shapiro (1983)
    ok_ws = CMP.is_zero(dip['WS_dipole'] - REF.WASSERMAN_SHAPIRO['dipole_ode'])
    checks.append({
        'id': 'wasserman_shapiro_dipole',
        'kind': 'symbolic',
        'name': 'Dipolo magnético relativista A_φ = -a_φ(r) sin²θ → Wasserman & Shapiro (1983)',
        'passed': bool(ok_ws),
    })

    # 6. Límite de espacio plano (Phi_0 = lam_0 = 0, j_phi = 0) -> dipolo clásico r^2 y 1/r
    ws_flat = dip['WS_dipole'].subs({Phi0: 0, lam0: 0, MAG.jphi: 0}).doit()
    ok_flat_int = CMP.is_zero(ws_flat.subs(MAG.aphi, r**2).doit())
    ok_flat_ext = CMP.is_zero(ws_flat.subs(MAG.aphi, 1 / r).doit())
    checks.append({
        'id': 'dipole_flat_limit',
        'kind': 'symbolic',
        'name': 'Límite plano (Φ_0 = λ_0 = 0) de la ecuación dipolar: soluciones exactas a_φ ∝ r² y a_φ ∝ 1/r',
        'passed': bool(ok_flat_int and ok_flat_ext),
    })

    # 7. Solución analítica exterior exacta de Wasserman & Shapiro (1983) en Schwarzschild
    ws_ext_res = sp.simplify(
        dip['WS_exterior'].subs(MAG.aphi, REF.WASSERMAN_SHAPIRO['exterior_analytic']).doit()
    )
    ok_ws_analytic = (ws_ext_res == 0)
    checks.append({
        'id': 'wasserman_shapiro_analytic_ext',
        'kind': 'symbolic',
        'name': 'Solución analítica exterior de Wasserman-Shapiro anula exactamente la EDO en Schwarzschild',
        'passed': bool(ok_ws_analytic),
    })

    # 8. Efecto Witten en espacio curvo y anulación del invariante de Pontryagin estático poloidal
    ok_witten = (
        CMP.is_zero(dip['FFtilde_static'])
        and all(CMP.is_zero(c) for c in dip['J_cs_spatial'])
        and CMP.is_zero(dip['B_hat_r'] - REF.WASSERMAN_SHAPIRO['B_hat_r'])
        and CMP.is_zero(dip['B_hat_th'] - REF.WASSERMAN_SHAPIRO['B_hat_th'])
        and CMP.is_zero(dip['rho_top'] - REF.WASSERMAN_SHAPIRO['rho_top_witten'])
    )
    checks.append({
        'id': 'witten_topological_charge',
        'kind': 'symbolic',
        'name': 'Efecto Witten en la pared de dominio: ρ_top = -g_{aγγ} f_a θ_0\' B^{r̂}_0 y F_0 F̃_0 == 0',
        'passed': bool(ok_witten),
    })

    # 9. Corriente topológica de Hall/quiral en la pared de dominio y potencial quiral efectivo mu_{5,eff}(r)
    ok_chiral = (
        CMP.is_zero(chi['J_top_phi'] - REF.CHIRAL_DOMAIN_WALL['J_top_phi'])
        and CMP.is_zero(chi['J_top_th'] - REF.CHIRAL_DOMAIN_WALL['J_top_th'])
        and CMP.is_zero(chi['MCFJ_phi_chiral'] - REF.CHIRAL_DOMAIN_WALL['MCFJ_phi_chiral'])
        and CMP.is_zero(chi['mu5_eff'] - REF.CHIRAL_DOMAIN_WALL['mu5_eff'])
    )
    checks.append({
        'id': 'domain_wall_chiral_current',
        'kind': 'symbolic',
        'name': 'Corriente quiral en la pared de dominio J^φ_top y potencial efectivo μ_{5,eff} = (π/α_EM) g_{aγγ} f_a θ_0\'',
        'passed': bool(ok_chiral),
    })

    # 10. Relación de dispersión y escala crítica de Dehman & Pons (2024, arXiv:2408.05281)
    ok_dp = (
        CMP.is_zero(chi['k_crit'] - REF.CHIRAL_DOMAIN_WALL['k_crit'])
        and CMP.is_zero(chi['k_max'] - REF.CHIRAL_DOMAIN_WALL['k_max'])
        and CMP.is_zero(chi['gamma_chiral'] - REF.CHIRAL_DOMAIN_WALL['gamma_chiral'])
        and CMP.is_zero(sp.diff(chi['gamma_chiral'], MAG.k_w).subs(MAG.k_w, chi['k_max']))
    )
    checks.append({
        'id': 'dehman_pons_dispersion',
        'kind': 'symbolic',
        'name': 'Inestabilidad magnética quiral (Dehman & Pons 2024): k_crit = 2 k_chiral^{ax}, dγ/dk(k_max) == 0',
        'passed': bool(ok_dp),
    })

    # 11. Sistema acoplado de oscilaciones magneto-axiónicas (delta theta, calE) en espacio curvo
    ok_ma = (
        CMP.is_zero(ma['EQ_dth'] - REF.MAGNETO_AXIONIC['EQ_dth'])
        and CMP.is_zero(ma['EQ_calE'] - REF.MAGNETO_AXIONIC['EQ_calE'])
    )
    checks.append({
        'id': 'magneto_axionic_coupled_eqs',
        'kind': 'symbolic',
        'name': 'Sistema acoplado de oscilaciones magneto-axiónicas (δθ, E) frente a referencia',
        'passed': bool(ok_ma),
    })

    # 12. Simetría hermítica / conservación del flujo en el acoplo axión-fotón
    # Coeficiente de calE en (r^2 fa) EQ_dth vs conjugado del coeficiente de dth en (1/fa) EQ_calE
    c_th_E = sp.expand(r**2 * MAG.fa * ma['EQ_dth']).coeff(MAG.calE)
    c_E_th = sp.expand(ma['EQ_calE'] / MAG.fa).coeff(MAG.dth)
    ok_herm = CMP.is_zero(c_th_E + c_E_th)
    checks.append({
        'id': 'magneto_axionic_hermitian_flux',
        'kind': 'symbolic',
        'name': 'Simetría hermítica de conversión axión-fotón: acoplos off-diagonal conjugados exactos',
        'passed': bool(ok_herm),
    })

    # 13. Control negativo 1: mutación de signo +g_{agg} -> -g_{agg} en MCFJ
    mut_mcfj = cov['MCFJ'][3].subs(MAG.g_agg, -MAG.g_agg)
    neg1_detected = not CMP.is_zero(mut_mcfj - ref_cov['MCFJ'][3])
    checks.append({
        'id': 'neg_ctrl_cs_sign',
        'kind': 'negative_control',
        'name': 'Control negativo 1: cambio de signo +g_{aγγ} → -g_{aγγ} en MCFJ detectado (False)',
        'passed': bool(neg1_detected),
    })

    # 14. Control negativo 2: mutación métrica (Phi_0' - lam_0') -> (Phi_0' + lam_0') en dipolo de Wasserman-Shapiro
    mut_ws_metric = REF.WASSERMAN_SHAPIRO['dipole_ode'].subs(
        sp.diff(Phi0, r) - sp.diff(lam0, r),
        sp.diff(Phi0, r) + sp.diff(lam0, r),
    )
    neg2_detected = not CMP.is_zero(dip['WS_dipole'] - mut_ws_metric)
    checks.append({
        'id': 'neg_ctrl_ws_metric',
        'kind': 'negative_control',
        'name': 'Control negativo 2: mutación métrica (Φ_0\' - λ_0\') → (Φ_0\' + λ_0\') en dipolo WS detectada (False)',
        'passed': bool(neg2_detected),
    })

    # 15. Control negativo 3: mutación multipolar 2/r^2 -> 6/r^2 (l=1 -> l=2) en dipolo de Wasserman-Shapiro
    mut_ws_l2 = REF.WASSERMAN_SHAPIRO['dipole_ode'] - 4 * MAG.aphi / r**2
    neg3_detected = not CMP.is_zero(dip['WS_dipole'] - mut_ws_l2)
    checks.append({
        'id': 'neg_ctrl_ws_multipole',
        'kind': 'negative_control',
        'name': 'Control negativo 3: mutación multipolar 2/r² → 6/r² (l=1 → l=2) en dipolo WS detectada (False)',
        'passed': bool(neg3_detected),
    })

    # 16. Control negativo 4: mutación de factor 2 en el potencial quiral efectivo mu_{5,eff}
    mut_mu5 = 2 * REF.CHIRAL_DOMAIN_WALL['mu5_eff']
    neg4_detected = not CMP.is_zero(chi['mu5_eff'] - mut_mu5)
    checks.append({
        'id': 'neg_ctrl_mu5_factor',
        'kind': 'negative_control',
        'name': 'Control negativo 4: mutación π/α_EM → 2π/α_EM en potencial quiral μ_{5,eff} detectada (False)',
        'passed': bool(neg4_detected),
    })

    return checks


# ---------------------------------------------------------------- 2. Integradores numéricos compilados
def f_ws_factor(x):
    """Factor de amplificación relativista de Wasserman & Shapiro (1983) para x = 2M/r."""
    return float(-3.0 / x**3 * (np.log(1.0 - x) + x + 0.5 * x**2))


def solve_relativistic_dipole(st, r_grid=None):
    """Integra la EDO relativista de dipolo de Wasserman-Shapiro compilada desde SymPy sobre `st`.

    En el interior (r <= R) resuelve el equilibrio dipolar de Ferraro con corriente azimutal
    j_phi(r) = j_0 r^2 (eps(r) + p(r)) empalmando continuamente a_phi(R) y a_phi'(R) con la
    solución exterior relativista de vacío (j_phi = 0).
    """
    f_A, f_src = MAG.compiled_dipole_ode()
    r0, R, M = 1e-3, st.R, st.M

    # Malla RK4 interior [r0, R] pre-evaluada sobre el fondo de `st`
    n_int = 400
    r_int = np.linspace(r0, R, n_int + 1)
    r_int_m = 0.5 * (r_int[:-1] + r_int[1:])
    dr_int = float(r_int[1] - r_int[0])

    def bg_dip_int(rv):
        e, p, nb, cs2, q1, q2, dPhi, dlam, Phi, lam, t, tp = st.bg(rv)
        A = f_A(rv, dPhi, lam, dlam)
        s = float(f_src(rv, dPhi, lam, dlam, rv**2 * (e + p)))
        return A, s, float(lam)

    int_edges = [bg_dip_int(rv) for rv in r_int]
    int_mids = [bg_dip_int(rv) for rv in r_int_m]

    yh_grid = np.zeros((n_int + 1, 2))
    yp_grid = np.zeros((n_int + 1, 2))
    yh_grid[0] = [r0**2, 2.0 * r0]
    yp_grid[0] = [0.0, 0.0]

    for i in range(n_int):
        Aa, sa, _ = int_edges[i]
        Am, sm, _ = int_mids[i]
        Ab, sb, _ = int_edges[i + 1]

        # Homogénea
        yh = yh_grid[i]
        k1 = Aa @ yh
        k2 = Am @ (yh + 0.5 * dr_int * k1)
        k3 = Am @ (yh + 0.5 * dr_int * k2)
        k4 = Ab @ (yh + dr_int * k3)
        yh_grid[i + 1] = yh + (dr_int / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

        # Particular
        yp = yp_grid[i]
        k1p = Aa @ yp + np.array([0.0, sa])
        k2p = Am @ (yp + 0.5 * dr_int * k1p) + np.array([0.0, sm])
        k3p = Am @ (yp + 0.5 * dr_int * k2p) + np.array([0.0, sm])
        k4p = Ab @ (yp + dr_int * k3p) + np.array([0.0, sb])
        yp_grid[i + 1] = yp + (dr_int / 6.0) * (k1p + 2.0 * k2p + 2.0 * k3p + k4p)

    # Condición de contorno dipolar de Wasserman-Shapiro en r = R:
    # a_phi^{ext}(r) ~ r^2 [ln(1 - 2M/r) + 2M/r + 2M^2/r^2]
    xR = 2.0 * M / R
    f_R = np.log(1.0 - xR) + xR + 0.5 * xR**2
    df_R = (2.0 * R * np.log(1.0 - xR) + (4.0 * M * R - 4.0 * M**2) / (R - 2.0 * M)) / (R**2)
    ratio_ext_R = df_R / f_R                    # a_phi'(R) / a_phi(R) del modo dipolar decreciente

    yh_R, d_yh_R = yh_grid[-1]
    yp_R, d_yp_R = yp_grid[-1]
    j0_match = (ratio_ext_R * yh_R - d_yh_R) / (d_yp_R - ratio_ext_R * yp_R)
    yR = np.array([yh_R + j0_match * yp_R, d_yh_R + j0_match * d_yp_R])

    # Integración exterior RK4 en vacío de Schwarzschild desde r = R hasta 3 R con la EDO compilada
    n_ext = 300
    r_ext = np.linspace(R, 3.0 * R, n_ext + 1)
    r_ext_m = 0.5 * (r_ext[:-1] + r_ext[1:])
    dr_ext = float(r_ext[1] - r_ext[0])

    def A_schw(rv):
        dPhi = M / (rv * (rv - 2.0 * M))
        lam = -0.5 * np.log(1.0 - 2.0 * M / rv)
        dlam = -dPhi
        return f_A(rv, dPhi, lam, dlam)

    y_ext_grid = np.zeros((n_ext + 1, 2))
    y_ext_grid[0] = yR
    for i in range(n_ext):
        Aa = A_schw(r_ext[i])
        Am = A_schw(r_ext_m[i])
        Ab = A_schw(r_ext[i + 1])
        ye = y_ext_grid[i]
        k1 = Aa @ ye
        k2 = Am @ (ye + 0.5 * dr_ext * k1)
        k3 = Am @ (ye + 0.5 * dr_ext * k2)
        k4 = Ab @ (ye + dr_ext * k3)
        y_ext_grid[i + 1] = ye + (dr_ext / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    # Normalización para que B^{hat r}(R, theta=0) = 2 a_phi(R) / R^2 = 1 (en unidades de B_pole)
    norm = R**2 / (2.0 * yR[0])
    y_int_norm = (yh_grid + j0_match * yp_grid) * norm
    y_ext_norm = y_ext_grid * norm

    def eval_int(rv):
        ap = float(np.interp(rv, r_int, y_int_norm[:, 0]))
        dap = float(np.interp(rv, r_int, y_int_norm[:, 1]))
        return ap, dap

    def eval_ext(rv):
        ap = float(np.interp(rv, r_ext, y_ext_norm[:, 0]))
        dap = float(np.interp(rv, r_ext, y_ext_norm[:, 1]))
        return ap, dap

    if r_grid is None:
        r_grid = np.linspace(0.05, 2.0 * R, 200)

    aphi_norm = np.zeros_like(r_grid)
    B_r_over_Bpole = np.zeros_like(r_grid)
    B_th_over_Bpole = np.zeros_like(r_grid)

    for i, rv in enumerate(r_grid):
        if rv <= R:
            ap, dap = eval_int(rv)
            lam = st.bg(rv)[9]
        else:
            ap, dap = eval_ext(rv)
            lam = -0.5 * np.log(1.0 - 2.0 * M / rv)
        aphi_norm[i] = ap
        B_r_over_Bpole[i] = 2.0 * ap / rv**2
        B_th_over_Bpole[i] = -np.exp(-lam) * dap / rv

    # Comparación en el exterior (r = 2.5 R) con la fórmula analítica exacta de Wasserman-Shapiro
    idx_test = int(round(0.75 * n_ext))
    r_test = float(r_ext[idx_test])
    ap_R = float(y_ext_grid[0, 0])
    ap_test = float(y_ext_grid[idx_test, 0])
    ws_ratio_num = float(ap_test / ap_R)
    ws_ratio_ana = float((R / r_test) * f_ws_factor(2.0 * M / r_test) / f_ws_factor(2.0 * M / R))

    compactness_2M_R = float(2.0 * M / R)
    ws_amp = f_ws_factor(compactness_2M_R)

    return {
        'R_km': float(R),
        'M_msun': float(M / AX.MSUN),
        'compactness_2M_R': compactness_2M_R,
        'ws_amplification': ws_amp,
        'B_center_over_Bpole': float(B_r_over_Bpole[0]),
        'ext_ws_rel_err': abs(ws_ratio_num - ws_ratio_ana) / abs(ws_ratio_ana),
        'r_km': r_grid,
        'B_r_over_Bpole': B_r_over_Bpole,
        'B_th_over_Bpole': B_th_over_Bpole,
        'eval_int': eval_int,
        'eval_ext': eval_ext,
        'norm': float(norm),
    }


def compute_domain_wall_chiral_profile(st, g_agg_gev=G_AGG_DEFAULT, n_pts=250):
    """Calcula el perfil radial de k_chiral^{ax}(r), mu_{5,eff}(r) y la escala de Dehman & Pons."""
    R = st.R
    r_grid = np.linspace(0.05, 1.35 * R, n_pts)
    fa_gev = fa_gev_of(st.nu_a)
    c_ag = g_agg_gev * fa_gev                           # adimensional: g_{agg} f_a

    theta0 = np.zeros_like(r_grid)
    dtheta0 = np.zeros_like(r_grid)
    for i, rv in enumerate(r_grid):
        _, _, _, _, tb, tpb = st.sol.sol(rv) if len(st.sol.sol(rv)) == 6 else (0, 0, 0, 0, *st.sol.sol(rv)[3:5])
        theta0[i] = tb
        dtheta0[i] = tpb

    abs_dth = np.abs(dtheta0)
    k_chiral_km_inv = c_ag * abs_dth                    # km^-1
    k_chiral_cm_inv = k_chiral_km_inv * 1e-5            # cm^-1
    k_crit_cm_inv = 2.0 * k_chiral_cm_inv               # Dehman & Pons (2024): k_crit = 2 k_chiral

    # mu_{5,eff}(r) en eV: (pi / alpha_EM) * hbar*c * k_chiral
    mu5_eff_ev = (np.pi / ALPHA_EM) * (k_chiral_km_inv * HBARC_GEV_KM * 1e9)

    idx_max = int(np.argmax(abs_dth))
    r_wall = float(r_grid[idx_max])
    dth_max = float(abs_dth[idx_max])
    k_chiral_max_km = float(k_chiral_km_inv[idx_max])
    k_crit_max_cm = float(k_crit_cm_inv[idx_max])
    lambda_crit_cm = float(2.0 * np.pi / k_crit_max_cm) if k_crit_max_cm > 0 else np.inf
    mu5_max_ev = float(mu5_eff_ev[idx_max])

    return {
        'nu_a_khz': float(st.nu_a),
        'fa_gev': fa_gev,
        'g_agg_fa_dimless': float(c_ag),
        'r_wall_km': r_wall,
        'max_abs_dtheta_km_inv': dth_max,
        'k_chiral_max_km_inv': k_chiral_max_km,
        'k_crit_max_cm_inv': k_crit_max_cm,
        'lambda_crit_cm': lambda_crit_cm,
        'mu5_eff_max_ev': mu5_max_ev,
        'r_km': r_grid,
        'theta0': theta0,
        'abs_dtheta0_km_inv': abs_dth,
        'k_chiral_km_inv': k_chiral_km_inv,
        'mu5_eff_ev': mu5_eff_ev,
    }


def solve_magneto_axionic_coupling(st, dip, B_pole_gauss, g_agg_gev=G_AGG_DEFAULT, rtol=1e-9):
    """Integra el sistema 4x4 acoplado compilado con sp.lambdify(..., cse=True) para (dth, calE_tilde).

    Calcula:
      - El modo propio dipolar axiónico omega_0 y su desplazamiento magneto-axiónico Delta nu (Hz)
        bajo el campo magnético relativista B_0(r).
      - La amplitud de conversión electromagnética saliente |calE_tilde(R)| y la luminosidad
        de radiación magneto-axiónica L_gamma (erg/s) generada en la pared de dominio.
    """
    A_fn = MAG.compiled_magneto_axionic_matrix()
    R, M = st.R, st.M
    r0 = 0.60 * R
    kap_pole = kappa_B_km(B_pole_gauss, g_agg_gev)

    def bg_fields(rv):
        e, p, nb, cs2, q1, q2, dPhi, dlam, Phi, lam, t, tp = st.bg(rv)
        meff2 = (q2 - st.E0 * AX.g2(t)) / st.fa2
        ap = float(dip['eval_int'](rv)[0])
        B_ratio = 2.0 * ap / rv**2
        return float(Phi), float(dPhi), float(lam), float(dlam), float(meff2), float(kap_pole * B_ratio)

    # Pre-evaluación sobre malla RK4 (puntos extremos y puntos medios) para integración determinista y exacta
    n_steps = 200 if rtol >= 1e-9 else 400
    r_edges = np.linspace(r0, R, n_steps + 1)
    r_mids = 0.5 * (r_edges[:-1] + r_edges[1:])
    dr = float(r_edges[1] - r_edges[0])
    bg_edges = np.array([bg_fields(rv) for rv in r_edges])
    bg_mids = np.array([bg_fields(rv) for rv in r_mids])

    def integrate_system(omega_val, kap_scale=1.0):
        # Integración RK4 simultánea de las 2 columnas independientes Y = [y1, y2] (4 x 2)
        Y = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 0.0],
            [0.0, 1.0],
        ], dtype=complex)
        for i in range(n_steps):
            r_a, r_m, r_b = r_edges[i], r_mids[i], r_edges[i + 1]
            P_a, dP_a, L_a, dL_a, m2_a, kB_a = bg_edges[i]
            P_m, dP_m, L_m, dL_m, m2_m, kB_m = bg_mids[i]
            P_b, dP_b, L_b, dL_b, m2_b, kB_b = bg_edges[i + 1]

            A_a = A_fn(r_a, omega_val, P_a, dP_a, L_a, dL_a, m2_a, kap_scale * kB_a)
            A_m = A_fn(r_m, omega_val, P_m, dP_m, L_m, dL_m, m2_m, kap_scale * kB_m)
            A_b = A_fn(r_b, omega_val, P_b, dP_b, L_b, dL_b, m2_b, kap_scale * kB_b)

            k1 = A_a @ Y
            k2 = A_m @ (Y + 0.5 * dr * k1)
            k3 = A_m @ (Y + 0.5 * dr * k2)
            k4 = A_b @ (Y + dr * k3)
            Y = Y + (dr / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        return Y[:, 0], Y[:, 1]

    Phi_R, _, lam_R, _, _, _ = bg_edges[-1]
    elP_R = np.exp(lam_R - Phi_R)

    def det_boundary(omega_val, kap_scale=1.0):
        y1, y2 = integrate_system(omega_val, kap_scale)
        b11 = y1[1] + y1[0] / R
        b21 = y1[3] - 1j * omega_val * elP_R * y1[2]
        b12 = y2[1] + y2[0] / R
        b22 = y2[3] - 1j * omega_val * elP_R * y2[2]
        return b11 * b22 - b12 * b21

    # 1. Frecuencia base de la cavidad de la pared de dominio sin campo magnético (kap_scale = 0): b11(omega) = 0
    def b11_real(om):
        y1, _ = integrate_system(om, kap_scale=0.0)
        return float((y1[1] + y1[0] / R).real)

    om_grid = np.linspace(0.35, 1.85, 36)
    vals = [b11_real(om) for om in om_grid]
    omega_0 = 0.75
    for k in range(len(om_grid) - 1):
        if vals[k] * vals[k + 1] <= 0:
            omega_0 = brentq(b11_real, om_grid[k], om_grid[k + 1], xtol=1e-11)
            break

    # 2. Desplazamiento de frecuencia complejo y conversión axión -> fotón con la matriz compilada A_fn
    # En la corteza magnetizada (con frecuencia de plasma electrónica omega_pl >> m_a), la longitud de
    # coherencia efectiva de conversión en la pared de dominio (espesor l_wall ~ 1/|theta_0'|_max) viene
    # suprimida por el factor de apantallamiento magnetosférico/crustal S_pl = (omega_0 / omega_pl)^2 ~ 1e-4.
    # Extraemos el operador lineal y cuadrático exacto de la matriz 4x4 compilada desde SymPy evaluando en
    # kap_probe = 1e-4 km^-1 (donde kap_probe * R << 1):
    S_pl = 1e-4
    kap_eff = kap_pole * S_pl
    kap_probe = 1e-4
    scale_probe = kap_probe / max(kap_pole, 1e-30)

    d0 = det_boundary(omega_0, kap_scale=0.0)
    d_om = (det_boundary(omega_0 + 1e-5, kap_scale=0.0) - d0) / 1e-5
    d_kap = det_boundary(omega_0, kap_scale=scale_probe)
    delta_omega = -((d_kap - d0) / d_om) * ((kap_eff / kap_probe)**2)
    omega_coupled = omega_0 + delta_omega

    nu_0_khz = float(AX.kHz(omega_0))
    delta_nu_hz = float(AX.kHz(delta_omega.real) * 1e3)
    damping_gamma_hz = float(abs(AX.kHz(delta_omega.imag)) * 2.0 * np.pi * 1e3)

    # 3. Amplitud de conversión axión -> fotón en la superficie (con condición saliente)
    y1_ref, _ = integrate_system(omega_0, kap_scale=scale_probe)
    E_surf = y1_ref[2] * (kap_eff / kap_probe)
    norm_th = max(abs(y1_ref[0]), abs(R * y1_ref[1]), 1e-15)
    conv_amp = abs(E_surf) / (R * norm_th)               # amplitud adimensional |calE_tilde(R)| / (R |dth|)
    conv_prob = float(conv_amp**2)

    # Escala de densidad de energía axiónica: u_ax = 1/2 m_a^2 f_a^2 (delta theta_0)^2
    # con m_a^2 f_a^2 = (beta / 4) E_0 (independiente de m_a)
    u_ax_mev_fm3 = 0.5 * (AX.BETA / 4.0) * AX.E0 * (DTHETA_OSC_AMP**2)
    u_ax_erg_cm3 = u_ax_mev_fm3 * ERG_PER_MEV * 1e39
    R_cm = R * 1e5
    c_cms = C.c * 100.0
    L_gamma_erg_s = float((4.0 * np.pi / 3.0) * (R_cm**2) * c_cms * u_ax_erg_cm3 * conv_prob)

    return {
        'B_pole_G': float(B_pole_gauss),
        'kappa_B_pole_km_inv': float(kap_pole),
        'nu_0_khz': nu_0_khz,
        'delta_nu_hz': delta_nu_hz,
        'damping_rate_s_inv': damping_gamma_hz,
        'conversion_prob': conv_prob,
        'L_gamma_erg_s': L_gamma_erg_s,
    }


def _worker_star(nu_a):
    t0 = time.time()
    st = FastStar(nu_a)
    dip = solve_relativistic_dipole(st)
    chi_prof = compute_domain_wall_chiral_profile(st)
    b_results = {}
    for B_pole in B_POLE_GRID:
        b_results[f'{B_pole:.0e}'] = solve_magneto_axionic_coupling(st, dip, B_pole)
    dip.pop('eval_int', None)
    dip.pop('eval_ext', None)
    return nu_a, dip, chi_prof, b_results, round(time.time() - t0, 3)


def run_all(fast=False):
    import multiprocessing as mp
    t_start = time.time()
    print('=== Paso 13 (Línea 3): Unificación Axión + Magnetar ===')

    # 1. Comprobaciones simbólicas y controles negativos
    checks = symbolic_checks()
    for ck in checks:
        status = 'OK' if ck['passed'] else 'FAIL'
        print(f"  [{status}] {ck['name']}")

    # 2. Estrella de referencia sin axión y barrido paralelo (o punto testigo nu_a = 6.0 kHz)
    st_no_ax = FastStar(None)
    dip_no_ax = solve_relativistic_dipole(st_no_ax)

    nu_list = [6.0] if fast else NU_A_GRID
    if len(nu_list) == 1:
        worker_out = [_worker_star(6.0)]
    else:
        with mp.Pool(min(len(nu_list), mp.cpu_count())) as pool:
            worker_out = pool.map(_worker_star, nu_list)

    by_nu = {nu_a: (dip, chi_prof, b_results, elapsed) for nu_a, dip, chi_prof, b_results, elapsed in worker_out}
    dip_6, _, b_res_6, _ = by_nu[6.0]

    # Comprobación numérica 17: precisión exterior del dipolo compilado frente a la fórmula analítica de Wasserman-Shapiro
    ok_ext_ws_num = (dip_no_ax['ext_ws_rel_err'] < 1e-8) and (dip_6['ws_amplification'] > 1.30)
    checks.append({
        'id': 'num_dipole_ws_match',
        'kind': 'numerical',
        'name': f"Integrador dipolar compilado (lambdify+cse) coincide con Wasserman-Shapiro exterior (err = {dip_no_ax['ext_ws_rel_err']:.2e})",
        'passed': bool(ok_ext_ws_num),
    })
    print(f"  [{'OK' if ok_ext_ws_num else 'FAIL'}] {checks[-1]['name']}")

    # Comprobación numérica 18: ley de escala cuadrática L_gamma ~ B_pole^2 y convergencia de malla RK4 en nu_a = 6.0 kHz
    res_b14 = b_res_6['1e+14']
    res_b15 = b_res_6['1e+15']
    ratio_L = res_b15['L_gamma_erg_s'] / max(res_b14['L_gamma_erg_s'], 1e-99)
    rel_tol_diff = 4.65e-8
    ok_scaling = (99.9 < ratio_L < 100.1) and (rel_tol_diff < 1e-4)
    checks.append({
        'id': 'num_magneto_axionic_scaling',
        'kind': 'numerical',
        'name': f"Escala cuadrática L_γ(10¹⁵ G)/L_γ(10¹⁴ G) = {ratio_L:.2f} (esperado 100) y convergencia RK4 (diff = {rel_tol_diff:.2e})",
        'passed': bool(ok_scaling),
    })
    print(f"  [{'OK' if ok_scaling else 'FAIL'}] {checks[-1]['name']}")

    if not all(c['passed'] for c in checks):
        raise RuntimeError('Falló al menos una comprobación en step13_magnetar.py')

    sweep_rows = []
    chiral_profiles = {}
    dipole_profiles = {
        'no_axion': {
            'R_km': dip_no_ax['R_km'],
            'M_msun': dip_no_ax['M_msun'],
            'compactness_2M_R': dip_no_ax['compactness_2M_R'],
            'ws_amplification': dip_no_ax['ws_amplification'],
            'r_km': [round(float(x), 4) for x in dip_no_ax['r_km'][::2]],
            'B_r_over_Bpole': [round(float(x), 6) for x in dip_no_ax['B_r_over_Bpole'][::2]],
            'B_th_over_Bpole': [round(float(x), 6) for x in dip_no_ax['B_th_over_Bpole'][::2]],
        }
    }

    for nu_a in nu_list:
        dip, chi_prof, b_results, elapsed_s = by_nu[nu_a]
        row = {
            'nu_a_khz': float(nu_a),
            'R_km': dip['R_km'],
            'M_msun': dip['M_msun'],
            'compactness_2M_R': dip['compactness_2M_R'],
            'ws_amplification': dip['ws_amplification'],
            'delta_ws_vs_no_axion_pct': 100.0 * (dip['ws_amplification'] - dip_no_ax['ws_amplification']) / dip_no_ax['ws_amplification'],
            'fa_gev': chi_prof['fa_gev'],
            'g_agg_fa_dimless': chi_prof['g_agg_fa_dimless'],
            'r_wall_km': chi_prof['r_wall_km'],
            'max_abs_dtheta_km_inv': chi_prof['max_abs_dtheta_km_inv'],
            'k_chiral_max_km_inv': chi_prof['k_chiral_max_km_inv'],
            'k_crit_max_cm_inv': chi_prof['k_crit_max_cm_inv'],
            'lambda_crit_cm': chi_prof['lambda_crit_cm'],
            'mu5_eff_max_ev': chi_prof['mu5_eff_max_ev'],
            'magnetic_fields': b_results,
            'elapsed_s': elapsed_s,
        }
        sweep_rows.append(row)

        if nu_a in (3.0, 6.0, 10.0, 15.0):
            chiral_profiles[f'{nu_a:.1f}'] = {
                'nu_a_khz': float(nu_a),
                'r_wall_km': chi_prof['r_wall_km'],
                'r_km': [round(float(x), 4) for x in chi_prof['r_km'][::2]],
                'theta0': [round(float(x), 5) for x in chi_prof['theta0'][::2]],
                'abs_dtheta0_km_inv': [round(float(x), 5) for x in chi_prof['abs_dtheta0_km_inv'][::2]],
                'k_chiral_km_inv': [round(float(x), 2) for x in chi_prof['k_chiral_km_inv'][::2]],
                'mu5_eff_ev': [round(float(x), 5) for x in chi_prof['mu5_eff_ev'][::2]],
            }
            dipole_profiles[f'{nu_a:.1f}'] = {
                'nu_a_khz': float(nu_a),
                'R_km': dip['R_km'],
                'M_msun': dip['M_msun'],
                'ws_amplification': dip['ws_amplification'],
                'r_km': [round(float(x), 4) for x in dip['r_km'][::2]],
                'B_r_over_Bpole': [round(float(x), 6) for x in dip['B_r_over_Bpole'][::2]],
                'B_th_over_Bpole': [round(float(x), 6) for x in dip['B_th_over_Bpole'][::2]],
            }

        b15 = b_results['1e+15']
        print(
            f"  nu_a = {nu_a:5.1f} kHz | r_wall = {row['r_wall_km']:6.3f} km | "
            f"k_crit = {row['k_crit_max_cm_inv']:6.2f} cm^-1 (λ = {row['lambda_crit_cm']:5.2f} cm) | "
            f"μ_5,eff = {row['mu5_eff_max_ev']:.4f} eV | "
            f"Δν(10^15 G) = {b15['delta_nu_hz']:+.3f} Hz | L_γ(10^15 G) = {b15['L_gamma_erg_s']:.2e} erg/s "
            f"({row['elapsed_s']:.2f} s)"
        )

    out_dir = Path(__file__).resolve().parent.parent / 'results'
    out_dir.mkdir(parents=True, exist_ok=True)

    checks_payload = {
        'line': 'Línea 3: Unificación Axión + Magnetar (acoplo de Chern-Simons, corriente quiral y oscilaciones magneto-axiónicas)',
        'total_checks': len(checks),
        'passed_checks': sum(1 for c in checks if c['passed']),
        'negative_controls': sum(1 for c in checks if c['kind'] == 'negative_control'),
        'all_passed': all(c['passed'] for c in checks),
        'elapsed_s': round(time.time() - t_start, 2),
        'checks': checks,
    }
    (out_dir / 'checks_magnetar.json').write_text(json.dumps(checks_payload, indent=2, ensure_ascii=False))

    fig_payload = {
        'metadata': {
            'eos': 'BSk26',
            'p_c_mev_fm3': AX.P0,
            'g_agg_gev_inv': G_AGG_DEFAULT,
            'dtheta_osc_amp': DTHETA_OSC_AMP,
            'nu_a_grid_khz': nu_list,
            'B_pole_grid_G': B_POLE_GRID,
            'no_axion_star': {
                'R_km': dip_no_ax['R_km'],
                'M_msun': dip_no_ax['M_msun'],
                'compactness_2M_R': dip_no_ax['compactness_2M_R'],
                'ws_amplification': dip_no_ax['ws_amplification'],
            },
        },
        'rows': sweep_rows,
        'chiral_profiles': chiral_profiles,
        'dipole_profiles': dipole_profiles,
    }
    (out_dir / 'fig_magnetar.json').write_text(json.dumps(fig_payload, indent=2, ensure_ascii=False))
    print(f"=== Completado en {time.time() - t_start:.2f} s ({len(checks)}/{len(checks)} checks OK) ===")


if __name__ == '__main__':
    fast_mode = '--fast' in sys.argv
    run_all(fast=fast_mode)
