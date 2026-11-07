"""Paso 11 — Estratificación de composición, modos de gravedad (g-modes), modo de interfaz
en la pared de dominio del axión (i_a) y resonancia de marea en LIGO / Einstein Telescope."""
import json
import sys
from multiprocessing import get_context
from pathlib import Path

import numpy as np
import sympy as sp

from . import axstar as AX
from . import gmodes as GM
from . import gmodes_ref as REF
from . import nonradial as NR
from .compare import is_zero
from .model import (
    G, Phi0, V, a0, eps0, lam0, mN, mstar, n0, p_of, r,
)
from .nonradial import dp, dth, fa, l, w
from .pipeline import report

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / 'results'
NU_LIST = [3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0]


def _worker_star(args):
    nu_a, l_val = args
    return GM.compute_star_gmodes(nu_a, l_val=l_val, n_grid=400, include_wall=True)


def run(fast=False):
    cov = GM.derive_covariant()
    th_data = GM.derive_thermodynamics()
    eqs = GM.derive_system()
    paper_bar = NR.PAPER
    wkb = GM.wkb_dispersion()
    orig = GM.origin_series()
    f_mat, M_sym = GM.compiled_matrix()

    out = []

    # 1. Transporte advectivo de composición congelada u^mu d_mu Y_p = 0
    out.append((
        'Composición congelada: u^μ ∂_μ Y_p = 0  ⇒  ΔY_p = δY_p + ξ^r Y_p0\' = 0',
        is_zero(th_data['dYp'] + GM.xi_r * sp.diff(GM.Yp0, r)),
        'transporte lagrangiano exacto desde la cuadrivelocidad perturbada u^μ',
    ))

    # 2. Discriminante relativista de Schwarzschild con axión A(r) (definición desde Δp = c_s^2(Δε - n dm*/da Δa))
    out.append((
        'Discriminante de Schwarzschild con axión A(r) = (ε\' - p\'/c_s^2 - n_b ∂m*/∂a a_0\')/(ε+p)',
        is_zero(th_data['A'] - REF.A_def),
        'deducido simbólicamente desde la relación adiabática lagrangiana congelada',
    ))

    # 3. Equivalencia algebraica A(r) = (1/c_e^2 - 1/c_s^2) p' / (ε+p)
    out.append((
        'Identidad termodinámica A(r) = (1/c_e^2 - 1/c_s^2) p\'/(ε+p)',
        is_zero(th_data['A'] - REF.A_ce) and is_zero(th_data['A'] - th_data['A_ce']),
        'separa la compresibilidad de equilibrio c_e^2 de la adiabática congelada c_s^2',
    ))

    # 4. Aceleración gravitatoria efectiva g_eff(r) módulo equilibrio hidrostático TOV
    out.append((
        'Gravedad efectiva con fuerza escalar: g_eff = Φ_0\' + n_b(∂m*/∂a)a_0\'/(ε+p) = -p\'/(ε+p)',
        is_zero(th_data['g_eff'] - REF.g_eff) and is_zero(GM.reduce_cowling(th_data['g_eff'] - th_data['g_eff_p'])),
        'identidad exacta módulo la ecuación hidrostática de fondo con axión',
    ))

    # 5. Frecuencia de Brunt-Väisälä relativista N^2(r)
    out.append((
        'Frecuencia de Brunt-Väisälä relativista N^2(r) = -g_eff(r) A(r) e^{2(Φ_0 - λ_0)}',
        is_zero(th_data['N2'] - REF.N2_def) and is_zero(GM.reduce_cowling(th_data['N2'] - REF.N2_ce)),
        'generalización de Reisenegger & Goldreich (1992) acoplada al condensado axiónico',
    ))

    # 6. Perturbaciones eulerianas δε y δn en materia estratificada con axión
    out.append((
        'Relación euleriana δε = δp/c_s^2 - (ε+p)A(r)ξ^r + n_b(∂m*/∂a)δa y δn/n = δp/(c_s^2(ε+p)) - A(r)ξ^r',
        is_zero(th_data['de'] - REF.de_euler) and is_zero(th_data['dn'] - REF.dn_euler),
        'término de flotabilidad -(ε+p)A(r)ξ^r y acoplo escalar n_b(∂m*/∂a)δa exactos',
    ))

    # 7-10. Sistema completo de 4 ecuaciones no barotrópicas en aproximación de Cowling
    for key in ('(1) momento radial', '(2) bariones', '(3) momento angular', '(4) KG interior'):
        diff = GM.reduce_cowling(eqs[key] - REF.REF_EQS[key])
        out.append((
            f'{key} no barotrópico vs referencia de Reisenegger-Goldreich + axión',
            is_zero(diff),
            'residuo simbólico idénticamente nulo tras separación en Y_{lm}',
        ))

    # 11. Límite barotrópico c_s^2 -> c_e^2 (A(r) -> 0) recupera Eqs. (8a)-(8d) del Paper II
    key_map = [
        ('(1) momento radial', '(8a) momento radial'),
        ('(2) bariones', '(8b) bariones'),
        ('(3) momento angular', '(8c) momento angular'),
        ('(4) KG interior', '(8d) KG interior'),
    ]
    bar_ok = all(
        is_zero(GM.reduce_cowling(eqs[k].subs({GM.cs2: GM.ce2, GM.dp: NR.dp}) - paper_bar[pk]))
        for k, pk in key_map
    )
    out.append((
        'Límite barotrópico c_s^2 → c_e^2 (A(r) → 0): recupera exactamente Eqs. (8a)–(8d) del Paper II',
        bar_ok,
        'las 4 ecuaciones colapsan idénticamente a las del fluido barotrópico del Paper II',
    ))

    # 12. Control negativo 1: invertir el signo de la flotabilidad en δε (o N^2 -> -N^2 en momento radial)
    eq1_wrong_sign = GM.reduce_cowling(
        eqs['(1) momento radial'] - 2 * (GM.eps + GM.p) * REF.N2_def * sp.exp(2 * (lam0 - Phi0)) * GM.W / r
    )
    detected_1 = not is_zero(eq1_wrong_sign - REF.REF_EQS['(1) momento radial'])
    out.append((
        'Control negativo 1: invertir signo de flotabilidad A(r) → -A(r) en momento radial',
        detected_1,
        'detectado (comparador devuelve False; cambiaría modos oscilatorios g por inestabilidad convectiva)',
    ))

    # 13. Control negativo 2: omitir el término no barotrópico -A(r)(ε+p)W en continuidad bariónica
    eq2_no_buoy = GM.reduce_cowling(eqs['(2) bariones'] + (GM.eps + GM.p) * REF.A_def * GM.W)
    detected_2 = not is_zero(eq2_no_buoy - REF.REF_EQS['(2) bariones'])
    out.append((
        'Control negativo 2: omitir término de estratificación -A(r)(ε+p)W en continuidad bariónica',
        detected_2,
        'detectado (comparador devuelve False; violaría la conservación covariante de bariones)',
    ))

    # 14. Control negativo 3: omitir el acoplo de flotabilidad -A(r)ξ^r en la fuente de Klein-Gordon
    eq4_no_buoy = GM.reduce_cowling(
        eqs['(4) KG interior'] - sp.diff(mstar(a0), a0) / fa * n0 * REF.A_def * sp.exp(-lam0) * GM.W / r
    )
    detected_3 = not is_zero(eq4_no_buoy - REF.REF_EQS['(4) KG interior'])
    out.append((
        'Control negativo 3: omitir acoplo de flotabilidad -(n_b ∂m*/∂θ)A(r)ξ^r en Klein-Gordon',
        detected_3,
        'detectado (comparador devuelve False; la fuente escalar depende de δn y por tanto de A(r)ξ^r)',
    ))

    # 15. Control negativo 4: reemplazar l(l+1) por l^2 en la relación de dispersión WKB
    wkb_wrong = wkb * l / (l + 1)
    detected_4 = not is_zero(wkb_wrong - REF.WKB_OMEGA2)
    out.append((
        'Control negativo 4: mutar l(l+1) → l^2 en el numerador de la relación de dispersión WKB',
        detected_4,
        'detectado (comparador devuelve False)',
    ))

    # 16. Regularidad en el centro r -> 0 (A(0) = 0, N^2(0) = 0, (N^2)'(0) = 0, δθ'(0) = 0, δp(0))
    reg_ok = (
        orig['A_at_0'] == 0
        and orig['N2_at_0'] == 0
        and orig['N2_r1'] == 0
        and orig['dth_p_0'] == [0]
        and is_zero(orig['dp_0'] - orig['dp_0_expected'])
        and orig['bar0_compatible']
    )
    out.append((
        'Regularidad en el origen r → 0: A(0) = 0, N^2(0) = 0, (N^2)\'(0) = 0, δθ\'(0) = 0 y Eq. (10a)',
        reg_ok,
        'el término de flotabilidad se anula como O(r^2) en el centro y preserva las condiciones iniciales',
    ))

    # 17. Relación de dispersión WKB local ω^2 = l(l+1) N^2(r) / (k_r^2 r^2 e^{2λ_0} + l(l+1))
    out.append((
        'Relación de dispersión WKB local ω^2 = l(l+1)N^2(r) / [k_r^2 r^2 e^{2λ_0} + l(l+1)]',
        is_zero(wkb - REF.WKB_OMEGA2),
        'deducida del determinante principal en el régimen de baja frecuencia e incompresible',
    ))

    # 18. Puente numérico directo sp.lambdify (cse=True) frente a axstar._perturb_rhs en el límite barotrópico
    st6 = GM.get_star(6.0)
    r_test = 6.0
    _, _, _, _, _, _, _, _, bg_t = GM.stratification_at(st6, r_test, include_wall=False)
    cs2_ax = st6.bg(r_test)[3]
    bg_bar = list(bg_t)
    bg_bar[10] = cs2_ax
    bg_bar[11] = cs2_ax  # c_s^2 = c_e^2 = cs2_ax
    w_test = AX.omega_of(1.8)
    M_num = np.asarray(f_mat(r_test, w_test, 2, *bg_bar), dtype=float)
    y_vec = np.array([1.2e-5, 3.4e-4, 0.08, -0.02])
    dy_lamb = M_num @ y_vec
    dy_ax = np.asarray(AX._perturb_rhs(st6, 2, w_test, True)(r_test, y_vec), dtype=float)
    rel_err = float(np.linalg.norm(dy_lamb - dy_ax) / np.linalg.norm(dy_ax))
    out.append((
        'Puente numérico directo sp.lambdify (cse=True) vs axstar._perturb_rhs en c_s^2 = c_e^2',
        rel_err < 1e-11,
        f'error relativo = {rel_err:.2e} en r = 6.0 km sobre BSk26 (ν_a = 6.0 kHz, l = 2)',
    ))

    # 19. Punto testigo (ν_a = 6.0 kHz, l = 2) o barrido completo en BSk26 (p_c = 100 MeV/fm^3)
    if fast:
        res_6 = GM.compute_star_gmodes(6.0, l_val=2, n_grid=400, include_wall=True)
        res_0 = GM.compute_star_gmodes(None, l_val=2, n_grid=400, include_wall=False)
        all_results = [res_0, res_6]
    else:
        tasks = [(None, 2)] + [(nu, 2) for nu in NU_LIST]
        with get_context('fork').Pool(min(len(tasks), 9)) as pool:
            all_results = pool.map(_worker_star, tasks)

    ax_results = [r_ for r_ in all_results if r_['nu_a_khz'] is not None]
    noax_res = next(r_ for r_ in all_results if r_['nu_a_khz'] is None)

    n_bulk_min = min(r_['N_bulk_hz'] for r_ in all_results)
    n_bulk_max = max(r_['N_comp_max_hz'] for r_ in all_results)
    n_wall_min = min(r_['N_max_hz'] for r_ in ax_results)
    n_wall_max = max(r_['N_max_hz'] for r_ in ax_results)
    bv_ok = (150.0 <= n_bulk_min <= n_bulk_max <= 250.0) and (500.0 <= n_wall_min <= n_wall_max <= 900.0)
    out.append((
        'Perfil de Brunt-Väisälä N(r)/(2π) en BSk26: fondo 150–250 Hz y pico axiónico 500–900 Hz',
        bv_ok,
        f'fondo nuclear = [{n_bulk_min:.1f}, {n_bulk_max:.1f}] Hz; pico en r_wall = [{n_wall_min:.1f}, {n_wall_max:.1f}] Hz',
    ))

    # 20. Ordenación Sturm-Liouville ν(i_a) > ν(g_1) > ν(g_2) > ν(g_3), conteo de nodos (0,1,2,3) y resonancia GW
    spec_ok = True
    for r_ in ax_results:
        m_ = r_['modes']
        if not (m_['i_a']['nu_hz'] > m_['g_1']['nu_hz'] > m_['g_2']['nu_hz'] > m_['g_3']['nu_hz'] > 0):
            spec_ok = False
        if not (m_['i_a']['nodes'] == 0 and m_['g_1']['nodes'] == 1 and m_['g_2']['nodes'] == 2 and m_['g_3']['nodes'] == 3):
            spec_ok = False
        if not (m_['i_a']['Q_nl'] > max(m_['g_1']['Q_nl'], m_['g_2']['Q_nl'], m_['g_3']['Q_nl']) > 0):
            spec_ok = False

    r6 = next(r_ for r_ in ax_results if r_['nu_a_khz'] == 6.0)
    m6 = r6['modes']
    out.append((
        'Espectro (i_a, g_1, g_2, g_3), nodos (0, 1, 2, 3) y desfase resonante ΔΦ_GW en LIGO A+ / ET',
        spec_ok and m6['i_a']['snr_ligo_aplus'] > 1.0 and m6['g_1']['snr_et_d'] > 1.0,
        f"ν_a=6 kHz: ν(i_a)={m6['i_a']['nu_hz']:.1f} Hz (f_GW={m6['i_a']['f_gw_hz']:.1f} Hz, "
        f"Q={m6['i_a']['Q_nl']:.2e}, ΔΦ_GW={m6['i_a']['dPhi_gw_rad']:.3f} rad, "
        f"SNR_LIGO={m6['i_a']['snr_ligo_aplus']:.1f}, SNR_ET={m6['i_a']['snr_et_d']:.1f})",
    ))

    rep = report('Paso 11 — Modos de gravedad (g-modes), modo de interfaz axiónica y marea GW', out)

    if not fast:
        RESULTS.mkdir(exist_ok=True)
        checks_payload = {
            'step11_gmodes': [[name, bool(ok), note] for name, ok, note in rep],
        }
        (RESULTS / 'checks_gmodes.json').write_text(json.dumps(checks_payload, indent=2, ensure_ascii=False))

        f_grid = np.geomspace(100.0, 1000.0, 120)
        thr_ligo, thr_et = GM.detector_phase_threshold(f_grid)
        fig_payload = {
            'meta': {
                'eos': 'BSk26',
                'pc_mev_fm3': 100.0,
                'l': 2,
                'M_binary_msun': [GM.M_BIN_MSUN, GM.M_BIN_MSUN],
                'delta_comp': GM.DELTA_COMP,
                'chi_wall': GM.CHI_WALL,
                'nu_a_list_khz': NU_LIST,
            },
            'no_axion': noax_res,
            'stars': ax_results,
            'detector_curves': {
                'f_gw_hz': [round(float(x), 2) for x in f_grid],
                'dPhi_min_ligo_aplus_rad': [round(float(x), 5) for x in thr_ligo],
                'dPhi_min_et_d_rad': [round(float(x), 6) for x in thr_et],
            },
        }
        (RESULTS / 'fig_gmodes.json').write_text(json.dumps(fig_payload, indent=2, ensure_ascii=False))
        print(f'\nGuardados: {RESULTS / "checks_gmodes.json"} y {RESULTS / "fig_gmodes.json"}')

    return rep


if __name__ == '__main__':
    is_fast = '--fast' in sys.argv
    res = run(fast=is_fast)
    if not all(ok for _, ok, _ in res):
        sys.exit(1)
