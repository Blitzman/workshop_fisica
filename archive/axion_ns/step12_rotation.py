"""Paso 12 (Línea 2): Rotación lenta O(Omega) de Hartle-Thorne, arrastre de Lense-Thirring y modos r con axión.

Verifica simbólicamente desde la acción (rotation.py vs rotation_ref.py) las 14 identidades exactas
(incluyendo 3 controles negativos deliberados) y ejecuta el puente numérico compilado sobre la estrella
BSk26 (p_c = 100 MeV/fm^3) sin axión y para nu_a in {3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0} kHz.
Guarda results/checks_rotacion.json y results/fig_rotacion.json.
"""
import datetime
import json
import multiprocessing as mp
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy
import sympy as sp

from . import rotation as ROT
from . import rotation_ref as REF
from .compare import is_zero
from .pipeline import report

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / 'results'
NUS_ROT = [3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0]


def symbolic_checks():
    """Comprobaciones simbólicas exactas (= 0) desde la acción y controles negativos obligatorios."""
    ht = ROT.derive_hartle_thorne()
    ax = ROT.derive_axial_rmodes()
    res = []

    # 1. T^t_{phi,axion} == 0 idénticamente
    ok1 = (ht['Ttphi_ax_1'] == 0)
    res.append((
        r'Acción axiónica en rotación lenta: T^t_{phi,axión} == 0 idénticamente para el condensado estático a_0(r)',
        ok1,
        'el escalar axisimétrico en m=0 no arrastra directamente el marco local'
    ))

    # 2. Ecuación de Hartle (1967, Ec. 43) en forma (Phi0' + lam0')
    ok2 = is_zero(ht['hartle_eq'] - REF.HARTLE_EQ)
    res.append((
        r"Ecuación de Hartle (1967, Ec. 43) para varpi(r) = Omega - omega_LT(r) desde G^t_phi = 8 pi G T^t_phi",
        ok2,
        r"varpi'' + (4/r - Phi_0' - lambda_0') varpi' - 16 pi G (eps + p) e^{2 lambda_0} varpi = 0"
    ))

    # 3. Forma con j(r) = exp(-(Phi0 + lam0)) y corrección cinética axiónica exacta
    ok3 = is_zero(ht['hartle_j_eq'] - REF.HARTLE_J_EQ_AXION)
    res.append((
        r"Identidad de Hartle con j(r) = e^{-(Phi_0+lambda_0)}: término nuevo -16 pi G j(r) (a_0')^2 varpi(r)",
        ok3,
        r"al sustituir Phi_0' + lambda_0' = 4 pi G r e^{2 lambda_0}(eps + p + e^{-2 lambda_0} a_0'^2)"
    ))

    # 4. Límite sin axión de la ecuación en forma j(r)
    noax_j = ht['hartle_j_eq'].subs(ROT.a0, 0).doit()
    ok4 = is_zero(noax_j - REF.HARTLE_J_EQ_NOAX)
    res.append((
        r'Límite sin axión (a_0 -> 0): recupera (1/r^4)(r^4 j varpi\')\' + (4/r) j\' varpi = 0 exactamente',
        ok4,
        ''
    ))

    # 5. Serie en el origen r -> 0
    ok5 = is_zero(ht['varpi_dd0'] - REF.VARPI_DD0_REF)
    res.append((
        r"Regularidad en el centro (r -> 0): varpi'(0) = 0 y varpi''(0) = (16 pi G / 5)(eps_c + p_c) varpi_c",
        ok5,
        ''
    ))

    # 6. Solución analítica exterior en vacío de Schwarzschild
    ext_res = sp.expand(ht['ext_eq'].subs(ROT.varpi, REF.VARPI_EXT).doit())
    ok6 = is_zero(ext_res)
    res.append((
        r'Exterior de Schwarzschild: varpi_ext(r) = Omega - 2 G I / r^3 satisface exactamente la ecuación en vacío',
        ok6,
        ''
    ))

    # 7. Desacoplo escalar de perturbaciones axiales (paridad (-1)^{l+1})
    ok7 = (ax['trace_h_ax'] == 0 and ax['div_u_ax'] == 0)
    res.append((
        r'Paridad axial (-1)^{l+1}: traza métrica g^{mu nu} h_{mu nu} == 0 y divergencia nabla_A delta u^A == 0',
        ok7,
        'las corrientes toroidales son incompresibles (delta n_b = delta a = 0 a orden axial)'
    ))

    # 8. Ecuación de vorticidad relativista de los modos r a orden O(Omega)
    ok8 = is_zero(ax['vorticity_eq'] - REF.KOJIMA_VORTICITY)
    res.append((
        r'Vorticidad radial nabla x (nabla_mu T^mu_A) = 0 a orden O(Omega): reproduce Kojima (1998, Ec. 20)',
        ok8,
        r'l(l+1) sigma (U - h_0) - 2 m varpi(r) U + 1/2 (l-1)(l+2) m omega_LT(r) h_0 = 0'
    ))

    # 9. Frecuencia corrotante del modo r y corrección relativista de arrastre kappa_{lm}(r)
    ok9 = is_zero(ax['sigma_corot'] - REF.SIGMA_RMODE_REF) and is_zero(ax['kappa_lm'] - REF.KAPPA_LM_REF)
    res.append((
        r'Relación de dispersión corrotante sigma = 2 m varpi(r)/[l(l+1)] y factor kappa_{lm}(r) = 1 - varpi/Omega',
        ok9,
        ''
    ))

    # 10. Ecuación maestra radial de Kojima (1998, Ec. 25) con condensado axiónico
    ok10 = is_zero(ax['kojima_raw'] - REF.KOJIMA_MASTER_RAW) and is_zero(ax['kojima_tov'] - REF.KOJIMA_MASTER_TOV)
    res.append((
        r'Ecuación maestra axial de Kojima (1998, Ec. 25) para h_0(r) modificada por el condensado axiónico',
        ok10,
        r'incluye los términos explícitos -4 pi G r a_0\'^2 h_0\' y +8 pi G a_0\'^2 h_0'
    ))

    # 11. Forma con masa encerrada M(r) y consistencia l = 1 con Hartle-Thorne
    l1_sub = ax['kojima_raw'].subs({ROT.l: 1, ROT.h0: -ROT.r**2 * (ROT.Om - ROT.varpi), ROT.U: -ROT.r**2 * ROT.Om}).doit()
    ok11 = is_zero(ax['kojima_M'] - REF.KOJIMA_MASTER_M) and is_zero(sp.expand(l1_sub / ROT.r**2) - ht['hartle_eq'])
    res.append((
        r'Consistencia multipolar: la ecuación maestra de Kojima en l = 1 se reduce exactamente a r^2 por la de Hartle',
        ok11,
        ''
    ))

    # 12. Control negativo 1: mutar 16 pi G -> 8 pi G en la ecuación de Hartle
    mut1 = REF.HARTLE_EQ + 8 * sp.pi * ROT.G * (REF.e + REF.p) * sp.exp(2 * ROT.lam0) * ROT.varpi
    neg1_detected = not is_zero(ht['hartle_eq'] - mut1)
    res.append((
        r'Control negativo 1: mutar 16 pi G -> 8 pi G en la ecuación de Hartle es detectado (devuelve False)',
        neg1_detected,
        'is_zero(derivada - mutada) == False'
    ))

    # 13. Control negativo 2: omitir el término cinético axiónico -16 pi G j(r) a_0'^2 varpi en la identidad de j(r)
    mut2 = REF.HARTLE_J_EQ_NOAX
    neg2_detected = not is_zero(ht['hartle_j_eq'] - mut2)
    res.append((
        r'Control negativo 2: omitir el término axiónico -16 pi G j(r) (a_0\')^2 varpi en la identidad de j(r) se detecta',
        neg2_detected,
        'is_zero(derivada - sin_correccion_axionica) == False'
    ))

    # 14. Control negativo 3: mutar l(l+1) -> 2(l+1) o el signo de a_0'^2 en la ecuación maestra de Kojima
    mut3 = REF.KOJIMA_MASTER_TOV + 16 * sp.pi * ROT.G * ROT.d(ROT.a0)**2 * ROT.h0
    neg3_detected = not is_zero(ax['kojima_tov'] - mut3)
    res.append((
        r'Control negativo 3: cambiar el signo del término axiónico 8 pi G (a_0\')^2 h_0 en Kojima se detecta',
        neg3_detected,
        'is_zero(derivada - mutada) == False'
    ))

    return res


def _worker(nu_a):
    st = ROT.FastStar(nu_a)
    rot = ROT.solve_rotation(st, l_val=2, n_pts=250)
    T_grid, f_crit = ROT.cfs_instability_window(rot)
    rot['T_grid_K'] = T_grid
    rot['f_crit_Hz'] = f_crit
    return rot


def compute_numerics(fast=False):
    nus = [6.0] if fast else NUS_ROT
    ref = _worker(None)
    if fast:
        rows = {6.0: _worker(6.0)}
    else:
        with mp.get_context('fork').Pool(min(6, len(nus))) as pool:
            rows_list = pool.map(_worker, nus)
        rows = {r_['nu_a']: r_ for r_ in rows_list}
    return ref, rows


def numeric_checks(ref, rows):
    res = []
    all_models = [ref] + [rows[k] for k in sorted(rows)]
    max_rel_I = max(m_['I_rel_diff'] for m_ in all_models)
    max_surf_diff = max(abs(m_['wLT_surface'] - m_['wLT_surface_analytic']) for m_ in all_models)

    # 15. Puente numérico e identidad integral de momento de inercia
    ok15 = (max_rel_I < 1e-6 and max_surf_diff < 1e-12)
    res.append((
        f'Puente numérico (lambdify + cse=True): I por EDO frente a la integral relativista coincide a {max_rel_I:.2e}',
        ok15,
        f'empalme superficial omega_LT(R)/Omega = 2GI/R^3 exacto a {max_surf_diff:.1e}'
    ))

    # 16. Estrella BSk26 de referencia sin axión
    ok16 = (1.50 < ref['I_45'] < 1.65 and 0.35 < ref['I_bar'] < 0.41 and 0.13 < ref['wLT_surface'] < 0.16)
    res.append((
        f"Estrella BSk26 sin axión (1.536 M_sun, 11.72 km): I = {ref['I_45']:.3f}x10^45 g cm^2, "
        f"I/(MR^2) = {ref['I_bar']:.4f}, omega_LT(0)/Omega = {ref['wLT_center']:.3f}, omega_LT(R)/Omega = {ref['wLT_surface']:.3f}",
        ok16,
        f"modo r (l=m=2): sigma/Omega = {ref['sigma_over_Omega']:.4f}, kappa_2 = {ref['kappa2_shift']:.4f}, tau_GW(1 kHz) = {ref['tau_gw_1kHz_s']:.1f} s"
    ))

    # 17. Efecto del condensado axiónico en I, w_LT y pared de dominio
    r6 = rows[6.0]
    dI_6 = 100.0 * (r6['I_45'] / ref['I_45'] - 1.0)
    max_wall = max(m_['max_wall_diff_pct'] for m_ in rows.values())
    ok17 = (dI_6 < -5.0 and r6['I_bar'] > ref['I_bar'] and max_wall > 0.01)
    res.append((
        f"Efecto del condensado axiónico (nu_a = 6 kHz): I baja un {abs(dI_6):.2f} % ({r6['I_45']:.3f}x10^45 g cm^2) "
        f"y I/(MR^2) sube a {r6['I_bar']:.4f}",
        ok17,
        f"ignorar el término -16 pi G j (a_0')^2 varpi desvía el perfil hasta un {max_wall:.2f} % en la pared"
    ))

    # 18. Ventana de inestabilidad CFS y tiempo de crecimiento de modos r
    fmin_ref = min(ref['f_crit_Hz'])
    fmin_6 = min(r6['f_crit_Hz'])
    dtau_6 = 100.0 * (r6['tau_gw_1kHz_s'] / ref['tau_gw_1kHz_s'] - 1.0)
    ok18 = (r6['tau_gw_1kHz_s'] > ref['tau_gw_1kHz_s'] and fmin_6 > fmin_ref)
    res.append((
        f"Inestabilidad CFS del modo r (l=m=2): el axión a 6 kHz alarga |tau_GW| un +{dtau_6:.1f} % "
        f"({r6['tau_gw_1kHz_s']:.1f} s vs {ref['tau_gw_1kHz_s']:.1f} s a 1 kHz) y eleva el umbral mínimo de espín "
        f"de {fmin_ref:.0f} Hz a {fmin_6:.0f} Hz",
        ok18,
        'estabiliza LMXBs rápidas en la banda T ~ 10^8-10^9 K'
    ))

    return res


def _round_obj(x, n=6):
    if isinstance(x, float):
        return round(x, n)
    if isinstance(x, list):
        return [_round_obj(v, n) for v in x]
    if isinstance(x, dict):
        return {k: _round_obj(v, n) for k, v in x.items()}
    return x


def save_outputs(checks_list, ref, rows, elapsed):
    RES.mkdir(parents=True, exist_ok=True)
    checks_payload = {
        'step': 'Paso 12: Rotación lenta de Hartle-Thorne, arrastre de Lense-Thirring y modos r con axión',
        'module': 'axion_ns.step12_rotation',
        'seconds': round(elapsed, 2),
        'n_ok': sum(1 for _, ok, _ in checks_list if ok),
        'n_total': len(checks_list),
        'env': {
            'python': sys.version.split()[0],
            'sympy': sp.__version__,
            'numpy': np.__version__,
            'scipy': scipy.__version__,
            'platform': f'{platform.system()} {platform.release()} ({platform.machine()})',
            'date': datetime.date.today().isoformat(),
        },
        'checks': [{'name': name, 'ok': bool(ok), 'note': note} for name, ok, note in checks_list],
    }
    (RES / 'checks_rotacion.json').write_text(json.dumps(checks_payload, indent=2, ensure_ascii=False))

    # Datos de figuras para informe_rotacion.html
    nus = sorted(rows)
    summary_table = []
    for m_ in [ref] + [rows[k] for k in nus]:
        summary_table.append({
            'nu_a': m_['nu_a'],
            'M_msun': round(m_['M_msun'], 4),
            'R_km': round(m_['R_km'], 3),
            'compactness': round(m_['compactness'], 4),
            'I_45': round(m_['I_45'], 4),
            'I_naive_45': round(m_['I_naive_45'], 4),
            'I_int_45': round(m_['I_int_45'], 4),
            'I_rel_diff': m_['I_rel_diff'],
            'I_bar': round(m_['I_bar'], 4),
            'wLT_center': round(m_['wLT_center'], 4),
            'wLT_surface': round(m_['wLT_surface'], 4),
            'max_wall_diff_pct': round(m_['max_wall_diff_pct'], 3),
            'kappa2_shift': round(m_['kappa2_shift'], 4),
            'varpi_eff': round(m_['varpi_eff'], 4),
            'sigma_over_Omega': round(m_['sigma_over_Omega'], 4),
            'tau_gw_1kHz_s': round(m_['tau_gw_1kHz_s'], 2),
            'tau_gw_600Hz_s': round(m_['tau_gw_600Hz_s'], 1),
            'f_crit_min_Hz': round(min(m_['f_crit_Hz']), 1),
        })

    # Perfiles radiales seleccionados (sin axión, 3, 4, 6, 10, 15 kHz) para la Gráfica 1 y la animación
    prof_keys = [k for k in [3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0] if k in rows]
    profiles = {
        'no_axion': _round_obj({
            'R_km': ref['R_km'],
            'M_msun': ref['M_msun'],
            'I_45': ref['I_45'],
            'sigma_over_Omega': ref['sigma_over_Omega'],
            'kappa2_shift': ref['kappa2_shift'],
            'tau_gw_1kHz_s': ref['tau_gw_1kHz_s'],
            'x': ref['x_grid'][::2],
            'wLT': ref['wLT_prof'][::2],
            'h0_over_U': ref['h0_over_U_prof'][::2],
            'theta_pi': ref['th_prof'][::2],
            'axion_kin': ref['axion_kin_prof'][::2],
        }, 5)
    }
    for k in prof_keys:
        rk = rows[k]
        profiles[str(k)] = _round_obj({
            'R_km': rk['R_km'],
            'M_msun': rk['M_msun'],
            'I_45': rk['I_45'],
            'sigma_over_Omega': rk['sigma_over_Omega'],
            'kappa2_shift': rk['kappa2_shift'],
            'tau_gw_1kHz_s': rk['tau_gw_1kHz_s'],
            'x': rk['x_grid'][::2],
            'wLT': rk['wLT_prof'][::2],
            'wLT_naive': rk['wLT_naive_prof'][::2],
            'h0_over_U': rk['h0_over_U_prof'][::2],
            'theta_pi': rk['th_prof'][::2],
            'axion_kin': rk['axion_kin_prof'][::2],
        }, 5)

    # Datos de LMXBs / púlsares de milisegundo observados para superponer en la ventana de inestabilidad CFS
    lmxb_obs = [
        {'name': 'PSR J0952−0607', 'T_K': 1.2e8, 'f_Hz': 707.3},
        {'name': '4U 1608−522', 'T_K': 4.5e8, 'f_Hz': 619.0},
        {'name': 'SAX J1750.8−2900', 'T_K': 3.4e8, 'f_Hz': 601.0},
        {'name': '4U 1636−536', 'T_K': 2.8e8, 'f_Hz': 581.0},
        {'name': 'Aql X-1', 'T_K': 3.1e8, 'f_Hz': 550.0},
        {'name': 'SAX J1808.4−3658', 'T_K': 1.5e8, 'f_Hz': 401.0},
    ]

    cfs_curves = {
        'T_grid_K': _round_obj(ref['T_grid_K'], 3),
        'no_axion': _round_obj(ref['f_crit_Hz'], 2),
    }
    for k in prof_keys:
        cfs_curves[str(k)] = _round_obj(rows[k]['f_crit_Hz'], 2)

    fig_payload = {
        'nus': nus,
        'summary': summary_table,
        'profiles': profiles,
        'cfs': cfs_curves,
        'lmxb': lmxb_obs,
    }
    (RES / 'fig_rotacion.json').write_text(json.dumps(fig_payload, indent=2, ensure_ascii=False))
    build_html_reports(checks_payload, fig_payload)


def build_html_reports(checks_payload, fig_payload):
    """Compila report/informe_rotacion.html con los datos verificados y copia a results/informe_rotacion.html."""
    import html
    import re

    rpt_path = ROOT / 'report' / 'informe_rotacion.html'
    if not rpt_path.exists():
        return
    src = rpt_path.read_text(encoding='utf-8')

    items_html = []
    for c in checks_payload['checks']:
        chip = '<span class="chip chip-ok">coincide</span>' if c['ok'] else '<span class="chip chip-bad">fallo</span>'
        note = f' <span class="muted">— {html.escape(c["note"])}</span>' if c['note'] else ''
        items_html.append(f'          <li>{chip}<span class="ck">{html.escape(c["name"])}{note}</span></li>')
    checks_block = '\n'.join(items_html)

    env = checks_payload['env']
    env_rows = [
        ('Python', env['python']),
        ('SymPy', env['sympy']),
        ('NumPy', env['numpy']),
        ('SciPy', env['scipy']),
        ('Sistema', env['platform']),
        ('Tiempo de ejecución', f"{checks_payload['seconds']} s"),
        ('Fecha', env['date']),
    ]
    env_block = ''.join(f'<dt>{k}</dt><dd>{html.escape(v)}</dd>' for k, v in env_rows)
    fig_json = json.dumps(fig_payload, separators=(',', ':'), ensure_ascii=False).replace('</', '<\\/')

    src = re.sub(r'<!--CHECKS_HTML-->.*?</ul>', f'<!--CHECKS_HTML-->\n{checks_block}\n        </ul>', src, flags=re.S)
    src = re.sub(r'<dl class="env" id="env-dl">.*?</dl>', f'<dl class="env" id="env-dl">{env_block}</dl>', src, flags=re.S)
    src = re.sub(r'<script type="application/json" id="figdata">.*?</script>',
                 f'<script type="application/json" id="figdata">{fig_json}</script>', src, flags=re.S)

    rpt_path.write_text(src, encoding='utf-8')
    (RES / 'informe_rotacion.html').write_text(src, encoding='utf-8')


def run(fast=False):
    t0 = time.time()
    sym_res = symbolic_checks()
    ref, rows = compute_numerics(fast=fast)
    num_res = numeric_checks(ref, rows)
    all_res = sym_res + num_res
    elapsed = time.time() - t0
    if not fast:
        save_outputs(all_res, ref, rows, elapsed)
    return report(f'Paso 12: Rotación lenta de Hartle-Thorne y modos r con axión ({elapsed:.1f} s)', all_res)


if __name__ == '__main__':
    is_fast = '--fast' in sys.argv
    run(fast=is_fast)
