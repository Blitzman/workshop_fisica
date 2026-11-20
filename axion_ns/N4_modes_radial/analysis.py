"""N4: the radial analysis (old steps 18, 21, 23 and 24): post-processing and controls on the radial modes, as extra parts of N4 (`--part 18,21,23,24`).

Port of axion_ns/step18_rydberg.py, step21_deep_artifacts.py, step23_eos_consistency.py and step24_scattering_width.py with the same algorithms and check texts. They read what N4's own
steps wrote: 18, 23 and 24 take the modes of step 15 (radial_star.json), 21 the deep zeros of step 20 (radial_damped.json), both from the output directory. Differences with the legacy:
  * step 21: the alternative matching radii are 0.85 and 0.7 of the matching radius (the legacy used the first |theta_0| < 1e-5 and 1e-6, which with a tail valid to r_end may not exist);
    matching at the surface (r_m = R) means no tail;
  * step 23: the c_s^2 definition is a switch of the EOS evaluator (star_table._EOSVec.cs2_mode) and the coefficient tables are rebuilt for each definition (the legacy patched `bsk26.EOS.cs2`
    and re-integrated the stars, which gives the same background: the star integration does not use c_s^2); F (6 kHz) is found with the pole-rejecting `real_root`; the P(rho) shifted by 1% in
    lg rho of the negative control is built by interpolating the table's own P(rho) instead of calling the analytic fit;
  * step 24: the real-frequency scans use N4's fork pool (`--procs`).
"""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit, least_squares

from common.checks import ok_field

from . import modes as MD
from . import validation as V
from .star_table import _EOSVec
from .tables import RadialTables

omega_of = MD.omega_of
MEV_FM3_TO_GCC = 1.782662e12            # constant of bskfit18.f, as axion_ns/bsk26.py


def _need(out, name, step):
    p = Path(out) / name
    if not p.exists():
        raise FileNotFoundError(f'{p} not found: it is the output of step {step}; run `python -m N4_modes_radial.run --part {step}` first')
    return json.loads(p.read_text(encoding='utf-8'))


# =============================================================================================== step 18
def kappa_whittaker(alpha, x):
    q = np.sqrt(1 - x * x)
    return alpha * (1 - 2 * q * q) / q


def n_hydrogen(alpha, x):                      # the non-relativistic formula used until then
    return alpha / np.sqrt(2 * (1 - x))


def kappa_newton(alpha, x):                    # negative control: Newtonian Coulomb with the exact q, without the 2 omega^2
    return alpha / np.sqrt(1 - x * x)


def step18(out):
    S = _need(out, 'radial_star.json', 15)
    a, nua = S['nu6']['alpha_g'], 6.0
    ns = np.array([l['n'] for l in S['nu6']['ladder']], dtype=float)
    xs = np.array([l['nu'] for l in S['nu6']['ladder']]) / nua
    d_w = ns - kappa_whittaker(a, xs)
    d_h = ns - n_hydrogen(a, xs)
    d_n = ns - kappa_newton(a, xs)
    rng = lambda d: float(d.max() - d.min())
    fit = least_squares(lambda p: p[0] + p[1] / (ns - p[0])**2 - d_w, [0.05, 0.0])
    res = [
        (f'H20(a): con el kappa de Whittaker (Ec. A14), n - kappa_n = {", ".join(f"{v:.4f}" for v in d_w)} para '
         f'n = {int(ns[0])}…{int(ns[-1])}: constante dentro de {rng(d_w):.4f}', rng(d_w) < 0.005, 'sin ajuste'),
        (f'H20(b): Rydberg-Ritz delta_n = delta_0 + delta_2/(n - delta_0)^2 con delta_0 = {fit.x[0]:.4f}, '
         f'delta_2 = {fit.x[1]:.4f}; residuo máximo {abs(fit.fun).max():.0e}', abs(fit.fun).max() < 1e-4, '2 parámetros, 5 niveles'),
        (f'Control negativo 1: con la fórmula hidrogenoide no relativista, n - n_eff varía {rng(d_h):.3f} '
         f'(lo que parecía un defecto creciente)', rng(d_h) > 0.02, ''),
        (f'Control negativo 2: con el Coulomb newtoniano alpha/q (sin el factor 2 omega^2 - m_a^2), n - kappa varía '
         f'{rng(d_n):.3f}', rng(d_n) > 0.02, ''),
    ]
    o = {'alpha_g': a, 'n': ns.tolist(), 'x': xs.tolist(), 'n_minus_kappa': d_w.tolist(), 'n_minus_neff': d_h.tolist(), 'n_minus_kappa_newton': d_n.tolist(),
         'ritz': {'delta0': float(fit.x[0]), 'delta2': float(fit.x[1]), 'resid': float(abs(fit.fun).max())},
         'checks': [{'name': n_, 'ok': ok_field(o_), 'note': t_} for n_, o_, t_ in res]}
    return res, o


# =============================================================================================== step 21
_S21 = {}


def _F21(nu):
    return complex(MD.det_entire(_S21['st'], _S21['tabs'], omega_of(nu), rm=_S21.get('rm'), inward=False, sheet=_S21.get('sheet', 1)))


def _root21(z, rm=None):
    _S21['rm'] = rm
    try:
        return complex(MD.muller(_F21, z, tol=1e-9))
    except RuntimeError:
        return complex(np.nan, np.nan)
    finally:
        _S21.pop('rm')


def step21(out):
    res, o = [], {}
    D = _need(out, 'radial_damped.json', 20)
    st, tabs = V.get_star(6.0)
    _S21.update(st=st, tabs=tabs)
    r0 = float(MD.matching_radius(st))
    alt = {'alt_0.85rm': 0.85 * r0, 'alt_0.7rm': 0.7 * r0}
    o['radii'] = {'default': r0, **alt, 'R': float(st.R), 'theta_1e-5': V.r_theta_below(st, r0, 1e-5), 'theta_1e-6': V.r_theta_below(st, r0, 1e-6)}
    ka, kb = list(alt)
    # (i) control: the fluid modes
    fl = []
    for z in (complex(*v) for v in D['nu6']['fluid_modes']):
        z0 = _root21(z)
        za, zb, zR = _root21(z0, alt[ka]), _root21(z0, alt[kb]), _root21(z0, float(st.R))
        sh = lambda w: float(100 * abs(w / z0 - 1))
        fl.append({'nu': [z0.real, z0.imag], 'shift_pct': {ka: sh(za), kb: sh(zb), 'R': sh(zR)},
                   'dIm_Hz': {ka: float(1e3 * (za - z0).imag), kb: float(1e3 * (zb - z0).imag), 'R': float(1e3 * (zR - z0).imag)}})
    # (ii) deep zeros, without duplicates
    deep = []
    for r_ in D['deep']['roots']:
        z = complex(*r_['nu'])
        if all(abs(z - complex(*d['nu'])) > 0.05 for d in deep):
            z0 = _root21(z)
            za, zb = _root21(z0, alt[ka]), _root21(z0, alt[kb])
            deep.append({'nu': [z0.real, z0.imag], f'nu_{ka}': [za.real, za.imag], f'nu_{kb}': [zb.real, zb.imag],
                         'shift_pct': {ka: float(100 * abs(za / z0 - 1)), kb: float(100 * abs(zb / z0 - 1))}, 'fluid': r_['fluid']})
    o['fluid'], o['deep'] = fl, deep
    ok_i = all(max(f['shift_pct'][ka], f['shift_pct'][kb]) < 0.1 for f in fl)
    res.append(('H23b(i), control: los modos fluidos casi no se mueven al cambiar el radio de empalme: ' + '; '.join(
        f"(Re {f['nu'][0]:.4f} kHz, Im {f['nu'][1] * 1e3:.2f} Hz): {f['shift_pct'][ka]:.2g} % y {f['shift_pct'][kb]:.2g} % "
        f"(Im: {f['dIm_Hz'][ka]:+.2g} y {f['dIm_Hz'][kb]:+.2g} Hz)" for f in fl), ok_i, 'radios alternativos: 0.85 y 0.7 del de empalme'))
    ok_ii = bool(deep) and all(min(d['shift_pct'].values()) > 5 for d in deep)
    res.append(('H23b(ii): ' + ('ningún cero profundo que mover: Muller no converge en la región (paso 20); la prueba no '
                                'decide' if not deep else 'los ceros profundos se mueven con los dos radios alternativos: ') + '; '.join(
        f"{d['nu'][0]:.3f}{d['nu'][1]:+.3f}i kHz -> {d[f'nu_{ka}'][0]:.3f}{d[f'nu_{ka}'][1]:+.3f}i ({d['shift_pct'][ka]:.2g} %) "
        f"y {d[f'nu_{kb}'][0]:.3f}{d[f'nu_{kb}'][1]:+.3f}i ({d['shift_pct'][kb]:.2g} %)" for d in deep),
        ok_ii, 'si (i) y (ii) pasan: artefactos del empalme, no modos'))
    # negative control 1: ingoing wave
    h1 = complex(*fl[0]['nu'])
    _S21['sheet'] = -1
    ratio_in = abs(_F21(h1)) / abs(_F21(h1 + 0.05))
    _S21['sheet'] = 1
    ratio_out = abs(_F21(h1)) / abs(_F21(h1 + 0.05))
    o['neg1'] = {'ratio_ingoing': float(ratio_in), 'ratio_outgoing': float(ratio_out)}
    res.append((f'Control negativo 1: con onda entrante, H1 deja de ser cero (|F(H1)|/|F(H1 + 50 Hz)| = {ratio_in:.2g}, '
                f'frente a {ratio_out:.1e} con la saliente)', ratio_in > 1e-2 and ratio_out < 1e-4, ''))
    # negative control 2: matching at the surface (no tail of the field)
    sR = fl[0]['shift_pct']['R']
    s_alt = max(fl[0]['shift_pct'][ka], fl[0]['shift_pct'][kb])
    o['neg2'] = {'shift_R_pct': sR, 'shift_alt_pct': s_alt}
    res.append((f'Control negativo 2: empalmando en la superficie (sin la cola del campo), H1 se mueve un {sR:.2g} %, '
                f'más que con los radios alternativos ({s_alt:.2g} %)', sR > 3 * s_alt, ''))
    o['checks'] = [{'name': n_, 'ok': ok_field(o_), 'note': t_} for n_, o_, t_ in res]
    return res, o


# =============================================================================================== step 23
def _modes23(R0, mode):
    """F (6 kHz), H1 and F (3.4 kHz) with the c_s^2 definition `mode`; the coefficient tables are rebuilt."""
    st6, _ = V.get_star(6.0)
    st34, _ = V.get_star(3.4)
    for s_ in (st6, st34):
        s_.eosv.cs2_mode = mode
    try:
        t6, t34 = RadialTables(st6, V._W['fi'], V._W['fe']), RadialTables(st34, V._W['fi'], V._W['fe'])
        rm6, rm34 = MD.matching_radius(st6), MD.matching_radius(st34)
        F6 = V.real_root(st6, t6, R0['F6'] - 0.01, R0['F6'] + 0.01, rm6)
        H1 = complex(MD.muller(lambda nu: MD.det(st6, t6, omega_of(nu), rm6), R0['H1']))
        F34 = complex(MD.muller(lambda nu: MD.det(st34, t34, omega_of(nu), rm34), R0['F34']))
    finally:
        for s_ in (st6, st34):
            s_.eosv.cs2_mode = 'eps'
    return {'F6': F6, 'H1': H1, 'F34': F34}


def step23(out):
    res, o = [], {}
    S = _need(out, 'radial_star.json', 15)
    R0 = {'F6': S['nu6']['F'], 'H1': complex(*S['nu6']['H1']), 'F34': complex(*S['nu34']['F'])}
    st6, _ = V.get_star(6.0)
    E = st6.eosv
    pc6 = float(st6.p0)
    o['mismatch'] = {str(p): E.mismatch(p) for p in (1, 10, 100, 300, 500, 900)}
    base = _modes23(R0, 'eps')
    o['pc6_MeVfm3'] = pc6
    o['mismatch_center'] = E.mismatch(pc6)
    alt, mid = _modes23(R0, 'n'), _modes23(R0, 'mid')
    cmp_ = lambda a, b: {'dRe_rel': float(abs((a.real - b.real) / b.real)), 'dIm_rel': float(abs((a.imag - b.imag) / b.imag)) if b.imag else 0.0,
                         'dRe_Hz': float(1e3 * (a.real - b.real)), 'dIm_Hz': float(1e3 * (a.imag - b.imag))}
    o['base'] = {k: [complex(v).real, complex(v).imag] for k, v in base.items()}
    o['alt'] = {k: [complex(v).real, complex(v).imag] for k, v in alt.items()}
    o['shift'] = {k: cmp_(complex(alt[k]), complex(base[k])) for k in ('F6', 'H1', 'F34')}
    o['shift_mid'] = {k: cmp_(complex(mid[k]), complex(base[k])) for k in ('F6', 'H1', 'F34')}
    sh = o['shift']
    res.append(('H27: inconsistencia de los ajustes BSk26, n deps/dn/(eps + p) - 1 = ' + ', '.join(f'{v:+.1e} (p = {k})' for k, v in o['mismatch'].items())
                + f' MeV/fm3; en el centro de la estrella de 6 kHz (p_c = {pc6:.0f}): {o["mismatch_center"]:+.1e}', None,
                'medido; no es un fallo del código sino de los ajustes publicados'))
    big = max(sh[k]['dRe_rel'] for k in sh)
    res.append(('H27: con c_s^2 = n dp/dn/(eps + p) en vez de dp/deps, Re nu se mueve ' + ', '.join(f"{k} {sh[k]['dRe_rel']:.1e} ({sh[k]['dRe_Hz']:+.2f} Hz)" for k in sh)
                + f"; Im nu de H1 {sh['H1']['dIm_rel']:.1e} y de F (3,4 kHz) {sh['F34']['dIm_rel']:.1e} relativo",
                1e-5 < big < 1e-2, 'mucho más que el 1e-6 de los controles del integrador: es la incertidumbre dominante del sector fluido'))
    # negative controls
    n_, e_, p_ = E.n, E.e, E.p
    lgrho, lgp = np.log10(e_ * MEV_FM3_TO_GCC), np.log10(p_)
    p_sh = 10**np.interp(1.01 * lgrho, lgrho, lgp)        # P(rho) shifted by 1 % in lg rho, from the table's own P(rho)
    ok = np.concatenate([[True], np.diff(p_sh) > 0])
    Bsh = _EOSVec(n_[ok], e_[ok], p_sh[ok])
    m_bad = abs(Bsh.mismatch(100.0))
    res.append((f'Control negativo 1: con P(rho) desplazado un 1 % en lg rho, el desajuste a 100 MeV/fm3 sube a {m_bad:.1e} '
                f'(frente a {abs(o["mismatch"]["100"]):.1e})', m_bad > 10 * abs(o['mismatch']['100']), ''))
    sm = o['shift_mid']
    inter = all(sm[k]['dRe_rel'] < sh[k]['dRe_rel'] for k in sh)
    res.append(('Control negativo 2: con la media geométrica de las dos definiciones, el desplazamiento es menor: ' + ', '.join(f"{k} {sm[k]['dRe_rel']:.1e}" for k in sm),
                inter, 'el efecto viene de c_s^2'))
    o['checks'] = [{'name': n2, 'ok': ok_field(o2), 'note': t2} for n2, o2, t2 in res]
    return res, o


# =============================================================================================== step 24
_S24 = {}


def _D24(arg):
    nu, sheet = arg
    return complex(MD.det_entire(_S24['st'], _S24['tabs'], omega_of(nu), inward=False, sheet=sheet, raw='inf'))


def _bw(x, c0, c1, nr, g):
    return c0 + c1 * (x - nr) + 2 * np.arctan((x - nr) / g)


def _bwd(deg):
    """arg S with a polynomial background of degree deg: p = (c0..c_deg, nu_r, Gamma/2)."""
    def f(x, *p):
        return np.polyval(np.asarray(p[:deg + 1])[::-1], x - p[deg + 1]) + 2 * np.arctan((x - p[deg + 1]) / p[deg + 2])
    return f


DEGS, WINS = (1, 2, 3, 4), (10.0, 3.0, 2.0, 1.5)        # degree of the background; half-window in widths


def scan(st, tabs, nu_c, half, n=61):
    """arg S and |S| in nu_c +- 10 half, denser at the centre. S = -W[u, psi_in] / W[u, psi_out]; half only sets the window."""
    _S24.update(st=st, tabs=tabs)
    nus = nu_c + half * np.tan(np.linspace(-1.47, 1.47, n))
    pool = V._pool()
    try:
        vals = V._map(pool, _D24, [(x, s) for x in nus for s in (1, -1)])
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    Dout, Din = np.array(vals[0::2]), np.array(vals[1::2])
    S = -Din / Dout
    return nus, np.unwrap(np.angle(S)), np.abs(S)


def fit(nus, ph, nu0, g0):
    """Linear background over the whole window (the naive fit; used in the negative control 2)."""
    x = (nus - nu0) * 1e3                       # Hz relative to the seed
    p0 = [float(np.mean(ph)), 0.0, 0.0, 1.3 * g0 * 1e3]     # width seed deliberately 30 % off
    p, _ = curve_fit(_bw, x, ph, p0=p0, maxfev=20000)
    res = float(np.max(np.abs(_bw(x, *p) - ph)))
    return nu0 + p[2] * 1e-3, abs(p[3]) * 1e-3, res


def fit_ladder(nus, ph, nu0, g0):
    """Fit for every background degree and half-window. Returns the table and the result (highest degree, narrowest window) with its error: the width difference with the previous degree there."""
    x = (nus - nu0) * 1e3
    tab = []
    for deg in DEGS:
        for W in WINS:
            m = np.abs(x) < W * g0 * 1e3
            if m.sum() < deg + 8:
                continue
            p0 = [float(np.mean(ph[m]))] + [0.0] * deg + [0.0, 1.3 * g0 * 1e3]
            p, _ = curve_fit(_bwd(deg), x[m], ph[m], p0=p0, maxfev=50000)
            tab.append({'deg': deg, 'win': W, 'n': int(m.sum()), 'nu_r': nu0 + p[-2] * 1e-3, 'half_width': abs(p[-1]) * 1e-3,
                        'resid_rad': float(np.max(np.abs(_bwd(deg)(x[m], *p) - ph[m])))})
    best = [t for t in tab if t['deg'] == DEGS[-1] and t['win'] == WINS[-1]][0]
    prev = [t for t in tab if t['deg'] == DEGS[-2] and t['win'] == WINS[-1]][0]
    return tab, best, abs(best['half_width'] / prev['half_width'] - 1)


def step24(out):
    res, o = [], {}
    S15 = _need(out, 'radial_star.json', 15)
    cases = {'H1 (6 kHz)': (6.0, complex(*S15['nu6']['H1'])), 'F (3,4 kHz)': (3.4, complex(*S15['nu34']['F']))}
    stars = {}
    for lab, (nu_a, z) in cases.items():
        st, tabs = stars[nu_a] = V.get_star(nu_a)
        nus, ph, aS = scan(st, tabs, z.real, -z.imag)
        tab, best, conv = fit_ladder(nus, ph, z.real, -z.imag)
        nr, g, rf = best['nu_r'], best['half_width'], best['resid_rad']
        lin = [t for t in tab if t['deg'] == 1 and t['win'] == WINS[0]][0]
        o[lab] = {'muller': [z.real, z.imag], 'nu_r': nr, 'half_width': g, 'fit_resid_rad': rf, 'ladder': tab, 'conv': float(conv),
                  'lin_rel_width': float(abs(lin['half_width'] / -z.imag - 1)), 'lin_resid_rad': lin['resid_rad'], 'nu_a': nu_a, 'thr_widths': float((z.real - nu_a) / -z.imag),
                  'unitarity': float(np.max(np.abs(aS - 1))), 'rel_width': float(abs(g / -z.imag - 1)), 'dnu_over_width': float(abs(nr - z.real) / -z.imag),
                  'phase_turn': float(ph[-1] - ph[0])}
        # negative control 1: window shifted by 20 Gamma
        nus_o, ph_o, _ = scan(st, tabs, z.real + 40 * -z.imag, -z.imag, n=21)
        lin_ = np.polyval(np.polyfit(nus_o, ph_o, 1), nus_o)
        o[lab]['off_turn'] = float(np.max(np.abs(ph_o - lin_)))
    # negative control 2: the old lapse normalisation (as step 20) on H1
    st, _ = stars[6.0]
    st_o, tabs_o = V.old_lapse(st)
    z = cases['H1 (6 kHz)'][1]
    nus, ph, _ = scan(st_o, tabs_o, z.real, -z.imag * 6)
    try:
        nr_o, g_o, _ = fit(nus, ph, z.real, -z.imag)
    except RuntimeError:
        nr_o, g_o = np.nan, np.nan
    o['old_lapse_H1'] = {'nu_r': float(nr_o), 'half_width': float(g_o), 'dnu_over_width': float(abs(nr_o - z.real) / -z.imag)}
    for lab in cases:
        o_ = o[lab]
        res.append((f"H28(a) {lab}: la dispersión a frecuencia real da Gamma/2 = {o_['half_width'] * 1e3:.4f} Hz frente a "
                    f"-Im nu = {-o_['muller'][1] * 1e3:.4f} Hz de Muller ({100 * o_['rel_width']:.2g} %), y nu_r a "
                    f"{o_['dnu_over_width']:.1e} anchuras de Re nu (fondo de grado {DEGS[-1]} en +-{WINS[-1]:g} anchuras; residuo "
                    f"{o_['fit_resid_rad']:.0e} rad; grado {DEGS[-2]} cambia la anchura {o_['conv']:.0e}). Con fondo lineal en "
                    f"+-{WINS[0]:g} anchuras: {100 * o_['lin_rel_width']:.2g} % con residuo {o_['lin_resid_rad']:.0e} rad (umbral a "
                    f"{o_['thr_widths']:.0f} anchuras)",
                    o_['rel_width'] < 0.01 and o_['dnu_over_width'] < 0.1 and o_['conv'] < 1e-3,
                    'sin plano complejo ni búsqueda de raíces; cerca del umbral el fondo de fase se curva'))
    un = max(o[lab]['unitarity'] for lab in cases)
    res.append((f'H28(b): |S| = 1 a {un:.0e} en las dos ventanas', un < 1e-8, 'lo garantizan unas ecuaciones reales: prueba la numérica del exterior, no la física'))
    off = max(o[lab]['off_turn'] for lab in cases)
    res.append((f'Control negativo 1: desplazando la ventana 20 Gamma, arg S se aparta de una recta como mucho {off:.2g} rad (frente a 2 pi en la resonancia)', off < 0.5, ''))
    ol = o['old_lapse_H1']
    res.append((f"Control negativo 2: con la normalización antigua del lapso, la resonancia de H1 cae a {ol['dnu_over_width']:.2g} anchuras del modo corregido",
                not np.isfinite(ol['dnu_over_width']) or ol['dnu_over_width'] > 1, ''))
    o['checks'] = [{'name': n_, 'ok': ok_field(o_), 'note': t_} for n_, o_, t_ in res]
    return res, o
