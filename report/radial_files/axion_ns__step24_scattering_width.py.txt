"""Paso 24: amortiguamiento desde la dispersión a frecuencia real (hipótesis H28).

Sin plano complejo ni búsqueda de raíces: a nu real se lanza una onda de axión entrante y se mide
    S(nu) = -W[u, psi_in] / W[u, psi_out],
con u la solución regular del interior (Delta p(R) = 0) integrada hasta r_m y psi_in/out normalizadas a
e^{-+i(kr + eta ln r)}/r (radial.det_entire con raw='inf', sheet = -1 / +1). Cerca de un modo, arg S gira 2 pi:
    arg S = c0 + c1 (nu - nu_r) + 2 arctan((nu - nu_r) / (Gamma/2)),
y el ajuste da nu_r y Gamma/2. Se compara con Muller (pasos 15 y 20): H1 a 6 kHz y F a 3,4 kHz.
El fondo no es lineal cerca del umbral: F a 3,4 kHz está a solo ~16 anchuras de nu_a, donde la fase de fondo
(que va como sqrt(nu - nu_a)) se curva. En el primer recálculo (41 puntos, fondo lineal en +-10 anchuras) ese ajuste
no cerró: residuo de 3 rad y anchura un 26 % fuera (H28). Con fondo lineal el residuo sigue siendo de radianes. Por eso el fondo es un polinomio de grado d y el ajuste se repite en ventanas cada vez más estrechas:
el resultado es el de mayor grado en la ventana más estrecha, y la diferencia con el grado anterior mide su error.
|S| = 1 lo garantizan unas ecuaciones reales: comprobarlo prueba la numérica del exterior, no la física.
Controles negativos: (1) una ventana desplazada 20 Gamma no muestra resonancia; (2) con la normalización antigua del
lapso, la resonancia de dispersión deja de coincidir con el modo corregido.
"""
import json
import multiprocessing as mp
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

from . import axstar as AX
from . import radial as RD
from .pipeline import ok_field, report

OUT = Path(__file__).resolve().parent.parent / 'results'
_S = {}


def _D(arg):
    nu, sheet = arg
    return complex(RD.det_entire(_S['st'], AX.omega_of(nu), inward=False, sheet=sheet, raw='inf'))


def _bw(x, c0, c1, nr, g):
    return c0 + c1 * (x - nr) + 2 * np.arctan((x - nr) / g)


def _bwd(deg):
    """arg S con fondo polinómico de grado deg: p = (c0..c_deg, nu_r, Gamma/2)."""
    def f(x, *p):
        return np.polyval(np.asarray(p[:deg + 1])[::-1], x - p[deg + 1]) + 2 * np.arctan((x - p[deg + 1]) / p[deg + 2])
    return f


DEGS, WINS = (1, 2, 3, 4), (10.0, 3.0, 2.0, 1.5)        # grado del fondo; semiventana en anchuras


def scan(st, nu_c, half, n=61):
    """arg S y |S| en nu_c +- 10 half, más denso en el centro. half solo fija la ventana."""
    _S['st'] = st                                         # antes de crear el pool: los procesos hijos lo heredan
    nus = nu_c + half * np.tan(np.linspace(-1.47, 1.47, n))
    with mp.get_context('fork').Pool(6) as pool:
        vals = pool.map(_D, [(x, s) for x in nus for s in (1, -1)])
    Dout, Din = np.array(vals[0::2]), np.array(vals[1::2])
    S = -Din / Dout
    return nus, np.unwrap(np.angle(S)), np.abs(S)


def fit(nus, ph, nu0, g0):
    """Fondo lineal en toda la ventana (el ajuste ingenuo; se usa en el control negativo 2)."""
    x = (nus - nu0) * 1e3                       # Hz respecto a la semilla
    p0 = [float(np.mean(ph)), 0.0, 0.0, 1.3 * g0 * 1e3]     # semilla de anchura deliberadamente desviada un 30 %
    p, _ = curve_fit(_bw, x, ph, p0=p0, maxfev=20000)
    res = float(np.max(np.abs(_bw(x, *p) - ph)))
    return nu0 + p[2] * 1e-3, abs(p[3]) * 1e-3, res


def fit_ladder(nus, ph, nu0, g0):
    """Ajuste para cada grado de fondo y semiventana. Devuelve la tabla y el resultado (mayor grado, ventana más
    estrecha) con su error: la diferencia de anchura con el grado anterior en esa ventana."""
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


def run():
    res, out = [], {}
    S15 = json.loads((OUT / 'radial_star.json').read_text())
    cases = {'H1 (6 kHz)': (6.0, complex(*S15['nu6']['H1'])), 'F (3,4 kHz)': (3.4, complex(*S15['nu34']['F']))}
    stars = {}
    if True:
        for lab, (nu_a, z) in cases.items():
            st = stars[nu_a] = AX.Star(nu_a)
            nus, ph, aS = scan(st, z.real, -z.imag)
            tab, best, conv = fit_ladder(nus, ph, z.real, -z.imag)
            nr, g, rf = best['nu_r'], best['half_width'], best['resid_rad']
            lin = [t for t in tab if t['deg'] == 1 and t['win'] == WINS[0]][0]
            out[lab] = {'muller': [z.real, z.imag], 'nu_r': nr, 'half_width': g, 'fit_resid_rad': rf, 'ladder': tab,
                        'conv': float(conv), 'lin_rel_width': float(abs(lin['half_width'] / -z.imag - 1)),
                        'lin_resid_rad': lin['resid_rad'],
                        'nu_a': nu_a, 'thr_widths': float((z.real - nu_a) / -z.imag),
                        'unitarity': float(np.max(np.abs(aS - 1))), 'rel_width': float(abs(g / -z.imag - 1)),
                        'dnu_over_width': float(abs(nr - z.real) / -z.imag), 'phase_turn': float(ph[-1] - ph[0])}
            # control negativo 1: ventana desplazada 20 Gamma
            nus_o, ph_o, _ = scan(st, z.real + 40 * -z.imag, -z.imag, n=21)
            lin = np.polyval(np.polyfit(nus_o, ph_o, 1), nus_o)
            out[lab]['off_turn'] = float(np.max(np.abs(ph_o - lin)))
        # control negativo 2: normalización antigua del lapso (como el paso 20) sobre H1
        st = stars[6.0]
        m_R, Phi_R = st.sol.sol(st.R)[0], st.sol.sol(st.R)[2]
        Ps_new = st.Phi_shift
        st._Ps = 0.5 * np.log(1 - 2 * m_R / st.R) - Phi_R
        z = cases['H1 (6 kHz)'][1]
        nus, ph, _ = scan(st, z.real, -z.imag * 6)
        st._Ps = Ps_new
        try:
            nr_o, g_o, _ = fit(nus, ph, z.real, -z.imag)
        except RuntimeError:
            nr_o, g_o = np.nan, np.nan
        out['old_lapse_H1'] = {'nu_r': float(nr_o), 'half_width': float(g_o),
                               'dnu_over_width': float(abs(nr_o - z.real) / -z.imag)}
    for lab in cases:
        o = out[lab]
        res.append((f"H28(a) {lab}: la dispersión a frecuencia real da Gamma/2 = {o['half_width'] * 1e3:.4f} Hz frente a "
                    f"-Im nu = {-o['muller'][1] * 1e3:.4f} Hz de Muller ({100 * o['rel_width']:.2g} %), y nu_r a "
                    f"{o['dnu_over_width']:.1e} anchuras de Re nu (fondo de grado {DEGS[-1]} en +-{WINS[-1]:g} anchuras; residuo "
                    f"{o['fit_resid_rad']:.0e} rad; grado {DEGS[-2]} cambia la anchura {o['conv']:.0e}). Con fondo lineal en "
                    f"+-{WINS[0]:g} anchuras: {100 * o['lin_rel_width']:.2g} % con residuo {o['lin_resid_rad']:.0e} rad (umbral a "
                    f"{o['thr_widths']:.0f} anchuras)",
                    o['rel_width'] < 0.01 and o['dnu_over_width'] < 0.1 and o['conv'] < 1e-3,
                    'sin plano complejo ni búsqueda de raíces; cerca del umbral el fondo de fase se curva'))
    un = max(out[lab]['unitarity'] for lab in cases)
    res.append((f'H28(b): |S| = 1 a {un:.0e} en las dos ventanas', un < 1e-8,
                'lo garantizan unas ecuaciones reales: prueba la numérica del exterior, no la física'))
    off = max(out[lab]['off_turn'] for lab in cases)
    res.append((f'Control negativo 1: desplazando la ventana 20 Gamma, arg S se aparta de una recta como mucho {off:.2g} rad '
                f'(frente a 2 pi en la resonancia)', off < 0.5, ''))
    ol = out['old_lapse_H1']
    res.append((f"Control negativo 2: con la normalización antigua del lapso, la resonancia de H1 cae a {ol['dnu_over_width']:.2g} "
                f"anchuras del modo corregido", not np.isfinite(ol['dnu_over_width']) or ol['dnu_over_width'] > 1, ''))
    out['checks'] = [{'name': n_, 'ok': ok_field(o_), 'note': t_} for n_, o_, t_ in res]
    OUT.mkdir(exist_ok=True)
    (OUT / 'radial_scattering.json').write_text(json.dumps(out, indent=1))
    return report('Paso 24: anchura de los modos desde la dispersión a frecuencia real', res)


if __name__ == '__main__':
    run()
