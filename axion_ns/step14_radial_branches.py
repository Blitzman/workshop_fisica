"""Paso 14: modos radiales del modelo simplificado (Paper I, Ec. 4.13). Transición cuasi-normal -> evanescente
y ramas del espectro (pregunta de José, 2026-11-07; hipótesis H8-H10 en LAB_NOTEBOOK.md).

Idea clave: en el modelo simplificado la función de dispersión solo depende de w^2 y del número de onda exterior
k = k_out (w = omega/m_a, unidades de m_a). Si se usa k como variable (w^2 = 1 + k^2), la función es univaluada y no
hay que elegir hoja de Riemann:
    Im k > 0  -> estado ligado (onda exterior evanescente, frecuencia real si Im k > 0 y Re k = 0)
    Im k < 0, Re k != 0 -> modo cuasi-normal (onda saliente, amortiguado)
    Re k = 0, Im k < 0  -> estado virtual (frecuencia real pero la onda crece en el infinito; no es un modo físico)

Limitación: f es la Ec. (4.13) tal como se usa en el paso 4 (`step4_numerics.f_toy`), no una compilación desde SymPy.
"""
import json
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.optimize import brentq
from scipy.optimize import minimize_scalar

from .pipeline import ok_field, report
from .step4_numerics import f_toy

OUT = Path(__file__).resolve().parent.parent / 'results'
mp.mp.dps = 30
ALPHA, CS2 = 10.0, 0.3                     # parámetros de la Fig. 3 del Paper I (lambda^2 = 1e-3, x0 = 0.6)
BETA = 4 * 0.48 / 1.48**2


def f_k(k, lam2, alpha, cs2, x0, mass2_in=None, sign_out=1):
    """Ec. (4.13) en la variable k = k_out. sign_out = -1 es la mutación «onda entrante» (control negativo)."""
    if mass2_in is None:
        mass2_in = 1 - alpha
    w2 = 1 + k**2
    ka2, ks2 = w2 - mass2_in, w2 / cs2
    A = alpha**2 * lam2
    sD = mp.sqrt((ka2 - ks2 + A)**2 + 4 * A * ks2)
    kp2, km2 = (ka2 + ks2 + A + sD) / 2, (ka2 + ks2 + A - sD) / 2
    kp, km = mp.sqrt(kp2), mp.sqrt(km2)
    cp, cm = mp.cos(kp * x0), mp.cos(km * x0)
    return ((ka2 - km2) * kp * mp.sin(kp * x0) * cm - (ka2 - kp2) * km * mp.sin(km * x0) * cp
            + sign_out * 1j * k * sD * cp * cm)


def w_of(k):
    w = complex(mp.sqrt(1 + mp.mpc(k)**2))
    return w if w.real >= 0 else -w


def root(k0, pars, **kw):
    k0 = mp.mpc(k0)
    return complex(mp.findroot(lambda k: f_k(k, *pars, **kw), (k0, k0 + 1e-3, k0 + 1e-3j), solver='muller', tol=1e-25))


def axis_roots(pars, kmax=1.0, n=1500, **kw):
    """Raíces sobre el eje imaginario k = i kappa, donde f es real (cambios de signo + refinado).
    Se descartan cambios de signo que no son ceros (|f| grande): ocurren si f deja de ser real."""
    g = lambda x: mp.re(f_k(mp.mpc(0, x), *pars, **kw))
    xs = np.linspace(-kmax, kmax, n)
    v = [float(g(x)) for x in xs]
    out = []
    for i in range(n - 1):
        if v[i] * v[i + 1] < 0:
            x = float(mp.findroot(g, (xs[i], xs[i + 1]), solver='anderson'))
            if abs(f_k(mp.mpc(0, x), *pars, **kw)) < 1e-8 * (1 + abs(v[i]) + abs(v[i + 1])):
                out.append(x)
    return out


def D_k(k, lam2, alpha, cs2, x0, mass2_in=None):
    """Discriminante D de la Ec. (4.13): f es impar en sqrt(D), así que f/sqrt(D) es la función entera de k
    y f tiene ceros espurios donde D = 0 (hay que comprobar que ninguna raíz cae ahí)."""
    if mass2_in is None:
        mass2_in = 1 - alpha
    w2 = 1 + k**2
    ka2, ks2 = w2 - mass2_in, w2 / cs2
    A = alpha**2 * lam2
    return (ka2 - ks2 + A)**2 + 4 * A * ks2


def g_k(k, *pars, **kw):
    return f_k(k, *pars, **kw) / mp.sqrt(D_k(k, *pars, **{x: v for x, v in kw.items() if x == 'mass2_in'}))


def winding(pars, re, im, n=600):
    """Número de ceros de f/sqrt(D) dentro del rectángulo (principio del argumento), sin buscador de raíces."""
    (a, b), (c, d) = re, im
    path = ([complex(a + (b - a) * t, c) for t in np.linspace(0, 1, n, endpoint=False)]
            + [complex(b, c + (d - c) * t) for t in np.linspace(0, 1, n, endpoint=False)]
            + [complex(b - (b - a) * t, d) for t in np.linspace(0, 1, n, endpoint=False)]
            + [complex(a, d - (d - c) * t) for t in np.linspace(0, 1, n, endpoint=False)])
    vals = [complex(g_k(mp.mpc(z), *pars)) for z in path] 
    vals.append(vals[0])
    dph = np.angle(np.array(vals[1:]) / np.array(vals[:-1]))
    return float(np.sum(dph) / (2 * np.pi)), float(np.max(np.abs(dph)))


# ------------------------------------------------------------------ H8 + H9: paso por el umbral
def threshold_track(lam2=1e-3, x_lo=0.80, x_hi=0.8525, dx=0.0025, **kw):
    k = root(0.37599 - 0.00786j, (lam2, ALPHA, CS2, x_lo), **kw)
    rows = []
    for x0 in np.arange(x_lo, x_hi + 1e-9, dx):
        k = root(k, (lam2, ALPHA, CS2, x0), **kw)
        rows.append((float(x0), k))
        if abs(k.real) < 1e-12:
            break
    return rows


def collision(lam2, k0, x00):
    """Colisión exacta del par cuasi-normal: raíz doble de g = f/sqrt(D) en (k, x0), con Newton complejo
    (g = 0 y dg/dk = 0). Devuelve (x0, k, residuo)."""
    P = lambda k, x: complex(g_k(k, lam2, ALPHA, CS2, x))

    def F(v):
        k, x = v
        h = 1e-7
        return np.array([P(k, x), (P(k + h, x) - P(k - h, x)) / (2 * h)], dtype=complex)
    v = np.array([k0, x00], dtype=complex)
    for _ in range(80):
        Fv, J = F(v), np.zeros((2, 2), complex)
        for j in range(2):
            dv = np.zeros(2, complex)
            dv[j] = 1e-8
            J[:, j] = (F(v + dv) - F(v - dv)) / 2e-8
        step = np.linalg.solve(J, -Fv)
        v = v + step
        if abs(step).max() < 1e-13:
            break
    return float(v[1].real), complex(v[0]), float(abs(F(v)).max())


def k0_crossing(lam2, x_from, width=0.05):
    """x0 en que un polo del eje imaginario pasa por k = 0 (w = 1): el virtual se convierte en ligado."""
    g0 = lambda xx: complex(g_k(1e-13j, lam2, ALPHA, CS2, xx)).real
    xs = np.linspace(x_from, x_from + width, 400)
    v = [g0(t) for t in xs]
    return next(brentq(g0, xs[i], xs[i + 1], xtol=1e-13) for i in range(len(xs) - 1) if v[i] * v[i + 1] < 0)


def damping_exponent(lam2=1e-4):
    k = root(complex(0.376, -0.0008), (lam2, ALPHA, CS2, 0.80))
    pts = []
    for x0 in np.arange(0.80, 0.86, 0.0002):
        k = root(k, (lam2, ALPHA, CS2, x0))
        if abs(k.real) < 1e-9:
            break
        w = w_of(k)
        pts.append((w.real - 1, -w.imag, k.real, -k.imag))
    p = np.array(pts)
    m = (p[:, 2] > 5 * p[:, 3]) & (p[:, 0] < 0.02)
    return float(np.polyfit(np.log(p[m, 0]), np.log(p[m, 1]), 1)[0]), p


# ------------------------------------------------------------------ H10: cruce evitado fluido-axión
def bound_w2(lam2, x0, **kw):
    return sorted(1 - x * x for x in axis_roots((lam2, ALPHA, CS2, x0), kmax=1.0, **kw) if x > 0)


def min_gap(lam2):
    def gap(x0):
        b = bound_w2(lam2, x0)
        return b[-1] - b[-2] if len(b) >= 2 else 9.0
    r = minimize_scalar(gap, bounds=(1.04, 1.10), method='bounded', options={'xatol': 2e-5})
    return float(r.fun), float(r.x)


# ------------------------------------------------------------------ inestabilidad taquiónica en theta = 0
def x0_unstable(alpha=ALPHA, mass2_in=None, lam2=0.0):
    """x0 mínimo con un modo w^2 = 0 (kappa = 1): a partir de ahí hay un modo que crece."""
    g = lambda x0: mp.re(f_k(mp.mpc(0, 1), lam2, alpha, CS2, x0, mass2_in=mass2_in))
    xs = np.linspace(0.02, 1.5, 600)
    v = [float(g(x)) for x in xs]
    for i in range(len(xs) - 1):
        if v[i] * v[i + 1] < 0:
            return float(mp.findroot(g, (xs[i], xs[i + 1]), solver='anderson'))
    return None


def run():
    res = []
    # 1. Coherencia con el paso 4: f_k(k) = f_toy(w) en la hoja física por encima del umbral
    w0 = mp.mpc(1.43, -0.003)
    k0 = mp.sqrt(w0**2 - 1)
    d = abs(f_k(k0, 1e-3, ALPHA, CS2, 0.6) - f_toy(w0, 1e-3, ALPHA, CS2, 0.6)) / abs(f_toy(w0, 1e-3, ALPHA, CS2, 0.6))
    res.append((f'f(k) en la variable k coincide con f_toy(w) del paso 4 (Re w > 1): diferencia relativa {float(d):.1e}',
                d < 1e-20, ''))
    kf = root(np.sqrt(1.43**2 - 1), (1e-3, ALPHA, CS2, 0.6))
    wf = w_of(kf)
    res.append((f'Modo fluido fundamental en x0 = 0.6: w = {wf.real:.5f}{wf.imag:+.2e}i (cuasi-normal, Im k = {kf.imag:+.4f} < 0)',
                wf.imag < 0 and kf.imag < 0 and kf.real > 0, ''))

    # 2. H9: secuencia resonancia -> colisión en el eje -> ligado + virtual
    tr = threshold_track()
    pair_ok = all(k.real > 0 and k.imag < 0 for _, k in tr[:-1])
    ax86 = axis_roots((1e-3, ALPHA, CS2, 0.86), kmax=0.4, n=800)
    bound = [x for x in ax86 if x > 0]
    virt = [x for x in ax86 if x < 0]
    seq_ok = pair_ok and len(bound) == 1 and len(virt) == 1
    wb = np.sqrt(1 - bound[0]**2) if bound else np.nan
    coll = {}
    for l2, k0, x00 in ((1e-3, -0.004j, tr[-1][0]), (1e-4, -0.0013j, 0.8594)):
        xc_, kc_, rc_ = collision(l2, k0, x00)
        coll[l2] = {'x0': xc_, 'k': [kc_.real, kc_.imag], 'w_minus_1': float((w_of(kc_) - 1).real),
                    'im_w': float(w_of(kc_).imag), 'resid': rc_, 'x0_k0': float(k0_crossing(l2, xc_))}
    c3 = coll[1e-3]
    res.append((f'H9: el par cuasi-normal k = ±a - ib (b ≈ {-tr[0][1].imag:.4f}...{-tr[-2][1].imag:.4f}) colisiona en el eje '
                f'imaginario negativo en x0 = {c3["x0"]:.7f} (raíz doble, residuo {c3["resid"]:.0e}; λ² = 1e-3), con '
                f'k = {c3["k"][1]:.7f}i y w - 1 = {c3["w_minus_1"]:.3e}: ya por debajo del umbral. Quedan dos virtuales; uno '
                f'cruza k = 0 (w = 1) en x0 = {c3["x0_k0"]:.7f} y pasa a ligado. En x0 = 0.86: un ligado (kappa = '
                f'+{bound[0] if bound else 0:.4f}, w = {wb:.5f}) y un virtual (kappa = {virt[0] if virt else 0:+.4f})',
                seq_ok and c3['w_minus_1'] < 0 and c3['x0_k0'] > c3['x0'] and c3['resid'] < 1e-10, ''))

    # 3. H8: |Im w| ∝ (Re w - 1)^s cerca del umbral
    s, damp_pts = damping_exponent()
    fit = (damp_pts[:, 2] > 5 * damp_pts[:, 3]) & (damp_pts[:, 0] < 0.02)
    imk = damp_pts[fit, 3]
    c4 = coll[1e-4]
    res.append((f'H8: |Im w| ~ (Re w - 1)^{s:.3f} en el tramo intermedio Re k > 5|Im k| (λ² = 1e-4; ½ si Im k '
                f'fuera constante, que varía un {100 * (imk.max() / imk.min() - 1):.0f} % en el tramo)', 0.40 < s < 0.52,
                f'ajuste en un tramo finito, no una ley asintótica: más cerca del umbral domina la colisión de polos. '
                f'El amortiguamiento solo se anula en la colisión, x0 = {c4["x0"]:.10f}, con w - 1 = {c4["w_minus_1"]:.4e} '
                f'(por debajo del umbral, = k²/2 a este orden, k = {c4["k"][1]:.7f}i); el virtual cruza k = 0 en '
                f'x0 = {c4["x0_k0"]:.10f}'))

    # 4. H10: separación mínima del cruce evitado lineal en lambda
    gaps = [(l2, *min_gap(l2)) for l2 in (1e-6, 1e-5, 1e-4, 1e-3)]
    ratios = [g / np.sqrt(l2) for l2, g, _ in gaps]
    lin_ok = max(ratios) / min(ratios) < 1.02
    res.append(('H10: cruce evitado fluido-axión, Δw²_min/λ = ' + ', '.join(f'{r:.3f}' for r in ratios)
                + f' para λ² = 1e-6…1e-3 (en x0 ≈ {gaps[0][2]:.3f}): separación ∝ λ', lin_ok, ''))

    # 5. Inestabilidad del modelo linealizado en theta = 0 (fuera de la ventana Im w <= 0.01 del paso 4)
    xc = x0_unstable()
    xc3 = x0_unstable(lam2=1e-3)
    xc_an = np.arctan(1 / np.sqrt(ALPHA - 1)) / np.sqrt(ALPHA - 1)
    un = [x for x in axis_roots((1e-3, ALPHA, CS2, 0.6), kmax=4.0, n=2000) if x > 1]
    res.append((f'Modo inestable del modelo simplificado (m²_in = m_a²(1 - α)): sin acoplo aparece para x0 > {xc:.8f} '
                f'(analítico arctan(1/√(α-1))/√(α-1) = {xc_an:.8f}); con λ² = 1e-3, para x0 > {xc3:.8f}; '
                f'en los parámetros de la figura (x0 = 0.6) w² = {1 - un[0]**2:.3f} (crece como e^({np.sqrt(un[0]**2 - 1):.2f} m_a t))',
                abs(xc - xc_an) < 1e-6 and len(un) == 1 and xc3 < xc,
                'el paso 4 no lo veía: buscaba solo Im w ≤ 0.01; el acoplo apenas mueve el umbral'))
    m2pi = (ALPHA - 1) / np.sqrt(1 - BETA)
    st_pi = [x for x in axis_roots((0.0, ALPHA, CS2, 0.6), kmax=4.0, n=2000, mass2_in=m2pi) if x > 0]
    res.append((f'Variante theta = pi (m²_in = {m2pi:.1f} m_a², sin acoplo): sin estados ligados ni inestables', not st_pi, ''))

    # 5b. Revisión externa (2026-11-08): identidad Im w con el signo correcto, ceros espurios D = 0, método independiente
    id_plus = abs(wf.imag - kf.real * kf.imag / wf.real)
    id_minus = abs(wf.imag + kf.real * kf.imag / wf.real)
    res.append((f'Identidad Im w = +Re k·Im k/Re w (de w² = 1 + k²): residuo {id_plus:.1e}; con el signo − que escribimos '
                f'por error en el informe: {id_minus:.1e}', id_plus < 1e-12 and id_minus > 1e-4, 'errata señalada por la revisión'))
    roots_all = [(kf, 0.6), *[(k, x0) for x0, k in tr], (1j * bound[0], 0.86), (1j * virt[0], 0.86), (1j * un[0], 0.6)]
    dmin = min(abs(complex(D_k(mp.mpc(k), 1e-3, ALPHA, CS2, x0))) for k, x0 in roots_all)
    res.append((f'Ninguna raíz reportada es un cero espurio de D (min |D| = {dmin:.1f})', dmin > 1, ''))
    boxes = [('par cuasi-normal, x0 = 0.80, Re k > 0', (1e-3, ALPHA, CS2, 0.80), (0.05, 0.6), (-0.05, 0.02), 1),
             ('par completo, x0 = 0.80', (1e-3, ALPHA, CS2, 0.80), (-0.6, 0.6), (-0.05, -0.002), 2),
             ('ligado, x0 = 0.86', (1e-3, ALPHA, CS2, 0.86), (-0.05, 0.05), (0.05, 0.25), 1),
             ('virtual, x0 = 0.86', (1e-3, ALPHA, CS2, 0.86), (-0.05, 0.05), (-0.25, -0.05), 1),
             ('inestable, x0 = 0.6', (1e-3, ALPHA, CS2, 0.6), (-0.1, 0.1), (2.4, 2.9), 1)]
    wres = [(lab, *winding(pp, re_, im_), exp) for lab, pp, re_, im_, exp in boxes]
    ok_w = all(abs(nw - exp) < 0.05 and mx < 1.0 for _, nw, mx, exp in wres)
    res.append(('Método independiente (principio del argumento sobre f/√D, sin Muller): ' +
                '; '.join(f'{lab}: {nw:.2f} (esperado {exp})' for lab, nw, _, exp in wres), ok_w, 'H14'))

    # 6. Controles negativos
    ax_mut = axis_roots((1e-3, ALPHA, CS2, 0.86), kmax=0.4, n=800, sign_out=-1)
    neg1 = not (len([x for x in ax_mut if x > 0]) == 1 and abs([x for x in ax_mut if x > 0][0] - bound[0]) < 1e-6)
    res.append(('Control negativo 1 (onda entrante, i k -> -i k): el estado ligado de x0 = 0.86 deja de estar en kappa > 0 '
                'en la misma posición', neg1, ''))
    xc_mut = x0_unstable(mass2_in=ALPHA - 1)
    res.append((f'Control negativo 2 (m²_in = 1 - α -> α - 1): desaparece el umbral de inestabilidad (x0_c = {xc_mut})',
                xc_mut is None, ''))

    # Datos para las figuras del informe (report/build_radial.py): nada se teclea a mano
    after = [(float(x0), axis_roots((1e-3, ALPHA, CS2, x0), kmax=0.5, n=600))
             for x0 in np.arange(0.8525, 0.9001, 0.0025)]
    branches = {}
    for l2 in (0.0, 1e-3):
        branches[str(l2)] = [(round(float(x0), 3), [1 - x * x for x in axis_roots((l2, ALPHA, CS2, x0), kmax=3.2, n=1600) if x > 0])
                             for x0 in np.arange(0.80, 1.3001, 0.01)]
    OUT.mkdir(exist_ok=True)
    (OUT / 'radial_branches.json').write_text(json.dumps({
        'params': {'alpha': ALPHA, 'cs2': CS2, 'lam2_track': 1e-3, 'lam2_damping': 1e-4},
        'checks': [{'name': n, 'ok': ok_field(o), 'note': nt} for n, o, nt in res],
        'fluid_x0_0p6': [wf.real, wf.imag, kf.real, kf.imag],
        'threshold_track': [[x0, k.real, k.imag, w_of(k).real, w_of(k).imag] for x0, k in tr],
        'axis_after': [{'x0': x0, 'kappa': ks} for x0, ks in after],
        'branches': branches,
        'damping_points': [[float(a_), float(b_)] for a_, b_, c_, d_ in damp_pts if c_ > 5 * d_ and a_ < 0.02],
        'unstable_x0_0p6_w2': float(1 - un[0]**2),
        'axis_x0_0p86': {'bound': bound, 'virtual': virt},
        'damping_exponent': s,
        'collision': {str(k_): v_ for k_, v_ in coll.items()},
        'gaps': [{'lam2': l2, 'gap_w2': g, 'x0': x} for l2, g, x in gaps],
        'x0_unstable': xc, 'x0_unstable_analytic': float(xc_an), 'x0_unstable_lam2_1e-3': xc3,
        'winding': [{'box': lab, 'count': nw, 'expected': exp} for lab, nw, _, exp in wres],
    }, indent=1))
    return report('Paso 14: modos radiales, umbral cuasi-normal/evanescente y ramas (modelo simplificado)', res)


if __name__ == '__main__':
    run()
