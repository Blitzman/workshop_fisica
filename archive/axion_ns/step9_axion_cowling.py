"""Paso 9: notas de la estrella con axión en Cowling, frente a las Tablas I, II y IV del Paper II.

Etapa A antes de integrar el sistema extendido con axión: comprueba el fondo realista (BSk26 +
condensado, axstar.py) y la parte axiónica de las ecuaciones con números publicados.

- Sin axión: M, R (Paper I: 1.53 M_sun, 11.7 km) y f, p1, p2 (Paper II: 2.43, 7.15, 10.18 kHz).
- Tabla IV: masa de las estrellas con axión de 8 a 15 kHz.
- Tabla I: notas con el fondo axiónico y delta theta apagado.
- Tabla II: notas con delta theta, parte real y amortiguamiento por emisión de axiones.
Como nuestro BSk26 es el ajuste analítico de Potekhin y no la tabla del grupo, las notas sin axión
difieren hasta un 0.3 %; las Tablas I y II se comparan como cociente nu/nu_0 con la referencia propia.
"""
import multiprocessing as mp
import time

import numpy as np
from scipy.optimize import brentq

from . import axstar as A

REF = {'M': 1.53, 'R': 11.7, 'modes': (2.43, 7.15, 10.18)}
# Paper II, Tabla I (solo fondo): nu_a -> (f, p1, p2) en kHz
T1 = {2.8: (2.44, 7.19, 10.22), 3.0: (2.49, 7.35, 10.42), 4.0: (2.64, 7.80, 11.17), 5.0: (2.70, 7.95, 11.49),
      6.0: (2.73, 8.01, 11.64), 8.0: (2.78, 8.08, 11.82), 12.0: (2.84, 8.15, 12.01), 15.0: (2.87, 8.19, 12.11)}
# Paper II, Tabla II (acoplado): nu_a -> ((Re kHz, Im Hz) de f, p1, p2)
T2 = {2.8: ((2.44, 0), (7.18, -0.28), (10.22, -0.74)), 3.0: ((2.47, 0), (7.32, -1.4), (10.40, -2.7)),
      4.0: ((2.55, 0), (7.72, -1.1), (11.14, -0.75)), 5.0: ((2.56, 0), (7.88, -0.074), (11.44, -12)),
      6.0: ((2.55, 0), (7.95, -1.4), (11.58, -13)), 7.0: ((2.54, 0), (7.98, -1.6), (11.68, -7.1)),
      8.0: ((2.54, 0), (7.99, 0), (11.75, -1.7)), 10.0: ((2.52, 0), (7.99, 0), (11.86, -1.5)),
      11.0: ((2.51, 0), (7.98, 0), (11.90, -3.8)), 12.0: ((2.51, 0), (7.98, 0), (11.91, 0)),
      13.0: ((2.50, 0), (7.97, 0), (11.95, 0)), 15.0: ((2.50, 0), (7.96, 0), (12.00, 0))}
# Paper II, Tabla IV: masa (M_sun)
T4M = {8.0: 1.494, 9.0: 1.501, 10.0: 1.507, 11.0: 1.512, 12.0: 1.515, 13.0: 1.519, 14.0: 1.521, 15.0: 1.524}
NUS = sorted(set(T1) | set(T2))
BRACKETS = [(2.2, 3.1), (6.9, 8.6), (9.9, 12.6)]
L = 2


def _real_root(fn, lo, hi, n=24):
    xs = np.linspace(lo, hi, n)
    v = [fn(x) for x in xs]
    for i in range(n - 1):
        if np.sign(v[i]) != np.sign(v[i + 1]):
            return brentq(fn, xs[i], xs[i + 1], xtol=1e-7)
    return np.nan


def _background_modes(st):
    return [_real_root(lambda x: A.surface_residual(st, L, A.omega_of(x)), *b) for b in BRACKETS]


def _coupled(st, guess):
    """Modo acoplado cerca de `guess` (kHz). Por debajo del umbral nu_a la raíz es real (estado ligado)."""
    det = lambda z: A.coupled_det(st, L, A.omega_of(z))
    if guess < st.nu_a:
        z = _real_root(lambda x: det(x).real, guess - 0.35, min(guess + 0.15, st.nu_a - 1e-4), n=16)
        if np.isfinite(z):
            return complex(z)
    # H29: con el arranque del interior integrado de verdad, el ruido de det deja z quieto solo hasta ~1e-9 |z|;
    # con tol = 1e-11 Muller agotaba sus 40 pasos y la fila quedaba en nan. 1e-9 |z| son < 1e-5 Hz.
    return A.muller(det, guess * 0.997, guess, guess * 1.003 - 2e-4j, tol=1e-9)


def row(nu):
    t0 = time.time()
    st = A.Star(nu)
    out = {'nu': nu, 'M': st.M / A.MSUN, 'R': st.R, 'thc': st.thc}
    xs = np.linspace(0.0, 3.0, 61)                       # perfil theta_0/pi para la animación de la página
    out['prof'] = [float(st.sol.sol(min(max(x * st.R, 1e-4), st.sol.t[-1]))[3] / np.pi) if x * st.R < st.sol.t[-1] else 0.0 for x in xs]
    if nu in T1 or nu in T2:
        bm = _background_modes(st)
        out['T1'] = bm
    if nu in T2:
        out['T2'] = []
        for g0 in out['T1']:
            try:
                out['T2'].append(_coupled(st, g0))
            except Exception:
                out['T2'].append(complex(np.nan, np.nan))
    out['secs'] = time.time() - t0
    return out


def no_condensate(nu=2.5):
    try:
        A.Star(nu)
        return False
    except ValueError:
        return True


def compute():
    """Referencia sin axión, una fila por nu_a (en paralelo) y la comprobación de que a 2.5 kHz no hay condensado."""
    st = A.Star(None)
    ref = {'M': st.M / A.MSUN, 'R': st.R, 'modes': _background_modes(st)}
    with mp.get_context('fork').Pool(6) as pool:
        none25 = pool.apply_async(no_condensate)
        rows = {r['nu']: r for r in pool.imap_unordered(row, sorted(NUS, key=lambda x: -x))}
        none25 = none25.get()
    return ref, rows, none25


# ---------------------------------------------------------------- comprobaciones
def _shift(v, v0):
    return 100 * (np.asarray(v) / np.asarray(v0) - 1)


def checks(ref, rows, none25):
    res = []
    r0, p0 = np.array(ref['modes']), np.array(REF['modes'])
    res.append((f"Fondo BSk26 sin axión, p_c = 100 MeV/fm^3: M = {ref['M']:.3f} M_sun, R = {ref['R']:.2f} km "
                f"(Paper I: 1.53 M_sun, 11.7 km)", abs(ref['M'] - REF['M']) < 0.01 and abs(ref['R'] - REF['R']) < 0.05, ''))
    dev = 100 * np.abs(r0 / p0 - 1)
    res.append((f"Cowling sin axión: f, p1, p2 = {r0[0]:.3f}, {r0[1]:.3f}, {r0[2]:.3f} kHz (Paper II: 2.43, 7.15, 10.18); "
                f"máx. dif. {dev.max():.2f} %", dev.max() <= 0.3,
                'ajuste analítico de BSk26 frente a la tabla del grupo; el paper da un 0.3 % de acuerdo con otros cálculos'))
    ok_thr = none25 and 2.8 in rows
    res.append(('Umbral del condensado: no hay a 2.5 kHz y sí a 2.8 kHz (Paper I: desaparece por debajo de ~2.8 kHz)', ok_thr, ''))
    dm = {nu: 100 * (rows[nu]['M'] / m - 1) for nu, m in T4M.items() if nu in rows}
    res.append((f"Tabla IV: masas de las estrellas con axión ({len(dm)} valores, 8-15 kHz) dentro del 0.5 % "
                f"(máx. {max(abs(v) for v in dm.values()):.2f} %)", max(abs(v) for v in dm.values()) < 0.5, ''))
    d1 = [_shift(rows[nu]['T1'], r0) - _shift(pv, p0) for nu, pv in T1.items()]
    m1 = np.abs(d1).max()
    res.append((f"Tabla I (solo fondo): los 24 desplazamientos nu/nu_0 - 1 coinciden con el paper dentro de 0.6 puntos "
                f"(máx. {m1:.2f})", m1 <= 0.6, 'nu_0: la nota sin axión de cada cálculo'))
    d2 = [_shift([z.real for z in rows[nu]['T2']], r0) - _shift([v[0] for v in pv], p0) for nu, pv in T2.items()]
    m2 = np.abs(d2).max()
    df = max(abs(d[0]) for d in d2)
    res.append((f"Tabla II (acoplado), parte real: los 36 desplazamientos dentro de 0.6 puntos (máx. {m2:.2f}; "
                f"modo f, {df:.2f})", m2 <= 0.6, ''))
    # patrón de amortiguamiento: Im = 0 por debajo del umbral nu_a
    mism, near = [], []
    for nu, pv in T2.items():
        for i, z in enumerate(rows[nu]['T2']):
            ours, theirs = abs(z.imag) * 1e3 > 1e-3, pv[i][1] != 0
            if ours != theirs:
                (near if abs(z.real / nu - 1) < 0.003 else mism).append((nu, ['f', 'p1', 'p2'][i]))
    res.append(("Tabla II: los modos amortiguados son los mismos (Im nu = 0 por debajo de nu_a)"
                + (f" salvo {len(near)} justo en el umbral" if near else ", sin excepciones"), not mism,
                ', '.join(f'{m} a {nu} kHz: Re nu a menos del 0.3 % de nu_a' for nu, m in near)))
    ratios = [abs(z.imag) * 1e3 / abs(pv[i][1]) for nu, pv in T2.items() for i, z in enumerate(rows[nu]['T2'])
              if abs(pv[i][1]) >= 1]
    im_max = max(abs(z.imag) for r in rows.values() for z in r.get('T2', []))
    tau = 1 / (2 * np.pi * im_max * 1e3) * 1e3
    res.append((f"Tabla II, amortiguamiento: los {len(ratios)} |Im nu| >= 1 Hz del paper, dentro de un factor 2 "
                f"(cociente {min(ratios):.2f}-{max(ratios):.2f}); tau_a mínimo = {tau:.1f} ms (paper: 12 ms)",
                0.5 <= min(ratios) and max(ratios) <= 2 and abs(tau - 12.2) < 1.5,
                'el amortiguamiento depende de delta theta(R), muy sensible al fondo (Fig. 6 del paper)'))
    return res


def save_figdata(ref, rows):
    from .figdata import save, r as rnd
    r0, p0 = np.array(ref['modes']), np.array(REF['modes'])
    nus = sorted(rows)
    out = {'nus': nus, 'modes': [], 'prof': [rows[nu].get('prof') for nu in nus]}
    for i, name in enumerate(('f', 'p1', 'p2')):
        out['modes'].append({
            'name': name,
            'bg': [_shift(rows[nu]['T1'][i], r0[i]) for nu in nus],
            'cp': [_shift(rows[nu]['T2'][i].real, r0[i]) for nu in nus],
            'im': [abs(rows[nu]['T2'][i].imag) * 1e3 for nu in nus],
            'pbg': [[nu, _shift(T1[nu][i], p0[i])] for nu in T1],
            'pcp': [[nu, _shift(T2[nu][i][0], p0[i])] for nu in T2],
            'pim': [[nu, abs(T2[nu][i][1])] for nu in T2 if T2[nu][i][1] != 0]})
    save('axcow', rnd(out))


def run():
    from .pipeline import report
    ref, rows, none25 = compute()
    save_figdata(ref, rows)
    return report('Paso 9: notas con axión en Cowling (Tablas I, II y IV del Paper II)', checks(ref, rows, none25))


if __name__ == '__main__':
    run()
