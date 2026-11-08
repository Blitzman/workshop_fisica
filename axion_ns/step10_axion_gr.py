"""Paso 10 (etapa B): notas con axión y con la métrica dinámica, sin la aproximación de Cowling.

No hay resultado publicado con el que comparar. Lo que se comprueba:
- que la conversión de fullgr.EXT a números es la misma que la del paso 8 con el axión apagado (a 1e-12)
  y que las sustituciones del axión dan las cantidades del Paper II;
- que este mismo programa, de principio a fin, reproduce el modo f del politropo de Lindblom, Mendell e
  Ipser (1997) con sus cinco cifras;
- que los resultados no dependen de la tolerancia, del radio de empalme ni de la distancia a la que se
  separa la onda entrante de la saliente.
Luego calcula f, p1 y p2 para las 12 masas del axión de la Tabla II del Paper II, con delta theta apagado
(solo el fondo) y encendido, y compara con Cowling (paso 9) y con la suma de amortiguamientos que
propone el Paper II, 1/tau = 1/tau_GW + 1/tau_a.
"""
import multiprocessing as mp
import time

import numpy as np

from . import axgr as AG
from . import axstar as AX
from . import step9_axion_cowling as S9

L = 2
MODES = ('f', 'p1', 'p2')


def _mode(st, guess, dyn=True, **kw):
    """Cero complejo de la amplitud entrante cerca de `guess` (kHz)."""
    f = lambda z: AG.incoming(st, AX.omega_of(z), L, dyn=dyn, **kw)
    try:
        return complex(AX.muller(f, guess * 0.997, guess, guess * 1.003 - 1e-5j, tol=1e-11))
    except RuntimeError:                 # suelo de ruido por encima de 1e-11: se acepta una tolerancia menor
        return complex(AX.muller(f, guess * 0.997, guess, guess * 1.003 - 1e-5j, tol=1e-9))


def reference():
    """Estrella sin axión: notas en Cowling (paso 9) y con la métrica dinámica."""
    st = AX.Star(None)
    cow = S9._background_modes(st)
    gr = []
    for i, g in enumerate((1.96, 6.59, 9.87)):
        gr.append(_mode(st, g))
    return {'cow': cow, 'gr': gr}


def row(args):
    nu, ref = args
    t0 = time.time()
    st = AX.Star(nu)
    p0 = np.array(S9.REF['modes'])
    out = {'nu': nu, 'bg': [], 'cp': []}
    for i in range(3):
        g0 = ref['gr'][i].real
        # puntos de partida: la nota sin axión, desplazada como en las Tablas I y II (Cowling)
        ks = sorted(S9.T1)
        gb = g0 * np.interp(nu, ks, [S9.T1[k][i] for k in ks]) / p0[i]
        gc = g0 * S9.T2[nu][i][0] / p0[i]
        for key, guess, dyn in (('bg', gb, False), ('cp', gc, True)):
            try:
                out[key].append(_mode(st, guess, dyn=dyn))
            except Exception:
                out[key].append(complex(np.nan, np.nan))
    out['secs'] = time.time() - t0
    return out


def convergence(nu=6.0, guess=None):
    """El modo f acoplado con otra tolerancia, otro radio de empalme y otra distancia de extracción."""
    st = AX.Star(nu)
    z0 = _mode(st, guess)
    rm = AG.r_match(st)
    alts = [_mode(st, z0.real, rtol=1e-11), _mode(st, z0.real, rm=0.8 * rm), _mode(st, z0.real, n_far=80)]
    return z0, alts


def lmi_polytrope():
    """El modo f l = 2 del politropo p = kappa eps^2 de 1.4 M_sun y 14.45 km, con este programa."""
    from scipy.optimize import brentq

    class Poly:
        kappa, p_min = 1e-3, 1e-14

        def eps(self, p):
            return np.sqrt(p / self.kappa)

        def nb(self, p):
            e = self.eps(p)
            return e / (1 + self.kappa * e)

        def cs2(self, p):
            return 2 * self.kappa * self.eps(p)

    keep, AX._eos = AX._eos, Poly()
    try:
        C = 1.4 * AX.MSUN / 14.45
        lp0 = brentq(lambda lp: (lambda s: s.M / s.R - C)(AX.Star(None, p0=10**lp)), -1, 3, xtol=1e-12)
        st = AX.Star(None, p0=10**lp0)
        unit = np.sqrt(st.M / st.R**3)
        z = AX.muller(lambda s: AG.incoming(st, s * unit, L), 1.19, 1.20, 1.21 - 1e-4j, tol=1e-12)
    finally:
        AX._eos = keep
    return z


def compute():
    ref = reference()
    with mp.get_context('fork').Pool(6) as pool:
        conv = pool.apply_async(convergence, (6.0, ref['gr'][0].real * S9.T2[6.0][0][0] / S9.REF['modes'][0]))
        lmi = pool.apply_async(lmi_polytrope)
        rows = {r['nu']: r for r in pool.imap_unordered(row, [(nu, ref) for nu in sorted(S9.T2, key=lambda x: -x)])}
        conv, lmi = conv.get(), lmi.get()
    return ref, rows, conv, lmi


# ---------------------------------------------------------------- comprobaciones
def matrix_checks():
    """La conversión a números, sin axión, frente a la del paso 8 en puntos al azar (politropo)."""
    from . import fmodes as FM
    Ain, Aref = AG.interior_matrix(), FM.full_system('rest')
    rng, worst = np.random.default_rng(1), 0.0
    for _ in range(5):
        rr, om = rng.uniform(0.1, 1), rng.uniform(0.5, 2)
        n, n1, n2 = rng.uniform(0.05, 0.3), -rng.uniform(0.01, 0.5), rng.uniform(-1, 1)
        P, P1, P2, Lm, L1 = rng.uniform(-1, 0), rng.uniform(0.01, 1), rng.uniform(-1, 1), rng.uniform(0, 0.5), rng.uniform(0.01, 1)
        ref = np.array(Aref(rr, om, L, n, n1, n2, P, P1, P2, Lm, L1), dtype=float)
        new = np.array(Ain(rr, om, L, P, P1, P2, Lm, L1, n, n1, 0, 0, n + n * n, n * n, 2 * n / (1 + 2 * n),
                           1.0, 0, 0, 0, 0, 1.0), dtype=float)
        worst = max(worst, np.abs(new[:4, :4] - ref).max() / np.abs(ref).max())
    return worst


def checks(ref, rows, conv, lmi, worst, subs):
    res = []
    res.append((f'Conversión de fullgr.EXT a números: con el axión apagado, igual que la del paso 8 '
                f'(máx. dif. relativa {worst:.0e})', worst < 1e-12, ''))
    res.append((f'Sustituciones del axión: {sum(ok for _, ok in subs)}/{len(subs)} cantidades del Paper II '
                f'(n dm*/dtheta, m_eff^2, eps(n, theta), c_s^2...)', all(ok for _, ok in subs), ''))
    res.append((f'Todo el programa, sin axión, con el politropo de LMI: Re w = {lmi.real:.5f} (LMI: 1.2016), '
                f'Im w = {-lmi.imag:.3e} (LMI: 3.95e-4)', abs(lmi.real - 1.2016) < 5e-5 and abs(-lmi.imag / 3.95e-4 - 1) < 0.01, ''))
    z0, alts = conv
    dre = max(abs(a.real / z0.real - 1) for a in alts)
    dim = max(abs(a.imag / z0.imag - 1) for a in alts)
    res.append((f'Convergencia (modo f, 6 kHz): tolerancia, radio de empalme y extracción mueven Re nu en '
                f'{dre:.0e} e Im nu en {100 * dim:.1f} %', dre < 1e-5 and dim < 0.01, ''))
    gr0, cow0 = np.array([z.real for z in ref['gr']]), np.array(ref['cow'])
    res.append((f"Sin axión, la métrica dinámica baja las notas de Cowling: f, p1, p2 = "
                f"{', '.join(f'{z.real:.3f}' for z in ref['gr'])} kHz frente a {', '.join(f'{c:.3f}' for c in cow0)} "
                f"({', '.join(f'{100 * (g / c - 1):+.0f} %' for g, c in zip(gr0, cow0))}); "
                f"f amortiguado por ondas gravitacionales en {1 / (2 * np.pi * -ref['gr'][0].imag * 1e3):.3f} s",
                all(gr0 < cow0) and ref['gr'][0].imag < 0, ''))
    cp = [z for r in rows.values() for z in r['cp']]
    res.append((f'Las {len(cp)} notas acopladas (f, p1, p2 para 12 masas del axión) encontradas y amortiguadas',
                all(np.isfinite(z.real) and z.imag < 0 for z in cp), ''))
    above = [(nu, i, z) for nu, r in rows.items() for i, z in enumerate(r['cp']) if i > 0]
    ok = all((abs(z.imag) * 1e3 > 0.05) == (z.real > nu) for nu, i, z in above)
    res.append(('Umbral: los modos p solo radian axiones por encima de nu_a (|Im nu| > 0.05 Hz si Re nu > nu_a; '
                'por debajo, solo ondas gravitacionales, < 0.05 Hz)', ok, ''))
    sf = [100 * (rows[nu]['cp'][0].real / gr0[0] - 1) for nu in sorted(rows)]
    res.append((f'Con la métrica dinámica, el axión sube el modo f hasta un {max(sf):.1f} % (Cowling, paso 9: 5.4 %) '
                f'y su amortiguamiento por ondas gravitacionales pasa de {1 / (2 * np.pi * -ref["gr"][0].imag * 1e3):.2f} '
                f'a {1 / (2 * np.pi * -rows[15.0]["cp"][0].imag * 1e3):.2f} s', True, 'resultado nuevo'))
    p2 = {nu: -rows[nu]['cp'][2].imag * 1e3 for nu in rows}
    num = max(p2, key=p2.get)
    res.append((f'Amortiguamiento máximo por axiones: p2 a {num:.0f} kHz, {p2[num]:.1f} Hz '
                f'(tau_a = {1e3 / (2 * np.pi * p2[num]):.1f} ms), como en Cowling; la suma 1/tau_GW + 1/tau_a con el '
                f'tau_a de Cowling solo acierta en orden de magnitud (p2 a 8 kHz: {p2[8.0]:.2f} Hz frente a 1.0)',
                True, 'resultado nuevo'))
    return res


def save_figdata(ref, rows):
    from .figdata import save, r as rnd
    nus = sorted(rows)
    gr0 = [z.real for z in ref['gr']]
    out = {'nus': nus, 'gr0': [[z.real, -z.imag * 1e3] for z in ref['gr']], 'cow0': ref['cow'], 'modes': []}
    for i, name in enumerate(MODES):
        out['modes'].append({'name': name,
                             'cp': [100 * (rows[nu]['cp'][i].real / gr0[i] - 1) for nu in nus],
                             'bg': [100 * (rows[nu]['bg'][i].real / gr0[i] - 1) for nu in nus],
                             'im': [-rows[nu]['cp'][i].imag * 1e3 for nu in nus],
                             're': [rows[nu]['cp'][i].real for nu in nus]})
    save('axgr', rnd(out, 6))


def run():
    from .pipeline import report
    worst, subs = matrix_checks(), AG.substitution_checks()
    ref, rows, conv, lmi = compute()
    save_figdata(ref, rows)
    return report('Paso 10: con axión y métrica dinámica (etapa B, sin publicar)',
                  checks(ref, rows, conv, lmi, worst, subs))


if __name__ == '__main__':
    run()
