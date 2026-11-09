"""Paso 8: comprobaciones numéricas del cálculo nuevo, y figuras de apoyo.

El paso 7 valida el sistema extendido con álgebra (límites y consistencia). Aquí se integra:

  * sin axión, el sistema extendido de la página (fullgr.EXT) con las condiciones de LD en el centro y
    Zerilli fuera, frente a las frecuencias "exactas" de Lindblom, Mendell e Ipser (1997, Tabla I) y de
    Yoshida y Kojima (1997, Tablas 1-2); y las (8a)-(8c) del paso 5 frente a sus valores en Cowling;
  * las mismas ecuaciones de LD tal como están impresas en las dos transcripciones con erratas;
  * con axión, el fondo de Klein-Gordon sobre un politropo (sin retroacción), para ver dónde viven los
    términos nuevos, y la masa efectiva del axión en función de theta (observación del Paper I).
"""
import multiprocessing as mp
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp, solve_bvp

from . import fmodes as FM

RES = Path(__file__).resolve().parent.parent / 'results'
CSS = {'ink': '#1f2328', 'muted': '#6b7280', 'cow': '#c2410c', 'full': '#1d4ed8', 'err': '#9333ea'}

# Lindblom, Mendell e Ipser (1997, PRD 56, 2118; gr-qc/9704046), Tabla I, columna "Exact":
# p = kappa rho^2 (rho: densidad de energía), kappa = 1.3346e5 cm^5/(g s^2), 1.4 M_sun, R = 14.45 km.
LMI97 = {2: (1.2016, 3.95e-4), 3: (1.5862, 3.96e-6), 4: (1.8738, 4.01e-8), 5: (2.1129, 3.63e-10)}
LMI_R_KM = 14.45
_G, _C, _MSUN = 6.6743e-8, 2.99792458e10, 1.98841e33
_LUNIT = np.sqrt(1.3346e5 / _G)                       # K^{1/2} en cm
MODES = ('f', 'p1', 'p2')
LS = (2, 3, 4, 5)

_W = {}                                               # matrices y estrellas, construidas antes del fork


def _setup():
    if not _W:
        _W['cow'] = FM.cowling_system('energy')
        _W['full'] = FM.full_system('energy')
        _W['ZL'] = FM.full_system('energy', 'ZL')
        _W['stars'] = {MR: FM.Star('energy', MR) for MR in (0.1, 0.2)}
        _W['stars']['LMI'] = FM.Star('energy', mass=1.4 * _G * _MSUN / _C**2 / _LUNIT)
    return _W


def _unit(st):
    return np.sqrt(st.M / st.R**3)


def _cowling(key):
    st = _W['stars'][key]
    u = _unit(st)
    out = {}
    for ll in LS:
        roots = FM.cowling_modes(st, _W['cow'], ll, np.linspace(0.8, 7.0, 64) * u)
        out[ll] = [x / u for x in roots[:3]]
    return key, out


def _full(job):
    key, ll, guess, which = job
    s0, g = FM.full_mode(_W['stars'][key], _W[which], ll, guess, width=0.06)
    return job, (s0, g)


def _converged(_):
    """El modo f l = 2 de LMI con otro punto de partida, otro corte en la superficie y extracción más lejana."""
    st = _W['stars']['LMI']
    orig = st.integrate
    st.integrate = lambda A, y0, o, ll, r0=None, r1=None: orig(A, y0, o, ll, r0=1e-5 * st.R, r1=st.R * (1 - 1e-9))
    u = _unit(st)
    from scipy.optimize import minimize_scalar
    m = minimize_scalar(lambda s: FM.outgoing_amplitude(st, _W['full'], s * u, 2, n_far=120),
                        bracket=(1.18, 1.2016, 1.22), tol=1e-10)
    st.integrate = orig
    return m.x


# |A_out|^2 frente a w (l = 2, modelo de LMI): rejilla más fina cerca de cada resonancia
_near = lambda w0: w0 * (1 + np.concatenate([-np.geomspace(1e-4, 3e-2, 12), [0], np.geomspace(1e-4, 3e-2, 12)]))
ZL_GUESS = 0.833                 # el mínimo con la errata de Zhao y Lattimer (rejilla 0.3-2)
SCAN = np.unique(np.concatenate([np.linspace(0.6, 2.0, 57), _near(1.2016), _near(ZL_GUESS)]))


def _scan(job):
    which, s = job
    st = _W['stars']['LMI']
    return job, FM.outgoing_amplitude(st, _W[which], s * _unit(st), 2)


def compute():
    _setup()
    jobs = [('LMI', ll, LMI97[ll][0], 'full') for ll in LS]
    jobs += [(MR, ll, FM.YK97[MR][ll][m][0], 'full') for MR in (0.1, 0.2) for ll in LS for m in MODES]
    jobs += [('LMI', 2, ZL_GUESS, 'ZL')]
    with mp.get_context('fork').Pool(min(6, mp.cpu_count())) as P:
        cow = dict(P.map(_cowling, [0.1, 0.2, 'LMI']))
        full = dict(P.map(_full, jobs))
        conv = P.map(_converged, [0])[0]
        sc = dict(P.map(_scan, [(w, s) for w in ('full', 'ZL') for s in SCAN]))
    scan = {w: np.array([sc[(w, s)] for s in SCAN]) for w in ('full', 'ZL')}
    return cow, full, conv, scan


# ---------------------------------------------------------------- fondo con axión (esquemático)
BETA = 0.88


def axion_background(maR, alpha_c=19.0, compactness=0.172):
    """theta_0(r) sobre un politropo (rho en reposo), sin retroacción: Klein-Gordon estática
        theta'' + (2/r + Phi' - lambda') theta' = -e^{2 lambda} (m_a R)^2 (alpha - 1) sin(theta)/g(theta),
    con alpha(r) = alpha_c n(r)/n(0), theta'(0) = 0 y caída de Yukawa fuera. Unidades r/R."""
    st = FM.Star('rest', compactness)
    R = st.R
    r0 = 1e-6
    b0 = st._start(st.nc, r0)
    b0[2] = st.Phic
    bgs = solve_ivp(lambda rr, y: st.bg(rr, *y)[1], (r0, R * (1 - 1e-9)), b0, dense_output=True,
                    rtol=1e-11, atol=1e-30)
    if bgs.status < 0:                                  # H32
        raise RuntimeError(f'paso 8, fondo: {bgs.message}')

    def bgx(x):
        if x >= 1:
            m, rr = st.M, x * R
            P1 = m / (rr**2 * (1 - 2 * m / rr)) * R
            return 0.0, 1 / (1 - 2 * m / rr), P1, -P1
        B, _ = st.bg(x * R, *bgs.sol(x * R))
        return B[0] / st.nc, np.exp(2 * B[6]), B[4] * R, B[7] * R

    xmax = 1 + 12 / maR

    def f(x, y):
        out = np.empty_like(y)
        for i, (xi, th, dth) in enumerate(zip(x, y[0], y[1])):
            nr, e2l, P1, L1 = bgx(xi)
            g = np.sqrt(1 - BETA * np.sin(th / 2)**2)
            out[:, i] = [dth, -(2 / xi + P1 - L1) * dth - e2l * maR**2 * (alpha_c * nr - 1) * np.sin(th) / g]
        return out

    bc = lambda ya, yb: np.array([ya[1], yb[1] + (maR + 1 / xmax) * yb[0]])
    x = np.concatenate([np.linspace(1e-3, 1, 150), np.linspace(1.002, xmax, 100)])
    s = 1 / (1 + np.exp((x - 0.8) / 0.05))
    y = np.vstack([np.pi * s, -np.pi / 0.05 * s * (1 - s)])
    sol = solve_bvp(f, bc, x, y, tol=1e-6, max_nodes=50000)
    return sol, (lambda xx: alpha_c * np.array([bgx(v)[0] for v in np.atleast_1d(xx)]))


def m_eff2(theta, alpha):
    """(V'' + n m*'')/m_a^2 con las funciones completas de la Sec. II del Paper I."""
    g = np.sqrt(1 - BETA * np.sin(theta / 2)**2)
    s, c = np.sin(theta), np.cos(theta)
    gpp = -BETA / 4 * c / g - BETA**2 / 16 * s**2 / g**3            # d^2 g / d theta^2
    return 4 * (alpha - 1) * gpp / BETA


# ---------------------------------------------------------------- comprobaciones
def checks(cow, full, conv, scan):
    res = []
    st = _W['stars']['LMI']
    Rkm = st.R * _LUNIT / 1e5
    res.append((f'Fondo: el modelo de 1.4 M_sun con kappa = 1.3346e5 cgs tiene R = {Rkm:.3f} km (LMI: 14.45 km)',
                abs(Rkm - LMI_R_KM) < 0.005, ''))
    pairs = [(cow[MR][ll][i], FM.YK97[MR][ll][m][1], f'{m} l={ll} M/R={MR}')
             for MR in (0.1, 0.2) for ll in LS for i, m in enumerate(MODES)]
    out = [f'{n}: {a:.4f} vs {b}' for a, b, n in pairs if abs(a - b) > 0.005]
    dmax = max(abs(a - b) for a, b, _ in pairs)
    res.append((f'Cowling, Ecs. (8a)-(8c) del paso 5: los 24 modos f, p1, p2 (l = 2-5, M/R = 0.1 y 0.2) de Yoshida '
                f'y Kojima (1997): {24 - len(out)} dentro de media unidad de su última cifra, máx. dif. {dmax:.4f}',
                dmax < 0.006, '; '.join(out)))
    dre = [abs(full[('LMI', ll, LMI97[ll][0], 'full')][0] - LMI97[ll][0]) for ll in LS]
    res.append((f'Métrica dinámica, sistema extendido sin axión: Re w de los modos f l = 2-5 igual a las 5 cifras de '
                f'Lindblom, Mendell e Ipser (1997) (máx. dif. {max(dre):.5f})', max(dre) <= 5e-5, ''))
    dim = [abs(full[('LMI', ll, LMI97[ll][0], 'full')][1] / LMI97[ll][1] - 1) for ll in LS]
    res.append((f'Im w (emisión de ondas gravitacionales) dentro del 2 % de LMI (máx. {100 * max(dim):.1f} %)',
                max(dim) < 0.02, ''))
    dyk = [abs(full[(MR, ll, FM.YK97[MR][ll][m][0], 'full')][0] / FM.YK97[MR][ll][m][0] - 1)
           for MR in (0.1, 0.2) for ll in LS for m in MODES]
    res.append((f'Métrica dinámica: los 24 modos "exactos" de Yoshida y Kojima dentro del {100 * max(dyk):.1f} % '
                '(a 3 cifras; LMI, a 5, coincide entero)', max(dyk) < 0.007, ''))
    s0 = full[('LMI', 2, LMI97[2][0], 'full')][0]
    res.append((f'Convergencia: otro punto de partida, corte en la superficie y extracción mueven Re w en '
                f'{abs(conv - s0):.0e}', abs(conv - s0) < 1e-5, ''))
    zl = full[('LMI', 2, ZL_GUESS, 'ZL')][0]
    res.append((f'Errata de signo de Zhao y Lattimer: con su X\' tal como está impreso, el modo f l = 2 sale en '
                f'{zl:.4f} en vez de 1.2016 ({100 * (zl / LMI97[2][0] - 1):+.0f} %)', abs(zl / LMI97[2][0] - 1) > 0.1, ''))
    return res


def background_checks(sols):
    res = []
    for maR, (sol, _) in sols.items():
        th = sol.sol(np.array([1e-3, 1 + 10 / maR]))[0]
        res.append((f'Fondo con axión (politropo + KG, m_a R = {maR}): converge, theta_0 = {th[0]:.2f} en el centro '
                    f'y {th[1]:.0e} fuera', sol.status == 0 and th[0] > 2.5 and abs(th[1]) < 1e-3, ''))
    return res


# ---------------------------------------------------------------- figuras
def plot_modes(cow, full, scan, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    a = scan['full']
    ax1.semilogy(SCAN, a / a.max(), color=CSS['full'], lw=2, label='métrica dinámica (sistema extendido, sin axión)')
    z = scan['ZL']
    ax1.semilogy(SCAN, z / z.max(), color=CSS['err'], lw=1.6, ls='--', label="ídem con el signo impreso por Zhao–Lattimer")
    zl = full[('LMI', 2, ZL_GUESS, 'ZL')][0]
    ax1.text(zl + 0.02, 3e-7, f'errata:\n{zl:.3f}', fontsize=9, color=CSS['err'])
    ax1.axvline(LMI97[2][0], color=CSS['ink'], lw=1, ls=':')
    ax1.text(LMI97[2][0] + 0.02, 2e-6, 'LMI (1997):\n1.2016', fontsize=9, color=CSS['ink'])
    wc = cow['LMI'][2][0]
    ax1.axvline(wc, color=CSS['cow'], lw=2)
    ax1.text(wc + 0.02, 0.3, f'Cowling (paso 5):\n{wc:.3f}', fontsize=9, color=CSS['cow'])
    ax1.set_xlabel(r'$\omega$ / $\sqrt{M/R^3}$ (real)')
    ax1.set_ylabel(r'$|A_\mathrm{out}|^2$ de la onda de Zerilli (normalizada)')
    ax1.set_title(r'Modo $f$, $l=2$, politropo de 1.4 $M_\odot$: resonancia', fontsize=11)
    ax1.legend(fontsize=8.5, frameon=False, loc='lower right')

    groups = [('Cowling vs. Yoshida–Kojima', CSS['cow'],
               [(cow[MR][ll][i], FM.YK97[MR][ll][m][1], 0.005) for MR in (0.1, 0.2) for ll in LS
                for i, m in enumerate(MODES)]),
              ('métrica dinámica vs. Yoshida–Kojima', CSS['full'],
               [(full[(MR, ll, FM.YK97[MR][ll][m][0], 'full')][0], FM.YK97[MR][ll][m][0], 0.005)
                for MR in (0.1, 0.2) for ll in LS for m in MODES]),
              ('métrica dinámica vs. Lindblom–Mendell–Ipser', CSS['ink'],
               [(full[('LMI', ll, LMI97[ll][0], 'full')][0], LMI97[ll][0], 5e-5) for ll in LS])]
    from .figdata import save, r as rnd
    tags = [f'M/R = {MR}, l = {ll}, {m}' for MR in (0.1, 0.2) for ll in LS for m in MODES]
    save('modes', rnd({'scan': SCAN, 'full': a / a.max(), 'zl': z / z.max(), 'lmi': LMI97[2][0], 'cow': wc, 'zlmin': zl,
                       'groups': [{'name': name, 'pts': [[100 * (p[0] / p[1] - 1), 100 * p[2] / p[1], p[0], p[1]]
                                                         for p in pts],
                                   'tags': tags if len(pts) == 24 else [f'l = {ll}, f' for ll in LS]}
                                  for name, _, pts in groups]}, 6))
    x0 = 0
    for name, col, pts in groups:
        xs = np.arange(x0, x0 + len(pts))
        rel = [100 * (a / b - 1) for a, b, _ in pts]
        err = [100 * e / b for _, b, e in pts]
        ax2.errorbar(xs, rel, yerr=err, fmt='o', color=col, ms=4.5, elinewidth=1, capsize=0, label=name)
        x0 += len(pts) + 3
    ax2.axhline(0, color=CSS['muted'], lw=1)
    ax2.set_xticks([])
    ax2.set_ylabel('este cálculo − referencia (%)')
    ax2.set_title('52 frecuencias publicadas', fontsize=11)
    ax2.legend(fontsize=8.5, frameon=False, loc='upper right')
    for ax in (ax1, ax2):
        ax.grid(alpha=0.25, lw=0.6)
        ax.spines[['top', 'right']].set_visible(False)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _axes(figsize):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    ax.grid(alpha=0.25, lw=0.6)
    ax.spines[['top', 'right']].set_visible(False)
    return plt, fig, ax


def plot_meff(path):
    plt, fig, ax = _axes((7.2, 3.9))
    th = np.linspace(0, np.pi, 400)
    for al, col, lab in [(10, CSS['full'], r'$\alpha = n_b/n_c = 10$ (interior)'),
                         (0.5, CSS['cow'], r'$\alpha = 0.5$ (exterior de la capa)')]:
        ax.plot(th, m_eff2(th, al), color=col, lw=2, label=lab)
    ax.axhline(0, color=CSS['muted'], lw=1)
    for t, txt, dx in [(0, r'$\theta=0$: $1-\alpha = -9$ (el del paper)', 10),
                       (np.pi, r'$\theta=\pi$: $(\alpha-1)/\sqrt{1-\beta} \approx +26$', -10)]:
        v = m_eff2(np.array([t]), 10)[0]
        ax.plot([t], [v], 'o', ms=8, mfc='white', mec=CSS['full'], mew=2)
        ax.annotate(txt, (t, v), xytext=(dx, 8 if t == 0 else -16), textcoords='offset points', fontsize=9.5,
                    ha='left' if t == 0 else 'right', color=CSS['ink'])
    from .figdata import save, r as rnd
    save('meff', rnd({'th': th[::2], 'a10': m_eff2(th, 10)[::2], 'a05': m_eff2(th, 0.5)[::2]}))
    ax.set_xticks([0, np.pi / 2, np.pi])
    ax.set_xticklabels(['0', r'$\pi/2$', r'$\pi$'])
    ax.set_xlim(-0.15, np.pi + 0.15)
    ax.set_ylim(-12, 30)
    ax.set_xlabel(r'$\theta = a/f_a$ (fondo en el que se linealiza)')
    ax.set_ylabel(r'$m_\mathrm{eff}^2/m_a^2$')
    ax.set_title(r'Masa efectiva del axión con las funciones completas ($\beta = 0.88$)', fontsize=11)
    ax.legend(fontsize=9, frameon=False, loc='upper left')
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_background(sols, path):
    plt, fig, ax = _axes((7.2, 3.9))
    x = np.linspace(1e-3, 1.6, 600)
    data = {'x': x[::2], 'bg': []}
    for (maR, (sol, alpha)), col, nu in zip(sols.items(), [CSS['cow'], CSS['full']], (5, 15)):
        t = sol.sol(x)[0]
        dm = np.abs(np.sin(t)) / np.sqrt(1 - BETA * np.sin(t / 2)**2)       # |dm*/dtheta|, sin constantes
        nb = np.where(x < 1, alpha(np.minimum(x, 1)), 0) / alpha(1e-3)[0]
        carrier = nb * dm
        ax.fill_between(x, 0, carrier / carrier.max(), color=col, alpha=0.18, lw=0)
        ax.plot(x, t / np.pi, color=col, lw=2, label=fr'$\nu_a$ = {nu} kHz ($m_a R$ = {maR})')
        data['bg'].append({'nu': nu, 'maR': maR, 'theta': (t / np.pi)[::2], 'carrier': (carrier / carrier.max())[::2]})
    from .figdata import save, r as rnd
    save('background', rnd(data))
    ax.axvline(1, color=CSS['muted'], ls=':', lw=1)
    ax.text(1.015, 0.95, 'superficie', color=CSS['muted'], fontsize=9)
    ax.set_xlabel('$r/R$')
    ax.set_ylabel(r'$\theta_0/\pi$  (línea);  $n_b\,|\partial m^*/\partial\theta|$ normalizado (sombra)')
    ax.yaxis.label.set_size(9.5)
    ax.set_ylim(0, 1.05)
    ax.set_title(r'Fondo con axión: dónde actúan los términos nuevos', fontsize=11)
    ax.legend(fontsize=9, frameon=False, loc='lower left')
    fig.savefig(path, dpi=150)
    plt.close(fig)


BG_MAR = (1.26, 3.77)            # nu_a = 5 y 15 kHz con R = 12 km


def run():
    from .pipeline import report
    cow, full, conv, scan = compute()
    sols = {maR: axion_background(maR) for maR in BG_MAR}
    plot_modes(cow, full, scan, RES / 'paso8_modos.png')
    plot_background(sols, RES / 'paso8_fondo_axion.png')
    plot_meff(RES / 'meff_theta.png')
    out = report('Paso 8: frecuencias de un politropo (sin axión) frente a la literatura', checks(cow, full, conv, scan))
    out += report('Paso 8: fondo con axión', background_checks(sols))
    return out


if __name__ == '__main__':
    run()
