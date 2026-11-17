"""N1, part `atom` (old step 6): cifras y tablas, y cálculo independiente de los niveles de 'átomo gravitatorio' (Tabla IV).


Moved from axion_ns/step6_nonradial_numerics.py. Differences: the tables T1-T4 come from R's data file via
init() instead of being typed here, the figure data go to outputs/levels.json (4-digit rounding kept as in figdata), and the
levels are also written as a long-format table.

Los niveles se calculan resolviendo la Ec. (13) en Schwarzschild sin axión de fondo (m_eff = m_a), con
la condición saliente/decreciente de las soluciones de Whittaker (A13)-(A16) en r~_inf, integrando hacia
dentro hasta R. El interior solo entra como condición en R; se acota con dos extremos:
  'plano'     dtheta = const dentro (dtheta_ext ~ r^l, L = l/R): masa efectiva interior pequeña
  'Dirichlet' dtheta_ext(R) = 0: masa efectiva interior muy grande
y con R entre 10.5 km (Fig. 5 del paper) y 12 km.
"""
import json
from pathlib import Path

import matplotlib
import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

OUT = Path(__file__).resolve().parent / 'outputs'      # overridden by init()
GMSUN_KM = 1.4766250614       # G M_sun / c^2 en km
C_KMS = 299792.458
HBARC_EVM = 1.973269804e-7    # hbar c en eV m
H_EVS = 4.135667696e-15       # h en eV s
L_ = 2

T1 = T2 = T3 = T4 = None      # set by init() from R's paper2_tables.ref.json


def alpha_of(nua, M):
    return GMSUN_KM * M * 2 * np.pi * nua * 1e3 / C_KMS


def _rnd(x, d=4):
    """Round numbers, arrays and nested lists to d significant digits (as axion_ns/figdata.r)."""
    if isinstance(x, dict):
        return {k: _rnd(v, d) for k, v in x.items()}
    if isinstance(x, (list, tuple, np.ndarray)):
        return [_rnd(v, d) for v in x]
    if isinstance(x, (int, np.integer)) or x is None or isinstance(x, str):
        return x if not isinstance(x, np.integer) else int(x)
    x = float(x)
    return 0.0 if x == 0 else float(f'{x:.{d}g}')


def init(tables, out):
    """Set the Paper II tables (T1..T4) and the output directory."""
    global T1, T2, T3, T4, OUT
    T1, T2, T3, T4 = tables
    OUT = Path(out)
    OUT.mkdir(parents=True, exist_ok=True)


def levels_columns(rows):
    """Long-format table of the bound levels: one row per (nu_a, n, interior variant)."""
    cols = {k: [] for k in ('nu_a_kHz', 'n', 'variant', 'R_km', 'nu_n_kHz', 'ratio_k2n2_over_alpha2', 'q')}
    for nua in sorted(rows):
        for n in (3, 4):
            for (bc, R), (nu_n, ratio, q) in zip(VARIANTS, rows[nua][n]):
                for k, v in zip(cols, (nua, n, bc, R, nu_n, ratio, q)):
                    cols[k].append(v)
    return {k: np.array(v) for k, v in cols.items()}


# ---------------------------------------------------------------- niveles ligados
def _seed(q, al, x):
    """L_ext = d ln dtheta_ext / dr~ en r~ = x para la solución de Whittaker decreciente, k~ = i q."""
    kap = al * (1 - 2 * q * q) / q
    mu = mp.sqrt((L_ + 0.5)**2 - 4 * al**2 * (2 - 3 * q * q))
    z = 2 * q * x
    W0, W1 = mp.whitw(kap, mu, z), mp.whitw(kap + 1, mu, z)
    Lp = float(2 * q * ((0.5 - kap / z) - W1 / (z * W0)))   # W' = (1/2 - kappa/z) W - W_{kappa+1}/z
    return Lp - 1 / x - al / (x * x * (1 - 2 * al / x))   # Ec. (A10)


def _rhs(q, al):
    w2 = 1 - q * q

    def f(x, y):
        u, up = y
        F = 1 - 2 * al / x
        return [up, -(2 * al / x**2 / F + 2 / x) * up - (w2 / F - 1 - L_ * (L_ + 1) / x**2) / F * u]
    return f


def _integrate(q, al, xR, dense=False):
    xinf = xR + 30 / q
    s = solve_ivp(_rhs(q, al), [xinf, xR], [1.0, _seed(q, al, xinf)], rtol=1e-10, atol=1e-14,
                  method='DOP853', dense_output=dense)
    if s.status < 0:                                    # H32
        raise RuntimeError(f'paso 6: {s.message}')
    return s


def _mismatch(q, al, xR, bc):
    u, up = _integrate(q, al, xR).y[:, -1]
    if bc == 'Dirichlet':
        return u / np.hypot(u, up)
    return (up - L_ / xR * u) / np.hypot(u, up)


def level(nua, M, R, n, bc):
    """Devuelve (nu_n en kHz, -k~^2 n^2/alpha^2, q)."""
    m = 2 * np.pi * nua * 1e3 / C_KMS
    al, xR = alpha_of(nua, M), m * R
    lo, hi = al / n * 0.95, al / n * 1.3
    q = brentq(_mismatch, lo, hi, args=(al, xR, bc), xtol=1e-13)
    return nua * np.sqrt(1 - q * q), q * q * n * n / al**2, q


VARIANTS = [('plano', 10.5), ('plano', 12.0), ('Dirichlet', 10.5), ('Dirichlet', 12.0)]


def levels_table():
    rows = {}
    for nua, (M, *_rest) in T4.items():
        rows[nua] = {n: [level(nua, M, R, n, bc) for bc, R in VARIANTS] for n in (3, 4)}
    return rows


# ---------------------------------------------------------------- comprobaciones
def numbers_checks():
    res = []
    lam_km = HBARC_EVM / 1e-11 / 1e3
    nu = 1e-11 / H_EVS / 1e3
    res.append((f'lambda_a = hbar/(m_a c) = {lam_km:.1f} km para m_a = 1e-11 eV (paper: ~20 km)',
                abs(lam_km - 20) < 0.5, ''))
    res.append((f'nu_a = m_a/2pi = {nu:.2f} kHz para m_a = 1e-11 eV (paper: ~2.4 kHz)', abs(nu - 2.4) < 0.05, ''))

    sh = [(T1[15][i] - T1[None][i]) / T1[None][i] * 100 for i in range(3)]
    res.append((f'Tabla I: desplazamientos f, p1, p2 a 15 kHz = {sh[0]:.1f}%, {sh[1]:.1f}%, {sh[2]:.1f}% '
                '(texto: 18%, 15%, 19%; resumen: hasta ~20%)',
                all(abs(a - b) < 0.6 for a, b in zip(sh, (18, 15, 19))) and max(sh) < 20.5, ''))
    txt = [(T1[None][0], T1[15][0], T2[15][0]), (T1[None][1], T1[15][1], T2[15][2]),
           (T1[None][2], T1[15][2], T2[15][4])]
    ok = txt == [(2.43, 2.87, 2.50), (7.15, 8.19, 7.96), (10.18, 12.11, 12.00)] and \
        max(T2.values(), key=lambda v: v[0])[0] == 2.56 and T2[5][0] == 2.56
    comp = (T1[15][0] - T2[15][0]) / (T1[15][0] - T1[None][0]) * 100
    res.append(('Texto vs tablas: f 2.43->2.87->2.50, p1 7.15->8.19->7.96, p2 10.18->12.11->12.00, '
                f'máximo de f = 2.56 kHz en 5-6 kHz', ok, f'el acoplo dinámico compensa el {comp:.0f}% '
                'del desplazamiento de f ("almost entirely")'))

    bad = []
    for nua, (fr, fi, p1r, p1i, p2r, p2i) in T2.items():
        for name, re_, im_ in (('f', fr, fi), ('p1', p1r, p1i), ('p2', p2r, p2i)):
            if (im_ != 0) != (re_ > nua):
                bad.append(f'{name}@{nua}')
    res.append(('Tabla II: Im nu != 0 exactamente cuando Re nu > nu_a (36 entradas)', not bad, ', '.join(bad)))
    imax = max(abs(v[i]) for v in T2.values() for i in (3, 5))
    tau = 1 / (2 * np.pi * imax) * 1e3
    res.append((f'tau_a mínimo = 1/(2 pi x {imax:g} Hz) = {tau:.0f} ms (texto: ~10 ms)', 8 < tau < 15, ''))
    ast = [nua for nua, v in T3.items() if v[4] != (v[0] < nua)]
    res.append(('Tabla III: el asterisco marca exactamente Re nu < nu_a; Fig. 4 = fila de 4 kHz '
                '(10.2-3.8i, 22.4-5.7i)', not ast and T3[4][:2] == (10.2, -3.8) and
                round(T3[4][2], 1) == 22.4 and round(T3[4][3], 1) == -5.7, ''))

    da = max(abs(alpha_of(nua, v[0]) - v[1]) for nua, v in T4.items())
    res.append((f'Tabla IV: alpha = G M m_a reproduce la columna alpha (máx. dif. {da:.4f})', da < 0.0015, ''))
    dn = max(abs(nua * np.sqrt(1 - v[2 + 2 * j + 1] * v[1]**2 / n**2) - v[2 + 2 * j])
             for nua, v in T4.items() for j, n in enumerate((3, 4)))
    res.append((f'Tabla IV: nu_n coherente con la columna -k~^2 n^2/alpha^2 (máx. dif. {dn:.3f} kHz)', dn < 0.006, ''))
    return res


def levels_checks():
    rows = levels_table()
    res, worst_low, outside = [], 0.0, []
    for nua, by_n in rows.items():
        for j, n in enumerate((3, 4)):
            nus = [v[0] for v in by_n[n]]
            paper = T4[nua][2 + 2 * j]
            if nua <= 10:
                worst_low = max(worst_low, max(abs(x - paper) for x in nus))
            if not (min(nus) - 0.006 <= paper <= max(nus) + 0.006):
                outside.append(f'n={n}@{nua}')
    res.append((f'Tabla IV reproducida desde el exterior de Schwarzschild: para nu_a <= 10 kHz, '
                f'máx. dif. {worst_low * 1e3:.0f} Hz (precisión de la tabla: 10 Hz)', worst_low < 0.006,
                'independiente del interior (plano/Dirichlet, R = 10.5-12 km)'))
    spread = max(max(v[0] for v in rows[15][n]) - min(v[0] for v in rows[15][n]) for n in (3, 4))
    res.append(('Tabla IV, nu_a = 11-15 kHz: valores del paper dentro de la banda debida al interior',
                not outside, f'la banda crece con m_a R: {spread * 1e3:.0f} Hz a 15 kHz; ' +
                (', '.join(outside) if outside else 'todas dentro')))

    # cruces evitados (Sec. III E): niveles n=3 a 12.19 kHz y n=4 a 12.06 kHz frente a p2 = 11.93 kHz
    Mi = lambda nua: np.interp(nua, [12, 13], [T4[12][0], T4[13][0]])
    c3 = [level(12.19, Mi(12.19), R, 3, bc)[0] for bc, R in VARIANTS]
    c4 = [level(12.06, Mi(12.06), R, 4, bc)[0] for bc, R in VARIANTS]
    ok = all(min(c) - 0.01 <= 11.93 <= max(c) + 0.01 for c in (c3, c4))
    res.append((f'Cruces evitados: nivel n=3 a 12.19 kHz = {min(c3):.2f}-{max(c3):.2f} kHz y n=4 a 12.06 kHz = '
                f'{min(c4):.2f}-{max(c4):.2f} kHz, frente a p2 = 11.93 kHz', ok, ''))
    return res, rows


# ---------------------------------------------------------------- figura
def figure(rows):
    css = {'ink': '#1f2328', 'muted': '#6b7280', 'n3': '#c2410c', 'n4': '#1d4ed8'}
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    nus = np.array(sorted(rows))
    data = {'nu': nus.tolist(), 'levels': {}, 'r': [], 'dV': [], 'psi': {}, 'R': 10.5}
    for j, n in enumerate((3, 4)):
        col = css[f'n{n}']
        vals = np.array([[v[1] for v in rows[nu][n]] for nu in nus])
        ax1.fill_between(nus, vals.min(1), vals.max(1), color=col, alpha=0.18, lw=0)
        data['levels'][n] = {'lo': vals.min(1), 'hi': vals.max(1), 'mid': vals[:, 0],
                             'table': [T4[nu][3 + 2 * j] for nu in nus]}
        ax1.plot(nus, vals[:, 0], color=col, lw=2, label=f'n = {n}: este cálculo (banda: interior)')
        ax1.plot(nus, [T4[nu][3 + 2 * j] for nu in nus], 'o', color=col, mfc='white', mew=2, ms=8,
                 label=f'n = {n}: Tabla IV del paper')
    ax1.axhline(1, color=css['muted'], lw=1, ls='--')
    ax1.text(8.1, 1.008, 'hidrogénico, Ec. (23)', color=css['muted'], fontsize=9)
    ax1.set_xlabel(r'$\nu_a$ (kHz)')
    ax1.set_ylabel(r'$-\tilde k_n^2\, n^2/\alpha^2$')
    ax1.set_title('Niveles ligados l = 2 (Tabla IV)', fontsize=11)
    ax1.legend(fontsize=8, frameon=False, loc='upper left')

    # eigenfunciones a 12 kHz, como la Fig. 5 del paper (exterior; interior 'plano', R = 10.5 km)
    nua, R = 12.0, 10.5
    M = T4[nua][0]
    m = 2 * np.pi * nua * 1e3 / C_KMS
    al = alpha_of(nua, M)
    r = np.linspace(R, 190, 800)
    x = m * r
    F = 1 - 2 * al / x
    dV = -2 * al / x + F * (L_ * (L_ + 1) / x**2 + 2 * al / x**3)
    ax2.plot(r, dV, color=css['ink'], lw=1.5, label=r'$\Delta V/m_a^2$, Ec. (22)')
    for n in (3, 4):
        _, _, q = level(nua, M, R, n, 'plano')
        sol = _integrate(q, al, m * R, dense=True)
        psi = x * sol.sol(x)[0]
        psi = psi / np.max(np.abs(psi)) * 0.022
        if psi[np.searchsorted(r, 35.0)] < 0:   # mismo signo que la Fig. 5 del paper
            psi = -psi
        k2 = -q * q
        ax2.axhline(k2, color=css[f'n{n}'], lw=0.8)
        ax2.plot(r, k2 + psi, color=css[f'n{n}'], lw=2, label=fr'$\psi$, n = {n}')
        data['psi'][n] = {'k2': k2, 'y': (k2 + psi)[::4]}
    ax2.axvline(R, color=css['muted'], ls=':', lw=1)
    ax2.set_ylim(-0.12, 0.06)
    ax2.set_xlim(0, 190)
    ax2.set_xlabel('r (km)')
    ax2.set_title(r'Eigenfunciones a $\nu_a$ = 12 kHz (cf. Fig. 5)', fontsize=11)
    ax2.legend(fontsize=8, frameon=False, loc='lower right')
    for ax in (ax1, ax2):
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(alpha=0.25, lw=0.6)
    fig.tight_layout()
    path = OUT / 'paper2_atom_levels.png'
    fig.savefig(path, dpi=150)
    data['r'], data['dV'] = r[::4], dV[::4]
    (OUT / 'levels.json').write_text(json.dumps(_rnd(data), ensure_ascii=False, separators=(',', ':')))
    return path


def run():
    res = numbers_checks()
    lv, rows = levels_checks()
    res += lv
    path = figure(rows)
    res.append((f'Figura guardada en results/{path.name}', path.exists(), ''))
    return res

