"""N4: radial (l = 0) modes of the star with the axion condensate in full GR, table-based (old axion_ns/radial.py).

`kHz`, `omega_of`, `_ok` and `exterior_L` are N1's (extracted verbatim: modules never import each other; `exterior_L` is the real-line exterior used by step 15 to show
that it fails for very damped modes); `winding`, `winding_converged` and `muller` are the legacy functions verbatim (argument principle with adaptive sampling; Muller with a
normalised-residual check). The rest is the port of the legacy `interior`, `outward`, `det`, `det_entire`, `exterior_psi` and `fluid_fraction` with the right-hand sides taken from
the coefficient tables (tables.py). Differences with the legacy: the background of the exterior wave uses -beta/4 for g''(0) (exact) instead of `axstar.g2(0.0)`;
the matching radius is `R + 0.7 (r_end - R)` (see matching_radius) instead of the first |theta| < 1e-7.

det(omega) = determinant of Delta p(R) = 0 and the matching with the exterior wave (outgoing, sheet = 1, or bound) at r_m; det_entire multiplies by the exterior
solution normalised at infinity (analytic, no poles) for the argument principle in Im omega > 0.
"""
import numpy as np
from scipy import constants as C
from scipy.integrate import solve_ivp

from .background import bgvals

CKMS = C.c / 1e3            # km/s (same value as axion_ns/axstar.py)


def kHz(omega):
    """omega in km^-1 -> nu in kHz."""
    return omega * CKMS / (2 * np.pi) / 1e3


def omega_of(nu_khz):
    return 2 * np.pi * nu_khz * 1e3 / CKMS


# ---------------------------------------------------------------- modos en Cowling


def _ok(sol):
    """H29: si solve_ivp no llega al final, se lanza (antes se tomaba s.y[:, -1] sin mirar status)."""
    if sol.status < 0:
        raise RuntimeError(f'solve_ivp: {sol.message}')
    return sol


def exterior_L(st, l, omega, r_stop=None, full=False, branch='auto', n_decay=40):
    """psi'/psi en R para psi = delta theta_ext, Ec. (13), con la condición de onda saliente
    (o ligada si Im k > 0) en r grande: psi ~ exp[i(k r + eta ln r)]/r, eta = M (2 w^2 - m_a^2)/k.
    branch: 'auto' (la física, H31), 'principal' (la de antes) o 'flipped' (la opuesta a la principal: controles).
    n_decay: el arranque está a n_decay/min(|k|, m_a) del radio final."""
    m2 = st.ma**2
    k = np.sqrt(complex(omega**2 - m2))
    # H31: por debajo del umbral el canal axiónico está cerrado y la solución debe decaer (Im k > 0). La rama principal
    # da Im k < 0 si el modo está amortiguado; antes solo se corregía con Re k ~ 0. Por encima del umbral, la principal
    # (Re k > 0) es la saliente. Con la rama anterior L salía igual (< 1e-13): la integración hacia dentro recorre
    # ~40 longitudes de decaimiento y filtra la solución creciente (paso 25)
    if branch == 'flipped':
        k = -k
    elif branch == 'auto' and (omega**2 - m2).real < 0 and k.imag < 0:
        k = -k
    Mt, R, rend = st.Mt, st.R, st.rend
    r_stop = r_stop or R
    rinf = min(r_stop + n_decay / min(abs(k), st.ma), 6000.0)
    eta = Mt * (2 * omega**2 - m2) / k

    def bgx(r):
        if r < rend:
            m, p, Phi, t, tp = st.y_at(r)
            dm = st.rhs(r, [m, p, Phi, t, tp])[0]
            dPhi = st.rhs(r, [m, p, Phi, t, tp])[2]
            Phi = Phi + st.Phi_shift
            meff2 = -st.E0 * st.g2(t) / st.fa2
        else:
            m, dm, t = Mt, 0.0, 0.0
            dPhi = m / (r * (r - 2 * m))
            Phi = 0.5 * np.log(1 - 2 * m / r)
            meff2 = m2
        dlam = (dm * r - m) / (r * (r - 2 * m))
        return dPhi, dlam, Phi, -0.5 * np.log(1 - 2 * m / r), meff2

    def f(r, y):
        dPhi, dlam, Phi, lam, meff2 = bgx(r)
        return [y[1], -(dPhi - dlam + 2 / r) * y[1]
                + np.exp(2 * lam) * (meff2 + l * (l + 1) / r**2 - omega**2 * np.exp(-2 * Phi)) * y[0]]

    psi = np.exp(1j * (k * rinf + eta * np.log(rinf))) / rinf
    dpsi = psi * (1j * (k + eta / rinf) - 1 / rinf)
    s = solve_ivp(f, (rinf, r_stop), [psi / abs(psi), dpsi / abs(psi)], method='DOP853', rtol=1e-10, atol=1e-60)
    if s.status < 0:                                    # H32: no se usa una integración que no terminó
        raise RuntimeError(f'exterior_L: {s.message}')
    if full:                      # (psi, psi') en r_stop, sin dividir: evita los polos de psi'/psi
        return s.y[0, -1], s.y[1, -1]
    return s.y[1, -1] / s.y[0, -1]


def winding(f, corners, max_step=np.pi / 3, n0=48, max_pts=6000, pool=None, cache=None):
    """Índice de f alrededor del polígono `corners` (sentido antihorario), con muestreo adaptativo: se
    subdivide todo tramo en que la fase salta más de max_step. Devuelve (índice, nº de puntos, |f| mínimo,
    salto de fase máximo). `cache` (dict) reutiliza evaluaciones entre llamadas.
    Ojo: el refinado no detecta saltos de una vuelta entera entre dos muestras (aliasing); para eso está
    winding_converged."""
    cache = {} if cache is None else cache
    _map = (pool.map if pool else map)

    def ev(f_, zs):
        todo = [z for z in zs if z not in cache]
        for z, v in zip(todo, _map(f_, todo)):
            cache[z] = v
        return [cache[z] for z in zs]
    pts = []
    for z0, z1 in zip(corners, corners[1:] + corners[:1]):
        pts += [complex(z) for z in z0 + (z1 - z0) * np.arange(n0) / n0]
    vals = list(ev(f, pts))
    while True:
        ph = [np.angle(vals[(i + 1) % len(vals)] / vals[i]) for i in range(len(vals))]
        bad = [i for i, d in enumerate(ph) if abs(d) > max_step]
        if not bad or len(pts) > max_pts:
            break
        new = [(i, complex(0.5 * (pts[i] + pts[(i + 1) % len(pts)]))) for i in bad]
        nv = list(ev(f, [z for _, z in new]))
        for (i, z), v in sorted(zip(new, nv), key=lambda t: -t[0][0]):
            pts.insert(i + 1, z)
            vals.insert(i + 1, v)
    ph = [np.angle(vals[(i + 1) % len(vals)] / vals[i]) for i in range(len(vals))]
    return float(sum(ph) / (2 * np.pi)), len(pts), float(min(abs(v) for v in vals)), max(abs(d) for d in ph)


def winding_converged(f, corners, n0=48, n_max=384, pool=None):
    """winding con n0, 2 n0, 4 n0...: se acepta el índice cuando dos muestreos seguidos coinciden (los puntos
    del grueso se reutilizan). Devuelve (índice, n0 final, puntos, salto máximo, lista de índices probados)."""
    cache, tried = {}, []
    n = n0
    while True:
        w = winding(f, corners, n0=n, pool=pool, cache=cache)
        tried.append(round(w[0]))
        if len(tried) >= 2 and tried[-1] == tried[-2] and w[3] <= np.pi / 3:   # H25: sin saltos sin resolver
            return tried[-1], n, len(cache), w[3], tried
        if n >= n_max:
            return None, n, len(cache), w[3], tried
        n *= 2


def muller(f, z, tol=1e-9, maxit=30, max_resid=0.1):
    """Raíz compleja de f por Muller desde z. Lanza RuntimeError si no converge en maxit pasos o si el residuo
    normalizado |f(z)|/|f(z (1 + 1e-6))| supera max_resid. Con un cero simple ese cociente es el paso de Newton
    que falta en unidades de 1e-6 |z|: max_resid = 0,1 exige un error < 1e-7 |z| (< 1 mHz a 8 kHz). Separa una raíz
    de un punto que no lo es (cociente ~1); no es la cota de precisión. Con 1e-2 fallaba en el empalme lejano del
    paso 15, donde el ruido de det limita el cociente a ~1e-2 (H25). Antes devolvía la última aproximación sin
    comprobar nada."""
    z0, z1, z2 = z, z + 0.01, z - 0.002j
    f0, f1, f2 = f(z0), f(z1), f(z2)
    for _ in range(maxit):
        h1, h2 = z1 - z0, z2 - z1
        d1, d2 = (f1 - f0) / h1, (f2 - f1) / h2
        a = (d2 - d1) / (h2 + h1)
        b = a * h2 + d2
        D = np.sqrt(b * b - 4 * f2 * a)
        den = b + D if abs(b + D) > abs(b - D) else b - D
        dz = -2 * f2 / den
        z0, z1, z2 = z1, z2, z2 + dz
        f0, f1, f2 = f1, f2, f(z2)
        if abs(dz) < tol * abs(z2):
            break
    else:
        raise RuntimeError(f'muller: sin converger en {maxit} pasos desde {z}')
    ref = abs(f(z2 * (1 + 1e-6)))
    if not np.isfinite(z2) or abs(f2) > max_resid * ref:
        raise RuntimeError(f'muller: residuo normalizado {abs(f2) / ref:.1e} en {z2}')
    return z2


def matching_radius(st):
    """70% of the valid axion tail (the tail of theta_0 is only trustworthy up to r_end, where the integration leaves the stable branch; see N3)."""
    return st.R + 0.7 * (st.rend - st.R)


def _mat(flat, n, rr, omega, B):
    v = flat(rr, omega, *B)
    return np.array([complex(x) for x in v], dtype=complex).reshape(n, n)


def exterior_psi(flat_e, omega, M, r_m, ma, fa_v, E0, beta, n_decay=40.0, rtol=1e-11, sheet=1, dphi=0.0, raw=False):
    """(psi, psi') of the axion at r_m with the outgoing (or bound) condition, integrating the derived exterior equation (exterior matrix with theta_0 = 0: Schwarzschild +
    massive Klein-Gordon = Eq. 3.7 of Paper I) along the complex ray r = r_m + s e^{i phi}, phi = pi/2 - arg(k), where the outgoing wave decays as e^{-|k| s}. On the real line
    (axstar.exterior_L) the integration inwards is unstable for very damped modes: the outgoing wave grows with r and any initial error mixes in the amplified ingoing one."""
    k = np.sqrt(complex(omega**2 - ma**2))
    if k.imag < 0 and abs(k.real) < 1e-14:
        k = -k
    if sheet == 0:                     # upper half plane: the branch with Im k > 0 (decays), also for Re omega < 0
        k = k if k.imag > 0 else -k
    else:
        k = sheet * k                  # sheet = -1: ingoing wave (negative control)
    phi = np.pi / 2 - np.angle(k) + dphi    # dphi != 0: test of the independence of the contour
    e = np.exp(1j * phi)
    S = n_decay / (abs(k) * np.cos(dphi))
    eta = M * (2 * omega**2 - ma**2) / k
    G2 = -beta / 4                     # g''(0) = -beta/4

    def bg(rr):
        f = 1 - 2 * M / rr
        P1 = M / (rr * (rr - 2 * M))
        P2 = -2 * (rr - M) * M / (rr * (rr - 2 * M))**2
        return (0.5 * np.log(f), P1, P2, -0.5 * np.log(f), -P1, 0, 0, 0.0, 0.0, 0, 0, 1.0, fa_v, 0.0, 0.0, E0 * G2, 0.0, -E0 * G2)

    rinf = r_m + S * e
    psi = np.exp(1j * (k * rinf + eta * np.log(rinf))) / rinf
    dpsi = psi * (1j * (k + eta / rinf) - 1 / rinf)
    sc = psi
    rhs = lambda s, y: e * (_mat(flat_e, 2, r_m + s * e, omega, bg(r_m + s * e)) @ y)
    sol = _ok(solve_ivp(rhs, (S, 0.0), [psi / sc, dpsi / sc], method='DOP853', rtol=rtol, atol=1e-60))
    if raw == 'inf':   # only the normalisation at infinity (negative control of step 16: the contour loses turns)
        return sol.y[0, -1] * sc, sol.y[1, -1] * sc
    if raw:     # normalised to the asymptotic form e^{i(kr + eta ln r)}/r: independent of the path, analytic in omega. Also divided by that same form at r_m (analytic,
        # no zeros in Im omega > 0): the number of zeros does not change and the fast phase e^{i k r_m}, which the contour sampling could not resolve, disappears.
        f_rm = np.exp(1j * (k * (rinf - r_m) + eta * np.log(rinf / r_m))) * r_m / rinf
        return sol.y[0, -1] * f_rm, sol.y[1, -1] * f_rm
    return sol.y[0, -1], sol.y[1, -1]


def interior(st, tabs, omega, seed, r0=1e-3, r1=0.5, rtol=1e-10, r_eval=None):
    """Regular solution from the centre to R. seed = 'fluid' (xi ~ r) or 'axion' (dth(0) = 1). The conditions at r0 come from the xi' row itself (regularity)."""
    M = tabs.mat_int(r0, omega)
    xi, th, dxi = (r0, 0.0, 1.0) if seed == 'fluid' else (0.0, 1.0, 0.0)
    y0 = np.array([xi, (dxi - M[0, 0] * xi - M[0, 2] * th) / M[0, 1], th, 0.0], dtype=complex)
    f = tabs.rhs_int(omega)
    s = _ok(solve_ivp(f, (r0, r1), y0, method='DOP853', rtol=rtol, atol=1e-40))
    return _ok(solve_ivp(f, (r1, st.R), s.y[:, -1], method='DOP853', rtol=rtol, atol=1e-40, t_eval=r_eval))


def outward(st, tabs, omega, y_R, rm, rtol=1e-10, r_eval=None):
    """(dth, dth') from R to r_m with the tail of the condensate (exterior matrix)."""
    return _ok(solve_ivp(tabs.rhs_ext(omega), (st.R, rm), y_R, method='DOP853', rtol=rtol, atol=1e-40, t_eval=r_eval))


def det(st, tabs, omega, rm=None, sheet=1, rtol=1e-10, r0=1e-3, n_decay=40.0):
    """Determinant of Delta p(R) = 0 and the matching with the exterior wave (outgoing or bound) at r_m. sheet = -1 is the 'ingoing wave' mutation (negative control)."""
    rm = rm or matching_radius(st)
    cols = []
    for seed in ('fluid', 'axion'):
        s = interior(st, tabs, omega, seed, r0=r0, rtol=rtol)
        cols.append((s.y[1, -1], outward(st, tabs, omega, s.y[2:, -1], rm, rtol=rtol).y[:, -1]))
    ps, dps = exterior_psi(tabs.flat_e, omega, st.Mt, rm, st.ma, np.sqrt(st.fa2), st.E0, st.beta, sheet=sheet, n_decay=n_decay, rtol=min(rtol, 1e-11))
    L = dps / ps
    (a, ya), (c, yc) = cols
    return a * (yc[1] - L * yc[0]) - c * (ya[1] - L * ya[0])


def det_entire(st, tabs, omega, rm=None, rtol=1e-10, n_decay=40.0, inward=True, sheet=0, raw=True, r0=1e-3):
    """Pole-free version of det for the argument principle in Im omega > 0: multiplied by Psi, normalised at infinity (raw = True), branch Im k > 0. Same zeros as det off the real axis.
    inward = True: Psi (which decays outwards) is integrated from r_m to R and matched at the surface; this is what is stable in Im omega > 0 (inwards Psi is the dominant solution).
    Carrying the seeds from R to r_m (inward = False, as det does) mixes the growing and the decaying solutions and loses the phase when Im omega is large."""
    rm = rm or matching_radius(st)
    ps, dps = exterior_psi(tabs.flat_e, omega, st.Mt, rm, st.ma, np.sqrt(st.fa2), st.E0, st.beta, sheet=sheet, n_decay=n_decay, rtol=min(rtol, 1e-11), raw=raw)
    sf, sa = (interior(st, tabs, omega, seed, r0=r0, rtol=rtol).y[:, -1] for seed in ('fluid', 'axion'))
    a, c = sf[1], sa[1]
    y = c * sf - a * sa                                   # Delta p(R) = 0
    if inward:
        u = _ok(solve_ivp(tabs.rhs_ext(omega), (rm, st.R), [ps, dps], method='DOP853', rtol=rtol, atol=1e-300)).y[:, -1]
        return y[3] * u[0] - y[2] * u[1]
    yo = outward(st, tabs, omega, y[2:], rm, rtol=rtol).y[:, -1]
    return yo[1] * ps - yo[0] * dps


def fluid_fraction(st, tabs, omega, rm=None):
    """Fraction of the energy (simple kinetic weights) carried by the fluid in the mode of frequency omega."""
    rm = rm or matching_radius(st)
    ri, re = np.linspace(0.5, st.R, 600), np.linspace(st.R, rm, 1500)
    sols = []
    for seed in ('fluid', 'axion'):
        s = interior(st, tabs, omega, seed, r_eval=ri)
        sols.append((s.y, outward(st, tabs, omega, s.y[2:, -1], rm, r_eval=re).y))
    (yf, ef), (ya, ea) = sols
    a, c = yf[1, -1], ya[1, -1]
    y, e = c * yf - a * ya, c * ef - a * ea
    Bi = np.array([bgvals(st, x) for x in ri])
    Be = np.array([bgvals(st, x) for x in re])
    eps = Bi[:, 9] + Bi[:, 5] * Bi[:, 13]
    Ef = np.trapezoid((eps + Bi[:, 10]) * np.exp(3 * Bi[:, 3] - Bi[:, 0]) * abs(y[0])**2 * ri**2, ri)
    Ea = (np.trapezoid(st.fa2 * np.exp(Bi[:, 3] - Bi[:, 0]) * abs(y[2])**2 * ri**2, ri)
          + np.trapezoid(st.fa2 * np.exp(Be[:, 3] - Be[:, 0]) * abs(e[0])**2 * re**2, re))
    return float(Ef / (Ef + Ea))
