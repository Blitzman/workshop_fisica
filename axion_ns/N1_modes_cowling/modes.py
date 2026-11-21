"""N1: Cowling modes of the star with an axion condensate (Paper II, Eqs. 8a-8d, centre 10, surface, exterior 13).

Moved from axion_ns/axstar.py (the mode half: perturbation RHS, interior, surface residual, exterior, coupled determinant,
Muller). The star is a StarTable (common/star_table.py) built from S1's tables instead of an in-memory Star. The ODEs are
typed by hand (verified symbol by symbol in D2 for the equations; the numerical implementation is not yet compiled
from the artifact: see README, "Known limitations").
"""
import numpy as np
from scipy import constants as C
from scipy.integrate import solve_ivp

CKMS = C.c / 1e3            # km/s (same value as axion_ns/axstar.py; also recorded in the S1 table metadata)


def kHz(omega):
    """omega in km^-1 -> nu in kHz."""
    return omega * CKMS / (2 * np.pi) / 1e3


def omega_of(nu_khz):
    return 2 * np.pi * nu_khz * 1e3 / CKMS


# ---------------------------------------------------------------- modos en Cowling
def _perturb_rhs(st, l, omega, dyn):
    """(8a)-(8d) con y = [dp, U = (eps+p) W, dtheta, dtheta']; dyn=False apaga delta theta."""
    def f(r, y):
        e, p, nb, cs2, q1, q2, dPhi, dlam, Phi, lam, t, tp = st.bg(r)
        dp, U, dt, dtp = y
        dnn = dp / (cs2 * (e + p))                        # delta n / n
        de = dp / cs2 + q1 * dt
        ddp = (-(dPhi + l / r) * dp + omega**2 * np.exp(lam - 2 * Phi) / r * U - de * dPhi
               - q1 * (dtp + l / r * dt) - tp * (dnn * q1 + q2 * dt))
        EV = -np.exp(2 * Phi) / omega**2 * (dp + q1 * dt)
        dU = -U * (dPhi + (l + 1) / r) - r * np.exp(lam) * (dp / cs2 + EV * l * (l + 1) / r**2)
        if not dyn:
            return [ddp, dU, 0 * dt, 0 * dtp]
        meff2 = (q2 - st.E0 * st.g2(t)) / st.fa2
        el = np.exp(-2 * lam)
        dtpp = (-el * (dPhi - dlam + 2 * (l + 1) / r) * dtp
                - (omega**2 * np.exp(-2 * Phi) - meff2 + l / r * el * (dPhi - dlam)
                   + l * (l + 1) / r**2 * (el - 1)) * dt + dnn * q1 / st.fa2) / el
        return [ddp, dU, dtp, dtpp]
    return f


def _ok(sol):
    """H29: si solve_ivp no llega al final, se lanza (antes se tomaba s.y[:, -1] sin mirar status)."""
    if sol.status < 0:
        raise RuntimeError(f'solve_ivp: {sol.message}')
    return sol


def interior(st, l, omega, U0, t0, dyn=True, r0=1e-3, frac=1.0, rtol=1e-9):
    """Integra desde el centro con las condiciones (10); devuelve el estado en la superficie."""
    e, p, nb, cs2, q1, q2, dPhi, dlam, Phi, lam, t, tp = st.bg(r0)
    dp0 = np.exp(-2 * Phi) * omega**2 / l * U0 - q1 * t0
    y0 = np.array([dp0, U0, t0, 0.0], dtype=complex)
    Rs = st.R * frac
    f = _perturb_rhs(st, l, omega, dyn)
    # H29: aquí había BDF «por rigidez» hasta r1, pero BDF fallaba en su primer paso y s.y[:, -1] era el valor
    # inicial en r0, usado como si fuera el de r1. El tramo no es rígido: lo que ahogaba al integrador era
    # atol = 1e-40 sobre componentes que arrancan en 0 (con la semilla axiónica, DOP853 necesitaba ~1e5 pasos).
    # Con atol = 1e-18, absoluto sobre semillas O(1), el estado en r1 coincide con el de 1e-40 a 1e-10 relativo.
    r1 = 0.5
    s = _ok(solve_ivp(f, (r0, r1), y0, method='DOP853', rtol=rtol, atol=1e-18))
    s = _ok(solve_ivp(f, (r1, Rs), s.y[:, -1], method='DOP853', rtol=rtol, atol=1e-40))
    return s.y[:, -1], Rs


def delta_p(st, Rs, y):
    """Perturbación lagrangiana de la presión en Rs: dp + p' xi^r."""
    e, p, nb, cs2, q1, q2, dPhi, dlam, Phi, lam, t, tp = st.bg(Rs)
    return y[0] - np.exp(-lam) * y[1] / Rs * (dPhi + q1 * tp / (e + p))


def surface_residual(st, l, omega, dyn=False):
    """|Delta p(R)| a frecuencia dada, sin delta theta (Tabla I y referencia)."""
    y, Rs = interior(st, l, omega, 1.0, 0.0, dyn=False)
    return delta_p(st, Rs, y).real


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


def coupled_det(st, l, omega):
    """Determinante de las dos condiciones de superficie (Delta p = 0 y empalme de delta theta)
    para las dos soluciones regulares en el centro."""
    L = exterior_L(st, l, omega) - l / st.R
    cols = []
    for U0, t0 in ((1.0, 0.0), (0.0, 1.0)):
        y, Rs = interior(st, l, omega, U0, t0, dyn=True)
        cols.append((delta_p(st, Rs, y), y[3] - L * y[2]))
    (a, b), (c, d) = cols
    return a * d - b * c


def muller(f, z0, z1, z2, tol=1e-10, maxit=40):
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
            return z2
    raise RuntimeError('Muller no converge')
