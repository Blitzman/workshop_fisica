"""Estrella de neutrones BSk26 con condensado axiónico y sus modos no radiales en Cowling.

Etapa A del trabajo con axión: reproducir las Tablas I y II del Paper II (arXiv:2610.37308) con las
ecuaciones derivadas desde la acción (pasos 1 y 5), integradas sobre un fondo realista:

  fondo     TOV + Klein-Gordon del Paper I, Ec. (2.6), con eps(n, theta) = eps0(n) - n sigma_N (1 - g),
            p(n) sin cambios (Ec. 2.7-2.8) y eps0, p de BSk26; presión central 100 MeV/fm^3.
  modos     (8a)-(8d) del Paper II (comprobadas símbolo a símbolo en el paso 5), centro (10), superficie
            Delta p(R) = 0 y, para delta theta, empalme con la Ec. (13) fuera (onda saliente o ligada).

Unidades geométricas G = c = 1 con longitudes en km. Todas las escalas del axión salen de dos números:
E0 = eps m_pi^2 f_pi^2 (= n_c sigma_N) y m_a; f_a^2 = beta E0 / (4 m_a^2) y x = n_b / n_c.
"""
import numpy as np
from scipy import constants as C
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from .bsk26 import EOS

KGEO = C.e * 1e6 / 1e-45 * C.G / C.c**4 * 1e6          # MeV/fm^3 -> km^-2
MSUN = C.G * 1.98847e30 / C.c**2 / 1e3                 # km
CKMS = C.c / 1e3
HBARC = C.hbar * C.c / C.e * 1e9                       # MeV fm

BETA = 0.88
SIGMA_N = 50.0                                         # MeV
EPS = 0.1
MPI, FPI = 134.98, 92.1                                # MeV (paso 4: con estos, eps m_pi^2 f_pi^2 = 2.0 MeV/fm^3)
E0 = EPS * (MPI * FPI)**2 / HBARC**3                   # MeV/fm^3
NC = E0 / SIGMA_N                                      # fm^-3
P0 = 100.0                                             # MeV/fm^3, presión central del modelo de referencia

_eos = EOS()


def kHz(omega):
    """omega en km^-1 -> nu en kHz."""
    return omega * CKMS / (2 * np.pi) / 1e3


def omega_of(nu_khz):
    return 2 * np.pi * nu_khz * 1e3 / CKMS


def g(t):
    return np.sqrt(1 - BETA * np.sin(t / 2)**2)


def g1(t):
    return -BETA * np.sin(t) / (4 * g(t))


def g2(t):
    gg = g(t)
    return -BETA * np.cos(t) / (4 * gg) - BETA**2 * np.sin(t)**2 / (16 * gg**3)


class Star:
    """Fondo estático. nu_a = None: sin axión. axion=False con nu_a: mismo fondo, pero en los modos
    se apaga delta theta (Tabla I). Guarda funciones r -> fondo dentro (y fuera, para el axión)."""

    def __init__(self, nu_a=None, p0=P0, e0=E0):
        self.nu_a, self.p0, self.E0 = nu_a, p0, e0 * KGEO
        self.nc = e0 / SIGMA_N
        self.ma = omega_of(nu_a) if nu_a else 0.0
        self.fa2 = BETA * self.E0 / (4 * self.ma**2) if nu_a else np.inf
        if nu_a:
            self._shoot()
        else:
            self.thc = 0.0
            self._integrate(0.0, keep=True)

    # ------------------------------------------------------------ materia
    def matter(self, p, t):
        """(eps, n_b, x, c_s^2 con el eps axiónico, q1 = n dm*/dtheta, q2 = n d2m*/dtheta2) en km^-2."""
        pm = p / KGEO
        if pm <= _eos.p_min * 1.0001:
            return 0.0, 0.0, 0.0, 1.0, 0.0, 0.0
        e0, nb, cs0 = _eos.eps(pm) * KGEO, _eos.nb(pm), _eos.cs2(pm)
        x = nb / self.nc
        e = e0 - self.E0 * x * (1 - g(t)) if self.nu_a else e0
        cs2 = (e0 + p) * cs0 / (e + p)                # dp/dn n / (eps + p), con eps(n, theta)
        q1 = self.E0 * x * g1(t) if self.nu_a else 0.0
        q2 = self.E0 * x * g2(t) if self.nu_a else 0.0
        return e, nb, x, cs2, q1, q2

    def rhs(self, r, y):
        m, p, Phi, t, tp = y
        inside = p / KGEO > _eos.p_min * 1.0001
        e, nb, x, cs2, q1, q2 = self.matter(p, t) if inside else (0, 0, 0, 1, 0, 0)
        pr = p if inside else 0.0
        F = 1 - 2 * m / r
        V = self.E0 * (1 - g(t)) if self.nu_a else 0.0
        kin = self.fa2 * tp**2 / 2 * F if self.nu_a else 0.0
        dPhi = (m + 4 * np.pi * r**3 * (pr - V + kin)) / (r**2 * F)
        dm = 4 * np.pi * r**2 * (e + V + kin)
        dp = (-(e + pr) * dPhi - q1 * tp) if inside else 0.0
        if self.nu_a:
            src = 4 * self.ma**2 / BETA * (x - 1) * g1(t)
            tpp = (src - 2 * tp / r * (1 - m / r - 2 * np.pi * r**2 * (e - pr + 2 * V))) / F
        else:
            tpp = 0.0
        return [dm, dp, dPhi, tp, tpp]

    def _start(self, r0, thc):
        p = self.p0 * KGEO
        e, nb, x, cs2, q1, q2 = self.matter(p, thc)
        V = self.E0 * (1 - g(thc)) if self.nu_a else 0.0
        tpp0 = (4 * self.ma**2 / BETA * (x - 1) * g1(thc)) / 3 if self.nu_a else 0.0
        m = 4 * np.pi / 3 * (e + V) * r0**3
        return [m, p - 2 * np.pi / 3 * (e + p) * (e + 3 * p) * r0**2, 0.0, thc + tpp0 * r0**2 / 2, tpp0 * r0]

    def _integrate(self, thc, keep=False, rmax=None):
        """Integra desde el centro. Devuelve +1 si theta se dispara hacia arriba, -1 si cruza cero."""
        r0 = 1e-5
        rmax = rmax or (60.0 + (40 / self.ma if self.nu_a else 0.0))
        up = lambda r, y: (y[3] - np.pi - 1e-3) if self.nu_a else 1.0
        down = lambda r, y: y[3] + 1e-3 if self.nu_a else 1.0
        turn = lambda r, y: y[4] if (self.nu_a and y[1] / KGEO < _eos.p_min * 1.0001) else -1.0
        for ev in (up, down):
            ev.terminal = True
        turn.terminal, turn.direction = True, 1
        s = solve_ivp(self.rhs, (r0, rmax), self._start(r0, thc), events=[up, down, turn],
                      method='DOP853', rtol=1e-11, atol=1e-30, dense_output=keep, max_step=0.25)
        if keep:
            self.sol = s
        if not self.nu_a:
            return 0
        if s.t_events[1].size:
            return -1
        return +1

    def _shoot(self):
        """theta_c = pi - delta, bisección en ln delta: el límite estable fuera es theta -> 0."""
        lo, hi = np.log(1e-300), np.log(np.pi - 0.3)    # delta pequeño -> se dispara; theta_c = 0.3 -> oscila y cruza
        if self._integrate(np.pi - np.exp(hi)) > 0:
            raise ValueError('sin condensado')
        for _ in range(64):
            mid = 0.5 * (lo + hi)
            if self._integrate(np.pi - np.exp(mid)) > 0:
                lo = mid
            else:
                hi = mid
        dlo, dhi = np.exp(lo), np.exp(hi)
        for _ in range(60):                         # afinado lineal
            mid = 0.5 * (dlo + dhi)
            if self._integrate(np.pi - mid) > 0:
                dlo = mid
            else:
                dhi = mid
            if dhi - dlo < 1e-16 * max(dhi, 1e-300) or mid in (dlo, dhi):
                break
        self.thc = np.pi - dlo
        self._integrate(self.thc, keep=True)

    # ------------------------------------------------------------ fondo como funciones de r
    @property
    def R(self):
        """Radio: donde p llega a la presión mínima de la tabla."""
        if not hasattr(self, '_R'):
            f = lambda r: self.sol.sol(r)[1] / KGEO - _eos.p_min * 1.001
            rr = np.linspace(self.sol.t[0], self.sol.t[-1], 4000)
            pv = np.array([f(x) for x in rr])
            i = np.argmax(pv < 0)
            self._R = brentq(f, rr[i - 1], rr[i], xtol=1e-12) if i else self.sol.t[-1]
        return self._R

    @property
    def M(self):
        return self.sol.sol(self.R)[0]

    def bg(self, r):
        """(eps, p, n, c_s^2, q1, q2, Phi', lambda', Phi, lambda, theta, theta') en r (dentro)."""
        m, p, Phi, t, tp = self.sol.sol(r)
        e, nb, x, cs2, q1, q2 = self.matter(p, t)
        F = 1 - 2 * m / r
        dPhi, dm = self.rhs(r, [m, p, Phi, t, tp])[2], self.rhs(r, [m, p, Phi, t, tp])[0]
        lam = -0.5 * np.log(F)
        dlam = (dm * r - m) / (r * (r - 2 * m))
        return e, p, nb, cs2, q1, q2, dPhi, dlam, Phi + self.Phi_shift, lam, t, tp

    @property
    def Phi_shift(self):
        """Constante que hace Phi -> 0 en el infinito. Se ajusta a Schwarzschild con la masa total al final de la
        integración del fondo, no en R: con condensado, la cola del axión fuera de la estrella lleva ~0,5 % de la
        masa y ajustar en R con m(R) dejaba Phi(inf) = +1,6e-3 (corregido el 2026-11-08, LAB_NOTEBOOK H16).
        Sin axión, m es constante fuera de R y ambos ajustes coinciden."""
        if not hasattr(self, '_Ps'):
            rend = self.sol.t[-1]
            m, p, Phi = self.sol.sol(rend)[:3]
            self._Ps = 0.5 * np.log(1 - 2 * m / rend) - Phi
        return self._Ps


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
        meff2 = (q2 - st.E0 * g2(t)) / st.fa2
        el = np.exp(-2 * lam)
        dtpp = (-el * (dPhi - dlam + 2 * (l + 1) / r) * dtp
                - (omega**2 * np.exp(-2 * Phi) - meff2 + l / r * el * (dPhi - dlam)
                   + l * (l + 1) / r**2 * (el - 1)) * dt + dnn * q1 / st.fa2) / el
        return [ddp, dU, dtp, dtpp]
    return f


def interior(st, l, omega, U0, t0, dyn=True, r0=1e-3, frac=1.0):
    """Integra desde el centro con las condiciones (10); devuelve el estado en la superficie."""
    e, p, nb, cs2, q1, q2, dPhi, dlam, Phi, lam, t, tp = st.bg(r0)
    dp0 = np.exp(-2 * Phi) * omega**2 / l * U0 - q1 * t0
    y0 = np.array([dp0, U0, t0, 0.0], dtype=complex)
    Rs = st.R * frac
    f = _perturb_rhs(st, l, omega, dyn)
    # cerca del centro los términos l/r y l(l+1)/r^2 hacen el sistema rígido: implícito hasta r1
    r1 = 0.5
    s = solve_ivp(f, (r0, r1), y0, method='BDF', rtol=1e-9, atol=1e-40)
    s = solve_ivp(f, (r1, Rs), s.y[:, -1], method='DOP853', rtol=1e-9, atol=1e-40)
    return s.y[:, -1], Rs


def delta_p(st, Rs, y):
    """Perturbación lagrangiana de la presión en Rs: dp + p' xi^r."""
    e, p, nb, cs2, q1, q2, dPhi, dlam, Phi, lam, t, tp = st.bg(Rs)
    return y[0] - np.exp(-lam) * y[1] / Rs * (dPhi + q1 * tp / (e + p))


def surface_residual(st, l, omega, dyn=False):
    """|Delta p(R)| a frecuencia dada, sin delta theta (Tabla I y referencia)."""
    y, Rs = interior(st, l, omega, 1.0, 0.0, dyn=False)
    return delta_p(st, Rs, y).real


def exterior_L(st, l, omega, r_stop=None, full=False):
    """psi'/psi en R para psi = delta theta_ext, Ec. (13), con la condición de onda saliente
    (o ligada si Im k > 0) en r grande: psi ~ exp[i(k r + eta ln r)]/r, eta = M (2 w^2 - m_a^2)/k."""
    m2 = st.ma**2
    k = np.sqrt(complex(omega**2 - m2))
    if k.imag < 0 and abs(k.real) < 1e-14:
        k = -k
    Mt, R, rend = st.sol.y[0, -1], st.R, st.sol.t[-1]
    r_stop = r_stop or R
    rinf = min(r_stop + 40 / min(abs(k), st.ma), 6000.0)
    eta = Mt * (2 * omega**2 - m2) / k

    def bgx(r):
        if r < rend:
            m, p, Phi, t, tp = st.sol.sol(r)
            dm = st.rhs(r, [m, p, Phi, t, tp])[0]
            dPhi = st.rhs(r, [m, p, Phi, t, tp])[2]
            Phi = Phi + st.Phi_shift
            meff2 = -st.E0 * g2(t) / st.fa2
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
