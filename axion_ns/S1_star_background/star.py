"""S1: neutron star BSk26 with an axion condensate: background (TOV + Klein-Gordon of Paper I, Eq. 2.6).

Moved from axion_ns/axstar.py (the background half: constants, EOS instance, Star class; unchanged). The mode half
(Cowling perturbations, exterior, Muller) moved to N1 and reads the tables written by export.py instead of importing this.

Geometric units G = c = 1 with lengths in km. All axion scales come from two numbers: E0 = eps m_pi^2 f_pi^2
(= n_c sigma_N) and m_a; f_a^2 = beta E0 / (4 m_a^2) and x = n_b / n_c.
Hand-coded ODEs: they are NOT yet compiled from D0's equations (see README, section "Known limitations").

Tolerances: the original code integrated with rtol = 1e-11 and located theta_c to 1e-16 (about 66 integrations of ~1.4 s
each). This module now defaults to rtol = 1e-6 and a shooting precision of 1e-6, which was measured to be 22x faster with
changes of ~1e-8 in M, R and the mode frequencies (MIGRATION.md, "Deliberate differences" 5). Pass rtol=1e-11,
shoot_tol=None to reproduce the original precision.
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

    def __init__(self, nu_a=None, p0=P0, e0=E0, rtol=1e-6, atol=1e-30, max_step=0.25, shoot_tol=1e-6):
        # Integration controls. DEFAULTS (consolidated, see MIGRATION.md): rtol = 1e-6 and shoot_tol = 1e-6 (the bisection of
        # delta = pi - theta_c stops when its bracket is 1e-6 relative). atol (1e-30) and max_step (0.25 km) are unchanged.
        # LEGACY precision of axion_ns: rtol=1e-11, shoot_tol=None (64 bisections in ln delta + linear refinement to 1e-16).
        self.rtol, self.atol, self.max_step, self.shoot_tol = rtol, atol, max_step, shoot_tol
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
    def _background_matter(self, p, shape, slope):
        """EOS and first-order scalar coupling used by hydrostatic equilibrium."""
        pm = p / KGEO
        if pm <= _eos.p_min * 1.0001:
            return 0.0, 0.0, 0.0, 0.0, 0.0
        e0, nb = _eos.eps(pm) * KGEO, _eos.nb(pm)
        x = nb / self.nc
        return e0 - self.E0 * x * shape, e0, nb, x, self.E0 * x * slope

    def matter(self, p, t):
        """(eps, n_b, x, c_s^2, q1, q2) for background and perturbation output."""
        pm = p / KGEO
        if pm <= _eos.p_min * 1.0001:
            return 0.0, 0.0, 0.0, 1.0, 0.0, 0.0
        e, e0, nb, x, q1 = self._background_matter(
            p, (1 - g(t)) if self.nu_a else 0.0,
            g1(t) if self.nu_a else 0.0)
        cs2 = (e0 + p) * _eos.cs2(pm) / (e + p)
        q2 = self.E0 * x * g2(t) if self.nu_a else 0.0
        return e, nb, x, cs2, q1, q2

    def _hydro_rhs(self, r, m, p, t, tp):
        """TOV derivatives and coefficients shared by both scalar equations."""
        inside = p / KGEO > _eos.p_min * 1.0001
        pr = p if inside else 0.0
        F = 1 - 2 * m / r
        if self.nu_a:
            gg = g(t)
            shape = 1 - gg
            slope = -BETA * np.sin(t) / (4 * gg)
        else:
            shape = slope = 0.0
        V = self.E0 * shape
        e, _, _, x, q1 = self._background_matter(p, shape, slope)
        kin = self.fa2 * tp**2 / 2 * F if self.nu_a else 0.0
        dPhi = (m + 4 * np.pi * r**3 * (pr - V + kin)) / (r**2 * F)
        dm = 4 * np.pi * r**2 * (e + V + kin)
        dp = (-(e + pr) * dPhi - q1 * tp) if inside else 0.0
        return dm, dp, dPhi, e, x, F, V, pr, slope

    def rhs(self, r, y):
        m, p, _, t, tp = y
        dm, dp, dPhi, e, x, F, V, pr, slope = self._hydro_rhs(r, m, p, t, tp)
        if self.nu_a:
            src = 4 * self.ma**2 / BETA * (x - 1) * slope
            tpp = (src - 2 * tp / r * (1 - m / r - 2 * np.pi * r**2 * (e - pr + 2 * V))) / F
        else:
            tpp = 0.0
        return [dm, dp, dPhi, tp, tpp]

    def _start(self, r0, thc):
        p = self.p0 * KGEO
        shape = 1 - g(thc) if self.nu_a else 0.0
        V = self.E0 * shape
        e, _, _, x, _ = self._background_matter(p, shape, 0.0)
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
                      method='DOP853', rtol=self.rtol, atol=self.atol, dense_output=keep, max_step=self.max_step)
        # H32: un disparo que no termina no se clasifica en silencio. Excepción: el integrador se para justo en la
        # superficie (p en el umbral p_min de la EOS, donde el lado derecho salta de la materia al vacío). Pasa con el
        # politropo del paso 17: es la llegada a la superficie y R sale de ahí igual que antes (14,155 km)
        if s.status < 0 and not s.y[1, -1] / KGEO <= _eos.p_min * 1.0002:
            raise RuntimeError(f'Star._integrate: {s.message} (r = {s.t[-1]:.6g}, p = {s.y[1, -1] / KGEO:.3g})')
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
            raise ValueError('el disparo con theta_c = 0.3 no cruza cero: o no hay condensado, o tiene theta_c < 0.3 '
                             '(fuera del intervalo de búsqueda)')
        if self.shoot_tol is None:
            n_bis, tol_lin = 64, 1e-16
        else:
            n_bis, tol_lin = 200, self.shoot_tol
        for _ in range(n_bis):
            mid = 0.5 * (lo + hi)
            if self._integrate(np.pi - np.exp(mid)) > 0:
                lo = mid
            else:
                hi = mid
            if self.shoot_tol is not None and hi - lo < self.shoot_tol:     # bracket of ln(delta) = relative bracket of delta
                break
        dlo, dhi = np.exp(lo), np.exp(hi)
        for _ in range(60):                         # linear refinement
            if self.shoot_tol is not None:
                break                               # already within shoot_tol
            mid = 0.5 * (dlo + dhi)
            if self._integrate(np.pi - mid) > 0:
                dlo = mid
            else:
                dhi = mid
            if dhi - dlo < tol_lin * max(dhi, 1e-300) or mid in (dlo, dhi):
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
