"""N3: modes of the star with the axion and the dynamic metric (Paper II extended, stage B), table-based.

`kHz`, `omega_of`, `_ok`, `exterior_L` and `muller` are N1's (extracted verbatim: modules never import each other). The rest is the port of
axion_ns/axgr.py (center, basis, zerilli_amplitudes, incoming) with the right-hand sides taken from the Laurent tables (tables.py), plus `find_mode`
(the Muller search of axion_ns/step10._mode) and `matching_radius`.

Mode = complex zero of the INGOING gravitational amplitude of the solution that satisfies Delta p(R) = 0 and the axion condition outside; the amplitude is analytic in
omega (the zeros are the quasi-normal modes). Seeds: (K0, W0) and, with the axion, delta theta(0). dyn=False: delta theta switched off (the axion background
only, as Table I of Paper II).
"""
import numpy as np
from scipy import constants as C
from scipy.integrate import solve_ivp

from .background import bgvals

CKMS = C.c / 1e3            # km/s (same value as axion_ns/axstar.py)
ATOL0 = 1e-15               # absolute tolerance on seeds of O(1) (H29: with 1e-40 the integrators choked on components that start at 0)


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


def matching_radius(st):
    """Radius where the interior solution is joined to the decoupled exterior: 70% of the valid axion tail. The tail of theta_0 is only trustworthy up to r_end
    (the integration leaves the stable branch there, at a distance set by the shooting tolerance); the legacy matched where |theta| < 1e-7 or at the end, which with
    a 1e-6 shooting tolerance would be inside the corrupted part. Test: for the 6 kHz f-mode the root moves 1.4e-8 (Re) and 3e-6 (Im) between 50 and 60 km."""
    return st.R + 0.7 * (st.rend - st.R)


def center(st, omega, ll, K0, W0, T0, r0):
    """Conditions at the centre of the extended system: LD plus the axion corrections in X(0) and delta theta''(0) (leading order)."""
    b = bgvals(st, r0)
    P, xv, e0v, p, M1v, Vvv = b[0], b[5], b[9], b[10], b[14], b[16]
    eps = e0v + xv * b[13]
    if st.nu_a:
        xc = st.matter(st.p0 * st.KGEO, st.thc)[2]
        tpp0 = (4 * st.ma**2 / st.beta * (xc - 1) * st.g1(st.thc)) / 3
    else:
        tpp0 = 0.0
    H1 = (2 * ll * K0 + 16 * np.pi * (eps + p) * W0) / (ll * (ll + 1))
    X = (eps + p) * np.exp(P) * ((4 * np.pi / 3 * (eps + 3 * p) - omega**2 * np.exp(-2 * P) / ll) * W0 + K0 / 2)
    if st.nu_a:
        X += np.exp(P) * (xv * M1v * (T0 + tpp0 * W0) - 8 * np.pi / 3 * (eps + p) * Vvv * W0)
    return np.array([H1, K0, W0, X, T0, 0.0], dtype=complex)


def basis(st, tabs, omega, ll, dyn=True, r0=1e-3, r1=0.5, rtol=1e-10, rm=None):
    """Regular solutions from the centre to R and their continuation outside to r_m. Returns [(X(R), y_exterior(r_m)) per seed]."""
    if ll != tabs.ell:
        raise ValueError(f'Mode ell={ll} does not match coefficient tables with ell={tabs.ell}')
    dyn = dyn and bool(st.nu_a)
    seeds = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0)] + ([(0.0, 0.0, 1.0)] if dyn else [])
    ki, ke = (None, None) if dyn else ((0, 1, 2, 3), (0, 1))
    fi, fo = tabs.rhs_int(omega, ki), tabs.rhs_ext(omega, ke)
    R, rm = st.R, rm or matching_radius(st)
    out = []
    for K0, W0, T0 in seeds:
        y0 = center(st, omega, ll, K0, W0, T0, r0)
        if ki is not None:
            y0 = y0[list(ki)]
        s = _ok(solve_ivp(fi, (r0, r1), y0, method='DOP853', rtol=rtol, atol=ATOL0))
        y1 = s.y[:, -1]
        s = _ok(solve_ivp(fi, (r1, R), y1, method='DOP853', rtol=rtol, atol=ATOL0 * np.max(np.abs(y1))))
        yR = s.y[:, -1]
        yR6 = yR if ki is None else np.concatenate([yR, np.zeros(2, dtype=complex)])
        ye = np.array([yR6[0], yR6[1], yR6[4], yR6[5]]) if ke is None else np.array([yR6[0], yR6[1]])
        if rm > R:
            s = _ok(solve_ivp(fo, (R, rm), ye, method='DOP853', rtol=rtol, atol=1e-40))
            ye = s.y[:, -1]
        if ke is not None:
            ye = np.concatenate([ye, np.zeros(2, dtype=complex)])
        out.append((yR6[3], ye))
    return out


def zerilli_amplitudes(st, tabs, omega, ll, ye, rm, n_far=40):
    """(A_out, A_in) of the Zerilli wave built with (H1, K) at r_m (outside, no background axion)."""
    dy = tabs.mat_ext(rm, omega) @ ye
    M = st.y_at(rm)[0]
    nn = (ll - 1) * (ll + 2) / 2
    f = lambda x_: 1 - 2 * M / x_
    h1, k, dh1, dk = ye[0], ye[1], dy[0], dy[1]
    pref = rm**(ll + 2) / (nn * rm + 3 * M)
    dpref = pref * ((ll + 2) / rm - nn / (nn * rm + 3 * M))
    Z = pref * (k - f(rm) * h1)
    dZ = dpref * (k - f(rm) * h1) + pref * (dk - 2 * M / rm**2 * h1 - f(rm) * dh1)
    VZ = lambda x_: f(x_) * (2 * nn**2 * (nn + 1) * x_**3 + 6 * nn**2 * M * x_**2 + 18 * nn * M**2 * x_ + 18 * M**3) / (x_**3 * (nn * x_ + 3 * M)**2)
    r_far = rm + n_far * 2 * np.pi / abs(omega.real)
    s = _ok(solve_ivp(lambda x_, u: [u[1] / f(x_), -(omega**2 - VZ(x_)) * u[0] / f(x_)], (rm, r_far), np.array([Z, f(rm) * dZ], dtype=complex),
                      method='DOP853', rtol=1e-11, atol=1e-60))
    Zf, Pf = s.y[0, -1], s.y[1, -1]
    L2 = ll * (ll + 1)
    rs = r_far + 2 * M * np.log(r_far / (2 * M) - 1)
    out = []
    for sg in (1, -1):        # u = outgoing, v = ingoing (analytic continuations in omega)
        e = np.exp(sg * 1j * omega * rs)
        u = e * (1 + sg * 1j * L2 / (2 * omega * r_far))
        du = sg * 1j * omega * u + f(r_far) * e * (-sg * 1j * L2 / (2 * omega * r_far**2))
        out.append((u, du))
    (u, du), (v, dv) = out
    return np.linalg.solve(np.array([[u, v], [du, dv]]), np.array([Zf, Pf]))


def incoming(st, tabs, omega, ll, dyn=True, rtol=1e-10, rm=None, n_far=40):
    """Ingoing amplitude of the gravitational wave for the solution with Delta p(R) = 0 and the axion condition outside. Analytic in omega."""
    rm = rm or matching_radius(st)
    B = basis(st, tabs, omega, ll, dyn=dyn, rtol=rtol, rm=rm)
    c1 = np.array([b[0] for b in B])                                  # X(R)
    if len(B) == 3:
        L = exterior_L(st, ll, omega, r_stop=rm) - ll / rm
        c2 = np.array([b[1][3] - L * b[1][2] for b in B])             # delta theta' - L delta theta at r_m
        c = np.cross(c1, c2)
    else:
        c = np.array([c1[1], -c1[0]])
    ye = sum(ci * b[1] for ci, b in zip(c, B))
    return zerilli_amplitudes(st, tabs, omega, ll, ye, rm, n_far)[1]


def find_mode(st, tabs, guess, ll=None, dyn=True, **kw):
    """Complex zero of the ingoing amplitude near `guess` (kHz), as axion_ns/step10._mode."""
    ll = tabs.ell if ll is None else ll
    f = lambda z: incoming(st, tabs, omega_of(z), ll, dyn=dyn, **kw)
    return complex(muller(f, guess * 0.997, guess, guess * 1.003 - 1e-5j, tol=1e-9))
