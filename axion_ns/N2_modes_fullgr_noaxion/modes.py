"""N2: Cowling and dynamic-metric frequencies of a polytrope without axion (moved from axion_ns/fmodes.py).

Same algorithms as the original; every integration takes an explicit `rtol` (default 1e-11 = the original). `cowling_X` and
`outgoing_amplitude` are compared with the original in archive/axion_ns_migration_scripts/compare_n2_legacy.py; `cowling_modes` and `full_mode` are ported verbatim
(plus the rtol argument) and not yet exercised in this layout.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar


def cowling_X(star, system, omega, ll, rtol=1e-11, atol=1e-14):
    y0 = star.center(omega, ll, 0.0, 1.0)[2:]
    return star.integrate(system, y0, omega, ll, rtol=rtol, atol=atol)[-1]


def cowling_modes(star, system, ll, grid, rtol=1e-11, xtol=1e-12, atol=1e-14):
    vals = [cowling_X(star, system, o, ll, rtol, atol) for o in grid]
    roots = []
    for a, b, fa_, fb in zip(grid[:-1], grid[1:], vals[:-1], vals[1:]):
        if fa_ * fb < 0:
            roots.append(brentq(lambda o: cowling_X(star, system, o, ll, rtol, atol), a, b, xtol=xtol))
    return roots


def outgoing_amplitude(star, system, omega, ll, n_far=40, rtol=1e-11, atol=1e-14):
    """|A_out|^2 of the exterior Zerilli wave for the interior solution with X(R) = 0 and W(R) = 1.
    At real frequency the exterior is a standing wave; its amplitude is minimal at the resonance."""
    ends = []
    for K0, W0 in [(1.0, 0.0), (0.0, 1.0)]:
        ends.append(star.integrate(system, star.center(omega, ll, K0, W0), omega, ll, rtol=rtol, atol=atol))
    y1, y2 = ends[0][3:], ends[1][3:]
    c = -y1[3] / y2[3]                              # X(R) = 0
    y = y1 + c * y2
    y = y / y[2]                                     # W(R) = 1
    M, R = star.M, star.R
    rr = R * (1 - 1e-7)
    B, _ = star.bg(rr, *ends[0][:3])
    dy = np.asarray(system.mat(rr, omega, ll, *B), dtype=float) @ y
    nn = (ll - 1) * (ll + 2) / 2
    f = lambda x: 1 - 2 * M / x
    h1, k, dh1, dk = y[0], y[1], dy[0], dy[1]
    pref = rr**(ll + 2) / (nn * rr + 3 * M)
    dpref = pref * ((ll + 2) / rr - nn / (nn * rr + 3 * M))
    Z = pref * (k - f(rr) * h1)
    dZ = dpref * (k - f(rr) * h1) + pref * (dk - 2 * M / rr**2 * h1 - f(rr) * dh1)
    VZ = lambda x: f(x) * (2 * nn**2 * (nn + 1) * x**3 + 6 * nn**2 * M * x**2 + 18 * nn * M**2 * x + 18 * M**3) \
        / (x**3 * (nn * x + 3 * M)**2)
    r_far = R + n_far * 2 * np.pi / omega
    s = solve_ivp(lambda x, u: [u[1] / f(x), -(omega**2 - VZ(x)) * u[0] / f(x)], (R, r_far), [Z, f(R) * dZ],
                  method='DOP853', rtol=rtol, atol=atol)
    if s.status < 0:                                # H32
        raise RuntimeError(f'modes: {s.message}')
    Zf, Pf = s.y[0, -1], s.y[1, -1]                 # Z and dZ/dr*
    L2 = ll * (ll + 1)
    rs = r_far + 2 * M * np.log(r_far / (2 * M) - 1)
    u = np.exp(1j * omega * rs) * (1 + 1j * L2 / (2 * omega * r_far))
    du = 1j * omega * u + f(r_far) * np.exp(1j * omega * rs) * (-1j * L2 / (2 * omega * r_far**2))
    # Z = 2 Re(a u), dZ/dr* = 2 Re(a du): two real equations for the complex a
    Mx = np.array([[2 * u.real, -2 * u.imag], [2 * du.real, -2 * du.imag]])
    ar, ai = np.linalg.solve(Mx, [Zf, Pf])
    return ar**2 + ai**2


def full_mode(star, system, ll, guess, width=0.04, rtol=1e-11, atol=1e-14):
    """Real part (minimum of |A_out|^2) and imaginary part (parabolic fit |A|^2 ~ (w-w0)^2 + g^2)."""
    unit = np.sqrt(star.M / star.R**3)
    amp = lambda s: outgoing_amplitude(star, system, s * unit, ll, rtol=rtol, atol=atol)
    m = minimize_scalar(amp, bracket=(guess * (1 - width), guess, guess * (1 + width)), tol=1e-10)
    s0 = m.x
    h = 1e-3 * s0
    a0_, ap, am = amp(s0), amp(s0 + h), amp(s0 - h)
    curv = (ap + am - 2 * a0_) / (2 * h**2)
    gamma = np.sqrt(max(a0_ / curv, 0.0))
    return s0, gamma
