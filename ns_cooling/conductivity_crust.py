"""Thermal electron conductivity of the crust, B = 0 (Potekhin's fit, 2019 version).

Strong-degeneracy, point-ion approximation. NOT included: thermal averaging,
magnetic branches, ion conduction, nuclear form factor, Umklapp freeze-out.
Qimp is the impurity parameter <(Z-<Z>)^2>.
Inputs: local T [K], total rho [g/cm^3], nuclear A, Z and mass fraction xh.
Only n_i = rho*xh/(m_u*A) enters: free-neutron mass is not counted as ions.
"""

from math import sqrt, exp, log, log1p, isfinite, pi
from constants import (
    ALPHA_FS,
    BOLTZMANN_MEV_PER_K,
    C,
    ELECTRON_MASS_MEV,
    K_B,
    LAMBDA_C,
    M_ELECTRON,
    M_U,
)

# Physical constants of the relativistic units of the fit, from constants.py (they had been rounded
# to single precision; decision D-020). UM1 is a parameter of the fit.
AUM = M_U / M_ELECTRON  # m_u / m_e
AUD = M_ELECTRON / LAMBDA_C**3  # m_e / lambda_C^3 [g cm^-3]
BOHR = 1.0 / ALPHA_FS  # inverse fine-structure constant
PI = pi
PI_F = pi
THIRD = 1.0 / 3.0
# UM1 = 2.8 is a parameter of the fit; the source of the value is not recorded in the original implementation.
UM1 = 2.8
CONDUCTIVITY_UNIT = K_B * C / LAMBDA_C**2  # k_B c / lambda_C^2 [erg cm^-1 s^-1 K^-1]


# The numerical coefficients of the functions below (Coulomb logarithm, collision integrals, thermal and impurity terms) are
# those of Potekhin's electron-conductivity fit (2019 version, COULAN 1999, TAUEESY 2016) as in the original implementation;
# they have not been checked against the papers (item RF-1).
def coulomb_log(xs, xw, v):
    """Coulomb logarithm, 1999 point-nucleus formulation (no finite nuclear-size correction)."""
    if xs <= 0 or xw < 0 or not 0 <= v < 1:
        raise ValueError("Coulomb-logarithm arguments out of domain")
    ea = exp(-xw)
    e1 = 1 - ea
    b = xs * xw
    if xw <= 0.01:
        cl0 = log((xs + 1) / xs)
        cl1 = 0.5 * xw * (2 - 1 / (xs + 1) - 2 * xs * cl0)
        cl2 = 0.5 * xw * (1.5 - 3 * xs - 1 / (xs + 1) + 3 * xs * xs * cl0)
    elif xw > 100 and b > 100:
        cl0 = log1p(1 / xs)
        cl1 = (cl0 - 1 / (1 + xs)) / 2
        cl2 = (2 * xs + 1) / (2 * xs + 2) - xs * cl0
    elif xs < 0.01 and b < 0.01:
        cl1 = 0.5 * (ea * expint_scaled(xw) + log(xw) + 0.5772156649)
        cl2 = 0.5 * (xw - e1) / xw
    else:
        cl0 = log((xs + 1) / xs)
        el = expint_scaled(b) - expint_scaled(b + xw) * ea
        cl1 = 0.5 * (cl0 + xs / (xs + 1) * e1 - (1 + b) * el)
        cl2 = 0.5 * ((xw - e1) / xw - xs * xs / (1 + xs) * e1 - 2 * xs * cl0 + xs * (2 + b) * el)
    return cl1 - v * v * cl2


def electron_electron_time(x, temp):
    """TAUEESY 2016 fit with partial-degeneracy correction.

    log1p avoids cancellation; Y>1e8 uses the asymptotic branch of the
    published fit, with a variable-name typo of the source (CITL -> CILT)
    corrected.
    """
    e = sqrt(1 + x * x)
    v = x / e
    y = 0.0963913 / temp * x * sqrt(v)
    if y > 1e8:
        cil = 20.4013123 / (v * y**3)
        cit = 2.404 * v / y**2
        # The factored power avoids forming y**8.
        cilt = 18.52 * v ** (1.0 / 3.0) * y ** (-8.0 / 3.0)
    else:
        c1 = 0.123636 + 0.016234 * v * v
        c2 = 0.0762 + 0.05714 * v**4
        a = 12.2 + 25.2 * v**3
        c = a * exp(c1 / c2)
        yv = y * v
        cil = (
            log1p(128.56 / (37.1 * y + 10.83 * y * y + y**3))
            * (0.1587 - 0.02538 / (1 + 0.0435 * y))
            / v
        )
        cit = v**3 * log1p(c / (a * yv + yv * yv)) * (2.404 / c + (c2 - 2.404 / c) / (1 + 0.1 * yv))
        cilt = (
            v
            * log1p(c / (a * y + 10.83 * yv * yv + yv ** (8 / 3)))
            * (18.52 * v * v / c + (c2 - 18.52 * v * v / c) / (1 + 0.1558 * y ** (1 - 0.75 * v)))
        )
    freq = 0.00021381 * x * y * sqrt(v) * (cil + cit + cilt)
    theta = temp / (e - 1) if x > 0.001 else 2 * temp / x**2
    t = 25 * theta
    return (1 + t + 0.4342 * sqrt(theta) * t * t) / (1 + t * t) / freq


# Qimp = 0.1 is a free default for the impurity parameter; its source is not recorded in the original implementation.
def crust_electron_conductivity(T_K, rho, A, Z, xh, Qimp=0.1):
    """Scalar isotropic electron conductivity kappa_e [erg cm^-1 s^-1 K^-1]."""
    T_K, rho, A, Z, xh, Qimp = map(float, (T_K, rho, A, Z, xh, Qimp))
    if not all(map(isfinite, (T_K, rho, A, Z, xh, Qimp))):
        raise ValueError("Non-finite crust conductivity input")
    if T_K <= 0 or rho <= 0 or A < 1 or Z < 0.5 or not 0 < xh <= 1 or Qimp < 0:
        raise ValueError("Invalid crust conductivity input")
    temp = T_K * BOLTZMANN_MEV_PER_K / ELECTRON_MASS_MEV  # k_B T / (m_e c^2)
    ni = rho * xh / (AUD * AUM * A)
    # Validity limit of the fit in the ion density (in the relativistic units of the fit); source not recorded.
    if ni > 1e4:
        raise ValueError("Ion density exceeds the validity limit of the conductivity fit")
    sphere = (0.75 / PI / ni) ** THIRD
    gamma = Z**2 / BOHR / temp / sphere
    x = (3 * PI**2 * ni * Z) ** THIRD
    e = sqrt(1 + x * x)
    v = x / e
    qe = 4 / PI / BOHR * x * e
    # The Coulomb-logarithm part uses a single-precision pi.
    sphere_log = (0.75 / PI_F / (x**3 / 3 / PI_F**2 / Z)) ** THIRD
    qi = 3 * gamma / sphere_log**2
    trp = Z / gamma * sqrt(A * AUM * sphere_log / 3 / BOHR)
    born = v * Z * PI_F / BOHR
    c = (1 + 0.06 * gamma) * exp(-sqrt(gamma))
    xs = (qi * c + qe) * exp(-born) / (2 * x) ** 2
    xw = 13 / qi * (1 + 0.3333 * born) * (2 * x) ** 2
    cl = coulomb_log(xs, xw, v)
    a0 = 1.683 * sqrt(x / A / Z)
    vib = exp(-a0 / 4 * UM1 * exp(-9.1 * trp))
    t0 = 0.19 / Z**0.16667
    g0 = trp / sqrt(trp * trp + t0 * t0) * (1 + (Z / 125) ** 2)
    cleff = cl * g0 * vib
    g2 = trp / sqrt(0.0081 + trp * trp) ** 3
    thermal = 1 + g2 / g0 * (1 + born * v**3) * 0.0105 * (1 - 1 / Z)
    impurity = coulomb_log(qe / (2 * x) ** 2, 1e99, v) if Qimp else 0.0
    tau = x**3 / e / 4 / PI / ni / (Z / BOHR) ** 2 / cleff
    clz = cleff * Z * Z
    tau = tau * clz / (clz * thermal + impurity * Qimp)
    ee = electron_electron_time(x, temp)
    tau = tau * ee / (tau + ee)
    cel = x**3 / 3 / PI**2 / e / BOHR
    k = cel * PI**2 * temp / 3 * BOHR * tau * CONDUCTIVITY_UNIT
    if not isfinite(k) or k <= 0:
        raise FloatingPointError("Invalid crust conductivity result")
    return k


def expint_scaled(xi, order=0):
    """Scaled exponential integral exp(x) E_(order+1)(x), 21-step approximation."""
    if xi < 0 or order < 0:
        raise ValueError("x >= 0 and order >= 0 required")
    if xi >= 1:
        c = 0.0
        for i in range(21, 0, -1):
            c = i / (xi + c)
            c = (order + i) / (1 + c)
        return 1 / (xi + c)
    psi = -0.5772156649 + sum(1 / k for k in range(1, order + 1))
    q, cmx = 0.0, 1.0
    for m in range(22):
        if m:
            cmx = -cmx * xi / m
        dq = cmx / (m - order) if m != order else cmx * (log(xi + 1e-20) - psi)
        q -= dq
    return exp(xi) * q
