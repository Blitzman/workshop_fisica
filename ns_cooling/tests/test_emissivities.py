"""Tests of the emissivities area (decision D-028). Run them with `make test-emissivities` when this part of the code changes, or on request.

Code covered: emissivity.py, emissivity_factors.py, emissivity_updated.py.
Every test states a hypothesis, a reference and a criterion, in the style of the validation steps of the GitHub code.
"""

import emissivity_updated as eu
import numpy as np
import pytest
from background import Points, fermi_momentum
from constants import MUON_MASS_MEV
from emissivity import direct_urca_emissivity, neutrino_emissivity
from scipy.integrate import quad
from state import LocalState
from synthetic import points_from_table


T9 = eu.K_B * 1e9 / eu.MEV_TO_ERG  # 1e9 K in MeV


RHO = eu.C_RHO  # the coupling of the rho meson (tests that switch it off restore it)


MN, MP = eu.NEUTRON_MASS_MEV, eu.PROTON_MASS_MEV


def core_points(kf_n, kf_p, kf_mu=0.0, mstar_n=1.0, mstar_p=1.0):
    """One homogeneous core point with the given Fermi momenta [fm^-1]."""
    n_n, n_p, n_mu = (k**3 / (3 * np.pi**2) for k in (kf_n, kf_p, kf_mu))
    nb = n_n + n_p
    rho = nb * 1e39 * 1.66e-24
    yn, yp, ymu = n_n / nb, n_p / nb, n_mu / nb
    a = lambda v: np.array([float(v)])  # noqa: E731
    return Points(
        rho=a(rho), ye=a(yp - ymu), yn=a(yn), yp=a(yp), xh=a(0), aa=a(1), zz=a(yp - ymu),
        ephi=a(1), tc_crust=a(0), tc_n=a(0), tc_p=a(0),
        kf_n=fermi_momentum(a(rho), a(yn)), kf_e=fermi_momentum(a(rho), a(yp - ymu)),
        mstar_n=a(mstar_n), mstar_p=a(mstar_p), transport_mstar_n=a(1), transport_mstar_p=a(1),
        ymu=a(ymu), kf_p=fermi_momentum(a(rho), a(yp)), kf_mu=fermi_momentum(a(rho), a(ymu)),
    )  # fmt: skip


def at_reference(T=1e9, **kwargs):
    """The reference state of the paper's figures: p_n = 1.75 fm^-1, p_p = 85 (p_n / 340 MeV)^2 MeV."""
    kf_n = 1.75
    kf_p = 85.0 * (kf_n * eu.HBARC_MEV_FM / 340.0) ** 2 / eu.HBARC_MEV_FM
    return LocalState(T, core_points(kf_n, kf_p, **kwargs), superfluid=False)


# ------------------------------------------------------------------ absolute values of the paper's figures
def erg(q_mev5):
    return q_mev5 * eu.MEV5_TO_CGS


def nn_figure_value(kf, c_rho, monkeypatch):
    monkeypatch.setattr(eu, "C_RHO", c_rho)
    p = kf * eu.HBARC_MEV_FM
    return erg(eu.N_FLAVORS * eu.BETA["nn"] * eu.identical_nucleon_bremsstrahlung(p, MN, T9))


# ------------------------------------------------------------------ the functions and integrals of the paper
@pytest.mark.parametrize("a, b", [(0.2, 0.9), (0.35, 1.8)])
def test_the_elementary_functions_are_the_integrals_they_stand_for(a, b):
    """Hypothesis: the closed forms of the updated rates are the integrals they stand for.
    Reference: numerical quadrature of their defining integrals (Bottaro, Caputo and Fiorillo 2024).
    Criterion: 1e-12 (phi, psi), 1e-10 (the two-argument functions) and the symmetry of the last one
        to 1e-14.
    """
    assert eu.phi(a) == pytest.approx(quad(lambda x: x**4 / (x**2 + a**2) ** 2, 0, 1)[0], rel=1e-12)
    assert eu.psi(a) == pytest.approx(quad(lambda x: x**2 / (x**2 + a**2), 0, 1)[0], rel=1e-12)
    cap_phi = quad(lambda x: x**4 / ((x**2 + a**2) * (x**2 + b**2)), 0, 1)[0]
    assert eu.cap_phi(a, b) == pytest.approx(cap_phi, rel=1e-10)

    # (2/pi) int int x^2 y^2 / ((x^2 + a^2)(y^2 + b^2) sqrt(1 - x^2 - y^2)); y = sqrt(1 - x^2) sin t
    def inner(x):
        return quad(
            lambda t: (1 - x * x) * np.sin(t) ** 2 / ((1 - x * x) * np.sin(t) ** 2 + b**2),
            0,
            np.pi / 2,
        )[0]

    cap_psi = 2 / np.pi * quad(lambda x: x**2 / (x**2 + a**2) * inner(x), 0, 1)[0]
    assert eu.cap_psi(a, b) == pytest.approx(cap_psi, rel=1e-10)
    assert eu.cap_psi(a, b) == pytest.approx(eu.cap_psi(b, a), rel=1e-14)


@pytest.mark.parametrize("p_f", [300.0])
@pytest.mark.parametrize("with_rho", [False, True])
def test_the_identical_nucleon_J_is_the_integral_of_the_amplitude(p_f, with_rho, monkeypatch):
    """Hypothesis: the closed form of J for identical nucleons is the integral of the amplitude (Eqs.
        3.3 and 3.6).
    Reference: a numerical double integral of the amplitude.
    Criterion: relative 1e-5, with and without the rho meson.
    """
    c_rho = eu.C_RHO if with_rho else 0.0
    monkeypatch.setattr(eu, "C_RHO", c_rho)

    def amplitude(x, mass_pi=eu.M_PI):  # A(x) of Eq. 3.3
        return x * x / (x * x + mass_pi**2) - c_rho * x * x / (x * x + eu.M_RHO**2)

    def m2(k, l):
        return (eu.F_PI / eu.M_PI) ** 4 * (
            amplitude(l) ** 2 + amplitude(k) ** 2 + amplitude(l) * amplitude(k)
        )

    def inner(k):  # l = sqrt(4 p_F^2 - k^2) sin t removes the endpoint singularity of Eq. 3.6
        return quad(lambda t: m2(k, np.sqrt(4 * p_f**2 - k * k) * np.sin(t)), 0, np.pi / 2)[0]

    numeric = quad(inner, 0, 2 * p_f)[0]
    assert eu.identical_nucleon_J(p_f) == pytest.approx(numeric, rel=1e-5)


@pytest.mark.parametrize("c_rho, landau", [(0.0, False), (0.0, True), (1.4, True)])
@pytest.mark.parametrize("p_n", [345.0])
def test_the_np_J_is_the_integral_of_the_amplitude(p_n, c_rho, landau, monkeypatch):
    """Hypothesis: the leading J^np of Eq. 3.16 is the integral of the amplitude of Eq. 3.13.
    Reference: a numerical integral over the exchanged momentum.
    Criterion: relative 1e-9: without rho and Landau parameters, with Landau parameters only, and
        with both.
    """
    monkeypatch.setattr(eu, "C_RHO", c_rho)
    if not landau:
        monkeypatch.setattr(eu, "F0_PRIME", 0.0)
        monkeypatch.setattr(eu, "G0", 0.0)
    p_p = 85.0 * (p_n / 340.0) ** 2
    ratio = (eu.F_PI / eu.M_PI) ** 2
    a_pi = lambda x: x * x / (x * x + eu.M_PI**2)  # noqa: E731
    a_rho = lambda x: x * x / (x * x + eu.M_RHO**2)  # noqa: E731
    h = lambda x: -ratio * (a_pi(x) - c_rho * a_rho(x))  # noqa: E731

    def m2(k):
        g_diff = eu.landau_f_minus_g(p_n) - c_rho * ratio * (
            a_rho(k) - a_rho(p_n)
        )  # g'_k - g'_l + f' - g
        hk, hl = h(k), h(p_n)
        return hk**2 + 2 * hl**2 + 4 * (hk - hl) * g_diff + 6 * g_diff**2 - 2 * hl * hk

    numeric = np.pi / 2 * quad(m2, 0, 2 * p_p)[0]
    assert eu.neutron_proton_J(p_n, p_p) == pytest.approx(numeric, rel=1e-9)


def test_nn_bremsstrahlung_reproduces_figure_1(monkeypatch):
    """Hypothesis: the nn bremsstrahlung reproduces Fig. 1 of Bottaro et al. (2024).
    Reference: values read from the figure at T = 1e9 K, 3 flavours, beta_nn = 0.56 and bare mass.
    Criterion: within 4 to 6 % of 4.35e19, 10.2e19, 2.5e19 and 4.55e19 erg cm^-3 s^-1.
    """
    assert nn_figure_value(1.75, RHO, monkeypatch) == pytest.approx(4.35e19, rel=0.04)
    assert nn_figure_value(1.75, 0.0, monkeypatch) == pytest.approx(
        10.2e19, rel=0.06
    )  # no rho meson
    assert nn_figure_value(1.0, RHO, monkeypatch) == pytest.approx(2.5e19, rel=0.06)
    assert nn_figure_value(2.0, RHO, monkeypatch) == pytest.approx(4.55e19, rel=0.05)


def test_np_bremsstrahlung_reproduces_figure_3(monkeypatch):
    """Hypothesis: the np bremsstrahlung reproduces Fig. 3 of Bottaro et al. (2024), and equals the
        standard rate without Landau parameters and rho meson.
    Reference: values read from the figure at p_n = 1.75 fm^-1, and the standard rate 1.5e20
        (Yakovlev et al. 2001, Eq. 166).
    Criterion: within 3 % of 5.6e19, 6.8e19 and 8.3e19 erg cm^-3 s^-1 and of the standard value.
    """
    p_n = 1.75 * eu.HBARC_MEV_FM
    p_p = 85.0 * (p_n / 340.0) ** 2

    def value(c_rho, landau):
        monkeypatch.setattr(eu, "C_RHO", c_rho)
        monkeypatch.setattr(eu, "F0_PRIME", 0.7 if landau else 0.0)
        monkeypatch.setattr(eu, "G0", 1.1 if landau else 0.0)
        return erg(
            eu.N_FLAVORS * eu.BETA["np"] * eu.neutron_proton_bremsstrahlung(p_n, p_p, MN, MP, T9)
        )

    assert value(1.4, True) == pytest.approx(5.6e19, rel=0.03)  # complete
    assert value(0.0, True) == pytest.approx(6.8e19, rel=0.03)  # Landau parameters, no rho
    assert value(1.4, False) == pytest.approx(8.3e19, rel=0.03)  # no Landau parameters
    standard_value = 3 * 1.5e20 * (p_p / eu.HBARC_MEV_FM / 1.68) * 1.06 * 0.66
    assert value(1.4, False) == pytest.approx(
        standard_value, rel=0.03
    )  # equals the standard rate there


# ------------------------------------------------------------------ scalings and regions
def test_the_rates_scale_as_t8_and_with_the_effective_masses():
    """Hypothesis: every updated rate scales as T^8, and with the effective masses as the formulas say.
    Reference: the analytic scaling of the formulas.
    Criterion: relative 1e-12.
    """
    processes = eu.make_processes()
    base = at_reference()
    hot = at_reference(T=2e9)
    for name, function in processes.items():
        assert function(hot)[0] == pytest.approx(2**8 * function(base)[0], rel=1e-12), name
    heavy = at_reference(mstar_n=0.8, mstar_p=0.6)
    expected = dict(
        modified_urca=0.8**3 * 0.6, brems_nn=0.8**4, brems_np=0.8**2 * 0.6**2, brems_pp=0.6**4
    )
    for name, factor in expected.items():
        if (
            name == "modified_urca"
        ):  # the proton branch scales as m*_n m*_p^3: compare the neutron branch only
            continue
        assert processes[name](heavy)[0] == pytest.approx(
            factor * processes[name](base)[0], rel=1e-12
        ), name


# ------------------------------------------------------------------ 1. the matrix element, with explicit spinors
S0 = np.eye(2, dtype=complex)


SIG = [
    np.array([[0, 1], [1, 0]], complex),
    np.array([[0, -1j], [1j, 0]]),
    np.array([[1, 0], [0, -1]], complex),
]


Z2 = np.zeros((2, 2), complex)


GAMMA = [np.block([[S0, Z2], [Z2, -S0]])] + [np.block([[Z2, s], [-s, Z2]]) for s in SIG]


GAMMA5 = 1j * GAMMA[0] @ GAMMA[1] @ GAMMA[2] @ GAMMA[3]


BASIS = [np.array([1, 0], complex), np.array([0, 1], complex)]


def sigma_dot(p):
    return sum(SIG[i] * p[i] for i in range(3))


def spin_summed_amplitude_squared(p, mass, k, g_a):
    """sum over all spins of |A|^2, A = L^0 <p|n> - g_A sum_i L^i <p|sigma_i|n>, L^mu = ubar gamma^mu (1 - g5) v.

    The lepton has momentum p and mass ``mass`` (Dirac spinor u), the antineutrino momentum k (massless v);
    the nucleons are non-relativistic Pauli spinors.
    """
    energy = np.sqrt(p @ p + mass**2)
    total = 0.0
    for chi in BASIS:
        u = np.concatenate(
            [np.sqrt(energy + mass) * chi, sigma_dot(p) @ chi / np.sqrt(energy + mass)]
        )
        ubar = u.conj() @ GAMMA[0]
        for eta in BASIS:
            e_nu = np.linalg.norm(k)
            v = np.concatenate([sigma_dot(k) @ eta / np.sqrt(e_nu), np.sqrt(e_nu) * eta])
            current = [ubar @ GAMMA[mu] @ (np.eye(4) - GAMMA5) @ v for mu in range(4)]
            for n_spin in BASIS:
                for p_spin in BASIS:
                    h0 = p_spin.conj() @ n_spin
                    hi = [p_spin.conj() @ SIG[i] @ n_spin for i in range(3)]
                    amplitude = current[0] * h0 - g_a * sum(
                        current[1 + i] * hi[i] for i in range(3)
                    )
                    total += abs(amplitude) ** 2
    return total


# ------------------------------------------------------------------ 2. the phase-space integral, exact kinematics
def fermi(x):
    return 0.5 * (1.0 - np.tanh(0.5 * x))


def phase_space_integral(T, p_fn, p_fp, p_fl, mass, m_n=500.0, m_p=500.0, limit=28.0, n=100):
    """Q / (G^2 cos^2 (1 + 3 g_A^2) m*_n m*_p T^6) * 2 pi^5 for one lepton, with exact kinematics.

    With p_nu neglected in the momentum conservation, nucleons with parabolic quasi-particle
    energies and the neutrino-angle averaged matrix element of the test above:
    Q = 2 G^2 c^2 (1 + 3 g_A^2) (2 pi)^-8 4 pi 4 pi 2 pi * 2 (decay and capture) * m*_n m*_p T^6 K,
    K = int dx_n dx_p dy E_l(y) (x_n - x_p - y)^3 f_n (1 - f_p)(1 - f_l) theta(|cos| <= 1), where
    cos is the angle between p_n and p_p fixed by p_l = |p_n - p_p| (the Jacobian is E_l dE_l / (p_n p_p)).
    """
    mu = np.sqrt(p_fl**2 + mass**2)
    x, w = np.polynomial.legendre.leggauss(n)
    x, w = limit * x, limit * w
    xn, xp, y = np.meshgrid(x, x, x, indexing="ij")
    weight = w[:, None, None] * w[None, :, None] * w[None, None, :]
    pn2, pp2 = p_fn**2 + 2 * m_n * T * xn, p_fp**2 + 2 * m_p * T * xp
    energy = mu + T * y
    pl2 = energy**2 - mass**2
    eps = xn - xp - y
    ok = (pn2 > 0) & (pp2 > 0) & (pl2 > 0) & (eps > 0)
    pn, pp = np.sqrt(np.where(ok, pn2, 1.0)), np.sqrt(np.where(ok, pp2, 1.0))
    cosine = (pn2 + pp2 - np.where(ok, pl2, 1.0)) / (2 * pn * pp)
    ok &= np.abs(cosine) <= 1
    integrand = np.where(ok, energy * eps**3 * fermi(xn) * fermi(-xp) * fermi(-y), 0.0)
    return float(np.sum(weight * integrand)), mu


KNOWN_COEFFICIENT = (
    457 * np.pi / 10080
)  # Q = coefficient G^2 cos^2 (1 + 3 g_A^2) m*_n m*_p mu_l T^6


@pytest.fixture(scope="module")
def integrals():
    cases = dict(
        electron=(0.511, 190.0), muon=(MUON_MASS_MEV, 190.0), muon_hot=(MUON_MASS_MEV, 210.0)
    )
    result = {}
    for name, (mass, mu) in cases.items():
        p_f = np.sqrt(mu**2 - mass**2)
        k, _ = phase_space_integral(
            0.1, 300.0, 150.0, p_f, mass
        )  # open: 150 + p_F > 300 with margin
        result[name] = (k / (2 * np.pi**5), mu)
    return result


# ------------------------------------------------------------------ 3. the implementation
KP = (1.3**3 + 0.8**3) ** (
    1 / 3
)  # proton Fermi momentum [fm^-1] for k_e = 1.3 and k_mu = 0.8 (charge neutrality)


@pytest.mark.parametrize("mass", [MUON_MASS_MEV])
def test_the_lepton_mass_enters_the_matrix_element_only_through_the_energy(mass):
    """Hypothesis: the lepton mass enters the squared matrix element of the direct Urca only through
        the lepton energy.
    Reference: sum|M|^2 = 8 G^2 [(1 + 3 g_A^2) E eps + (1 - g_A^2) p.k], derived analytically,
        against an explicit spinor calculation.
    Criterion: relative 1e-12 for the muon mass and four random kinematics.
    """
    rng = np.random.default_rng(11)
    g_a = 1.27
    for _ in range(4):
        p = rng.normal(size=3)
        p *= rng.uniform(80, 300) / np.linalg.norm(p)
        k = rng.normal(size=3)
        k *= rng.uniform(0.2, 5.0) / np.linalg.norm(k)
        energy, eps = np.sqrt(p @ p + mass**2), np.linalg.norm(k)
        expected = 16 * ((1 + 3 * g_a**2) * energy * eps + (1 - g_a**2) * (p @ k))
        assert spin_summed_amplitude_squared(p, mass, k, g_a) == pytest.approx(expected, rel=1e-12)


def test_each_lepton_gives_the_known_coefficient_times_its_fermi_energy(integrals):
    """Hypothesis: the exact phase-space integral gives the known coefficient of the direct Urca for
        the electron and the muon alike.
    Reference: Q = (457 pi / 10080) G^2 cos^2(theta_C) (1 + 3 g_A^2) m*_n m*_p mu_l T^6 (Lattimer et
        al. 1991; Yakovlev et al. 2001). It validates the prefactor independently.
    Criterion: the ratio to the coefficient is 1 within 5e-3 (O((T/E_F)^2) corrections).
    """
    for name, (q, mu) in integrals.items():
        assert q / mu / KNOWN_COEFFICIENT == pytest.approx(
            1.0, abs=5e-3
        ), name  # O((T/E_F)^2) corrections


def test_each_lepton_has_its_own_threshold():
    """Hypothesis: the electron and the muon channels open at their own thresholds, k_n <= k_p + k_l.
    Reference: the triangle condition, on constructed charge-neutral points.
    Criterion: the open and closed pattern is exact in four configurations: both open, only the
        electron, none, and npe matter with Yp > 1/9.
    """
    # charge neutral points, k_p^3 = k_e^3 + k_mu^3 (k_e = 1.3, k_mu = 0.8 fm^-1): k_p + k_e = 2.69, k_p + k_mu = 2.19
    open_both = LocalState(1e9, core_points(2.1, KP, 0.8), superfluid=False)
    only_electron = LocalState(
        1e9, core_points(2.4, KP, 0.8), superfluid=False
    )  # k_p + k_mu < k_n < k_p + k_e
    closed = LocalState(1e9, core_points(2.8, KP, 0.8), superfluid=False)
    assert open_both.direct_urca_allowed[0] and open_both.direct_urca_muon_allowed[0]
    assert only_electron.direct_urca_allowed[0] and not only_electron.direct_urca_muon_allowed[0]
    assert not closed.direct_urca_allowed[0] and not closed.direct_urca_muon_allowed[0]
    no_muons = LocalState(1e9, core_points(1.9, 1.0, 0.0), superfluid=False)  # npe: Y_p = 1/8 > 1/9
    assert no_muons.direct_urca_allowed[0] and not no_muons.direct_urca_muon_allowed[0]


def npe_points(ye, rho=1e15):
    """Core-only npe points (Yp = Ye, Yn = 1 - Ye) with the given electron fractions."""
    ye = np.asarray(ye, dtype=float)
    n = len(ye)
    zeros = np.zeros(n)
    table = dict(
        radius_km=np.arange(n, dtype=float),
        rho=np.full(n, rho),
        ye=ye,
        yn=1 - ye,
        yp=ye,
        xh=zeros,
        aa=np.ones(n),
        zz=ye,
        ephi=np.full(n, 0.8),
        elam=np.full(n, 1.2),
        tc_crust=zeros,
        tc_n=zeros,
        tc_p=zeros,
    )
    return points_from_table(table)


# ------------------------------------------------------------------ Direct Urca threshold
YE_SCAN = np.array([0.10, 0.1105, 0.11110, 0.11112, 0.1112, 0.12, 0.2])  # 1/9 = 0.111111...


def test_the_threshold_does_not_depend_on_the_density():
    """Hypothesis: in npe matter the threshold is Yp = 1/9 at any density.
    Reference: the triangle condition with k_p = k_e.
    Criterion: the pattern of open points is the same at 3e14, 1e15 and 2e15 g/cm^3.
    """
    for rho in (3e14, 1e15, 2e15):
        state = LocalState(1e9, npe_points(YE_SCAN, rho=rho))
        assert list(state.direct_urca_allowed) == [False, False, False, True, True, True, True]


def test_direct_urca_emissivity_switches_on_at_the_threshold():
    """Hypothesis: the emissivity vanishes where the channel is closed, is positive where it is open,
        and scales as T^6.
    Reference: the threshold and the T^6 law of the direct Urca. The prefactor is checked by the
        exact-coefficient test, not here.
    Criterion: zero below the threshold, positive above, a ratio of 64 for a doubled temperature to
        1e-12, and the registered process equal to the function.
    """
    faces = npe_points(YE_SCAN)
    state = LocalState(1e9, faces, superfluid=False)
    allowed = state.direct_urca_allowed
    q = direct_urca_emissivity(state)
    assert np.all(q[~allowed] == 0) and np.all(q[allowed] > 0)
    q2 = direct_urca_emissivity(LocalState(2e9, faces, superfluid=False))
    np.testing.assert_allclose(q2[allowed] / q[allowed], 64.0, rtol=1e-12)  # T^6
    np.testing.assert_array_equal(
        neutrino_emissivity(state)["direct_urca"], q
    )  # it is the registered process
