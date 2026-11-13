"""Direct Urca with muons: verification of the lepton dependence and of the muon channel."""

import copy

import numpy as np
import pytest

import emissivity as em
from constants import HBARC_MEV_FM, MUON_MASS_MEV
from state import LocalState
from stellar_structure import StructureOptions, build_stellar_model
from test_emissivity_updated import core_points

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


@pytest.mark.parametrize("mass", [0.511, MUON_MASS_MEV, 400.0])
def test_the_lepton_mass_enters_the_matrix_element_only_through_the_energy(mass):
    """sum|M|^2 = 8 G^2 [(1 + 3 g_A^2) E eps + (1 - g_A^2) p.k]: no explicit mass term (G^2/2 * sum|A|^2)."""
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


def test_each_lepton_gives_the_known_coefficient_times_its_fermi_energy(integrals):
    """Q = (457 pi / 10080) G^2 cos^2 (1 + 3 g_A^2) m*_n m*_p mu_l T^6 for the electron and the muon alike."""
    for name, (q, mu) in integrals.items():
        assert q / mu / KNOWN_COEFFICIENT == pytest.approx(
            1.0, abs=5e-3
        ), name  # O((T/E_F)^2) corrections


def test_the_muon_to_electron_ratio_is_the_ratio_of_the_fermi_energies(integrals):
    q_e, mu_e = integrals["electron"]
    for name in ("muon", "muon_hot"):
        q_mu, mu_mu = integrals[name]
        assert q_mu / q_e == pytest.approx(
            mu_mu / mu_e, rel=5e-3
        )  # equal at mu_mu = mu_e (beta equilibrium)


def test_a_closed_channel_does_not_emit():
    closed, mu = phase_space_integral(0.1, 300.0, 150.0, 120.0, MUON_MASS_MEV)  # 150 + 120 < 300
    open_, _ = phase_space_integral(
        0.1, 300.0, 150.0, np.sqrt(190.0**2 - MUON_MASS_MEV**2), MUON_MASS_MEV
    )
    assert closed < 1e-6 * open_


# ------------------------------------------------------------------ 3. the implementation
KP = (1.3**3 + 0.8**3) ** (
    1 / 3
)  # proton Fermi momentum [fm^-1] for k_e = 1.3 and k_mu = 0.8 (charge neutrality)


def test_each_lepton_has_its_own_threshold():
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


def test_the_muon_channel_is_the_electron_formula_with_the_muon_fermi_energy():
    state = LocalState(1e9, core_points(2.1, KP, 0.8), superfluid=False)
    stripped = copy.copy(state)
    stripped.kf_mu = np.zeros_like(state.kf_mu)
    stripped.__dict__.pop("direct_urca_muon_allowed", None)
    electron_only = em.direct_urca_emissivity(stripped)[0]
    both = em.direct_urca_emissivity(state)[0]
    mu_muon = np.hypot(HBARC_MEV_FM * state.kf_mu[0], MUON_MASS_MEV)
    electron_factor = (state.ye[0] * state.rho[0] / 2.8e14) ** (1 / 3)
    expected = 1.0 + mu_muon / em.DIRECT_URCA_REFERENCE_ENERGY / electron_factor
    assert both / electron_only == pytest.approx(expected, rel=1e-12)
    assert em.DIRECT_URCA_REFERENCE_ENERGY == pytest.approx(337.295, abs=1e-3)


def test_the_pairing_factor_of_the_direct_urca_covers_the_muon_channel():
    paired = LocalState(2e8, core_points(2.1, KP, 0.8), superfluid=True)
    paired.tau_n[:] = 0.5
    paired.tau_p[:] = 0.3
    assert (
        em.process_factors(paired).Rdu[0] < 1.0
    )  # superfluid suppression applies where either channel is open


def test_in_a_star_the_muon_channel_opens_inside_the_electron_one(root, input_config):
    """BSk22 at 1.8 Msun: the muon region is a core subset of the electron region and the emission doubles there."""
    options = StructureOptions("BSk22", central_pressure=1.627427981e35)
    star = build_stellar_model(options, root / "eos", input_config)
    state = LocalState(1e9, star.background.cells, superfluid=False)
    electron, muon = state.direct_urca_allowed, state.direct_urca_muon_allowed
    assert muon.any() and electron.sum() > muon.sum() and np.all(electron[muon])
    stripped = copy.copy(state)
    stripped.kf_mu = np.zeros_like(state.kf_mu)
    stripped.__dict__.pop("direct_urca_muon_allowed", None)
    ratio = em.direct_urca_emissivity(state)[muon] / em.direct_urca_emissivity(stripped)[muon]
    np.testing.assert_allclose(ratio, 2.0, atol=0.01)  # mu_mu = mu_e in beta equilibrium


def test_the_default_stars_have_no_muon_direct_urca(background):
    """In the default BSk26 star both Direct Urca channels are closed everywhere."""
    state = LocalState(1e9, background.cells, superfluid=False)
    assert not state.direct_urca_muon_allowed.any()
