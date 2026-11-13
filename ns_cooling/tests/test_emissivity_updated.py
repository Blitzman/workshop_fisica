"""Updated modified Urca and nucleon bremsstrahlung (Bottaro, Caputo and Fiorillo, arXiv:2406.18640)."""

import copy

import numpy as np
import pytest
from scipy.integrate import quad

import compare_emissivities as ce
import emissivity as standard
import emissivity_updated as eu
from background import Points, fermi_momentum
from emissivity import NEUTRINO_PROCESSES
from state import LocalState

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


# ------------------------------------------------------------------ the functions and integrals of the paper
@pytest.mark.parametrize("a, b", [(0.2, 0.9), (0.35, 1.8)])
def test_the_elementary_functions_are_the_integrals_they_stand_for(a, b):
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


@pytest.mark.parametrize("p_f", [300.0, 400.0])
@pytest.mark.parametrize("with_rho", [False, True])
def test_the_identical_nucleon_J_is_the_integral_of_the_amplitude(p_f, with_rho, monkeypatch):
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


@pytest.mark.parametrize("c_rho, landau", [(0.0, False), (0.0, True), (1.4, False), (1.4, True)])
@pytest.mark.parametrize("p_n", [280.0, 345.0, 420.0])
def test_the_np_J_is_the_integral_of_the_amplitude(p_n, c_rho, landau, monkeypatch):
    """The leading J^np of Eq. 3.16 against the amplitude of Eq. 3.13 integrated over the exchanged momentum."""
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


# ------------------------------------------------------------------ absolute values of the paper's figures
def erg(q_mev5):
    return q_mev5 * eu.MEV5_TO_CGS


def nn_figure_value(kf, c_rho, monkeypatch):
    monkeypatch.setattr(eu, "C_RHO", c_rho)
    p = kf * eu.HBARC_MEV_FM
    return erg(eu.N_FLAVORS * eu.BETA["nn"] * eu.identical_nucleon_bremsstrahlung(p, MN, T9))


def test_nn_bremsstrahlung_reproduces_figure_1(monkeypatch):
    """T = 1e9 K, 3 flavors, beta_nn = 0.56, bare mass; values read from the figure at p_n = 1.75 fm^-1."""
    assert nn_figure_value(1.75, RHO, monkeypatch) == pytest.approx(4.35e19, rel=0.04)
    assert nn_figure_value(1.75, 0.0, monkeypatch) == pytest.approx(
        10.2e19, rel=0.06
    )  # no rho meson
    assert nn_figure_value(1.0, RHO, monkeypatch) == pytest.approx(2.5e19, rel=0.06)
    assert nn_figure_value(2.0, RHO, monkeypatch) == pytest.approx(4.55e19, rel=0.05)


def test_np_bremsstrahlung_reproduces_figure_3(monkeypatch):
    """p_p = 85 (p_n / 340 MeV)^2 MeV, 3 flavors, beta_np = 0.66, bare masses, p_n = 1.75 fm^-1."""
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


@pytest.mark.parametrize("kf_n, low, high", [(1.5, 0.55, 0.70), (2.2, 0.28, 0.40)])
def test_the_murca_neutron_branch_relative_to_the_standard_rate(kf_n, low, high):
    """In npe matter with p_e = p_p the ratio depends only on p_n; the paper's Fig. 6 gives ~0.66 and ~0.37."""
    p_n = kf_n * eu.HBARC_MEV_FM
    p_p = 85.0 * (p_n / 340.0) ** 2
    new = eu.BETA["murca"] * erg(
        eu.murca_neutron_branch(p_n, p_p, eu.lepton_velocity(p_p, eu.ELECTRON_MASS_MEV), MN, MP, T9)
    )
    old = 8.55e21 * (p_p / eu.HBARC_MEV_FM / 1.68) * (1.76 - 0.63 * (1.68 / kf_n) ** 2) * 0.68
    assert low < new / old < high


def test_without_the_rho_meson_the_murca_neutron_branch_is_close_to_the_standard_rate(monkeypatch):
    monkeypatch.setattr(eu, "C_RHO", 0.0)
    p_n = 1.5 * eu.HBARC_MEV_FM
    p_p = 85.0 * (p_n / 340.0) ** 2
    new = eu.BETA["murca"] * erg(eu.murca_neutron_branch(p_n, p_p, 1.0, MN, MP, T9))
    old = 8.55e21 * (p_p / eu.HBARC_MEV_FM / 1.68) * (1.76 - 0.63 * (1.68 / 1.5) ** 2) * 0.68
    assert new / old == pytest.approx(1.0, abs=0.06)


# ------------------------------------------------------------------ the muon channels and their thresholds
def test_the_neutron_branch_is_proportional_to_the_sum_of_the_lepton_velocities():
    p_n, p_p = 450.0, 150.0
    v_e = eu.lepton_velocity(150.0, eu.ELECTRON_MASS_MEV)
    v_mu = eu.lepton_velocity(130.0, eu.MUON_MASS_MEV)
    only_e = eu.murca_neutron_branch(p_n, p_p, v_e, MN, MP, T9)
    both = eu.murca_neutron_branch(p_n, p_p, v_e + v_mu, MN, MP, T9)
    assert both / only_e == pytest.approx((v_e + v_mu) / v_e, rel=1e-14)
    assert 0.5 < v_mu < 0.9 and v_e == pytest.approx(
        1.0, abs=1e-5
    )  # the muons are only mildly relativistic


def test_each_lepton_has_its_own_threshold_in_the_proton_branch():
    p_n, p_p = 450.0, 130.0
    mu_e = np.sqrt(130.0**2 + eu.ELECTRON_MASS_MEV**2)
    mu_mu = np.sqrt(50.0**2 + eu.MUON_MASS_MEV**2)
    assert eu.murca_proton_branch(p_n, p_p, 130.0, mu_e, MN, MP, T9) > 0  # 130 + 390 > 450
    assert (
        eu.murca_proton_branch(p_n, p_p, 50.0, mu_mu, MN, MP, T9) == 0.0
    )  # 50 + 390 < 450: closed
    assert eu.murca_proton_branch(p_n, p_p, 70.0, mu_mu, MN, MP, T9) > 0  # 70 + 390 > 450: open


def test_the_muon_channels_add_emissivity_only_where_there_are_muons(bsk22_model):
    cells = bsk22_model.background.cells
    state = LocalState(1e9, cells, superfluid=False)
    full = eu.make_processes()["modified_urca"](state)
    stripped = copy.copy(state)
    stripped.kf_mu = np.zeros_like(state.kf_mu)
    no_muons = eu.make_processes()["modified_urca"](stripped)
    muons = cells.ymu > 0
    assert muons.any() and np.all(full[muons] > no_muons[muons])
    np.testing.assert_array_equal(full[~muons], no_muons[~muons])
    gain = full[muons] / no_muons[muons] - 1.0
    assert (
        0.3 < gain.max() < 1.2
    )  # up to ~+87% where the muons are most abundant (centre of BSk22 1.4 Msun)
    assert np.all(np.diff(gain) <= 1e-9)  # the gain falls outwards with the muon fraction


# ------------------------------------------------------------------ scalings and regions
def test_the_rates_scale_as_t8_and_with_the_effective_masses():
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


def test_the_updated_rates_apply_to_the_core_only(bsk22_model):
    cells = bsk22_model.background.cells
    state = LocalState(1e9, cells, superfluid=False)
    crust = cells.xh > 0
    updated = eu.make_processes()
    for name in ("modified_urca", "brems_np", "brems_pp"):
        assert np.all(updated[name](state)[crust] == 0.0), name
    # the crust keeps the standard neutron bremsstrahlung
    np.testing.assert_array_equal(
        updated["brems_nn"](state)[crust], standard.neutron_neutron_bremsstrahlung(state)[crust]
    )
    assert np.all(updated["brems_nn"](state)[~crust] > 0)


def test_the_pairing_reduction_factors_are_those_of_the_standard_rates(bsk22_model):
    cells = bsk22_model.background.cells
    normal = LocalState(1e9, cells, superfluid=False)
    paired = LocalState(2e8, cells, superfluid=True)
    updated = eu.make_processes()
    f = standard.process_factors(paired)
    core = cells.xh == 0
    # a single factor multiplies the np bremsstrahlung: the ratio to the unpaired rate at the same T
    same_T = LocalState(2e8, cells, superfluid=False)
    ratio = updated["brems_np"](paired)[core] / updated["brems_np"](same_T)[core]
    np.testing.assert_allclose(ratio, f.Rnp[core], rtol=1e-12)
    urca = updated["modified_urca"](paired)[core] / updated["modified_urca"](same_T)[core]
    assert (
        np.all(urca <= 1.0 + 1e-12) and urca.min() < 0.1
    )  # the pairing suppresses the modified Urca strongly
    assert normal.core.any()


def test_in_the_star_the_updated_rates_are_below_the_standard_ones(bsk22_model):
    cells = bsk22_model.background.cells
    state = LocalState(1e9, cells, superfluid=False)
    core = state.core
    updated = {n: f(state) for n, f in eu.make_processes().items()}
    ratio = {
        n: updated[n][core] / getattr(standard, _STANDARD_NAMES[n])(state)[core] for n in updated
    }
    assert 0.4 < np.mean(ratio["modified_urca"]) < 0.9
    assert (
        0.05 < np.mean(ratio["brems_nn"]) < 0.35
    )  # includes the normalization of the standard nn rate
    assert 0.4 < np.mean(ratio["brems_np"]) < 0.9


_STANDARD_NAMES = dict(
    modified_urca="modified_urca_emissivity",
    brems_nn="neutron_neutron_bremsstrahlung",
    brems_np="neutron_proton_bremsstrahlung",
    brems_pp="proton_proton_bremsstrahlung",
)


# ------------------------------------------------------------------ the registry
def test_the_registry_is_replaced_and_restored():
    before = list(NEUTRINO_PROCESSES.items())
    with eu.updated_rates():
        inside = dict(NEUTRINO_PROCESSES.items())
        assert NEUTRINO_PROCESSES.names() == tuple(name for name, _ in before)  # same order
        for name in ("modified_urca", "brems_nn", "brems_np", "brems_pp"):
            assert inside[name] is not dict(before)[name]
        assert inside["pbf_core"] is dict(before)["pbf_core"]  # the other processes are untouched
    assert list(NEUTRINO_PROCESSES.items()) == before


def test_without_the_suppression_factors_the_rates_are_larger_by_the_betas():
    state = at_reference()
    with_beta, without = eu.make_processes(True), eu.make_processes(False)
    for name, beta in (("brems_nn", 0.56), ("brems_np", 0.66), ("brems_pp", 0.70)):
        assert without[name](state)[0] == pytest.approx(with_beta[name](state)[0] / beta, rel=1e-12)
    assert without["modified_urca"](state)[0] == pytest.approx(
        with_beta["modified_urca"](state)[0] / 0.68, rel=1e-12
    )


# ------------------------------------------------------------------ the driver flag and the comparison script
def test_the_comparison_script_runs_variants_and_restores_the_registry(background):
    before = list(NEUTRINO_PROCESSES.items())
    runs = ce.compare(
        background, final_years=300.0, variants=["standard", "updated", "updated, no muon channels"]
    )
    assert list(runs) == ["standard", "updated", "updated, no muon channels"]
    assert list(NEUTRINO_PROCESSES.items()) == before
    at = np.array([100.0])
    assert (
        abs(
            ce.luminosity_at(runs["updated"], at)[0] / ce.luminosity_at(runs["standard"], at)[0]
            - 1.0
        )
        > 0.05
    )
    text = ce.summary(runs, ages=(1, 2))
    assert "updated [%]" in text and "L_inf standard" in text
    ages = np.array([0.001, 1.0, 1e6])
    values = ce.luminosity_at(runs["standard"], ages)
    assert np.isnan(values[0]) and np.isnan(values[2]) and np.isfinite(values[1])


def test_an_unknown_variant_is_refused():
    with pytest.raises(ValueError, match="unknown variant"):
        ce.configure("nonsense")
    eu.use_standard_rates()


def test_the_normal_matter_input_has_the_superfluidity_switched_off(root):
    text = (root / "input_bsk26").read_text()
    off = ce.switch_off_superfluidity(text).split("\n")
    original = text.split("\n")
    changed = [(a, b) for a, b in zip(original, off) if a != b]
    assert [a[0] for a, _ in changed] == ["1", "1", "1"] and [b[0] for _, b in changed] == [
        "0",
        "0",
        "0",
    ]
    assert [a[1:] for a, _ in changed] == [b[1:] for _, b in changed]  # the gap models are kept
    with ce.normal_matter_input(root / "input_bsk26") as temporary:
        assert temporary.exists() and temporary.parent == root
        assert (
            ce.switch_off_superfluidity((root / "input_bsk26").read_text()) == temporary.read_text()
        )
    assert not temporary.exists()
