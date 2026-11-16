"""Vectorized B=0 neutrino emissivities, local CGS units.

Each process is a function ``process(state)`` registered in
``NEUTRINO_PROCESSES``; the total emissivity is the sum over the registry.
``state`` is a ``state.LocalState`` (temperature ``state.T`` [K] and the
background fields at the points of interest). A process returns an array of the
same shape as ``state.T`` in erg cm^-3 s^-1, zero where it does not apply.
Processes work on masked arrays; there is no cell loop.

To test new physics, register an additional process (or replace one)::

    @NEUTRINO_PROCESSES.register("my_process")
    def my_process(state):
        return 1e20 * (state.T * 1e-9) ** 8 * (state.rho > 2.8e14)

Processes that need the pairing-reduction factors get them from
``process_factors(state)``. These factors depend on T/Tc, not on T,
so Q(T) and Q(1.01 T) at fixed pairing share them (``state.with_temperature``).
Several constants are single-precision values and some exponentials are floored at
exp(-150), as the formulas they implement require.
"""

import numpy as np
from numpy import sqrt, exp, log10, arctan, pi

import emissivity_factors as sf
from constants import (
    ALPHA_FS,
    BOLTZMANN_MEV_PER_K,
    ELECTRON_MASS_MEV,
    HBARC_MEV_FM,
    MUON_MASS_MEV,
    N0_FM3,
    N_NU,
    NEUTRON_MASS_MEV,
    RHO0,
)
from registry import Registry
from superfluid import dimensionless_pairing_gap

NEUTRINO_PROCESSES = Registry("neutrino process")

THIRD = 1.0 / 3.0
# hbar c (3 pi^2 n0)^(1/3) [MeV] for n0 = 0.16 fm^-3 (decision D-014): the energy that the
# electron factor (Ye rho / RHO0)^(1/3) of the Direct Urca stands for. The exact coefficient of the
# Direct Urca rate gives 4.0e27 for mu = 330 MeV, so this value (331.4 MeV) is consistent with it.
DIRECT_URCA_REFERENCE_ENERGY = HBARC_MEV_FM * (3.0 * np.pi**2 * N0_FM3) ** THIRD
# Exponents of the density in the electron-proton bremsstrahlung (EP_EXP) and in the pair process (PAIR_EXP).
EP_EXP = -2.0 / 3.0
PAIR_EXP = 1.0 / 3.0
# Coefficients of the fit of the pair-annihilation neutrino process: Itoh et al. (1996), ApJS 102, 411, Eqs. (2.1)-(2.10) and Table 1
# (checked in the audit: every coefficient agrees). The fit takes rho Z/A with Z/A the electron fraction: in the core cells of the star
# Z/A equals Y_e, but in the crust cells with free neutrons it is larger (1/x_h times Y_e), where the process is negligible. The process
# matters only in the first months: switching it off changes L by -77 % at 0.03 yr, -2 % at 0.3 yr and less than 0.02 % after 1 yr (item NU-14, kept in D-039).
PAIR_CONSTANTS = (
    5.5924,
    1.2383,
    -0.8141,
    4.9924,
    13.04,
    918.6,
    1.0050,
    0.3967,
    10.7480,
    -0.3,
    0.840766,
    0.090766,
)


class ProcessFactors:
    """Urca threshold factor, pairing-reduction factors and PBF functions of one state.

    ``urca_alpha``; ``Rmn``, ``Rmp`` (modified Urca), ``Rnp``, ``Rpp``, ``Rnn``
    (bremsstrahlung), ``Rdu`` (direct Urca); ``Fc``, ``Fn``, ``Fp`` (pair
    breaking and formation). Computed once per state, from T/Tc only.
    """

    def __init__(self, s):
        self.s = s
        taus = (s.tau_crust, s.tau_n, s.tau_p)
        self.vc, self.vn, self.vp = [self.gap(t, k) for t, k in zip(taus, (1, 2, 1))]
        # alpha(n) = 1.76 - 0.63 (n0/n)^(2/3): Yakovlev and Levenfish (1995), Eq. (12), checked in the paper by the literature audit. It is
        # the density factor of the legacy rate only (modified_urca_yl95_emissivity). The standard rate (D-038) uses the constant alpha_n =
        # 1.13 of Yakovlev et al. (2001), Eq. (140); here alpha > 0 only selects the cells where the modified Urca acts and where its pairing
        # factors are computed (the 53 core cells of the example star; alpha vanishes at n_n = 0.034 fm^-3, below the core).
        # The regulariser 1e-20 and the mask alpha > 0 are choices of this code.
        self.urca_alpha = 1.76 - 0.63 / (s.rho / RHO0 * s.yn + 1e-20) ** (2.0 / 3.0)
        m = (self.urca_alpha > 0) & (s.yp > 0)
        self.Rmn = self.combined(
            m,
            sf.modified_urca_n_neutron_triplet_factor,
            sf.modified_urca_n_proton_factor,
            sf.modified_urca_n_combined_factor,
        )
        self.Rmp = self.combined(
            m,
            sf.modified_urca_p_neutron_triplet_factor,
            sf.modified_urca_p_proton_factor,
            sf.modified_urca_p_combined_factor,
        )
        self.Rnp = self.combined(
            s.yp > 0,
            sf.brems_np_neutron_triplet_factor,
            sf.brems_np_proton_factor,
            sf.brems_np_combined_factor,
        )
        self.Rdu = self.combined(
            s.direct_urca_allowed | s.direct_urca_muon_allowed,
            sf.direct_urca_triplet_m0_factor,
            sf.direct_urca_singlet_factor,
            sf.direct_urca_combined_factor,
        )
        self.Rpp = self.single(s.yp > 0, s.tau_p, self.vp, sf.brems_pp_proton_factor)
        self.Rnn = self.single(
            (s.yn > 0) & (s.xh == 0), s.tau_n, self.vn, sf.brems_nn_neutron_triplet_factor
        )
        m = (s.yn > 0) & (s.xh > 0) & (s.tau_crust <= 1)
        self.Rnn[m] = sf.brems_nn_neutron_singlet_factor(self.vc[m])
        self.Fc = self.single(
            (s.yn > 0) & (s.xh > 0) & (s.xh < 1),
            s.tau_crust,
            self.vc,
            sf.pbf_singlet_function,
            0.0,
        )
        self.Fn = self.single(
            (s.yn > 0) & (s.xh == 0), s.tau_n, self.vn, sf.pbf_triplet_m0_function, 0.0
        )
        self.Fp = self.single(
            (s.yp > 0) & (s.xh == 0), s.tau_p, self.vp, sf.pbf_singlet_function, 0.0
        )

    @staticmethod
    def gap(tau, kind):
        v = np.zeros_like(tau)
        mask = tau <= 1.0
        v[mask] = dimensionless_pairing_gap(tau[mask], kind)
        return v

    def single(self, mask, tau, v, function, normal=1.0):
        result = np.full(self.s.shape, normal)
        active = mask & (tau <= 1.0)
        result[active] = function(v[active])
        return result

    def combined(self, mask, neutron, proton, both):
        result = np.ones(self.s.shape)
        n = self.s.tau_n <= 1.0
        p = self.s.tau_p <= 1.0
        active = mask & n & p
        result[active] = both(self.vn[active], self.vp[active])
        active = mask & n & ~p
        result[active] = neutron(self.vn[active])
        active = mask & ~n & p
        result[active] = proton(self.vp[active])
        return result


def process_factors(state):
    """The ``ProcessFactors`` of a state (computed on first use, shared at fixed T/Tc)."""
    return state.fixed_memo("neutrino_process_factors", lambda: ProcessFactors(state))


def modified_urca_yl95_emissivity(s):
    """Legacy modified Urca: the density-dependent rate of Yakovlev and Levenfish (1995), electrons only.

    Kept to reproduce the runs made before decision D-038 (``--murca yl95``, see ``use_yl95_modified_urca``). The standard rate is
    ``modified_urca_emissivity`` below.
    """
    T = s.T
    f = process_factors(s)
    out = np.zeros_like(T)
    m = (f.urca_alpha > 0) & (s.yp > 0)
    rho0 = s.rho[m] / RHO0
    ye = s.ye[m]
    yp = s.yp[m]
    yn = s.yn[m]
    t9 = T[m] * 1e-9
    # Prefactors 8.55e21 (neutron branch) and 8.53e21 (proton branch), beta = 0.68 and the factor 1 - (Ye/Yp)^(1/3)/4 of the proton
    # branch: Yakovlev and Levenfish (1995), Eqs. (14) and (16), after Friman and Maxwell (1979), Eqs. (61)-(65c) (checked in the audit
    # of the emissivities). The proton branch is the angular approximation of YL95 that Yakovlev et al. (2001) correct after their
    # Eq. (142). Y01 give 8.1e21 with alpha_n beta_n = 1.13 x 0.68, n_p in place of rho Ye, a proton branch equal to the neutron branch
    # times (m*_p/m*_n)^2 (p_Fl + 3 p_Fp - p_Fn)^2 / (8 p_Fl p_Fp), which opens smoothly, and both lepton channels: that is the standard
    # rate since D-038 (items NU-1 to NU-3); this legacy function is not changed.
    qn = (
        8.55e21
        * s.mstar_n[m] ** 3
        * s.mstar_p[m]
        * (rho0 * ye) ** THIRD
        * t9**8
        * f.urca_alpha[m]
        * 0.68
    )
    qp = np.zeros_like(qn)
    a = ye**THIRD + 3.0 * yp**THIRD > yn**THIRD
    qp[a] = (
        8.53e21
        * s.mstar_p[m][a] ** 3
        * s.mstar_n[m][a]
        * (rho0[a] * ye[a]) ** THIRD
        * t9[a] ** 8
        * f.urca_alpha[m][a]
        * 0.68
        * np.maximum(1.0 - 0.25 * (ye[a] / yp[a]) ** THIRD, 0.0)
    )
    out[m] = qn * f.Rmn[m] + qp * f.Rmp[m]
    return out


def use_yl95_modified_urca():
    """Register the legacy Yakovlev-Levenfish (1995) modified Urca as the process ``modified_urca``."""
    NEUTRINO_PROCESSES.register("modified_urca", replace=True)(modified_urca_yl95_emissivity)


@NEUTRINO_PROCESSES.register("modified_urca")
def modified_urca_emissivity(s):
    """Modified Urca process of Yakovlev et al. (2001), neutron and proton branches, electron and muon channels (decision D-038).

    Neutron branch, Eq. (140), for each lepton l = e, mu:
        Q_n(l) = 8.1e21 (m*_n/m_n)^3 (m*_p/m_p) (n_p/n_0)^(1/3) T_9^8 alpha_n beta_n (v_Fl/c),  alpha_n = 1.13, beta_n = 0.68,
    where v_Fl/c = p_Fl / sqrt(p_Fl^2 + m_l^2 c^2) is the lepton velocity (Y01, Sect. 3.4(d), after Eq. 142: v_Fe is close to c and
    v_Fmu/c is close to (n_mu/n_e)^(1/3) in beta equilibrium). Proton branch, Eq. (142), for each lepton with its own opening condition:
        Q_p(l) = Q_n(l) (m*_p/m*_n)^2 (p_Fl + 3 p_Fp - p_Fn)^2 / (8 p_Fl p_Fp)   where p_Fl + 3 p_Fp > p_Fn, zero otherwise,
    which vanishes quadratically at the threshold. The effective masses of the background are ratios to one nucleon mass, so
    (m*_p/m*_n)^2 is the ratio of the ratios. n_0 = 0.16 fm^-3 (D-014). The cells are those of the legacy rate (core cells with
    neutrons and protons). The pairing factors are those of ProcessFactors. Not implemented: the phase space above the direct-Urca
    threshold (Gusakov 2002, Sect. 2), which Y01 Eq. (141) does not cover.
    """
    T = s.T
    f = process_factors(s)
    out = np.zeros_like(T)
    m = (f.urca_alpha > 0) & (s.yp > 0)
    n_p = s.rho[m] / RHO0 * s.yp[m]  # n_p / n_0
    t9 = T[m] * 1e-9
    base = 8.1e21 * s.mstar_n[m] ** 3 * s.mstar_p[m] * n_p**THIRD * t9**8 * 1.13 * 0.68
    mass_ratio2 = (s.mstar_p[m] / s.mstar_n[m]) ** 2
    k_n, k_p = s.kf_n[m], s.kf_p[m]
    qn = np.zeros_like(base)
    qp = np.zeros_like(base)
    for kf_lepton, mass in ((s.kf_e[m], ELECTRON_MASS_MEV), (s.kf_mu[m], MUON_MASS_MEV)):
        present = kf_lepton > 0
        k_l = np.where(present, kf_lepton, 1.0)
        p_l = HBARC_MEV_FM * k_l
        q = base * np.where(present, p_l / np.hypot(p_l, mass), 0.0)  # times v_Fl / c
        excess = k_l + 3.0 * k_p - k_n
        qn += q
        qp += np.where(present & (excess > 0), q * mass_ratio2 * excess**2 / (8.0 * k_l * k_p), 0.0)
    out[m] = qn * f.Rmn[m] + qp * f.Rmp[m]
    return out


@NEUTRINO_PROCESSES.register("brems_nn")
def neutron_neutron_bremsstrahlung(s):
    """Neutron-neutron bremsstrahlung."""
    T = s.T
    f = process_factors(s)
    out = np.zeros_like(T)
    m = s.yn > 0
    rho, yn, xh, mn, t9 = s.rho[m], s.yn[m], s.xh[m], s.mstar_n[m], T[m] * 1e-9
    alpha = np.full_like(rho, 0.59)
    fmu = np.ones_like(rho)
    cr = (xh > 0) & (xh < 1)
    fmu[cr] = yn[cr]
    kfn = s.kf_n[m][cr] * HBARC_MEV_FM
    # u = m_pi c / (2 p_F); the crust shape follows Yakovlev et al. (2001), Eqs. (88)-(89), after Friman and Maxwell (1979), Eqs. (45) and
    # (52) (checked in the audit). m_pi c^2 = 139 MeV is the rounded charged-pion value; Y01 (165) specify the neutral pion (135 MeV),
    # which changes the crust factor by up to 11 % at small p_F and L by less than 1 % (item NU-14, kept in D-039). Y01 define the filling factor as a
    # volume fraction, while yn is a baryon fraction here (audit; to be checked).
    u = 139.0 / 2.0 / kfn
    alpha[cr] = 1.0 - 3.0 * u * arctan(1.0 / u) / 2.0 + u**2 / 2.0 / (1.0 + u**2)
    alpha[xh == 1] = 0.0
    # nn bremsstrahlung: Yakovlev et al. (2001), Eq. 165: 7.5e19, alpha_nn = 0.59 (core), beta_nn = 0.56, N_nu = 3 once (D-020).
    # In the crust alpha follows the pion-exchange form above (see u).
    out[m] = (
        7.5e19
        * mn**4
        * (rho * yn / RHO0) ** THIRD
        * alpha
        * 0.56
        * N_NU
        * t9**8
        * fmu
        * f.Rnn[m]
    )
    return out


@NEUTRINO_PROCESSES.register("brems_np")
def neutron_proton_bremsstrahlung(s):
    """Neutron-proton bremsstrahlung."""
    T = s.T
    f = process_factors(s)
    out = np.zeros_like(T)
    m = s.yp > 0
    # np bremsstrahlung: Yakovlev et al. (2001), Eq. 166: 1.5e20, alpha_np = 1.06, beta_np = 0.66, N_nu = 3 once.
    out[m] = (
        1.5e20
        * (s.mstar_n[m] * s.mstar_p[m]) ** 2
        * (s.rho[m] * s.yp[m] / RHO0) ** THIRD
        * 1.06
        * 0.66
        * N_NU
        * (T[m] * 1e-9) ** 8
        * f.Rnp[m]
    )
    return out


@NEUTRINO_PROCESSES.register("brems_pp")
def proton_proton_bremsstrahlung(s):
    """Proton-proton bremsstrahlung."""
    T = s.T
    f = process_factors(s)
    out = np.zeros_like(T)
    m = s.yp > 0
    # pp bremsstrahlung: Yakovlev et al. (2001), Eq. 167: 7.5e19, alpha_pp = 0.11, beta_pp = 0.70, N_nu = 3 once.
    out[m] = (
        7.5e19
        * s.mstar_p[m] ** 4
        * (s.rho[m] * s.yp[m] / RHO0) ** THIRD
        * 0.11
        * 0.7
        * N_NU
        * (T[m] * 1e-9) ** 8
        * f.Rpp[m]
    )
    return out


@NEUTRINO_PROCESSES.register("brems_ep")
def electron_proton_bremsstrahlung(s):
    """Electron-proton bremsstrahlung in the core.

    Not changed by D-014: its reference density 2.82e14 g cm^-3 and the hard density cut below are
    kept until the source of this rate is checked (the cut switches the process off between the
    crust-core interface and 2.82e14 g cm^-3).
    """
    T = s.T
    out = np.zeros_like(T)
    # Prefactor 2.4e17, the exponent EP_EXP, the reference density and the cut at 2.82e14 g/cm^3: source not identified. Kaminker and
    # Haensel (1999), Eqs. (40)-(41) and (45), agree with it only in the T^8 dependence and the factor m*_p^2 (audit): their low-T
    # coefficient is 3.69e14 (n0/n_p)^(1/3) R_ep / y_s with a screening parameter y_s, there is no hard density cut, and the pairing
    # factor of their Eq. (41) is absent here (item NU-14, kept in D-039). Switching the process off changes L by less than 0.6 %. Not changed by D-014.
    m = (s.rho >= 2.82e14) & (s.yp > 0)
    out[m] = (
        2.4e17 * (s.rho[m] * s.yp[m] / 2.82e14) ** EP_EXP * (T[m] * 1e-9) ** 8 * s.mstar_p[m] ** 2
    )
    return out


@NEUTRINO_PROCESSES.register("pbf_core")
def core_pbf_emissivity(s):
    """Neutron and proton pair breaking and formation (PBF) in the core."""
    T = s.T
    f = process_factors(s)
    out = np.zeros_like(T)
    m = (s.xh == 0) & (s.tau_n <= 1) & (s.yn > 0)
    pf = HBARC_MEV_FM * s.kf_n[m] / NEUTRON_MASS_MEV  # p_F / (m_N c) from the background (D-014)
    # Pair breaking and formation. Rate normalisation: Yakovlev et al. (2001), Eq. (236) (1.17e21, N_nu = 3). The weak coefficients are the
    # vector-suppressed ones of Page et al. (2009), Eqs. (9)-(11): triplet 2 g_A^2 = 3.175 (Y01's bare vertex gives 1 + 2 g_A^2 = 4.175;
    # Leinson 2010, Eq. 96, and Potekhin and Chabrier 2018, Sect. 3.4, give about g_A^2/2 = 0.794, the factor 0.19 of the bare rate);
    # singlet 4/81 (p_F/m*)^4 with the proton axial term of Kaminker, Haensel and Yakovlev (1999). g_A = 1.26 and sin^2(theta_W) = 0.23 are
    # the values of Page et al. (2009) (Bottaro et al. use g_A = 1.27, which would raise the axial terms by 1.6 %). Checked in the papers
    # by the literature audit (item NU-7: the default keeps the local value). The three triplet values are: 3.175, the bare axial term
    # only (the vector contribution 1 of Y01 is omitted); 4.175, Y01's bare vertex; 0.794, the anomalous-vertex result of Leinson (2010)
    # and Potekhin and Chabrier (2018). The singlet expression is a hybrid (vector term suppressed by 4 (v_F/c)^4 / 81, after Kolomeitsev
    # and Voskresensky 2010, Eqs. 140-143, plus the axial term with 11/42), not the exact result of PC18 (second audit); the weighting
    # of the crust rate by the free-neutron volume is to be checked.
    out[m] = 1.17e21 * s.mstar_n[m] * pf * N_NU * (2.0 * 1.26**2) * f.Fn[m] * (T[m] * 1e-9) ** 7
    m = (s.xh == 0) & (s.tau_p <= 1) & (s.yp > 0)
    pf = HBARC_MEV_FM * s.kf_p[m] / NEUTRON_MASS_MEV  # p_F / (m_N c) from the background (D-014)
    mp = s.mstar_p[m]
    S = 4.0 / 81.0 * (pf / mp) ** 4
    a = (4.0 * 0.23 - 1.0) ** 2 * S + 1.26**2 * pf**2 * (1.0 + 11.0 / 42.0 / mp**2)
    out[m] += 1.17e21 * mp * pf * N_NU * a * f.Fp[m] * (T[m] * 1e-9) ** 7
    return out


@NEUTRINO_PROCESSES.register("direct_urca")
def direct_urca_emissivity(s):
    """Direct Urca process, electron and muon channels (each above its own threshold).

    Q_l = 4.0e27 m*_n m*_p T9^6 (mu_l / E_ref) R_du, with mu_l the Fermi energy of the lepton and
    E_ref = hbar c (3 pi^2 n0)^(1/3), n0 = 0.16 fm^-3: for the electrons
    mu_e / E_ref = (Ye rho / RHO0)^(1/3). The lepton enters only through mu_l and its
    threshold (see README, REFERENCES).
    """
    T = s.T
    f = process_factors(s)
    out = np.zeros_like(T)
    m = s.direct_urca_allowed
    # Electron channel. Prefactor 4.0e27 with n0 = 0.16 fm^-3: consistent with the exact coefficient of the rate (README, D-014); Y01
    # Eqs. (120)-(121) and Sect. 3.3(c) (checked in the audit). The analytic coefficient gives 4.015e27 with g_A = 1.26 and
    # cos(theta_C) = 0.974, so the rounded value is 0.37 % lower. Threshold: Lattimer et al. (1991).
    out[m] = (
        4.0e27
        * s.mstar_n[m]
        * s.mstar_p[m]
        * (T[m] * 1e-9) ** 6
        * (s.ye[m] * s.rho[m] / RHO0) ** THIRD
        * f.Rdu[m]
    )
    mu = s.direct_urca_muon_allowed
    mu_muon = np.hypot(HBARC_MEV_FM * s.kf_mu[mu], MUON_MASS_MEV)
    # Muon channel: the same expression with the muon Fermi energy (README, DIRECT URCA THRESHOLD). The rule is published in Y01, Sect.
    # 3.3(c), after Eq. (121) (audit), and was also verified numerically by the tests. In beta equilibrium mu_e = mu_mu, so a second open
    # channel doubles the rate. Item NU-8 stays as it is (D-036).
    out[mu] += (
        4.0e27
        * s.mstar_n[mu]
        * s.mstar_p[mu]
        * (T[mu] * 1e-9) ** 6
        * (mu_muon / DIRECT_URCA_REFERENCE_ENERGY)
        * f.Rdu[mu]
    )
    return out


@NEUTRINO_PROCESSES.register("brems_e_ions")
def electron_ion_bremsstrahlung(s):
    """Electron-ion bremsstrahlung in the crust (fit range 5e7-2e9 K, 1e9-1.4e14 g/cm^3)."""
    T = s.T
    out = np.zeros_like(T)
    m = (T >= 5e7) & (T <= 2e9) & (s.rho >= 1e9) & (s.rho <= 1.4e14) & (s.xh > 0)
    tau = log10(T[m] * 1e-8)
    xr = log10(s.rho[m] * 1e-12)
    # Coefficients of the fit (T = 5e7-2e9 K, rho = 1e9-1.4e14 g/cm^3), including the term 0.228 rho/2.8e14: Kaminker, Pethick, Potekhin,
    # Thorsson and Yakovlev (1999), A&A 343, 1009, Eq. (40) (p. 1021), and Yakovlev et al. (2001), Eq. (76) (checked in the audit; Haensel,
    # Kaminker and Yakovlev 1996 is the earlier liquid-phase calculation). The last term uses log10 in the papers (audit); the code had the
    # natural logarithm, which lowered the rate by 0.93, 0.50 and 0.39 at 1e13, 1e14 and 1.4e14 g/cm^3. Corrected in D-038 (item NU-13):
    # L changed by about -1 % at 100 yr and by less than 0.2 % after 841 yr.
    out[m] = 10.0 ** (
        11.204
        + 7.304 * tau
        + 0.2976 * xr
        - 0.370 * tau**2
        + 0.188 * tau * xr
        - 0.103 * xr**2
        + 0.0547 * tau**2 * xr
        - 6.77 * log10(1.0 + 0.228 * s.rho[m] / 2.8e14)
    )
    return out


@NEUTRINO_PROCESSES.register("plasma")
def plasma_emissivity(s):
    """Plasma neutrino process (plasmon decay) in the crust."""
    T = s.T
    out = np.zeros_like(T)
    m = (s.xh > 0) & (s.ye > 0)
    xtr = BOLTZMANN_MEV_PER_K * T[m] / ELECTRON_MASS_MEV
    xr = s.kf_e[m] * HBARC_MEV_FM / ELECTRON_MASS_MEV
    alpha = ALPHA_FS
    fp = sqrt(4.0 * alpha * xr**3 / 3.0 / pi / sqrt(1.0 + xr**2)) / xtr
    # Coefficients 16.23, 4.604 and 0.9248: Yakovlev et al. (2001), Eqs. (38) and (42), after Itoh et al. (1996) (checked in the audit). The
    # normalisation Q_c = G_F^2 m_e^9 is 1.0239e23 erg cm^-3 s^-1 from the constants of constants.py and 1.023e23 in Y01 Eq. (23). The code
    # had 1.203e23, two digits swapped (+17.5 %), corrected in D-038 (item NU-12): L changed by +11 % at 0.1 yr, +4 % at 1 yr and less
    # than 0.5 % after 100 yr.
    xip = xtr**9 * (16.23 * fp**6 + 4.604 * fp**7.5) * sf.exp_with_floor(-fp)
    out[m] = 1.023e23 * xip * 0.9248 / 96.0 / pi**4 / alpha
    return out


@NEUTRINO_PROCESSES.register("pbf_crust")
def crust_pbf_emissivity(s):
    """Neutron pair breaking and formation (PBF) in the crust."""
    T = s.T
    f = process_factors(s)
    out = np.zeros_like(T)
    m = (s.tau_crust <= 1) & (s.xh > 0) & (s.xh < 1) & (s.yn > 0)
    pf = HBARC_MEV_FM * s.kf_n[m] / NEUTRON_MASS_MEV  # p_F / (m_N c) from the background (D-014)
    mn = s.mstar_n[m]
    S = 4.0 / 81.0 * (pf / mn) ** 4
    # Singlet neutron weak coefficient from Page et al. (2009), Eq. (9).
    a = S + 1.26**2 * pf**2 * (1.0 + 11.0 / 42.0 / mn**2)
    out[m] = 1.17e21 * mn * pf * N_NU * a * f.Fc[m] * (T[m] * 1e-9) ** 7
    return out


@NEUTRINO_PROCESSES.register("pair")
def pair_annihilation_emissivity(s):
    """Electron-positron pair annihilation (T >= 3e8 K)."""
    T = s.T
    out = np.zeros_like(T)
    m = T >= 3e8
    if np.any(s.aa[m] <= 0):
        raise ValueError("QPAIR requires A>0 for T>=3e8 K")
    t = T[m]
    x = t / 5.9302e9
    x2 = x * x
    xi = 1.0 / x
    rza = s.rho[m] * s.zz[m] / s.aa[m]
    chi = 1e-3 * rza**PAIR_EXP * xi
    c_lo, b1_hi, b2_hi, c_hi, g13, g918, q1, q4, q10, qexp, f84, f09 = PAIR_CONSTANTS
    # np.where selects constants only; no unsafe branch expressions evaluated.
    hot = t >= 1e10
    b1 = np.where(hot, b1_hi, 0.9383)
    b2 = np.where(hot, b2_hi, -0.4141)
    b3 = np.where(hot, 0.0, 0.05829)
    c = np.where(hot, c_hi, c_lo)
    f = (
        (6.002e19 + chi * (2.084e20 + chi * 1.872e21))
        * exp(-c * chi)
        / (chi**3 + xi * (b1 + xi * (b2 + xi * b3)))
    )
    g = 1.0 - x2 * (g13 - x2 * (133.5 + x2 * (1534.0 + x2 * g918)))
    q = (q1 + q4 * sqrt(x) + q10 * x2) ** (-1.0) * (
        1.0 + rza / (7.692e7 * x**3 + 9.715e6 * sqrt(x))
    ) ** qexp
    out[m] = (f84 + f09 * q) * g * exp(-2.0 * xi) * f
    return out


def neutrino_emissivity(state, *, processes=None):
    """Evaluate every registered process; return one array per process plus 'total'.

    ``processes`` optionally supplies an ordered snapshot of (name, callable) pairs.
    Emissivities are in erg cm^-3 s^-1, with the shape of ``state.T``.
    """
    selected = NEUTRINO_PROCESSES.items() if processes is None else processes
    result = {name: function(state) for name, function in selected}
    if any(not np.isfinite(v).all() or np.any(v < 0) for v in result.values()):
        raise FloatingPointError("Nonfinite or negative neutrino emissivity")
    result["total"] = sum(result.values())
    return result


def neutrino_emissivity_total(state, *, processes=None):
    """Total neutrino emissivity [erg cm^-3 s^-1]."""
    return neutrino_emissivity(state, processes=processes)["total"]


def linearized_neutrino_source(state):
    """Internal offset and slope of the neutrino source, holding T/Tc fixed in the 1% difference.

    ``state`` is the state at the old temperature of the cells. The source is
    linearized as offset + slope * u, with u = T_inf / 1e8 K.
    """
    phi = state.ephi
    q = neutrino_emissivity(state)["total"]
    q1 = neutrino_emissivity(state.with_temperature(1.01 * state.T))["total"]
    source_slope = (q1 - q) / (0.01 * (state.T / 1e8)) * 1e15 * 1e-40 * phi
    source_offset = (q - (q1 - q) / 0.01) * 1e15 * 1e-40 * phi**2
    return source_offset, source_slope
