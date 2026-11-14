"""Local thermodynamic state shared by every registered piece of physics.

A ``LocalState`` is the temperature at a set of points together with the
background fields at those points. Neutrino processes, conductivity channels and
heat-capacity terms are all functions ``f(state)``, so a new ingredient is
written the same way whatever it computes::

    @NEUTRINO_PROCESSES.register("my_process")        # emissivity [erg cm^-3 s^-1]
    def my_process(state):
        return 1e20 * (state.T * 1e-9) ** 8 * (state.rho > 2.8e14)

Attributes (arrays of the common shape ``state.shape``)
    T                  local temperature [K]
    rho, ye, yn, yp, xh, aa, zz, ephi, tc_crust, tc_n, tc_p, kf_n, kf_e,
    mstar_n, mstar_p, transport_mstar_n, transport_mstar_p
                       the background fields (see ``background.Points``)
    ymu, kf_p, kf_mu   muon fraction and the proton and muon Fermi momenta [fm^-1]
                       (charge neutrality: n_p = n_e + n_mu)
    n_n, n_p, n_e, n_mu  number densities [cm^-3] = nb * fraction
    direct_urca_allowed  mask of the Direct Urca threshold in the core: kf_n <= kf_p + kf_e
    direct_urca_muon_allowed  the same for the muon channel: kf_n <= kf_p + kf_mu
    nb                 baryon number density rho / m_u [cm^-3]
    crust, core        masks: crust = nuclei present (xh > 0), core = the rest
    tau_crust, tau_n, tau_p   T/Tc of crust neutrons, core neutrons and core
                       protons (10 everywhere if superfluidity is off or Tc = 0)
    superfluid         whether pairing is active

``state.T`` may have more dimensions than the background (for example a grid of
temperatures against the cells): everything is broadcast to a common shape.
"""

import copy
from functools import cached_property

import numpy as np

from background import POINT_FIELDS
from numerical_validation import checked


class LocalState:
    """Temperature plus background at the same points; validated once, used by all physics."""

    def __init__(self, T_K, points, superfluid=True):
        arrays = np.broadcast_arrays(
            checked(T_K, "T", True), *[getattr(points, name) for name in POINT_FIELDS]
        )
        self.shape = arrays[0].shape
        self.T = arrays[0]
        for name, array in zip(POINT_FIELDS, arrays[1:]):
            setattr(self, name, array)
        self.superfluid = bool(superfluid)
        self._memo = {}  # intermediate results that depend on T
        self._fixed = {}  # intermediate results that depend only on the background and T/Tc

    @cached_property
    def nb(self):
        # 1.66e-24 g: the mass unit of the microphysics (constants.MICROPHYSICS_MASS_UNIT).
        return self.rho / 1.66e-24

    @cached_property
    def n_n(self):
        return self.nb * self.yn

    @cached_property
    def n_p(self):
        return self.nb * self.yp

    @cached_property
    def n_e(self):
        return self.nb * self.ye

    @cached_property
    def n_mu(self):
        return self.nb * self.ymu

    @cached_property
    def direct_urca_allowed(self):
        """Direct Urca n -> p e nu is kinematically allowed: momentum conservation needs
        kf_n <= kf_p + kf_e (homogeneous matter only). In npe matter this is Yp >= 1/9."""
        return self.core & (self.yp > 0) & (self.ye > 0) & (self.kf_n <= self.kf_p + self.kf_e)

    @cached_property
    def direct_urca_muon_allowed(self):
        """Direct Urca n -> p mu nu is allowed: kf_n <= kf_p + kf_mu (homogeneous matter only).
        In beta equilibrium kf_mu < kf_e, so this channel opens at a higher density than the electron one.
        """
        return self.core & (self.yp > 0) & (self.ymu > 0) & (self.kf_n <= self.kf_p + self.kf_mu)

    @cached_property
    def crust(self):
        return self.xh > 0

    @cached_property
    def core(self):
        return ~self.crust

    def _tau(self, tc):
        if not self.superfluid:
            return np.full_like(self.T, 10.0)
        return self.T / np.maximum(tc, 0.1 * self.T)

    @cached_property
    def tau_crust(self):
        return self._tau(self.tc_crust)

    @cached_property
    def tau_n(self):
        return self._tau(self.tc_n)

    @cached_property
    def tau_p(self):
        return self._tau(self.tc_p)

    def with_temperature(self, T_K):
        """The same state at another temperature, keeping T/Tc fixed.

        Used for finite differences in T "at fixed pairing" (the neutrino-source
        slope). The temperature-independent intermediate results are shared.
        """
        self.tau_crust, self.tau_n, self.tau_p  # pin T/Tc before copying
        new = copy.copy(self)
        new.T = np.broadcast_to(checked(T_K, "T", True), self.shape)
        new._memo = {}
        return new

    def memo(self, key, function):
        """Per-state cache for expensive intermediates that depend on the temperature."""
        if key not in self._memo:
            self._memo[key] = function()
        return self._memo[key]

    def fixed_memo(self, key, function):
        """Cache shared with ``with_temperature`` copies: only for results that do not depend on T."""
        if key not in self._fixed:
            self._fixed[key] = function()
        return self._fixed[key]
