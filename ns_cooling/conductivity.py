"""Isotropic B=0 thermal conductivity; local temperatures in Kelvin.

Each conductivity channel is a function ``channel(state)`` registered in
``CONDUCTIVITY_CHANNELS``; the total is the sum over the registry. ``state`` is a
``state.LocalState`` (temperature and background at the points of interest). A
channel returns kappa in CGS [erg cm^-1 s^-1 K^-1] with the shape of
``state.T``, zero where it does not apply.

To test new physics, register an additional channel (or replace one)::

    @CONDUCTIVITY_CHANNELS.register("my_channel")
    def my_channel(state):
        return np.where(state.crust, 1e18 * state.T, 0.0)
"""

import numpy as np

from conductivity_core import core_conductivity
from conductivity_crust import crust_electron_conductivity
from registry import Registry

CONDUCTIVITY_CHANNELS = Registry("conductivity channel")


def core_components(state):
    """Electron, muon and neutron conductivities on the core points (None without a core).

    Computed once per state and shared by the channels that need them.
    """

    def compute():
        c = state.core
        if not np.any(c):
            return None
        return core_conductivity(
            state.T[c],
            state.rho[c],
            state.ye[c],
            state.tau_n[c],
            state.tau_p[c],
            effective_masses=(state.transport_mstar_n[c], state.transport_mstar_p[c]),
            Y_mu=state.ymu[c],
        )

    return state.memo("core_conductivity", compute)


@CONDUCTIVITY_CHANNELS.register("core_electrons")
def core_electron_channel(state):
    """Electron conductivity in the core."""
    out = np.zeros(state.shape)
    components = core_components(state)
    if components is not None:
        out[state.core] = components["electrons"]
    return out


@CONDUCTIVITY_CHANNELS.register("core_muons")
def core_muon_channel(state):
    """Muon conductivity in the core (zero where there are no muons)."""
    out = np.zeros(state.shape)
    components = core_components(state)
    if components is not None:
        out[state.core] = components["muons"]
    return out


@CONDUCTIVITY_CHANNELS.register("core_neutrons")
def core_neutron_channel(state):
    """Neutron conductivity in the core."""
    out = np.zeros(state.shape)
    components = core_components(state)
    if components is not None:
        out[state.core] = components["neutrons"]
    return out


@CONDUCTIVITY_CHANNELS.register("crust_electrons")
def crust_electron_channel(state):
    """Electron conductivity in the crust (scalar fit, evaluated point by point)."""
    out = np.zeros(state.shape)
    c = state.crust
    out[c] = np.fromiter(
        (
            crust_electron_conductivity(float(t), float(r), float(a), float(z), float(x))
            for t, r, a, z, x in zip(
                state.T[c], state.rho[c], state.aa[c], state.zz[c], state.xh[c]
            )
        ),
        dtype=float,
        count=np.count_nonzero(c),
    )
    return out


def thermal_transport(state):
    """CGS conductivity of each registered channel (key ``kappa_<name>``), the total
    (``kappa_total``) and the internal coefficient (``transport_internal``).

    Crust: Potekhin's electron fit with Qimp = 0.1.
    """
    result = {f"kappa_{name}": channel(state) for name, channel in CONDUCTIVITY_CHANNELS.items()}
    total = sum(result.values())
    result = {"kappa_total": total, **result}
    result["transport_internal"] = total * 1e13 * 1e-40
    return result
