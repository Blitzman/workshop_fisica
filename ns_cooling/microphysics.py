"""Convenience evaluation of all the microphysics of a state (diagnostics, notebooks)."""

from conductivity import thermal_transport
from emissivity import neutrino_emissivity
from heat_capacity import heat_capacity


def evaluate_microphysics(state):
    """Heat capacity, neutrino emissivity and B=0 conductivity of a ``state.LocalState``."""
    return dict(
        heat_capacity=heat_capacity(state),
        neutrino_emissivity=neutrino_emissivity(state),
        conductivity=thermal_transport(state),
    )
