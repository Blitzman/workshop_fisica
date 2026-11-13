"""Regression checks for run-bound diagnostics and muonic transport masses."""
from dataclasses import replace

import numpy as np
import pytest

import emissivity_updated as updated
from background import derived_fermi_quantities
from compare_emissivities import compare
from config import CoolingConfig
from cooling import run_cooling
from effective_mass import effective_mass_ratios
from emissivity import NEUTRINO_PROCESSES
from energy_balance import energy_budget
from heat_capacity import HEAT_CAPACITY_TERMS


def test_energy_budget_keeps_run_physics_after_registry_changes(background):
    cv_before = HEAT_CAPACITY_TERMS.items()
    q_before = NEUTRINO_PROCESSES.items()
    try:
        # A registered Cv contribution must be remembered too, not just Q.
        HEAT_CAPACITY_TERMS.register('review_extra')(lambda s: np.full_like(s.T, 1e18))
        with updated.updated_rates():
            run = run_cooling(background, replace(CoolingConfig(), final_years=300))
            before = energy_budget(run, background)
        HEAT_CAPACITY_TERMS.restore(cv_before)
        q_active = NEUTRINO_PROCESSES.items()
        after = energy_budget(run, background)
        np.testing.assert_array_equal(after.as_columns()[1], before.as_columns()[1])
        assert NEUTRINO_PROCESSES.items() == q_active
        assert HEAT_CAPACITY_TERMS.items() == cv_before
        with pytest.raises(ValueError, match='snapshot'):
            energy_budget(replace(run, neutrino_processes=None), background)
    finally:
        HEAT_CAPACITY_TERMS.restore(cv_before)
        NEUTRINO_PROCESSES.restore(q_before)


def test_nested_context_and_comparison_restore_custom_registry(background):
    before = NEUTRINO_PROCESSES.items()
    try:
        NEUTRINO_PROCESSES.register('review_extra')(lambda s: np.zeros_like(s.T))
        custom = NEUTRINO_PROCESSES.items()
        with updated.updated_rates(False):
            outer = NEUTRINO_PROCESSES.items()
            with pytest.raises(RuntimeError):
                with updated.updated_rates():
                    NEUTRINO_PROCESSES.unregister('review_extra')
                    raise RuntimeError('test restoration')
            assert NEUTRINO_PROCESSES.items() == outer
            compare(background, final_years=.03, variants=['standard', 'updated'])
            assert NEUTRINO_PROCESSES.items() == outer
            with pytest.raises(ValueError):
                compare(background, final_years=.03, variants=['unknown'])
            assert NEUTRINO_PROCESSES.items() == outer
        assert NEUTRINO_PROCESSES.items() == custom
    finally:
        NEUTRINO_PROCESSES.restore(before)


def test_transport_masses_use_actual_nucleon_composition():
    rho = np.array([3e14, 5e14])
    ye = np.array([.03, .05])
    yp = ye + np.array([.01, .03])
    yn = 1 - yp
    result = derived_fermi_quantities(rho, ye, yn, yp)
    expected = effective_mass_ratios(rho / 1.66e-24 * yn, rho / 1.66e-24 * yp)
    for key, value in zip(('transport_mstar_n', 'transport_mstar_p'), expected):
        np.testing.assert_array_equal(result[key], value)
    wrong_p = effective_mass_ratios(rho / 1.66e-24 * (1 - ye), rho / 1.66e-24 * ye)[1]
    assert not np.allclose(result['transport_mstar_p'], wrong_p)
