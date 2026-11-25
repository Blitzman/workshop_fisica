"""Fast controls for the extended D2 verifier and export gate."""
import pytest
import sympy as sp

from common.notation import V, a0, dth
from D2_derive_nonradial import checks_fullgr as checks


@pytest.mark.parametrize('order', (1, 2))
def test_potential_derivative_survives_residual_reduction(order):
    a = sp.Dummy('a')
    derivative = sp.diff(V(a0), a0, order)
    residual = dth * derivative
    via_subs = sp.Subs(sp.diff(V(a), a, order), a, a0)
    assert not checks.vanishes(residual)
    assert not checks.vanishes(dth * via_subs)
    assert checks.vanishes(dth * (via_subs - derivative))


def test_export_requires_all_equations_and_omission_controls(monkeypatch):
    from types import SimpleNamespace

    from D2_derive_nonradial.ext_system import NEW

    # Exercise the gate here without invoking the slow action-derived reduction.
    monkeypatch.setattr(checks, 'F', SimpleNamespace(EXT={k: sp.S.One for k in NEW}, NEW=NEW))
    monkeypatch.setattr(checks, 'implied', lambda e: e in (0, 1))
    results = checks.extended_checks()
    assert checks.extended_checks_passed(results)
    for i, (name, _, detail) in enumerate(results):
        assert not checks.extended_checks_passed(results[:i] + results[i + 1:])
        assert not checks.extended_checks_passed(results[:i] + [(name, False, detail)] + results[i + 1:])
