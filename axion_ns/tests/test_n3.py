"""N3: fast tests of the table-based pieces (no symbolic compile; the 12-star comparison with the legacy is `python -m tests.n3_gate`)."""
import numpy as np
import pytest

from tests.conftest import ROOT


@pytest.fixture(scope='module')
def tables(tmp_path_factory):
    """S1's ln p star at 6 kHz written to a table, read by N3."""
    from S1_star_background.export import write_eos, write_star
    from S1_star_background.star_logp import StarLogP
    d = tmp_path_factory.mktemp('s1')
    st = StarLogP(6.0)
    write_eos(d, 'test')
    meta, dense = write_star(d, st, 'star_nua6', 'test')
    from N3_modes_fullgr_axion.star_table import StarTable
    return st, StarTable(d / 'eos_bsk26.csv', d / 'star_nua6.csv'), dense


def test_s1_table_interpolation_error_is_small(tables):
    _, _, dense = tables
    assert dense['interp_max_err'] <= 1e-5, dense['interp_err']


def test_table_state_matches_the_exact_state(tables):
    st, tab, _ = tables
    rng = np.random.default_rng(2)
    r = np.sort(np.concatenate([rng.uniform(0.02, st.R_stop, 60), st.R_stop + rng.uniform(0.01, 30, 20)]))
    ex = st.state_at_r(r)
    scale = np.abs(ex).max(axis=1)
    for k, rr in enumerate(r):
        got = np.array(tab.y_at(rr))
        # relative to the larger of the local value and 1% of the column maximum (m, Phi and theta' vanish at the centre; p is 0 in the tail)
        assert np.all(np.abs(got - ex[:, k]) <= 1e-5 * np.maximum(np.abs(ex[:, k]), 1e-2 * scale)), (rr, got, ex[:, k])


def test_analytic_phi_second_derivative_matches_a_finite_difference(tables):
    from N3_modes_fullgr_axion.background import bg18
    st, tab, _ = tables
    h = 1e-4
    dphi = lambda x: tab.rhs(x, list(tab.y_at(x)))[2]
    for r in (0.3, 2.0, 5.0, 9.0, 11.0, 11.3):
        m, p, Phi, t, tp = tab.y_at(r)
        an = bg18(tab, r, m, p, Phi, t, tp)[2][0]
        fd = (dphi(r + h) - dphi(r - h)) / (2 * h)
        assert abs(an / fd - 1) < 2e-5, (r, an, fd)


def test_bg18_is_consistent_with_the_star_reader(tables):
    from N3_modes_fullgr_axion.background import bg18
    st, tab, _ = tables
    for r in (0.5, 4.0, 8.0, 11.0):
        m, p, Phi, t, tp = tab.y_at(r)
        b = bg18(tab, r, m, p, Phi, t, tp)
        dm, dp, dPhi, dt, dtp = tab.rhs(r, [m, p, Phi, t, tp])
        assert abs(b[1][0] / dPhi - 1) < 1e-12 and abs(b[8][0] / tp - 1) < 1e-12
        bg = tab.bg(r)                                   # (e, p, n, cs2, q1, q2, dPhi, dlam, Phi + shift, lam, t, tp)
        assert abs(b[4][0] / bg[7] - 1) < 1e-10 and abs(b[3][0] / bg[9] - 1) < 1e-10
