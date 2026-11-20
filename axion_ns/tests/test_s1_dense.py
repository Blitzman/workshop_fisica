"""S1's ln p star and dense state table, and N1's reader: nodes, exact state, interpolation accuracy, stacked evaluators vs SciPy, round trip
(cross-module on purpose)."""
import numpy as np
from scipy.interpolate import CubicSpline

from N1_modes_cowling.star_table import StarTable, _EOS, _FastSpline
from S1_star_background import star as S
from S1_star_background.dense import R_MIN, dense_state, table_nodes, validate_dense
from S1_star_background.export import write_eos, write_star
from S1_star_background.star_logp import StarLogP


def test_table_nodes_are_fine_near_the_surface_and_include_the_hand_over_node():
    st = StarLogP(None)
    r, k = table_nodes(st)
    assert r[0] == R_MIN and np.all(np.diff(r) > 0) and r[k] == st.R_stop
    d = np.diff(r[:k + 1])
    assert d[-1] < 1e-4 and d[-1] < np.median(d) / 5          # last interior spacing below 10 cm (the old grid: 19 cm) and far below the typical spacing
    assert r[-1] == st.R_stop + 5.0                            # no-axion star: a short vacuum tail


def test_fastspline_and_stacked_eos_match_scipy():
    rng = np.random.default_rng(1)
    x = np.sort(rng.uniform(0, 10, 300))
    Y = np.column_stack([np.sin(x), x**2, np.exp(-x)])
    fs, cs = _FastSpline(x, Y), CubicSpline(x, Y, axis=0)
    for r in rng.uniform(0, 10, 200):
        assert np.allclose(fs(r), cs(r), rtol=1e-12, atol=1e-13)
    e = S._eos
    eos = _EOS(e.n, e.e, e.p)
    for lp in rng.uniform(np.log(e.p[0]) + 1e-6, np.log(e.p[-1]) - 1e-6, 300):
        w = eos(lp)
        assert np.allclose(w, [e._le(lp), e._ln(lp), e._le(lp, 1)], rtol=1e-11, atol=1e-13)


def test_noaxion_star_state_at_r_is_consistent_and_the_table_round_trips(tmp_path):
    st = StarLogP(None)
    # the exact state at the hand-over radius is the end of the structure integration; the tail is Schwarzschild vacuum
    m, p, Phi, t, tp = st.state_at_r(np.array([st.R_stop - 1e-6, st.R_stop, st.R_stop + 1.0]))
    assert abs(m[1] / st.M_stop - 1) < 1e-12 and m[2] == st.M_stop and p[2] == 0.0 and p[1] > 0
    cols, k = dense_state(st)
    assert np.isneginf(cols['ln_p'][k + 1:]).all() and np.isfinite(cols['ln_p'][:k + 1]).all()
    err = validate_dense(st, cols, k)
    assert max(err.values()) <= 1e-5, err
    write_eos(tmp_path, 'test')
    meta, dense = write_star(tmp_path, st, 'star_noaxion', 'test')
    assert dense['interp_max_err'] <= 1e-5
    tab = StarTable(tmp_path / 'eos_bsk26.csv', tmp_path / 'star_noaxion.csv')
    assert tab.R == st.R_stop and abs(tab.M / st.M_stop - 1) < 1e-14
    rng = np.random.default_rng(4)
    for r in rng.uniform(0.05, st.R_stop - 0.05, 40):
        m, p, Phi, t, tp = tab.y_at(r)
        ex = st.state_at_r(r)
        assert abs(m / ex[0] - 1) < 1e-6 and abs(p / ex[1] - 1) < 1e-6, r
