"""S1's ln p star and dense state table, and the shared reader (common.star_table): nodes, exact state at the hand-over, stacked evaluators vs
SciPy, metadata round trip. The interpolation error of the table is checked only by S1 itself, when it writes the table (dense.validate_dense)."""
import numpy as np
from scipy.interpolate import CubicSpline

from common.star_table import _EOS, _FastSpline
from S1_star_background import star as S
from S1_star_background.dense import R_MIN, dense_state, table_nodes


def test_table_nodes_are_fine_near_the_surface_and_include_the_hand_over_node(star_noaxion):
    st = star_noaxion
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


def test_noaxion_star_hand_over_and_metadata_round_trip(star_noaxion, table_noaxion):
    st, (tab, _) = star_noaxion, table_noaxion
    # the exact state at the hand-over radius is the end of the structure integration; the tail is Schwarzschild vacuum
    m, p, Phi, t, tp = st.state_at_r(np.array([st.R_stop - 1e-6, st.R_stop, st.R_stop + 1.0]))
    assert abs(m[1] / st.M_stop - 1) < 1e-12 and m[2] == st.M_stop and p[2] == 0.0 and p[1] > 0
    cols, k = dense_state(st)
    assert np.isneginf(cols['ln_p'][k + 1:]).all() and np.isfinite(cols['ln_p'][:k + 1]).all()
    assert tab.R == st.R_stop and abs(tab.M / st.M_stop - 1) < 1e-14
