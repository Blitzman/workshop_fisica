"""S1's dense state table and N1's reader: grid, interpolation accuracy, stacked evaluators vs SciPy, round trip (cross-module on purpose)."""
import numpy as np
from scipy.interpolate import CubicSpline

from N1_modes_cowling.star_table import StarTable, _EOS, _FastSpline
from S1_star_background import star as S
from S1_star_background.dense import dense_state, surface_clustered, validate_dense
from S1_star_background.export import write_eos, write_star


def test_surface_clustered_grid():
    r = surface_clustered(1e-3, 11.0, 2000)
    assert r[0] == 1e-3 and abs(r[-1] - 11.0) < 1e-12 and np.all(np.diff(r) > 0)
    d = np.diff(r)
    assert d[-1] < d[0] / 50 and d[-1] < 1e-3      # strongly clustered toward the surface (last spacing below 1 m)


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


def test_noaxion_star_dense_table_accuracy_and_roundtrip(tmp_path):
    st = S.Star(None)                                      # default tolerances (rtol 1e-6)
    cols, rk = dense_state(st, 2000, 500)
    assert np.all(np.diff(cols['r_km']) > 0) and cols['r_km'][0] == 1e-3 and rk > st.R
    err = validate_dense(st, cols, rk)
    assert max(err.values()) <= 1e-5, err
    write_eos(tmp_path, 'test')
    meta, dense = write_star(tmp_path, st, 'star_noaxion', 'test')
    assert dense['interp_max_err'] <= 1e-5 and dense['n_rows'] == 2499
    tab = StarTable(tmp_path / 'eos_bsk26.csv', tmp_path / 'star_noaxion.csv')
    assert tab.R == st.R and tab.M == st.M
    rng = np.random.default_rng(4)
    for r in rng.uniform(1e-3, st.R - 0.5, 40):
        a, b = np.array(tab.bg(r)), np.array(st.bg(r))
        assert np.allclose(a, b, rtol=1e-4, atol=0), (r, np.abs(a / b - 1).max())
