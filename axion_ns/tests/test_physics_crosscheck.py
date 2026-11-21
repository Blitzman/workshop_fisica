"""The shared reader (common.star_table) re-implements the two short background functions (matter, rhs). They must agree with S1's Star
(cross-module on purpose)."""
import numpy as np
import pytest


@pytest.mark.parametrize('which', ['table_noaxion', 'table_6'])
def test_rhs_and_matter_agree_with_S1(which, request):
    from S1_star_background import star as S
    st, _ = request.getfixturevalue(which)
    ref = object.__new__(S.Star)
    ref.nu_a, ref.p0, ref.E0, ref.nc, ref.ma, ref.fa2 = st.nu_a, st.p0, st.E0, st.nc, st.ma, st.fa2
    rng = np.random.default_rng(3)
    for r in rng.uniform(0.05, 0.98 * st.R, 60):
        y = st.y_at(r)
        a, b = np.array(st.rhs(r, y)), np.array(ref.rhs(r, list(y)))
        # the reader evaluates the EOS with a stacked polynomial (S1: scipy splines); matter() includes a spline derivative (c_s^2) and
        # cancellations at the condensate wall, so last-bit differences are amplified (seen 1.2e-12): tolerance 1e-10
        assert np.allclose(a, b, rtol=1e-10, atol=0), (r, a, b)
        assert np.allclose(st.matter(y[1], y[3]), ref.matter(y[1], y[3]), rtol=1e-10, atol=0)
