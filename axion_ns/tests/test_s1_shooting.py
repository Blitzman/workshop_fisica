"""S1 shooting: log-offset interior + Brent match onto the decaying exterior branch (star_logp.StarLogP)."""
import math

import numpy as np
import pytest

from S1_star_background.star_logp import StarLogP


def test_matching_radius_does_not_move_the_solution(star_12):
    """theta_c is a property of the star, not of where the decaying branch is imposed (5 vs 7 reduced Compton wavelengths)."""
    a, b = star_12, StarLogP(12.0, match=7.0)
    assert abs(a.delta_c / b.delta_c - 1) < 1e-8
    assert abs(a.R - b.R) < 1e-7 and abs(a.M / b.M - 1) < 1e-9
    assert abs(a.decay_mismatch) < 1e-5 and abs(b.decay_mismatch) < 1e-5


def test_high_frequency_star_keeps_delta_c_and_a_decaying_tail():
    """At 50 kHz theta_c rounds to pi in double precision; delta_c = exp(u_c) is still resolved (lab notebook: 2.48e-27 rad)."""
    st = StarLogP(50.0)
    assert 1e-28 < st.delta_c < 1e-26
    r = np.linspace(st.R_stop * (1 + 1e-9), st.r_end, 400)
    m, p, Phi, t, tp = st.state_at_r(r)
    assert np.all(t > 0) and np.all(tp < 0)                     # positive, monotonically decaying to r_end
    assert t[-1] < 1e-6 * t[0]
    assert st.r_end == pytest.approx(st.R + st.tail / st.ma)
    # theta' is continuous across the matching radius (nonlinear solution inside, decaying branch outside)
    rm = st.r_match
    lo, hi = st.state_at_r(rm * (1 - 1e-9)), st.state_at_r(rm * (1 + 1e-9))
    assert abs(hi[4] / lo[4] - 1) < 1e-5 and abs(hi[3] / lo[3] - 1) < 1e-7


def test_no_condensate_below_threshold():
    with pytest.raises(ValueError):
        StarLogP(2.5)


def test_seed_only_narrows_the_bracket(star_6):
    ref = star_6
    seeded = StarLogP(6.0, seed=ref.u_c + 0.1)
    assert seeded.nint < ref.nint
    assert abs(seeded.delta_c / ref.delta_c - 1) < 1e-10
    assert math.isfinite(seeded.theta_match)
