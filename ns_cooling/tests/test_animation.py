"""Comparison animation: the common time grid, the interpolation of the profiles, the figure."""

from dataclasses import replace

import matplotlib

matplotlib.use("Agg")  # no display needed; must come before the animation module imports pyplot

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402

import cooling_animation as ca  # noqa: E402
from config import CoolingConfig  # noqa: E402
from cooling import run_cooling  # noqa: E402


@pytest.fixture(scope="module")
def eos_run(bsk22_model):
    return run_cooling(bsk22_model.background, replace(CoolingConfig(), final_years=100.0))


@pytest.fixture(scope="module")
def entries(background, short_run, bsk22_model, eos_run):
    """The default BSk26 star (to 1000 yr) and the BSk22 star (to 100 yr: it ends earlier)."""
    return [
        ("table", short_run, background.cells.rho),
        ("1.40 Msun", eos_run, bsk22_model.background.cells.rho),
    ]


def stored(run, step):
    return np.log10(run.profiles[step] * 1e8)


def test_interpolation_reproduces_the_stored_profiles_at_the_steps(short_run):
    steps = [1, 5, 40, short_run.completed_steps // 2, short_run.completed_steps]
    got = ca.profile_log_temperature(short_run, short_run.step_times_years[steps])
    for row, step in zip(got, steps):
        np.testing.assert_allclose(row, stored(short_run, step), rtol=1e-12)


def test_interpolation_is_linear_in_log_time_and_log_temperature(short_run):
    n = 30
    t0, t1 = short_run.step_times_years[[n, n + 1]]
    mid = np.sqrt(t0 * t1)  # the midpoint in log t
    got = ca.profile_log_temperature(short_run, [mid])[0]
    np.testing.assert_allclose(
        got, 0.5 * (stored(short_run, n) + stored(short_run, n + 1)), rtol=1e-12
    )


def test_times_beyond_the_end_of_a_run_give_its_final_state(short_run):
    end = short_run.step_times_years[short_run.completed_steps]
    rows = ca.profile_log_temperature(short_run, [end, 10 * end, 1e6])
    np.testing.assert_allclose(rows[1], rows[0], rtol=0, atol=0)
    np.testing.assert_allclose(rows[2], stored(short_run, short_run.completed_steps), rtol=1e-12)


def test_frame_times_cover_every_run_and_are_equally_spaced_in_log_time(
    entries, short_run, eos_run
):
    times, stars = ca.frame_data(entries, frames=50)
    first = min(r.step_times_years[r.step_times_years > 0].min() for _, r, _ in entries)
    assert times[0] == pytest.approx(first) and times[-1] == pytest.approx(
        short_run.step_times_years[-1]
    )
    np.testing.assert_allclose(np.diff(np.log(times)), np.diff(np.log(times))[0], rtol=1e-9)
    assert len(stars) == 2 and [s["label"] for s in stars] == ["table", "1.40 Msun"]
    for s, (_, run, density) in zip(stars, entries):
        assert s["log_T"].shape == (50, len(density)) and np.all(
            np.diff(s["log_rho"]) >= 0
        )  # cells sorted by density
        assert np.all(np.diff(s["count"]) >= 0) and s["count"][-1] == len(
            run.times_years
        )  # an ended run keeps its curve
    np.testing.assert_allclose(
        stars[1]["log_T"][-1],
        ca.profile_log_temperature(eos_run, [times[-1]])[0][np.argsort(entries[1][2])],
    )


def test_a_star_that_ends_earlier_stops_changing(entries, eos_run):
    times, stars = ca.frame_data(entries, frames=60)
    ended = times >= eos_run.step_times_years[eos_run.completed_steps]
    assert ended.sum() > 5
    np.testing.assert_array_equal(np.diff(stars[1]["log_T"][ended], axis=0), 0.0)


def test_invalid_requests_are_refused(entries):
    with pytest.raises(ValueError, match="nothing to animate"):
        ca.frame_data([])
    with pytest.raises(ValueError, match="two frames"):
        ca.frame_data(entries, frames=1)


def test_the_figure_shows_every_star(entries):
    fig, animation = ca.cooling_animation(entries, interval_ms=1, frames=30)
    fig.canvas.draw()  # starts the animation, as the window does
    draw = fig._cooling_controls[3]
    ax_l, ax_t = fig.axes[0], fig.axes[1]
    assert [t.get_text() for t in ax_l.get_legend().get_texts()] == ["table", "1.40 Msun"]
    draw(0)
    assert "Frame 1/30" in fig._suptitle.get_text()
    draw(29)
    assert "Frame 30/30" in fig._suptitle.get_text() and "Final" in fig._suptitle.get_text()
    profiles = [line for line in ax_t.lines]
    assert [len(p.get_xdata()) for p in profiles] == [len(e[2]) for e in entries]
    curves = [line for line in ax_l.lines if line.get_linestyle() == "-"]
    assert [len(c.get_xdata()) for c in curves] == [
        len(e[1].times_years) for e in entries
    ]  # complete at the last frame
    plt.close(fig)


def test_one_star_needs_no_legend(entries):
    fig, animation = ca.cooling_animation(entries[:1], interval_ms=1, frames=10)
    fig.canvas.draw()
    assert fig.axes[0].get_legend() is None
    fig._cooling_controls[3](9)
    plt.close(fig)
