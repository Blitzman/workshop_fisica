"""Tests of the output area (decision D-028). Run them with `make test-output` when this part of the code changes, or on request.

Code covered: cooling_driver.py, cooling_animation.py, compare_emissivities.py, neutrino_dominance.py.
Every test states a hypothesis, a reference and a criterion, in the style of the validation steps of the GitHub code.
"""

import cooling_animation as ca
import matplotlib.pyplot as plt
import pytest
import subprocess
import sys
from config import CoolingConfig
from cooling import run_cooling
from dataclasses import replace


# ------------------------------------------------------------------ command line
def run_driver(root, *args):
    """Run the driver; the default star (input_bsk26) is used unless --input is given."""
    given = "--input" in map(str, args)
    default = [] if given else ["--input", str(root / "input_bsk26")]
    return subprocess.run(
        [sys.executable, str(root / "cooling_driver.py"), "--no-animation", *default, *map(str, args)],
        capture_output=True, text=True, cwd=root,
    )  # fmt: skip


def test_driver_writes_all_outputs(root, tmp_path):
    """Hypothesis: one short run of the driver writes every output file (smoke test).
    Reference: the list of output files in the README.
    Criterion: all files are non-empty, the energy summary is printed, and a single star uses the
        flat layout.
    """
    result = run_driver(root, "--final-years", 30, "--out", tmp_path / "o", "--check-energy")
    assert result.returncode == 0, result.stderr
    for name in (
        "cooling_curve.csv",
        "temperature_profiles.csv",
        "profile_times.csv",
        "time_steps.csv",
        "energy_budget.csv",
    ):
        assert (tmp_path / "o" / name).stat().st_size > 0
    assert "Energy balance" in result.stdout and "Accepted steps" in result.stdout
    assert not (tmp_path / "o" / "cooling_curves.csv").exists()  # one star: the flat layout


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


def test_the_figure_shows_every_star(entries):
    """Hypothesis: the figure shows every star, and its complete curve at the last frame (smoke test).
    Reference: the structure of the figure: legend, frame titles and curves.
    Criterion: legend labels, frame titles and the length of every line.
    """
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
