"""Input file parsing and the command-line driver."""

import subprocess
import sys

import pytest

from input_model import read_parameters

LINES = [
    "structure file", "none", "initial T", "1.d2", "unused", "800",
    "sf", "a", "1 7", "b", "1 12", "c", "1 8", "envelope_model = accreted", "", "Options: ...",
]  # fmt: skip


def write(tmp_path, lines):
    path = tmp_path / "input_test"
    path.write_text("\n".join(lines) + "\n")
    return path


def test_included_input_file(root):
    p = read_parameters(root / "input_bsk26")
    assert p["structure"] == "none" and p["Tinit8"] == 100.0
    assert p["envelope_model"] == "gudmundsson" and p["eos"] == "BSk26"
    assert p["crust_n"] == {"enabled": True, "model": 7} and p["core_p"] == {
        "enabled": True,
        "model": 8,
    }


def test_the_unused_line_is_ignored(tmp_path):
    lines = list(LINES)
    lines[5] = "not a number"
    assert read_parameters(write(tmp_path, lines))["Tinit8"] == 100.0


def test_default_envelope_when_the_option_is_absent(tmp_path):
    assert read_parameters(write(tmp_path, LINES[:13]))["envelope_model"] == "gudmundsson"


@pytest.mark.parametrize(
    "edit, message",
    [
        (lambda L: L.__setitem__(8, "2 7"), "Invalid configuration"),
        (lambda L: L.__setitem__(10, "1 99"), "Invalid configuration"),
        (lambda L: L.__setitem__(3, "-5"), "Initial temperature"),
        (lambda L: L.__setitem__(13, "envelope_model = iron"), "envelope_model must be one of"),
        (lambda L: L.__setitem__(14, "envelope_model = fe"), "Duplicate"),
        (lambda L: L.__setitem__(14, "something = 1"), "Unknown input option"),
    ],
)
def test_invalid_input_files_are_rejected(tmp_path, edit, message):
    lines = list(LINES)
    edit(lines)
    with pytest.raises(ValueError, match=message):
        read_parameters(write(tmp_path, lines))


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


def test_driver_rejects_invalid_options(root, tmp_path):
    for args in (("--theta", 0.2), ("--final-years", 0.001), ("--interval-ms", -1)):
        result = run_driver(root, *args, "--out", tmp_path / "x")
        assert result.returncode != 0 and "error" in result.stderr


def test_scheme_flags_select_backward_euler(root, tmp_path):
    result = run_driver(
        root, "--final-years", 30, "--theta", 1, "--no-predictor-corrector", "--out", tmp_path / "o"
    )
    assert result.returncode == 0, result.stderr


def test_a_single_star_keeps_the_flat_layout(root, tmp_path):
    result = run_driver(root, "--final-years", 30, "--out", tmp_path / "o")
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "o" / "cooling_curve.csv").is_file() and not (
        tmp_path / "o" / "cooling_curves.csv"
    ).exists()


def test_a_face_table_of_an_eos_star_cannot_be_used_as_a_background(root, tmp_path, bsk22_model):
    path = tmp_path / "bsk22.csv"
    bsk22_model.background.to_csv(path)
    result = run_driver(root, "--background", path, "--final-years", 30, "--out", tmp_path / "o")
    assert result.returncode != 0 and "inspection only" in result.stderr


def test_folder_names_from_labels():
    from cooling_driver import folder_name

    assert (
        folder_name("1.40 Msun") == "1.40_Msun"
        and folder_name("a/b  c") == "a_b_c"
        and folder_name("???") == "star"
    )
