"""Shared fixtures. The tests import the modules of the parent folder."""

import sys
from dataclasses import replace
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from config import CoolingConfig  # noqa: E402
from cooling import run_cooling  # noqa: E402
from input_model import load_background  # noqa: E402


@pytest.fixture(scope="session")
def root():
    return ROOT


@pytest.fixture(scope="session")
def background():
    """The background of the default star: BSk26, 1.4 Msun (read-only, shared by all tests)."""
    return load_background(ROOT / "input_bsk26")


@pytest.fixture(scope="session")
def short_run(background):
    """Default scheme from 1e10 K to 1000 yr (about 2 s)."""
    return run_cooling(background, replace(CoolingConfig(), final_years=1000.0))


# ---------------------------------------------------------------- equation-of-state fixtures
CENTRAL_PRESSURE_14 = 7.94037e34  # dyn/cm^2: central pressure of the 1.4 Msun BSk22 star


@pytest.fixture(scope="session")
def eos_dir(root):
    return root / "eos" / "BSk22"


@pytest.fixture(scope="session")
def eos_table(eos_dir):
    from eos_table import CompOSETable

    return CompOSETable(eos_dir)


@pytest.fixture(scope="session")
def bsk22_star(eos_table):
    """The 1.4 Msun BSk22 structure (integrated once, shared by all tests)."""
    from tov import integrate_star

    return integrate_star(eos_table, CENTRAL_PRESSURE_14)


@pytest.fixture(scope="session")
def input_config(root):
    from input_model import read_parameters

    return read_parameters(root / "input_bsk22")


@pytest.fixture(scope="session")
def bsk22_model(root, input_config):
    """Stellar model (mesh A) built from the BSk22 table; the central pressure is given
    directly so that the fixture does not repeat the mass search tested elsewhere."""
    import warnings

    from stellar_structure import StructureOptions, build_stellar_model

    options = StructureOptions("BSk22", central_pressure=CENTRAL_PRESSURE_14)
    with warnings.catch_warnings():  # mesh A leaves the Direct Urca edge unresolved: tested on purpose elsewhere
        warnings.simplefilter("ignore", UserWarning)
        return build_stellar_model(options, root / "eos", input_config)


@pytest.fixture(scope="session")
def branch(eos_table):
    """Minimum and maximum mass of the stable branch of BSk22 (about 15 s, computed once)."""
    from tov import stable_branch

    return stable_branch(eos_table)
