import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
REPO = ROOT.parent
for p in (ROOT, REPO):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))


# ---- stars built ONCE per test session (each StarLogP costs 0.5-2 s) and shared by every test that needs them
@pytest.fixture(scope='session')
def star_noaxion():
    from S1_star_background.star_logp import StarLogP
    return StarLogP(None)


@pytest.fixture(scope='session')
def star_6():
    from S1_star_background.star_logp import StarLogP
    return StarLogP(6.0)


@pytest.fixture(scope='session')
def star_12():
    from S1_star_background.star_logp import StarLogP
    return StarLogP(12.0)


def _write(tmp_path_factory, st, name):
    """S1's table for `st` in a fresh directory, read back with the shared reader (common.star_table)."""
    from S1_star_background.export import write_eos, write_star
    from common.star_table import StarTable
    d = tmp_path_factory.mktemp(name)
    write_eos(d, 'test')
    write_star(d, st, name, 'test')
    return StarTable(d / 'eos_bsk26.csv', d / f'{name}.csv'), d


@pytest.fixture(scope='session')
def table_noaxion(tmp_path_factory, star_noaxion):
    """(StarTable, directory) of the no-axion star."""
    return _write(tmp_path_factory, star_noaxion, 'star_noaxion')


@pytest.fixture(scope='session')
def table_6(tmp_path_factory, star_6):
    """(StarTable, directory) of the 6 kHz star."""
    return _write(tmp_path_factory, star_6, 'star_nua6')
