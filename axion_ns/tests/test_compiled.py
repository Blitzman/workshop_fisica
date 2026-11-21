"""The compiled-equation cache (common/compiled.py): round trip, stale detection, and that the committed caches of N2, N3, N4 match their
equation artifacts and reduction code (they load without recompiling)."""
import numpy as np
import pytest
import sympy as sp

from common import compiled as C


def test_cache_round_trip_stale_detection_and_recompile(tmp_path):
    x, y = sp.symbols('x y')
    path = tmp_path / 'f.json'
    built = []

    def build():
        built.append(1)
        return {'f': ((x, y), [sp.exp(x) * sp.sin(y) + x**2, sp.sqrt(y)], 'numpy'), 'g': ((x,), sp.cos(x)**2, 'math')}
    f1 = C.compiled(path, {'equations': 'a', 'code': 'b'}, build)
    f2 = C.compiled(path, {'equations': 'a', 'code': 'b'}, build)
    assert len(built) == 1                                                   # second call loaded the cache
    assert np.array_equal(f1['f'](0.3, 1.7), f2['f'](0.3, 1.7)) and f1['g'](0.4) == f2['g'](0.4)
    with pytest.raises(C.StaleCache, match='equations'):
        C.compiled(path, {'equations': 'changed', 'code': 'b'}, build)
    C.RECOMPILE = True
    try:
        C.compiled(path, {'equations': 'changed', 'code': 'b'}, build)
    finally:
        C.RECOMPILE = False
    assert len(built) == 2
    C.compiled(path, {'equations': 'changed', 'code': 'b'}, build)       # the rebuilt cache now matches
    assert len(built) == 2


def test_committed_caches_match_the_equation_artifacts(monkeypatch):
    """Loads the committed caches; StaleCache here means an equation artifact or reduction code changed: rerun the module with --recompile."""
    from tests.conftest import ROOT
    from N2_modes_fullgr_noaxion import systems as S2
    from N3_modes_fullgr_axion import systems as S3
    from N4_modes_radial import systems as S4

    orig = C.compiled

    def strict(path, key, build):
        def missing():
            raise AssertionError(f'{path} is missing: run the module once (or with --recompile)')
        return orig(path, key, missing)
    monkeypatch.setattr(C, 'compiled', strict)
    ext, p2, _ = S2.load_inputs(ROOT / 'D2_derive_nonradial' / 'outputs', ROOT / 'R_paper_refs' / 'outputs')
    L, B, _ = S4.load_inputs(ROOT / 'D0_derive_action' / 'outputs', ROOT / 'D1_derive_radial' / 'outputs')
    S2.cowling_system(p2, 'energy')
    S2.full_systems(ext, 'energy')
    S2.full_systems(ext, 'rest')                 # used by tests/compare_n3_n2_matrix.py
    S3.compile_flat(ext)
    S4.compile_flat(L, B)

