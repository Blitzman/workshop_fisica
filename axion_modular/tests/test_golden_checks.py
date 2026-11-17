"""The refactored modules must report the same checks, with the same outcome, as the original steps (committed checks.json)."""
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
OLD = ROOT.parent / 'results' / 'checks.json'


def _old(*modules):
    d = json.loads(OLD.read_text(encoding='utf-8'))
    return {c['name']: c['ok'] for s in d['steps'] if s['module'] in modules for c in s['checks']}


def _new(*mods):
    out = {}
    for m in mods:
        p = ROOT / m / 'outputs' / 'manifest.json'
        if not p.exists():
            pytest.skip(f'{m} not run')
        out.update({c['name']: c['ok'] for c in json.loads(p.read_text(encoding='utf-8'))['checks']})
    return out


@pytest.mark.skipif(not OLD.exists(), reason='legacy results/checks.json not present')
def test_paper1_checks_match_old_steps_1_to_3():
    old = _old('step1_background', 'step2_perturbations', 'step3_toy_model')
    new = _new('D0_derive_action', 'D1_derive_radial')
    assert set(old) == set(new), (set(old) ^ set(new))
    assert {k: v for k, v in old.items()} == {k: new[k] for k in old}


@pytest.mark.skipif(not OLD.exists(), reason='legacy results/checks.json not present')
def test_paper2_checks_match_old_steps_5_and_7():
    old = _old('step5_nonradial', 'step7_fullgr')
    new = _new('D2_derive_nonradial')
    assert set(old) == set(new), (set(old) ^ set(new))
    assert old == {k: new[k] for k in old}


@pytest.mark.skipif(not OLD.exists(), reason='legacy results/checks.json not present')
def test_toy_checks_match_old_steps_4_and_14():
    old = _old('step4_numerics', 'step14_radial_branches')
    new = {k: v for k, v in _new('T_toy_model').items() if not k.startswith('Puente')}
    assert set(old) == set(new), (set(old) ^ set(new))      # names embed the computed numbers
    assert old == {k: new[k] for k in old}


@pytest.mark.skipif(not OLD.exists(), reason='legacy results/checks.json not present')
def test_atom_level_checks_match_old_step_6():
    old = _old('step6_nonradial_numerics')
    new = _new('N1_modes_cowling')
    missing = [k for k in old if k not in new]
    assert not missing, missing
    assert old == {k: new[k] for k in old}


@pytest.mark.skipif(not OLD.exists(), reason='legacy results/checks.json not present')
def test_n2_checks_match_old_step_8():
    """Full run only (`make N2`). Names embed computed numbers (noisy at 1e-5), so they are compared with the digits masked."""
    import re
    p = ROOT / 'N2_modes_fullgr_noaxion' / 'outputs' / 'manifest.json'
    if not p.exists():
        pytest.skip('N2 not run')
    man = json.loads(p.read_text(encoding='utf-8'))
    if man.get('extra', {}).get('quick'):
        pytest.skip('N2 quick run')
    norm = lambda n: re.sub(r'[-+]?\d[\d.eE+-]*', '#', n)
    d = json.loads(OLD.read_text(encoding='utf-8'))
    old = [(norm(c['name']), c['ok']) for s in d['steps'] if s['module'] == 'step8_numerics' for c in s['checks']
           if not c['name'].startswith('Fondo con axi')]
    new = [(norm(c['name']), c['ok']) for c in man['checks']]
    assert old == new, (old, new)
