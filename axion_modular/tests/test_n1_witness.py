"""Regression of N1 against the OLD committed results (fig_axcow.json, shifts in %, 4 significant digits) at nu_a = 6 kHz."""
import json

import numpy as np
import pytest

from common import artifacts as A
from tests.conftest import REPO, ROOT

OLD = REPO / 'results' / 'fig_axcow.json'
NEW = ROOT / 'N1_modes_cowling' / 'outputs' / 'modes_cowling.csv'


@pytest.mark.skipif(not (OLD.exists() and NEW.exists()), reason='needs legacy fig_axcow.json and an N1 run')
def test_witness_matches_legacy_percent_shifts():
    old = json.loads(OLD.read_text())
    i = old['nus'].index(6.0)
    cols, meta = A.load_table(NEW)
    sel = cols['nu_a_kHz'] == 6.0
    assert sel.sum() == 3
    nu0 = cols['nu_noaxion_kHz'][sel]
    import math

    def close(new, old):
        """Legacy values keep 4 significant digits: accept up to 0.6 unit of the last digit (rounding boundary)."""
        if old == 0:
            return abs(new) < 1e-9
        return abs(new - old) <= 0.6 * 10 ** (math.floor(math.log10(abs(old))) - 3)

    for k, mode in enumerate(old['modes']):
        bg = 100 * (cols['nu_bg_kHz'][sel][k] / nu0[k] - 1)
        cp = 100 * (cols['nu_cp_re_kHz'][sel][k] / nu0[k] - 1)
        im = abs(cols['nu_cp_im_kHz'][sel][k]) * 1e3
        assert close(bg, mode['bg'][i]), (mode['name'], 'bg', bg, mode['bg'][i])
        assert close(cp, mode['cp'][i]), (mode['name'], 'cp', cp, mode['cp'][i])
        assert close(im, mode['im'][i]), (mode['name'], 'im', im, mode['im'][i])
