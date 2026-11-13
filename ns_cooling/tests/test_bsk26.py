"""BSk26 input and monotonic interpolation of the supplied table."""

import numpy as np
import pytest
from eos_table import CompOSETable


def test_bsk26_small_pressure_reversal(root):
    with pytest.warns(UserWarning, match="276"):
        eos = CompOSETable(root / "eos" / "BSk26")
    assert eos.omitted_rows == [276]
    assert np.all(np.diff(eos.pressure) > 0)
    points = eos.at_pressure(eos.pressure)
    np.testing.assert_allclose(points.nb_fm3, eos.nb_fm3, rtol=1e-12)
    assert len(eos.mstar_n) == len(eos.pressure)


def test_large_reversal_is_rejected():
    eos = object.__new__(CompOSETable)
    eos.pressure = np.array([1.0, 2.0, 1.9, 3.0])
    with pytest.raises(ValueError, match="reversal"):
        eos._pressure_keep_mask()
