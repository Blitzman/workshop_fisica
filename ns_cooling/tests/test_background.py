"""The Background: validation of its arrays and scalars."""

from dataclasses import replace

import numpy as np
import pytest

from background import POINT_FIELDS, Points


def test_arrays_are_read_only(background):
    with pytest.raises(ValueError):
        background.cells.rho[0] = 1.0
    with pytest.raises(ValueError):
        background.radius_km[1] = 1.0


@pytest.mark.parametrize(
    "change, message",
    [
        (lambda f: f.update(rho=-f["rho"]), "rho"),
        (lambda f: f.update(ye=f["ye"] + 2.0), "ye must lie"),
        (lambda f: f.update(ephi=0 * f["ephi"]), "ephi"),
    ],
)
def test_invalid_points_are_rejected(background, change, message):
    fields = {k: np.array(getattr(background.cells, k)) for k in POINT_FIELDS}
    change(fields)
    with pytest.raises(ValueError, match=message):
        Points(**fields)


@pytest.mark.parametrize(
    "change, message",
    [
        (lambda b: {"radius_km": b.radius_km[::-1].copy()}, "radius"),
        (lambda b: {"elam": -b.elam}, "elam"),
    ],
)
def test_invalid_geometry_is_rejected(background, change, message):
    with pytest.raises(ValueError, match=message):
        replace(background, **change(background))


def test_invalid_scalars_are_rejected(background):
    with pytest.raises(ValueError, match="g14"):
        replace(background, g14=-1.0)
    with pytest.raises(ValueError, match="envelope_model"):
        replace(background, envelope_model="iron_maiden")


def test_points_require_equal_lengths(background):
    fields = {k: getattr(background.cells, k) for k in POINT_FIELDS}
    fields["tc_p"] = fields["tc_p"][:-1]
    with pytest.raises(ValueError):
        Points(**fields)
