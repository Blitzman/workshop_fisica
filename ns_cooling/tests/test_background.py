"""The Background: validation of its arrays, tables and scalars."""

from dataclasses import replace

import numpy as np
import pytest

from background import POINT_FIELDS, Background, Points


def test_arrays_are_read_only(background):
    with pytest.raises(ValueError):
        background.cells.rho[0] = 1.0
    with pytest.raises(ValueError):
        background.radius_km[1] = 1.0


@pytest.mark.parametrize(
    "change, message",
    [
        (lambda t: t.update(rho=-t["rho"]), "rho"),
        (lambda t: t.update(ye=t["ye"] + 2.0), "ye|particle_fraction"),
        (lambda t: t.update(ephi=0 * t["ephi"]), "ephi"),
        (lambda t: t.update(radius_km=t["radius_km"][::-1].copy()), "radius"),
        (lambda t: t.update(elam=-t["elam"]), "elam"),
        (lambda t: t.pop("tc_n"), "tc_n"),
    ],
)
def test_invalid_tables_are_rejected(background, change, message):
    table = {k: np.array(v) for k, v in background.face_table().items()}
    change(table)
    with pytest.raises(ValueError, match=message):
        Background.from_face_table(table, g14=background.g14)


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
