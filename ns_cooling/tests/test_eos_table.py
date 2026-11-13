"""CompOSE reader: parsing, validation, interpolation and the crust-core interface."""

import shutil

import numpy as np
import pytest

from eos_table import PHASE_CORE, PHASE_INNER_CRUST, PHASE_OUTER_CRUST, CompOSETable


def copy_package(eos_dir, tmp_path, skip=()):
    target = tmp_path / "pkg"
    shutil.copytree(eos_dir, target, ignore=shutil.ignore_patterns(*skip) if skip else None)
    return target


def test_grid_and_ranges(eos_table):
    assert len(eos_table.nb_fm3) == 489 and eos_table.name == "BSk22"
    assert eos_table.nb_fm3[0] == pytest.approx(4.67962e-10) and eos_table.nb_fm3[
        -1
    ] == pytest.approx(1.4922156)
    assert eos_table.pressure_min == pytest.approx(1.442e22, rel=1e-3)
    assert eos_table.pressure_max == pytest.approx(3.320e36, rel=1e-3)
    assert np.all(np.diff(eos_table.pressure) > 0) and np.all(np.diff(eos_table.rho_energy) > 0)


def test_thermodynamic_identity_from_the_raw_columns(eos_dir, eos_table):
    """e + p = n_b mu_b at T = 0, evaluated here independently of the reader."""
    rows = np.array(
        [
            [float(v) for v in line.split()[:10]]
            for line in (eos_dir / "eos.thermo").read_text().splitlines()[1:]
        ]
    )
    m_n = float((eos_dir / "eos.thermo").read_text().split()[0])
    nb = eos_table.nb_fm3
    p, mu_b, e = rows[:, 3] * nb, rows[:, 5], nb * m_n * (1 + rows[:, 9])
    assert np.max(np.abs((e + p) / (nb * m_n * (1 + mu_b)) - 1)) < 1e-6


def test_phases_are_ordered_and_counted(eos_table):
    phase = eos_table.phase
    assert (phase == PHASE_OUTER_CRUST).sum() == 217
    assert (phase == PHASE_INNER_CRUST).sum() == 83
    assert (phase == PHASE_CORE).sum() == 189
    assert np.all(np.diff((phase == PHASE_CORE).astype(int)) >= 0)  # crust rows first
    assert eos_table.pressure_core == eos_table.pressure[217 + 83]


def test_core_composition_is_charge_neutral(eos_table):
    core = eos_table.phase == PHASE_CORE
    assert np.max(np.abs(eos_table.yp[core] - eos_table.ye[core] - eos_table.ymu[core])) < 1e-6
    assert np.max(np.abs(eos_table.yn[core] + eos_table.yp[core] - 1)) < 1e-6
    assert eos_table.ymu[core][0] == 0 and eos_table.ymu[core][-1] > 0.2  # muons appear in the core


def test_crust_composition_follows_the_cooling_convention(eos_table):
    crust = eos_table.phase != PHASE_CORE
    rel = (
        np.abs(
            eos_table.ye[crust]
            - eos_table.z[crust] * (1 - eos_table.yn[crust]) / eos_table.a[crust]
        )
        / eos_table.ye[crust]
    )
    assert np.median(rel) < 1e-6  # Ye = Z (1 - Yn) / A holds in almost every row...
    assert rel.max() < 0.1  # ...and is off by at most a few percent near the crust-core interface


def test_rows_are_reproduced_exactly(eos_table):
    state = eos_table.at_pressure(eos_table.pressure)
    np.testing.assert_allclose(state.rho_energy, eos_table.rho_energy, rtol=1e-12)
    np.testing.assert_allclose(state.nb_fm3, eos_table.nb_fm3, rtol=1e-12)
    np.testing.assert_allclose(state.ye, eos_table.ye, rtol=1e-9, atol=1e-15)
    np.testing.assert_array_equal(state.core, eos_table.phase == PHASE_CORE)


def test_no_mixing_across_the_crust_core_interface(eos_table):
    n_crust = 217 + 83
    last_crust = eos_table.at_pressure(eos_table.pressure[n_crust - 1])
    below = eos_table.at_pressure(eos_table.pressure_core * (1 - 1e-9))
    first_core = eos_table.at_pressure(eos_table.pressure_core)
    assert not bool(below.core) and bool(first_core.core)
    assert bool(below.xh > 0) and float(first_core.xh) == 0.0
    # between the last crust row and the first core row the crust values are held constant
    assert float(below.aa) == pytest.approx(float(last_crust.aa), rel=1e-12)
    assert float(below.zz) == pytest.approx(float(last_crust.zz), rel=1e-12)
    assert float(first_core.yp) == pytest.approx(eos_table.yp[n_crust], rel=1e-12)
    assert float(first_core.aa) == 1.0 and float(below.yp) == 0.0


def test_conventions_in_the_crust_and_the_core(eos_table):
    outer = eos_table.at_pressure(1e25)  # outer crust: iron-group nuclei, no free neutrons
    assert (
        float(outer.xh) == pytest.approx(1.0) and float(outer.yn) == 0.0 and float(outer.yp) == 0.0
    )
    inner = eos_table.at_pressure(3e31)  # inner crust: free neutrons
    assert float(inner.yn) > 0.5 and float(inner.xh) == pytest.approx(1 - float(inner.yn))
    core = eos_table.at_pressure(1e35)
    assert float(core.xh) == 0.0 and float(core.zz) == pytest.approx(float(core.ye))
    assert float(core.yp) == pytest.approx(float(core.ye) + float(core.ymu), rel=1e-6)


def test_no_extrapolation(eos_table):
    for bad in (eos_table.pressure_min * 0.5, eos_table.pressure_max * 2, np.nan):
        with pytest.raises(ValueError):
            eos_table.at_pressure(bad)
    with pytest.raises(ValueError):
        eos_table.pressure_at_rho_baryon(1e5)  # below the table
    assert float(eos_table.at_pressure(eos_table.pressure_min).rho_energy) == pytest.approx(
        eos_table.rho_energy[0]
    )


def test_baryon_density_convention(eos_table):
    p = eos_table.pressure_at_rho_baryon(1e10)
    state = eos_table.at_pressure(p)
    assert float(state.rho_baryon) == pytest.approx(1e10, rel=1e-9)
    assert float(state.rho_baryon) == pytest.approx(
        float(state.nb_fm3) * 1e39 * 1.66e-24, rel=1e-12
    )
    assert float(eos_table.at_pressure(1e35).rho_baryon) < float(
        eos_table.at_pressure(1e35).rho_energy
    )


def test_micro_data_only_in_the_core(eos_table):
    assert eos_table.has_micro
    core = eos_table.phase == PHASE_CORE
    assert np.all(np.isfinite(eos_table.mstar_n[core])) and np.all(eos_table.mstar_n[core] > 0)
    assert np.all(np.isnan(eos_table.mstar_n[~core]))
    crust_point, core_point = eos_table.at_pressure(1e30), eos_table.at_pressure(1e35)
    assert np.isnan(float(crust_point.mstar_n)) and 0.3 < float(core_point.mstar_n) < 0.7
    assert np.nanmin(eos_table.gap_n) >= 0 and np.nanmax(eos_table.gap_p) < 0.1  # MeV


def test_mass_radius_curve_is_available(eos_table):
    mr = eos_table.mass_radius
    assert mr.shape[1] == 3 and mr[:, 1].max() == pytest.approx(2.265, abs=1e-3)


def test_optional_files_may_be_missing(eos_dir, tmp_path):
    pkg = copy_package(eos_dir, tmp_path, skip=("eos.micro", "eos.mr", "eos.pdf"))
    table = CompOSETable(pkg)
    assert not table.has_micro and table.mass_radius is None
    assert np.isnan(float(table.at_pressure(1e35).mstar_n))


def test_required_files_are_checked(eos_dir, tmp_path):
    pkg = copy_package(eos_dir, tmp_path, skip=("eos.compo",))
    with pytest.raises(FileNotFoundError, match="eos.compo"):
        CompOSETable(pkg)


def edit(path, function):
    path.write_text(function(path.read_text()))


def test_corrupted_tables_are_refused(eos_dir, tmp_path):
    # a grid that is not increasing
    pkg = copy_package(eos_dir, tmp_path)
    edit(pkg / "eos.nb", lambda t: t.replace("0.58496400E-09", "0.38496400E-09"))
    with pytest.raises(ValueError, match="strictly increasing"):
        CompOSETable(pkg)
    # a thermodynamic column that is not consistent (e + p != n_b mu_b)
    shutil.rmtree(pkg)
    pkg = copy_package(eos_dir, tmp_path)
    edit(pkg / "eos.thermo", lambda t: t.replace("0.13885259E+04", "0.23885259E+04"))
    with pytest.raises(ValueError, match="e \\+ p"):
        CompOSETable(pkg)
    # a header that announces a different number of points
    shutil.rmtree(pkg)
    pkg = copy_package(eos_dir, tmp_path)
    edit(pkg / "eos.nb", lambda t: t.replace("489", "488", 1))
    with pytest.raises(ValueError, match="announces"):
        CompOSETable(pkg)


def test_charge_neutrality_violation_is_refused(eos_dir, tmp_path):
    pkg = copy_package(eos_dir, tmp_path)
    lines = (pkg / "eos.compo").read_text().splitlines()
    t = lines[-1].split()
    t[6] = "0.30000000E+00"  # Ye of the last core row no longer matches Yp - Ymu
    lines[-1] = " ".join(t)
    (pkg / "eos.compo").write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError, match="charge neutrality"):
        CompOSETable(pkg)
