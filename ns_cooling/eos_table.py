"""Cold equation of state read from a CompOSE package (one folder per EOS, files ``eos.*``).

Required files
    eos.nb       baryon number density grid [fm^-3]
    eos.thermo   thermodynamics: m_n, m_p [MeV] and, per grid point, p/n_b [MeV],
                 mu_b/m_n - 1 and e/(n_b m_n) - 1 (e = energy density including the rest mass)
    eos.compo    composition and phase (0 core, 1 outer crust, 2 inner crust)
Optional files
    eos.micro    core only: Landau effective masses of n and p and the 1S0 pairing gaps [MeV]
    eos.mr       mass-radius curve (R [km], M [Msun], tidal deformability) of the authors,
                 used by the tests to validate the structure integration

The table must be a one-dimensional (T = 0, beta-equilibrium) EOS. Nothing is
extrapolated: asking for a pressure outside the table is an error.

Conventions of the returned ``EOSPoints`` (those of the cooling code)
    crust   xh = 1 - Yn (mass fraction in nuclei), yn = free-neutron fraction, yp = 0,
            aa = A and zz = Z of the nuclei, ye = Ye
    core    xh = 0, yn = Yn, yp = Yp (= Ye + Ymu), ye = Ye, aa = 1, zz = Ye
Linear interpolation in log P is used for every quantity. Crust and core are
interpolated separately, so values are never mixed across the crust-core interface:
below the first core row the crust values are held constant after the last crust row.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import warnings

from constants import FM3_PER_CM3, MEV_FM3_TO_DYN_CM2, MEV_TO_G, MICROPHYSICS_MASS_UNIT

PHASE_CORE, PHASE_OUTER_CRUST, PHASE_INNER_CRUST = 0, 1, 2

# CompOSE particle and quantity indices
_ELECTRON, _MUON, _NEUTRON, _PROTON = 0, 1, 10, 11
_LANDAU_MASS_N, _LANDAU_MASS_P, _GAP_N_1S0, _GAP_P_1S0 = 10040, 11040, 700060, 702060

_RANGE_TOLERANCE = 1e-12  # relative slack for pressures that sit exactly on a table edge


@dataclass(frozen=True)
class EOSPoints:
    """Equation of state and composition at a set of pressures.

    pressure [dyn cm^-2]; rho_energy = e/c^2 [g cm^-3]; nb_fm3 [fm^-3]; composition in
    the conventions of the cooling code (see the module docstring); ``core`` is True in
    homogeneous matter. mstar_n, mstar_p: Landau effective-mass ratios (NaN where the
    table has none: crust, or a table without eos.micro). gap_n, gap_p: 1S0 pairing
    gaps [MeV] (NaN where not provided).
    """

    pressure: np.ndarray
    rho_energy: np.ndarray
    nb_fm3: np.ndarray
    ye: np.ndarray
    ymu: np.ndarray
    yn: np.ndarray
    yp: np.ndarray
    xh: np.ndarray
    aa: np.ndarray
    zz: np.ndarray
    core: np.ndarray
    mstar_n: np.ndarray
    mstar_p: np.ndarray
    gap_n: np.ndarray
    gap_p: np.ndarray

    @property
    def rho_baryon(self):
        """Baryon mass density in the convention of the microphysics, 1.66e-24 g * n_b [g cm^-3]."""
        return self.nb_fm3 * FM3_PER_CM3 * MICROPHYSICS_MASS_UNIT


def _read_grid(path):
    lines = Path(path).read_text().split()
    if len(lines) < 3:
        raise ValueError(f"{path.name}: too short")
    n = int(lines[1])
    values = np.array([float(v) for v in lines[2:]])
    if values.size != n:
        raise ValueError(f"{path.name}: header announces {n} points, found {values.size}")
    return values


def _read_thermo(path, n):
    rows = Path(path).read_text().splitlines()
    first = rows[0].split()
    if len(first) < 2:
        raise ValueError(f"{path.name}: the first line must give m_n and m_p")
    m_n, m_p = float(first[0]), float(first[1])
    data = np.array([[float(v) for v in row.split()[:10]] for row in rows[1:] if row.strip()])
    if data.shape != (n, 10):
        raise ValueError(
            f"{path.name}: expected {n} rows of at least 10 columns, found {data.shape}"
        )
    if not (np.all(data[:, 0] == 0) and np.all(data[:, 2] == 1)):
        raise ValueError(
            f"{path.name}: only one-dimensional tables (T and Yq indices fixed) are supported"
        )
    if not np.array_equal(data[:, 1], np.arange(1, n + 1)):
        raise ValueError(f"{path.name}: the grid indices must run 1..{n} in order")
    return m_n, m_p, data[:, 3], data[:, 5], data[:, 9]  # p/n_b, mu_b/m_n - 1, e/(n_b m_n) - 1


def _read_compo(path, n):
    rows = Path(path).read_text().splitlines()
    rows = [r.split() for r in rows if r.strip()]
    if len(rows) != n:
        raise ValueError(f"{path.name}: expected {n} rows, found {len(rows)}")
    phase = np.empty(n, dtype=int)
    fractions = {k: np.zeros(n) for k in (_ELECTRON, _MUON, _NEUTRON, _PROTON)}
    a = np.zeros(n)
    z = np.zeros(n)
    for i, t in enumerate(rows):
        if int(t[1]) != i + 1:
            raise ValueError(f"{path.name}: row {i + 1} has grid index {t[1]}")
        phase[i] = int(t[3])
        n_pairs = int(t[4])
        for j in range(n_pairs):
            index, value = int(t[5 + 2 * j]), float(t[6 + 2 * j])
            if index in fractions:
                fractions[index][i] = value
        k = 5 + 2 * n_pairs
        n_nuclei = int(t[k])
        if n_nuclei > 1:
            raise ValueError(
                f"{path.name}: row {i + 1} has {n_nuclei} nuclear species; only one is supported"
            )
        if n_nuclei == 1:
            a[i], z[i] = float(t[k + 2]), float(t[k + 3])
    return phase, fractions, a, z


def _read_micro(path, n):
    columns = {
        k: np.full(n, np.nan) for k in (_LANDAU_MASS_N, _LANDAU_MASS_P, _GAP_N_1S0, _GAP_P_1S0)
    }
    rows = [r.split() for r in Path(path).read_text().splitlines() if r.strip()]
    if len(rows) != n:
        raise ValueError(f"{path.name}: expected {n} rows, found {len(rows)}")
    for i, t in enumerate(rows):
        for j in range(int(t[3])):
            index = int(t[4 + 2 * j])
            if index in columns:
                columns[index][i] = float(t[5 + 2 * j])
    return columns


class CompOSETable:
    """A cold EOS table; see the module docstring.

    Public arrays (one entry per table row): nb_fm3, pressure [dyn cm^-2],
    rho_energy [g cm^-3], phase, ye, ymu, yn, yp, a, z, and (NaN where absent)
    mstar_n, mstar_p, gap_n, gap_p. ``mass_radius`` is None or an array with the columns
    R [km], M [Msun], Lambda.
    """

    def __init__(self, folder):
        folder = Path(folder).resolve()
        for name in ("eos.nb", "eos.thermo", "eos.compo"):
            if not (folder / name).is_file():
                raise FileNotFoundError(f"{folder}: required CompOSE file {name} is missing")
        self.folder = folder
        self.name = folder.name
        self.nb_fm3 = _read_grid(folder / "eos.nb")
        n = self.nb_fm3.size
        m_n_mev, m_p_mev, p_over_nb, mu_b, e_over_nb = _read_thermo(folder / "eos.thermo", n)
        self.m_n_mev, self.m_p_mev = m_n_mev, m_p_mev
        pressure_mev = p_over_nb * self.nb_fm3
        energy_mev = self.nb_fm3 * m_n_mev * (1.0 + e_over_nb)
        self.pressure = pressure_mev * MEV_FM3_TO_DYN_CM2
        self.rho_energy = energy_mev * MEV_TO_G * FM3_PER_CM3
        self._validate_thermodynamics(pressure_mev, energy_mev, mu_b, m_n_mev)
        phase, fractions, self.a, self.z = _read_compo(folder / "eos.compo", n)
        self.phase = phase
        self.ye, self.ymu = fractions[_ELECTRON], fractions[_MUON]
        self.yn, self.yp = fractions[_NEUTRON], fractions[_PROTON]
        self._validate_composition()
        self._set_micro(folder, n)
        mr = folder / "eos.mr"
        self.mass_radius = np.loadtxt(mr) if mr.is_file() else None
        # Drop tiny pressure reversals consistently from all interpolated fields.
        # Keep source files unchanged; never sort P independently of composition.
        keep = self._pressure_keep_mask()
        self.omitted_rows = (np.flatnonzero(~keep) + 1).tolist()
        if self.omitted_rows:
            warnings.warn(
                f"{self.name}: omitting non-monotonic pressure rows {self.omitted_rows} "
                "from interpolation (relative reversal <= 1e-3); source files unchanged",
                UserWarning,
            )
            for name in (
                "nb_fm3",
                "pressure",
                "rho_energy",
                "phase",
                "ye",
                "ymu",
                "yn",
                "yp",
                "a",
                "z",
                "mstar_n",
                "mstar_p",
                "gap_n",
                "gap_p",
            ):
                setattr(self, name, getattr(self, name)[keep])
            self._validate_composition()
        self._logp = np.log(self.pressure)
        self._log_eps = np.log(self.rho_energy)
        self._log_nb = np.log(self.nb_fm3)

    # ------------------------------------------------------------------ validation
    def _validate_thermodynamics(self, pressure_mev, energy_mev, mu_b, m_n):
        if np.any(self.nb_fm3 <= 0) or np.any(np.diff(self.nb_fm3) <= 0):
            raise ValueError("eos.nb: the grid must be positive and strictly increasing")
        if np.any(pressure_mev <= 0):
            raise ValueError("eos.thermo: pressure must be positive and strictly increasing")
        if np.any(np.diff(energy_mev) <= 0):
            raise ValueError("eos.thermo: energy density must be strictly increasing")
        # T = 0 thermodynamic identity e + p = n_b mu_b guards against a misread column
        residual = np.max(
            np.abs((energy_mev + pressure_mev) / (self.nb_fm3 * m_n * (1.0 + mu_b)) - 1.0)
        )
        if residual > 1e-5:
            raise ValueError(
                f"eos.thermo: e + p = n_b mu_b violated (relative residual {residual:.1e})"
            )

    def _pressure_keep_mask(self):
        keep = np.ones(len(self.pressure), dtype=bool)
        last = 0
        for i in range(1, len(self.pressure)):
            if self.pressure[i] <= self.pressure[last]:
                if self.pressure[i] < self.pressure[last] * (1.0 - 1e-3):
                    raise ValueError("eos.thermo: pressure reversal exceeds 1e-3")
                keep[i] = False
            else:
                last = i
        return keep

    def _validate_composition(self):
        n = self.nb_fm3.size
        if not set(np.unique(self.phase)) <= {PHASE_CORE, PHASE_OUTER_CRUST, PHASE_INNER_CRUST}:
            raise ValueError("eos.compo: unknown phase index")
        crust = self.phase != PHASE_CORE
        self._n_crust = int(crust.sum())
        if not np.all(crust[: self._n_crust]) or np.any(crust[self._n_crust :]):
            raise ValueError("eos.compo: the crust rows must precede the core rows")
        if np.any(self.a[crust] <= 0) or np.any(self.z[crust] <= 0):
            raise ValueError("eos.compo: A and Z must be positive in the crust")
        for name in ("ye", "ymu", "yn", "yp"):
            v = getattr(self, name)
            if np.any(v < 0) or np.any(v > 1):
                raise ValueError(f"eos.compo: {name} outside [0, 1]")
        core = ~crust
        if core.any():
            if np.max(np.abs(self.yp[core] - self.ye[core] - self.ymu[core])) > 1e-5:
                raise ValueError("eos.compo: charge neutrality Yp = Ye + Ymu violated in the core")
            if np.max(np.abs(self.yn[core] + self.yp[core] - 1.0)) > 1e-5:
                raise ValueError("eos.compo: Yn + Yp = 1 violated in the core")
        if self._n_crust == n:
            self.pressure_core = np.inf
        else:
            self.pressure_core = float(self.pressure[self._n_crust])

    def _set_micro(self, folder, n):
        self.mstar_n, self.mstar_p = np.full(n, np.nan), np.full(n, np.nan)
        self.gap_n, self.gap_p = np.full(n, np.nan), np.full(n, np.nan)
        self.has_micro = False
        path = folder / "eos.micro"
        if not path.is_file():
            return
        c = _read_micro(path, n)
        core = slice(self._n_crust, n)
        masses = (c[_LANDAU_MASS_N], c[_LANDAU_MASS_P])
        if n > self._n_crust and all(
            np.all(np.isfinite(m[core])) and np.all(m[core] > 0) for m in masses
        ):
            self.mstar_n, self.mstar_p = masses
            self.gap_n, self.gap_p = c[_GAP_N_1S0], c[_GAP_P_1S0]
            self.has_micro = True

    # ------------------------------------------------------------------ ranges
    @property
    def pressure_min(self):
        return float(self.pressure[0])

    @property
    def pressure_max(self):
        return float(self.pressure[-1])

    def _log_pressure(self, pressure):
        p = np.asarray(pressure, dtype=float)
        lo, hi = self.pressure_min, self.pressure_max
        if (
            not np.all(np.isfinite(p))
            or np.any(p < lo * (1 - _RANGE_TOLERANCE))
            or np.any(p > hi * (1 + _RANGE_TOLERANCE))
        ):
            raise ValueError(
                f"pressure must lie in [{lo:.4g}, {hi:.4g}] dyn/cm^2 (no extrapolation)"
            )
        return np.log(np.clip(p, lo, hi))

    # ------------------------------------------------------------------ queries
    def energy_density(self, pressure):
        """e/c^2 [g cm^-3] at a pressure (fast path used by the TOV integration)."""
        return np.exp(np.interp(self._log_pressure(pressure), self._logp, self._log_eps))

    def pressure_at_rho_baryon(self, rho_baryon):
        """Pressure at a baryon mass density 1.66e-24 g * n_b [g cm^-3] (the convention of ``rho``)."""
        nb = np.asarray(rho_baryon, dtype=float) / (FM3_PER_CM3 * MICROPHYSICS_MASS_UNIT)
        lo, hi = self.nb_fm3[0], self.nb_fm3[-1]
        if (
            not np.all(np.isfinite(nb))
            or np.any(nb < lo * (1 - _RANGE_TOLERANCE))
            or np.any(nb > hi)
        ):
            raise ValueError(
                f"baryon density must lie in [{lo * FM3_PER_CM3 * MICROPHYSICS_MASS_UNIT:.4g}, "
                f"{hi * FM3_PER_CM3 * MICROPHYSICS_MASS_UNIT:.4g}] g/cm^3"
            )
        return np.exp(np.interp(np.log(np.maximum(nb, lo)), self._log_nb, self._logp))

    def at_pressure(self, pressure):
        """Equation of state and composition at the given pressures (see ``EOSPoints``)."""
        x = self._log_pressure(pressure)
        p = np.exp(x)
        n_c = self._n_crust
        is_core = p >= self.pressure_core

        def crust(values):
            return np.interp(x, self._logp[:n_c], values[:n_c]) if n_c else np.zeros_like(x)

        def core(values):
            return (
                np.interp(x, self._logp[n_c:], values[n_c:])
                if n_c < len(values)
                else np.zeros_like(x)
            )

        ye = np.where(is_core, core(self.ye), crust(self.ye))
        yn = np.where(is_core, core(self.yn), crust(self.yn))
        ymu = np.where(is_core, core(self.ymu), 0.0)
        yp = np.where(is_core, core(self.yp), 0.0)
        xh = np.where(is_core, 0.0, 1.0 - yn)
        aa = np.where(is_core, 1.0, crust(self.a))
        zz = np.where(is_core, ye, crust(self.z))
        nan = np.full_like(x, np.nan)
        masses = [
            np.where(is_core, core(m), nan) if self.has_micro else nan
            for m in (self.mstar_n, self.mstar_p)
        ]
        gaps = [
            np.where(is_core, core(g), nan) if self.has_micro else nan
            for g in (self.gap_n, self.gap_p)
        ]
        return EOSPoints(
            pressure=p,
            rho_energy=np.exp(np.interp(x, self._logp, self._log_eps)),
            nb_fm3=np.exp(np.interp(x, self._logp, self._log_nb)),
            ye=ye,
            ymu=ymu,
            yn=yn,
            yp=yp,
            xh=xh,
            aa=aa,
            zz=zz,
            core=is_core,
            mstar_n=masses[0],
            mstar_p=masses[1],
            gap_n=gaps[0],
            gap_p=gaps[1],
        )
