"""Time-independent background of the star.

``Background`` is the single object that holds everything the thermal evolution
needs and that does not change with time or temperature. The evolution, the
microphysics registries and the diagnostics only see this object: they never
read files and they never know where the numbers came from.

Layout
    faces   N points: conductivity is evaluated at the interior faces
    cells   N-1 points: heat capacity and neutrino emissivity are evaluated at
            the cell centers
    radius_km, elam   N face values: geometry and metric factor exp(Lambda)
    g14, envelope_model   surface gravity [1e14 cm s^-2] and Tb-Ts relation
    surface_radius_km, surface_ephi   optional radius and exp(Phi) of the outer edge of
            the star when it lies above the last face (the envelope is neglected in the
            thermal domain). They only enter the luminosity at infinity and the flux
            leaving the base, through ``surface_factor`` (1 if not given)

Fields of ``Points`` (one value per point; units in brackets)
    rho [g cm^-3], ye, yn, yp (electron, neutron, proton fractions per baryon), ymu (muon
    fraction: optional, 0 without muons),
    xh (mass fraction in heavy nuclei; 0 in the core), aa, zz (nucleus mass and
    charge numbers, used only where xh > 0), ephi (exp(Phi), redshift),
    tc_crust, tc_n, tc_p (critical temperatures [K] of crust neutrons, core
    neutrons and core protons; 0 = not superfluid),
    kf_n, kf_p, kf_e, kf_mu (Fermi momenta [fm^-1], each from its own number density
    rho*y/1.66e-24; in the core charge neutrality means kf_p^3 = kf_e^3 + kf_mu^3),
    mstar_n, mstar_p (Landau effective-mass ratios m*/m of the neutrons and protons),
    transport_mstar_n, transport_mstar_p (the ones used by the conductivity).

How a ``Background`` is produced (providers)
    The files of this package are one provider (``input_model.load_background``).
    A future provider only has to produce the *face table* below, either as a
    dictionary of arrays or as a CSV file, and call ``Background.from_face_table``
    or ``Background.from_csv``:

        radius_km, rho, ye, yn, yp, xh, aa, zz, ephi, elam, tc_crust, tc_n, tc_p [, ymu]

    ``ymu`` is optional (zeros if absent). Cell values are face averages. The Fermi momenta and effective masses are
    derived from the table by ``derived_fermi_quantities``; a different routine
    (for example one that returns effective masses from the equation of state)
    can be passed as ``derived=``. If every point of the background should come
    from an external routine, build the two ``Points`` objects directly.
"""

import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from effective_mass import effective_mass_ratios
from superfluid import fermi_quantities
from thermal_fluxes import MODELS as ENVELOPE_MODELS

POINT_FIELDS = (
    "rho",
    "ye",
    "yn",
    "yp",
    "xh",
    "aa",
    "zz",
    "ephi",
    "tc_crust",
    "tc_n",
    "tc_p",
    "kf_n",
    "kf_e",
    "mstar_n",
    "mstar_p",
    "transport_mstar_n",
    "transport_mstar_p",
    "ymu",
    "kf_p",
    "kf_mu",
)
# Primary columns of the face table (the other fields of Points are derived).
AVERAGED_FIELDS = (
    "rho", "ye", "yn", "yp", "xh", "aa", "zz", "ephi", "tc_crust", "tc_n", "tc_p", "ymu",
)  # fmt: skip
OPTIONAL_COLUMNS = ("ymu",)  # absent from older face tables: zeros
FACE_TABLE_COLUMNS = ("radius_km",) + AVERAGED_FIELDS[:8] + ("elam",) + AVERAGED_FIELDS[8:]


def _readonly(values, name, size=None):
    a = np.array(values, dtype=float, copy=True)
    if a.ndim != 1 or a.size == 0:
        raise ValueError(f"{name}: a non-empty one-dimensional array is required")
    if size is not None and a.size != size:
        raise ValueError(f"{name}: expected {size} values, got {a.size}")
    if not np.isfinite(a).all():
        raise ValueError(f"{name}: non-finite values")
    a.setflags(write=False)
    return a


@dataclass(frozen=True)
class Points:
    """Background fields at a set of locations (all arrays have the same length)."""

    rho: np.ndarray
    ye: np.ndarray
    yn: np.ndarray
    yp: np.ndarray
    xh: np.ndarray
    aa: np.ndarray
    zz: np.ndarray
    ephi: np.ndarray
    tc_crust: np.ndarray
    tc_n: np.ndarray
    tc_p: np.ndarray
    kf_n: np.ndarray
    kf_e: np.ndarray
    mstar_n: np.ndarray
    mstar_p: np.ndarray
    transport_mstar_n: np.ndarray
    transport_mstar_p: np.ndarray
    ymu: np.ndarray
    kf_p: np.ndarray
    kf_mu: np.ndarray

    def __post_init__(self):
        size = None
        for name in POINT_FIELDS:
            a = _readonly(getattr(self, name), name, size)
            size = a.size
            object.__setattr__(self, name, a)
        if np.any(self.rho <= 0):
            raise ValueError("rho must be positive")
        for name in ("ye", "yn", "yp", "xh", "ymu"):
            a = getattr(self, name)
            if np.any((a < 0) | (a > 1)):
                raise ValueError(f"{name} must lie in [0, 1]")
        if np.any(self.ymu > self.yp + 1e-9):
            raise ValueError("ymu cannot exceed yp: the muons are neutralized by protons")
        if np.any(self.ephi <= 0):
            raise ValueError("ephi must be positive")
        for name in ("tc_crust", "tc_n", "tc_p", "kf_n", "kf_p", "kf_e", "kf_mu"):
            if np.any(getattr(self, name) < 0):
                raise ValueError(f"{name} must be non-negative")
        for name in ("mstar_n", "mstar_p", "transport_mstar_n", "transport_mstar_p"):
            if np.any(getattr(self, name) <= 0):
                raise ValueError(f"{name} must be positive")
        nuclei = self.xh > 0
        if np.any(self.aa[nuclei] <= 0) or np.any(self.zz[nuclei] <= 0):
            raise ValueError("aa and zz must be positive where xh > 0")

    def __len__(self):
        return self.rho.size

    def interior(self):
        """The points without the first and the last one (interior faces)."""
        return Points(**{name: getattr(self, name)[1:-1] for name in POINT_FIELDS})


def derived_fermi_quantities(rho, ye, yn, yp):
    """Fermi momenta and effective-mass ratios derived from the composition.

    Returns a dictionary with ``kf_n``, ``kf_p``, ``kf_e`` [fm^-1], each from its own
    number density (the muon momentum ``kf_mu`` is added by the caller, which knows ``ymu``), ``mstar_n``,
    ``mstar_p`` (from the neutron and proton number densities),
    ``transport_mstar_n``, ``transport_mstar_p`` (the same composition, using
    the actual Yn and Yp also when muons are present). A replacement routine must accept the same
    four arrays and return ``kf_n``, ``kf_e`` and the four mass keys; ``kf_p`` is computed
    from ``yp`` if the routine does not return it.
    """

    def fermi(fraction):
        # Point-by-point evaluation, done once when the background is built.
        values = [fermi_quantities(float(r), float(y)) for r, y in zip(rho, fraction)]
        return (
            np.array([float(v[0]) for v in values]),
            np.array([float(v[1]) for v in values]),
        )

    n_n, kf_n = fermi(yn)
    n_p, kf_p = fermi(yp)
    _, kf_e = fermi(ye)
    nb = rho / 1.66e-24
    # Evaluated at each location, not by averaging effective masses of other locations.
    mstar_n, mstar_p = effective_mass_ratios(n_n, n_p)
    transport_n, transport_p = effective_mass_ratios(nb * yn, nb * yp)
    return dict(
        kf_n=kf_n,
        kf_p=kf_p,
        kf_e=kf_e,
        mstar_n=mstar_n,
        mstar_p=mstar_p,
        transport_mstar_n=transport_n,
        transport_mstar_p=transport_p,
    )


def fermi_momentum(rho, fraction):
    """Fermi momentum [fm^-1] of a species with ``fraction`` particles per baryon at mass density rho.

    The scalar path of ``superfluid.fermi_quantities`` is used, as for the nucleons and electrons.
    """
    return np.array([float(fermi_quantities(float(r), float(y))[1]) for r, y in zip(rho, fraction)])


def composition_residuals(points):
    """Charge-neutrality and baryon-number residuals at the points (dictionary of arrays).

    core (xh = 0):   charge = yp - ye - ymu       baryon = yn + yp - 1
    crust (xh > 0):  charge = ye - xh zz/aa       baryon = yn + xh - 1
    (the nuclei carry zz/aa charges per nucleon). ``charge_relative`` is |charge| / ye.
    """
    crust = points.xh > 0
    safe_a = np.where(crust, points.aa, 1.0)
    charge = np.where(
        crust, points.ye - points.xh * points.zz / safe_a, points.yp - points.ye - points.ymu
    )
    baryon = np.where(crust, points.yn + points.xh - 1.0, points.yn + points.yp - 1.0)
    relative = np.abs(charge) / np.maximum(points.ye, 1e-300)
    return dict(charge=charge, baryon=baryon, charge_relative=relative)


def check_composition(
    points, core_tolerance=1e-5, baryon_tolerance=1e-5, crust_charge_warning=0.25
):
    """Verify charge neutrality and baryon-number conservation; return the largest residuals.

    The core is strict: a charge residual above ``core_tolerance`` or a baryon-number residual
    above ``baryon_tolerance`` (anywhere) raises ValueError. In the crust the charge is
    compared relative to ye. Tabulated equations of state satisfy ye = xh zz/aa to
    rounding except in the last rows before the crust-core interface, where the nuclear
    description gives way to the core (up to ~15% for BSk22 and BSk26), so only a residual
    above ``crust_charge_warning`` produces a UserWarning. The electron fraction of the
    table is used in any case.

    Returns a dict with ``core_charge``, ``baryon`` and ``crust_charge_relative`` (the
    largest absolute residuals; 0 when the region is empty).
    """
    r = composition_residuals(points)
    crust = points.xh > 0
    worst = dict(
        core_charge=float(np.max(np.abs(r["charge"][~crust]), initial=0.0)),
        baryon=float(np.max(np.abs(r["baryon"]), initial=0.0)),
        crust_charge_relative=float(np.max(r["charge_relative"][crust], initial=0.0)),
    )
    problems = []
    for label, bad, value in (
        ("charge neutrality in the core", (~crust) & (np.abs(r["charge"]) > core_tolerance), r["charge"]),
        ("baryon number", np.abs(r["baryon"]) > baryon_tolerance, r["baryon"]),
    ):  # fmt: skip
        if bad.any():
            i = int(np.argmax(np.where(bad, np.abs(value), 0.0)))
            problems.append(f"{label}: residual {value[i]:.3e} at point {i}")
    if problems:
        raise ValueError("Inconsistent composition: " + "; ".join(problems))
    if worst["crust_charge_relative"] > crust_charge_warning:
        warnings.warn(
            f"Crust charge neutrality: ye differs from xh*zz/aa by {worst['crust_charge_relative']:.0%} "
            f"(more than {crust_charge_warning:.0%}); the electron fraction of the table is used",
            UserWarning,
            stacklevel=2,
        )
    return worst


@dataclass(frozen=True)
class Background:
    """Everything the evolution needs that does not depend on time (see module docstring)."""

    faces: Points
    cells: Points
    radius_km: np.ndarray
    elam: np.ndarray
    g14: float
    envelope_model: str = "gudmundsson"
    surface_radius_km: float = None
    surface_ephi: float = None
    description: str = ""
    face_table_exact: bool = True  # False: cells/masses cannot be rebuilt from the face table

    def __post_init__(self):
        n = len(self.faces)
        if n < 3 or len(self.cells) != n - 1:
            raise ValueError("A background needs N >= 3 faces and N-1 cells")
        r = _readonly(self.radius_km, "radius_km", n)
        el = _readonly(self.elam, "elam", n)
        if r[0] != 0 or np.any(np.diff(r) <= 0):
            raise ValueError("radius_km must start at 0 and increase")
        if np.any(el <= 0):
            raise ValueError("elam must be positive")
        if not (np.isfinite(self.g14) and self.g14 > 0):
            raise ValueError("g14 must be positive and finite")
        if self.envelope_model not in ENVELOPE_MODELS:
            raise ValueError(f"envelope_model must be one of {', '.join(ENVELOPE_MODELS)}")
        surface = (self.surface_radius_km, self.surface_ephi)
        if any(v is None for v in surface) and any(v is not None for v in surface):
            raise ValueError("surface_radius_km and surface_ephi must be given together")
        if surface[0] is not None:
            if not all(np.isfinite(v) and v > 0 for v in surface):
                raise ValueError("surface_radius_km and surface_ephi must be positive and finite")
            if surface[0] < r[-1]:
                raise ValueError("The surface cannot lie below the last face")
            object.__setattr__(self, "surface_radius_km", float(surface[0]))
            object.__setattr__(self, "surface_ephi", float(surface[1]))
        object.__setattr__(self, "radius_km", r)
        object.__setattr__(self, "elam", el)
        object.__setattr__(self, "g14", float(self.g14))

    @property
    def n_faces(self):
        return len(self.faces)

    @property
    def n_cells(self):
        return len(self.cells)

    @property
    def surface_factor(self):
        """(R_s exp(Phi_s))^2 / (R_b exp(Phi_b))^2: surface over last face (1 if no surface given).

        The luminosity at infinity is conserved across the static envelope between the
        last face (base b) and the outer edge (surface s): L = 4 pi R^2 exp(2 Phi) F is
        the same at both. The envelope relation gives the surface flux F_s, so the
        luminosity and the flux leaving the base are the base-geometry values times this factor.
        """
        if self.surface_radius_km is None:
            return 1.0
        surface = self.surface_radius_km * self.surface_ephi
        return float((surface / (self.radius_km[-1] * self.faces.ephi[-1])) ** 2)

    # ----------------------------------------------------------------- table
    @classmethod
    def from_face_table(
        cls,
        table,
        *,
        g14,
        envelope_model="gudmundsson",
        derived=derived_fermi_quantities,
        surface_radius_km=None,
        surface_ephi=None,
        description="",
    ):
        """Build the background from a face table (see the module docstring).

        ``table`` maps each column name of ``FACE_TABLE_COLUMNS`` to an array of
        N face values. Cell values are the averages of the two neighboring faces
        (critical temperatures included: they are averaged after being computed
        at the faces). ``derived`` supplies the Fermi momenta and effective masses.
        """
        missing = [c for c in FACE_TABLE_COLUMNS if c not in table and c not in OPTIONAL_COLUMNS]
        if missing:
            raise ValueError(f"Face table lacks the columns: {', '.join(missing)}")
        f = {
            c: (
                np.asarray(table[c], dtype=float)
                if c in table
                else np.zeros(len(np.asarray(table["rho"])))
            )
            for c in FACE_TABLE_COLUMNS
        }

        if np.any((f["ymu"] < 0) | (f["ymu"] > 1)):
            raise ValueError("ymu must lie in [0, 1]")

        def complete(values, rho, yp, ymu):
            values = dict(values)
            values.setdefault("kf_p", fermi_momentum(rho, yp))
            values["kf_mu"] = fermi_momentum(rho, ymu)
            return values

        faces = Points(
            **{k: f[k] for k in AVERAGED_FIELDS},
            **complete(derived(f["rho"], f["ye"], f["yn"], f["yp"]), f["rho"], f["yp"], f["ymu"]),
        )
        c = {k: 0.5 * (f[k][1:] + f[k][:-1]) for k in AVERAGED_FIELDS}
        cells = Points(
            **c,
            **complete(derived(c["rho"], c["ye"], c["yn"], c["yp"]), c["rho"], c["yp"], c["ymu"]),
        )
        return cls(
            faces,
            cells,
            f["radius_km"],
            f["elam"],
            g14,
            envelope_model,
            surface_radius_km,
            surface_ephi,
            description,
        )

    def face_table(self):
        """The primary face table: the complete description of this background."""
        table = {"radius_km": self.radius_km, "elam": self.elam}
        table.update({k: getattr(self.faces, k) for k in AVERAGED_FIELDS})
        return {c: table[c] for c in FACE_TABLE_COLUMNS}

    def to_csv(self, path):
        """Write the face table as CSV (exact round trip; header lines start with '#')."""
        table = self.face_table()
        lines = [f"g14 = {self.g14!r}", f"envelope_model = {self.envelope_model}"]
        if self.surface_radius_km is not None:
            lines += [
                f"surface_radius_km = {self.surface_radius_km!r}",
                f"surface_ephi = {self.surface_ephi!r}",
            ]
        if self.description:
            lines.append("description = " + " ".join(self.description.split()))
        lines += [f"exact = {self.face_table_exact}", f"columns: {','.join(FACE_TABLE_COLUMNS)}"]
        header = "\n".join(lines)
        data = np.column_stack([table[c] for c in FACE_TABLE_COLUMNS])
        np.savetxt(path, data, delimiter=",", fmt="%.17g", header=header, comments="# ")

    @classmethod
    def from_csv(cls, path, *, derived=derived_fermi_quantities):
        """Read a face table written by ``to_csv`` (or prepared by hand in the same format)."""
        meta = {}
        for line in Path(path).read_text().splitlines():
            if not line.startswith("#"):
                break
            key, _, value = line[1:].partition("=" if "=" in line else ":")
            meta[key.strip()] = value.strip()
        if "columns" not in meta or "g14" not in meta:
            raise ValueError("The CSV header must define 'g14 = ...' and 'columns: ...'")
        if meta.get("exact", "True") != "True":
            raise ValueError(
                "This face table was written for inspection only: its cell values and effective "
                "masses come from the equation of state and cannot be rebuilt from the faces"
            )
        columns = [c.strip() for c in meta["columns"].split(",")]
        data = np.loadtxt(path, delimiter=",", ndmin=2)
        if data.shape[1] != len(columns):
            raise ValueError("The number of columns does not match the header")
        return cls.from_face_table(
            dict(zip(columns, data.T)),
            g14=float(meta["g14"]),
            envelope_model=meta.get("envelope_model", "gudmundsson"),
            derived=derived,
            surface_radius_km=(
                float(meta["surface_radius_km"]) if "surface_radius_km" in meta else None
            ),
            surface_ephi=float(meta["surface_ephi"]) if "surface_ephi" in meta else None,
            description=meta.get("description", ""),
        )
