"""Stellar model from an equation-of-state table: builds the ``Background`` of the cooling code.

Pipeline
    1. read the CompOSE table (``eos_table``) from ``<eos root>/<eos name>/``;
    2. choose the central pressure, directly or by the target mass (``tov``);
    3. integrate the structure once, from the center to the outer edge of the table;
       exp(Phi) is normalized globally there (Schwarzschild matching);
    4. place the faces from the center to the base of the cooled domain (``mesh``);
    5. evaluate the metric and the EOS at the faces and, independently, at the middle
       (in radius) of every cell, so that no cell mixes crust and core values;
    6. attach what the evolution needs from the outer edge: the gravity g14 of the
       envelope relation and the radius and redshift of the edge.

``rho`` of the background is the baryon mass density 1.66e-24 g * n_b, so that the
number densities of the microphysics are exactly those of the table. The effective
masses are the Landau masses of the EOS in the core; in the crust, where the EOS
has none, the tables of ``effective_mass`` are used. The critical temperatures come
from the pairing models of the input file, evaluated at the cells and at the faces.

Input file options (``key = value`` lines, see ``input_model``)
    eos              name of the folder in the eos directory (required)
    mass | central_pressure   target gravitational mass [Msun] or central pressure [dyn/cm^2]
    base_density     baryon mass density of the base of the cooled domain [g/cm^3] (1e10)
    mesh             preset A, B, C or D (A); dlnp_core, dlnp_crust, blend, dr_core_max_m,
                     dr_crust_max_m, refine override single parameters (see ``mesh``)
"""

import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from background import (
    Background,
    Points,
    check_composition,
    derived_fermi_quantities,
    fermi_momentum,
)
from eos_table import CompOSETable
from mesh import MeshSpec, mesh_spec, place_faces
from state import LocalState
from superfluid import pairing_critical_temperatures
from tov import central_pressure_for_mass, integrate_star

MESH_KEYS = ("dlnp_core", "dlnp_crust", "blend", "dr_core_max_m", "dr_crust_max_m", "refine")
STRUCTURE_KEYS = ("eos", "mass", "central_pressure", "base_density", "mesh") + MESH_KEYS


def _as_tuple(value):
    """None, a number or a sequence of numbers as a tuple of floats (None stays None)."""
    if value is None:
        return None
    values = tuple(float(v) for v in np.atleast_1d(value))
    if not values:
        raise ValueError("An empty list of stars was given")
    return values


@dataclass(frozen=True)
class StructureOptions:
    """What to build: see the module docstring."""

    eos: str
    mass_solar: float = None
    central_pressure: float = None
    base_density: float = 1e10
    mesh: MeshSpec = None
    rtol: float = 1e-9
    max_step: float = 0.05

    def __post_init__(self):
        if (self.mass_solar is None) == (self.central_pressure is None):
            raise ValueError("Give exactly one of 'mass' and 'central_pressure'")
        for name in ("mass_solar", "central_pressure"):
            v = getattr(self, name)
            if v is not None and not (np.isfinite(v) and v > 0):
                raise ValueError(f"{name} must be positive and finite")
        if not (np.isfinite(self.base_density) and self.base_density > 0):
            raise ValueError("base_density must be positive and finite")
        if self.mesh is None:
            object.__setattr__(self, "mesh", mesh_spec("A"))

    @classmethod
    def list_from_config(cls, cfg):
        """One ``StructureOptions`` per star of the dictionary read by ``input_model.read_parameters``.

        ``mass`` or ``central_pressure`` (exactly one of them) may hold several values.
        """
        if "eos" not in cfg:
            raise ValueError("The input file has no 'eos' option")
        masses, pressures = _as_tuple(cfg.get("mass")), _as_tuple(cfg.get("central_pressure"))
        if (masses is None) == (pressures is None):
            raise ValueError("Give exactly one of 'mass' and 'central_pressure'")
        common = dict(
            eos=cfg["eos"],
            base_density=cfg.get("base_density", 1e10),
            mesh=mesh_spec(cfg.get("mesh", "A"), **{k: cfg.get(k) for k in MESH_KEYS}),
        )
        if masses is not None:
            return [cls(mass_solar=m, **common) for m in masses]
        return [cls(central_pressure=p, **common) for p in pressures]

    @classmethod
    def from_config(cls, cfg):
        """The options of a configuration that defines exactly one star."""
        options = cls.list_from_config(cfg)
        if len(options) != 1:
            raise ValueError(f"The input defines {len(options)} stars; use list_from_config")
        return options[0]


@dataclass
class StellarModel:
    """Result of the construction: the background and the objects it came from."""

    background: Background
    eos: CompOSETable
    star: object  # tov.StarSolution
    x_faces: np.ndarray
    x_interface: float
    x_base: float
    options: StructureOptions

    def summary(self):
        """One line describing the star and the mesh."""
        return self.background.description

    @property
    def label(self):
        """Short name of the star, from its integrated mass (for legends and folder names)."""
        return f"{self.star.mass_solar:.2f} Msun"


DIRECT_URCA_MAX_CELL_KM = 0.35


def direct_urca_boundary_widths(background):
    """Widths [km] of the cells on both sides of every edge of the Direct Urca region.

    The process is strong and switches on sharply at the threshold, so the cooling is
    only converged in space if the cells there are narrow. Empty if the region has no edge.
    """
    allowed = LocalState(1e9, background.cells).direct_urca_allowed
    widths = np.diff(background.radius_km)
    edges = np.nonzero(np.diff(allowed.astype(int)))[
        0
    ]  # cell i is allowed, cell i + 1 is not (or vice versa)
    return np.array([max(widths[i], widths[i + 1]) for i in edges])


def _warn_if_the_direct_urca_edge_is_unresolved(background):
    widths = direct_urca_boundary_widths(background)
    if widths.size and widths.max() > DIRECT_URCA_MAX_CELL_KM:
        warnings.warn(
            f"The Direct Urca threshold lies in cells up to {1e3 * widths.max():.0f} m wide: its volume is "
            "quantized by the mesh and the cooling may not be converged in space "
            "(a test star changed by ~16% in L_inf). Use a mesh with a cap on the core cells, "
            "for example mesh = C",
            UserWarning,
            stacklevel=3,
        )


def _describe(eos, star, options, g14, n_faces):
    return (
        f"EOS {eos.name}: M = {star.mass_solar:.4f} Msun, R = {star.radius_cm / 1e5:.3f} km, "
        f"Pc = {star.pressure_central:.4e} dyn/cm^2, g14 = {g14:.4f}, "
        f"{n_faces} faces (mesh {options.mesh.name}), base rho = {options.base_density:.3g} g/cm^3"
    )


def _points(state, ephi, cfg):
    """``Points`` of the background from EOS points (see the module docstring)."""
    rho = state.rho_baryon
    derived = derived_fermi_quantities(rho, state.ye, state.yn, state.yp)
    has_mass = np.isfinite(state.mstar_n) & np.isfinite(state.mstar_p)
    core = state.xh == 0.0
    if core.any() and not has_mass[core].all():  # decision D-006: no masses assumed for the core
        raise ValueError(
            "The EOS package provides no Landau effective masses for the core "
            "(eos.micro missing or incomplete); packages without them are not used"
        )
    if has_mass.any():  # the EOS provides Landau masses in the core
        for key, value in (
            ("mstar_n", state.mstar_n),
            ("transport_mstar_n", state.mstar_n),
            ("mstar_p", state.mstar_p),
            ("transport_mstar_p", state.mstar_p),
        ):
            derived[key] = np.where(has_mass, value, derived[key])
    tc = pairing_critical_temperatures(cfg, dict(rho=rho, yn=state.yn, yp=state.yp, xh=state.xh))
    return Points(
        rho=rho,
        ye=state.ye,
        yn=state.yn,
        yp=state.yp,
        xh=state.xh,
        aa=state.aa,
        zz=state.zz,
        ephi=ephi,
        ymu=state.ymu,
        kf_mu=fermi_momentum(rho, state.ymu),
        **tc,
        **derived,
    )


def build_stellar_model(options, eos_root, cfg):
    """Build the stellar model; ``cfg`` is the dictionary of ``input_model.read_parameters``
    (pairing switches and gap models, envelope model)."""
    eos = CompOSETable(Path(eos_root) / options.eos)
    if options.mass_solar is not None:
        pc = central_pressure_for_mass(eos, options.mass_solar)
    else:
        pc = options.central_pressure
    star = integrate_star(eos, pc, rtol=options.rtol, max_step=options.max_step)
    p_base = float(eos.pressure_at_rho_baryon(options.base_density))
    if not star.pressure_outer < p_base < pc:
        raise ValueError(
            "The base density must lie inside the star, above the outer edge of the table"
        )
    x_base = float(np.log(p_base / pc))
    has_core = np.isfinite(eos.pressure_core) and eos.pressure_core < pc
    x_interface = float(np.log(eos.pressure_core / pc)) if has_core else np.inf
    x_faces = place_faces(star, x_base, x_interface, options.mesh)

    # Faces: the interface and the base carry exactly the pressures of the table.
    p_faces = pc * np.exp(x_faces)
    if has_core:
        p_faces[np.isclose(x_faces, x_interface, rtol=0.0, atol=1e-13)] = eos.pressure_core
    p_faces[-1] = p_base
    face_profile = star.profile(x_faces)
    faces = _points(eos.at_pressure(p_faces), face_profile["ephi"], cfg)

    # Cells: EOS and metric at the middle (in radius) of each cell.
    r_mid = 0.5 * (face_profile["r_cm"][1:] + face_profile["r_cm"][:-1])
    x_cells = star.x_at_radius(r_mid)
    cell_profile = star.profile(x_cells)
    cells = _points(eos.at_pressure(star.pressure_at_x(x_cells)), cell_profile["ephi"], cfg)
    for points in (faces, cells):  # charge neutrality and baryon-number conservation of the EOS
        check_composition(points)

    g14 = star.surface_gravity_cm_s2 / 1e14
    background = Background(
        faces=faces,
        cells=cells,
        radius_km=face_profile["r_cm"] / 1e5,
        elam=face_profile["elam"],
        g14=g14,
        envelope_model=cfg["envelope_model"],
        surface_radius_km=star.radius_cm / 1e5,
        surface_ephi=star.surface_ephi,
        description=_describe(eos, star, options, g14, len(x_faces)),
        face_table_exact=False,
    )
    _warn_if_the_direct_urca_edge_is_unresolved(background)
    return StellarModel(background, eos, star, x_faces, x_interface, x_base, options)


def build_stellar_models(cfg, eos_root):
    """One ``StellarModel`` per star defined by the input options ``cfg`` (see ``list_from_config``)."""
    return [build_stellar_model(o, eos_root, cfg) for o in StructureOptions.list_from_config(cfg)]


def build_background(cfg, eos_root):
    """The ``Background`` of the stellar model described by the input options ``cfg`` (one star)."""
    return build_stellar_model(StructureOptions.from_config(cfg), eos_root, cfg).background
