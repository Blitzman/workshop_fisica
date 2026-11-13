"""Placement of the faces of the cooling mesh in the variable x = ln(P/Pc).

The faces are placed by marching outward from the center (x = 0, r = 0) to the
base of the cooled domain. The step in ln P is ``dlnp_core`` above the crust-core
interface and ``dlnp_crust`` below it; the step shrinks linearly from one to the
other over ``blend`` (in ln P) before reaching the interface. Optional caps limit
the radial size of the cells. Two faces are exact: the crust-core interface and the
base. When the regular step would leave a very short cell next to one of them, the
previous face is moved onto it instead, so there are no sliver cells.

``refine`` divides every spacing and cap, so ``refine = 2`` doubles the resolution
(used for convergence studies).

Presets: A uniform ln P with step 0.101;
B crust step 0.05; C = B with core cells limited to 300 m; D = twice as fine as C.
"""

from dataclasses import dataclass, replace

import numpy as np

_MAX_FACES = 200000


@dataclass(frozen=True)
class MeshSpec:
    """Mesh parameters (see the module docstring)."""

    name: str = "A"
    dlnp_core: float = 0.101
    dlnp_crust: float = 0.101
    blend: float = 0.0
    dr_core_max_m: float = None
    dr_crust_max_m: float = None
    refine: float = 1.0

    def __post_init__(self):
        for name in ("dlnp_core", "dlnp_crust", "refine"):
            v = getattr(self, name)
            if not (np.isfinite(v) and v > 0):
                raise ValueError(f"{name} must be positive and finite")
        if not (np.isfinite(self.blend) and self.blend >= 0):
            raise ValueError("blend must be non-negative and finite")
        for name in ("dr_core_max_m", "dr_crust_max_m"):
            v = getattr(self, name)
            if v is not None and not (np.isfinite(v) and v > 0):
                raise ValueError(f"{name} must be positive and finite (or omitted)")


PRESETS = {
    "A": MeshSpec("A", 0.101, 0.101),
    "B": MeshSpec("B", 0.10, 0.05, blend=0.5),
    "C": MeshSpec("C", 0.10, 0.05, blend=0.5, dr_core_max_m=300.0),
    "D": MeshSpec("D", 0.05, 0.025, blend=0.5, dr_core_max_m=150.0),
}


def mesh_spec(name="A", **overrides):
    """A preset with some parameters overridden (None values are ignored)."""
    if name not in PRESETS:
        raise ValueError(f"Unknown mesh '{name}'; presets: {', '.join(PRESETS)}")
    changes = {k: v for k, v in overrides.items() if v is not None}
    unknown = set(changes) - set(MeshSpec.__dataclass_fields__)
    if unknown:
        raise ValueError(f"Unknown mesh parameters: {', '.join(sorted(unknown))}")
    return replace(PRESETS[name], **changes) if changes else PRESETS[name]


def _step(x, x_interface, spec):
    core, crust = spec.dlnp_core / spec.refine, spec.dlnp_crust / spec.refine
    if x <= x_interface:
        return crust
    if spec.blend > 0:
        return crust + (core - crust) * min((x - x_interface) / spec.blend, 1.0)
    return core


def place_faces(star, x_base, x_interface, spec):
    """Values of x = ln(P/Pc) of the faces, from the center (0) to the base (x_base).

    ``star`` is a ``tov.StarSolution`` (it provides r(x) and its inverse for the radial
    caps). ``x_interface`` is the crust-core interface (+inf for a star without core).
    """
    if not x_base < 0:
        raise ValueError("The base must lie below the center (x_base < 0)")
    targets = sorted(
        {t for t in (x_interface, x_base) if np.isfinite(t) and x_base <= t < 0}, reverse=True
    )
    faces = [0.0]
    x = 0.0
    while x > x_base + 1e-12:
        h = _step(x, x_interface, spec)
        cap = spec.dr_core_max_m if x > x_interface else spec.dr_crust_max_m
        if cap is not None:
            r = float(star.profile(x)["r_cm"][0])
            h = min(h, x - float(star.x_at_radius(r + 100.0 * cap / spec.refine)))
        if h < 1e-13:
            raise RuntimeError("The mesh step collapsed; check the radial caps")
        x_next = x - h
        for t in targets:
            if x > t + 1e-12 and x_next <= t + 1e-12:  # the step reaches or crosses a target
                if x - t < 0.5 * h and len(faces) > 1:  # too close: move the last face onto it
                    faces[-1] = t
                    x = t
                    x_next = None
                else:
                    x_next = t
                break
        if x_next is None:
            continue
        faces.append(x_next)
        x = x_next
        if len(faces) > _MAX_FACES:
            raise RuntimeError("Too many faces; check the mesh parameters")
    faces[-1] = x_base
    return np.array(faces)


def describe_mesh(star, x_faces, x_interface=np.inf):
    """Statistics of a mesh: face and cell counts, radial cell sizes [m], worst neighbour ratio."""
    x = np.asarray(x_faces)
    dr = np.diff(star.profile(x)["r_cm"]) / 100.0
    core = 0.5 * (x[1:] + x[:-1]) > x_interface

    def stats(values):
        return None if values.size == 0 else (values.min(), float(np.median(values)), values.max())

    ratio = np.maximum(dr[1:] / dr[:-1], dr[:-1] / dr[1:])
    return dict(
        n_faces=len(x),
        n_core_cells=int(core.sum()),
        n_crust_cells=int((~core).sum()),
        core_dr_m=stats(dr[core]),
        crust_dr_m=stats(dr[~core]),
        worst_neighbour_ratio=float(ratio.max()) if ratio.size else 1.0,
    )
