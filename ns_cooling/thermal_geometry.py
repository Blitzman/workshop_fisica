"""Finite-volume geometry with general-relativistic metric factors.

N faces bound N-1 cells. Solver vectors contain only the physical cells.
Distances are in km.
"""

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class ThermalGeometry:
    radius_faces: np.ndarray
    radius_cells: np.ndarray
    widths: np.ndarray
    areas: np.ndarray
    volumes: np.ndarray
    face_left_weight: np.ndarray
    center_distances: np.ndarray
    ephi_faces: np.ndarray
    elam_faces: np.ndarray

    @property
    def n_cells(self):
        return len(self.widths)


def thermal_geometry(background):
    """Finite-volume geometry of a ``background.Background`` (radii, redshift and Lambda at the faces)."""
    r = np.array(background.radius_km, dtype=float, copy=True)
    ep = np.array(background.faces.ephi, dtype=float, copy=True)
    el = np.array(background.elam, dtype=float, copy=True)
    if r.ndim != 1 or len(r) < 3 or ep.shape != r.shape or el.shape != r.shape:
        raise ValueError("Geometry requires matching 1D face arrays, at least 3 faces")
    if not all(np.isfinite(v).all() for v in (r, ep, el)):
        raise ValueError("Non-finite geometry")
    if r[0] != 0 or np.any(np.diff(r) <= 0) or np.any(ep <= 0) or np.any(el <= 0):
        raise ValueError("Require r[0]=0, increasing radii and positive metric factors")
    h = np.diff(r)
    centers = 0.5 * (r[1:] + r[:-1])
    area = 4 * np.pi * r**2
    # Deliberately use OUTER-face elam, not its cell average.
    volume = el[1:] * np.pi * (r[1:] ** 3 - r[:-1] ** 3) * 4 / 3
    weight = h[1:] / (h[:-1] + h[1:])
    distances = 0.5 * (h[:-1] + h[1:])
    arrays = (r, centers, h, area, volume, weight, distances, ep, el)
    for a in arrays:
        a.setflags(write=False)
    return ThermalGeometry(*arrays)


def cell_temperatures(Tinf_8, geometry):
    u = np.asarray(Tinf_8, dtype=float)
    if u.shape != (geometry.n_cells,) or not np.isfinite(u).all() or np.any(u <= 0):
        raise ValueError("Tinf_8 requires one finite positive value per physical cell")
    return u


def interior_face_temperatures(Tinf_8, geometry):
    """Return local T/1e8 at the interior faces.

    Cell temperatures are linearly interpolated to the faces (weights from the
    cell widths) and divided by the redshift factor of the face.
    """
    u = cell_temperatures(Tinf_8, geometry)
    w = geometry.face_left_weight
    return ((1 - w) * u[1:] + w * u[:-1]) / geometry.ephi_faces[1:-1]
