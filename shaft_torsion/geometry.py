"""Step 1 - CAD import, shaft axis, splitting and section inertias.

Coordinates along the shaft are expressed by ``s``: the signed distance from
``axis.origin`` along ``axis.direction``.  In the GUI positions are shown as
distance from the shaft start (``s - s_min``).
"""
from __future__ import annotations

import os
import tempfile
from dataclasses import dataclass, field, asdict

import numpy as np

from .gmsh_session import gmsh_session


# --------------------------------------------------------------------------- #
# data classes
# --------------------------------------------------------------------------- #
@dataclass
class Axis:
    origin: np.ndarray            # point on the axis [m]
    direction: np.ndarray         # unit vector

    def __post_init__(self):
        self.origin = np.asarray(self.origin, float)
        d = np.asarray(self.direction, float)
        self.direction = d / np.linalg.norm(d)

    # --- helpers ---------------------------------------------------------- #
    def s_of(self, pts):
        """Axial coordinate of points (N,3)."""
        return (np.asarray(pts, float) - self.origin) @ self.direction

    def radial(self, pts):
        """Radial vectors (N,3) from the axis to the points."""
        p = np.asarray(pts, float) - self.origin
        return p - np.outer(p @ self.direction, self.direction)

    def perpendicular(self):
        """A unit vector perpendicular to the axis (reference for pointers)."""
        a = self.direction
        trial = np.array([0.0, 0.0, 1.0]) if abs(a[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
        v = np.cross(a, trial)
        return v / np.linalg.norm(v)

    def rotate(self, pts, angles):
        """Rotate points (N,3) about the axis by angles (N,) or scalar [rad]."""
        p = np.asarray(pts, float) - self.origin
        a = self.direction
        along = np.outer(p @ a, a)
        rad = p - along
        c = np.cos(angles)[..., None] if np.ndim(angles) else np.cos(angles)
        s = np.sin(angles)[..., None] if np.ndim(angles) else np.sin(angles)
        return self.origin + along + rad * c + np.cross(a, rad) * s

    def to_dict(self):
        return {"origin": self.origin.tolist(), "direction": self.direction.tolist()}


@dataclass
class Section:
    index: int
    s_start: float                # [m]
    s_end: float                  # [m]
    volume: float                 # [m^3]
    mass: float                   # [kg]
    s_com: float                  # axial position of the centre of mass [m]
    J: float                      # mass moment of inertia about the shaft axis [kg m^2]
    volume_tags: list = field(default_factory=list)

    @property
    def s_mid(self):
        return 0.5 * (self.s_start + self.s_end)

    @property
    def length(self):
        return self.s_end - self.s_start

    def to_dict(self):
        d = asdict(self)
        d["s_mid"] = self.s_mid
        return d


@dataclass
class ShaftGeometry:
    """Everything step 1 produces."""
    step_path: str
    axis: Axis
    s_min: float
    s_max: float
    density: float
    cuts: list                    # axial positions s of the cut planes [m]
    sections: list                # list[Section]
    split_brep: str               # path of the split model (BREP, metres)
    total_mass: float = 0.0
    total_J: float = 0.0
    surface: dict = field(default_factory=dict)   # section index -> (points, triangles)

    @property
    def length(self):
        return self.s_max - self.s_min

    def mid_positions(self, use_com=False):
        return np.array([sec.s_com if use_com else sec.s_mid for sec in self.sections])


# --------------------------------------------------------------------------- #
# mass properties
# --------------------------------------------------------------------------- #
def _volume_props(occ, tag):
    vol = occ.getMass(3, tag)
    com = np.array(occ.getCenterOfMass(3, tag))
    inertia = np.array(occ.getMatrixOfInertia(3, tag)).reshape(3, 3)  # about COM, per unit density
    return vol, com, inertia


def _inertia_about_axis(vol, com, i_com, axis: Axis):
    a = axis.direction
    d = axis.radial(com[None, :])[0]
    return float(a @ i_com @ a + vol * (d @ d))


def _combined_tensor(props):
    vols = np.array([p[0] for p in props])
    coms = np.array([p[1] for p in props])
    vtot = vols.sum()
    com = (vols[:, None] * coms).sum(0) / vtot
    tensor = np.zeros((3, 3))
    for v, c, i in props:
        d = c - com
        tensor += i + v * ((d @ d) * np.eye(3) - np.outer(d, d))
    return vtot, com, tensor


def principal_axis(com, tensor):
    """Axis of an (almost) axisymmetric body: principal direction whose moment
    differs most from the other two (the two others are equal for a body of
    revolution).  Works for long shafts and short discs alike."""
    w, v = np.linalg.eigh(tensor)
    best, score = 0, -1.0
    for k in range(3):
        j, l = [m for m in range(3) if m != k]
        spread = abs(w[j] - w[l]) / max(abs(w).max(), 1e-300)
        dist = abs(w[k] - 0.5 * (w[j] + w[l])) / max(abs(w).max(), 1e-300)
        sc = dist - 2 * spread
        if sc > score:
            best, score = k, sc
    d = v[:, best]
    # snap to a global axis when the CAD is (almost) aligned with it
    k = int(np.argmax(abs(d)))
    if abs(d[k]) > np.cos(np.radians(1.0)):
        d = np.eye(3)[k] * np.sign(d[k])
    # make the dominant component positive for reproducibility
    if d[np.argmax(abs(d))] < 0:
        d = -d
    return Axis(com, d)


# --------------------------------------------------------------------------- #
# STEP import
# --------------------------------------------------------------------------- #
@dataclass
class StepInfo:
    path: str
    n_volumes: int
    bbox: tuple
    volume: float
    com: np.ndarray
    tensor: np.ndarray            # inertia tensor about COM per unit density
    auto_axis: Axis
    surface: tuple                # (points, triangles) for display, metres


def _import(gm, path):
    occ = gm.model.occ
    ext = os.path.splitext(path)[1].lower()
    dimtags = occ.importShapes(path, highestDimOnly=True)
    occ.synchronize()
    vols = [t for d, t in dimtags if d == 3]
    if not vols:
        vols = [t for d, t in gm.model.getEntities(3)]
    if not vols:
        raise ValueError(f"No solid found in {os.path.basename(path)} ({ext}). "
                         "The STEP file must contain closed solids.")
    return vols


def _tessellate(gm, volume_groups, size):
    """Surface triangulation per group of volumes -> {key: (pts, tris)}."""
    gm.option.setNumber("Mesh.MeshSizeMax", size)
    gm.option.setNumber("Mesh.MeshSizeMin", size / 20)
    gm.option.setNumber("Mesh.MeshSizeFromCurvature", 24)
    gm.option.setNumber("Mesh.Algorithm", 6)
    gm.model.mesh.generate(2)
    ntags, ncoord, _ = gm.model.mesh.getNodes()
    coords = ncoord.reshape(-1, 3)
    index = np.zeros(int(ntags.max()) + 1, dtype=np.int64)
    index[ntags.astype(np.int64)] = np.arange(len(ntags))
    out = {}
    for key, vols in volume_groups.items():
        tris = []
        surfs = set()
        for v in vols:
            for d, s in gm.model.getBoundary([(3, v)], oriented=False, recursive=False):
                surfs.add(abs(s))
        for s in surfs:
            types, _, nodes = gm.model.mesh.getElements(2, s)
            for t, n in zip(types, nodes):
                if t == 2:   # 3-node triangle
                    tris.append(index[n.astype(np.int64)].reshape(-1, 3))
        if not tris:
            continue
        tris = np.vstack(tris)
        used, inv = np.unique(tris, return_inverse=True)
        out[key] = (coords[used].copy(), inv.reshape(-1, 3))
    return out


def load_step(path, display_size=None) -> StepInfo:
    """Read a STEP (or BREP) file, compute global mass properties and the
    automatic shaft axis."""
    with gmsh_session() as gm:
        vols = _import(gm, path)
        occ = gm.model.occ
        props = [_volume_props(occ, v) for v in vols]
        vtot, com, tensor = _combined_tensor(props)
        bb = gm.model.getBoundingBox(-1, -1)
        diag = np.linalg.norm(np.array(bb[3:]) - np.array(bb[:3]))
        surf = _tessellate(gm, {0: vols}, display_size or diag / 60).get(0)
    axis = principal_axis(com, tensor)
    return StepInfo(path, len(vols), tuple(bb), vtot, com, tensor, axis, surf)


def axial_range(points, axis: Axis):
    s = axis.s_of(points)
    return float(s.min()), float(s.max())


# --------------------------------------------------------------------------- #
# splitting
# --------------------------------------------------------------------------- #
def split_model(step_path, axis: Axis, cuts, density, s_range=None,
                out_brep=None, display_size=None, progress=None) -> ShaftGeometry:
    """Cut the CAD model with planes normal to ``axis`` at the axial positions
    ``cuts`` (absolute ``s`` in metres) and compute the inertia of every
    section about the shaft axis.

    All solids of the file are fragmented together, so touching parts of an
    assembly share their contact faces (bonded) in the later FEM mesh.
    """
    say = progress or (lambda msg: None)
    with gmsh_session() as gm:
        occ = gm.model.occ
        say("Importing CAD ...")
        vols = _import(gm, step_path)
        bb = gm.model.getBoundingBox(-1, -1)
        diag = float(np.linalg.norm(np.array(bb[3:]) - np.array(bb[:3])))
        corners = np.array([[x, y, z] for x in (bb[0], bb[3]) for y in (bb[1], bb[4])
                            for z in (bb[2], bb[5])])
        if s_range is None:
            s_range = axial_range(corners, axis)
        s_min, s_max = s_range
        cuts = sorted(float(c) for c in cuts)
        valid = [c for c in cuts if s_min + 1e-9 * diag < c < s_max - 1e-9 * diag]
        if len(valid) != len(cuts):
            raise ValueError("Cut positions must lie strictly inside the shaft "
                             f"({0.0:.1f} .. {(s_max - s_min) * 1e3:.1f} mm from the start).")
        tools = []
        xdir = axis.perpendicular()
        for c in cuts:
            p = axis.origin + c * axis.direction
            tag = occ.addDisk(*p, 2 * diag, 2 * diag, zAxis=axis.direction.tolist(),
                              xAxis=xdir.tolist())
            tools.append((2, tag))
        say("Splitting ...")
        if tools or len(vols) > 1:
            occ.fragment([(3, v) for v in vols], tools)
            occ.synchronize()
            # remove the parts of the cutting disks that lie outside the solid
            free = [(2, t) for d, t in gm.model.getEntities(2)
                    if len(gm.model.getAdjacencies(2, t)[0]) == 0]
            if free:
                occ.remove(free, recursive=True)
                occ.synchronize()
        vols = [t for d, t in gm.model.getEntities(3)]

        # ---- assign volumes to sections by the axial position of their COM
        bounds = [s_min] + cuts + [s_max]
        nsec = len(bounds) - 1
        groups = {i: [] for i in range(nsec)}
        acc = {i: [0.0, 0.0, 0.0] for i in range(nsec)}  # volume, J/rho, s_com*vol
        say("Mass properties ...")
        for v in vols:
            vol, com, icom = _volume_props(occ, v)
            sc = float(axis.s_of(com[None])[0])
            i = int(np.clip(np.searchsorted(cuts, sc), 0, nsec - 1))
            groups[i].append(v)
            acc[i][0] += vol
            acc[i][1] += _inertia_about_axis(vol, com, icom, axis)
            acc[i][2] += vol * sc
        sections = []
        for i in range(nsec):
            vol, jr, sv = acc[i]
            sections.append(Section(
                index=i, s_start=bounds[i], s_end=bounds[i + 1], volume=vol,
                mass=vol * density, s_com=(sv / vol) if vol > 0 else 0.5 * (bounds[i] + bounds[i + 1]),
                J=jr * density, volume_tags=groups[i]))
        empty = [s.index + 1 for s in sections if s.volume <= 0]
        if empty:
            raise ValueError(f"Section(s) {empty} contain no material - check the cut positions.")

        if out_brep is None:
            fd, out_brep = tempfile.mkstemp(suffix=".brep", prefix="shaft_split_")
            os.close(fd)
        gm.write(out_brep)
        say("Tessellating for display ...")
        surface = _tessellate(gm, groups, display_size or diag / 60)

    geo = ShaftGeometry(step_path=step_path, axis=axis, s_min=s_min, s_max=s_max,
                        density=density, cuts=cuts, sections=sections, split_brep=out_brep,
                        total_mass=sum(s.mass for s in sections),
                        total_J=sum(s.J for s in sections), surface=surface)
    return geo


def import_split_brep(gm, geo: ShaftGeometry):
    """Re-import the split model into an (open) gmsh session and return the
    volume tags grouped by section index (assignment by centre of mass)."""
    occ = gm.model.occ
    occ.importShapes(geo.split_brep, highestDimOnly=True)
    occ.synchronize()
    vols = [t for d, t in gm.model.getEntities(3)]
    if len(vols) > 1:
        # make sure neighbouring volumes share their faces (conformal mesh)
        occ.fragment([(3, v) for v in vols], [])
        occ.synchronize()
        vols = [t for d, t in gm.model.getEntities(3)]
    groups = {s.index: [] for s in geo.sections}
    for v in vols:
        com = np.array(occ.getCenterOfMass(3, v))
        sc = float(geo.axis.s_of(com[None])[0])
        i = int(np.clip(np.searchsorted(geo.cuts, sc), 0, len(geo.sections) - 1))
        groups[i].append(v)
    return groups
