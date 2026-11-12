"""Step 2 - static torsion FEM of the split shaft.

* 3-D solid mesh (linear or quadratic tetrahedra, gmsh) of the split model.
* One end face is clamped (all DOFs = 0), the other end face is rigidly
  coupled (RBE2-like) and twisted about the shaft axis.  The reaction torque
  of that coupling is scaled to the requested torque, i.e. the result is the
  same as applying the torque through a rigid end plate.
* The rotation of the shaft at any axial position is the least-squares rigid
  rotation of the nodes in a thin slice around that position.
* Spring stiffness between two lumped inertias:  k = T / |theta_j - theta_i|.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .geometry import ShaftGeometry, import_split_brep
from .gmsh_session import gmsh_session


@dataclass
class FemSettings:
    E: float = 210e9              # Young's modulus [Pa]
    nu: float = 0.3               # Poisson ratio
    torque: float = 1000.0        # applied torque [N m]
    mesh_size: float = 0.0        # [m]; 0 -> automatic (length/40, limited by diameter)
    order: int = 2                # 1 = linear tets, 2 = quadratic tets (recommended)
    fixed_end: str = "start"      # "start" (s_min) or "end" (s_max)
    end_band: float = 0.0         # [m] nodes within this distance of an end belong to it; 0 -> auto
    mid_from_com: bool = False    # lumped inertia position: section middle (False) or COM (True)

    @property
    def G(self):
        return self.E / (2 * (1 + self.nu))


@dataclass
class TorsionResult:
    settings: FemSettings
    nodes: np.ndarray             # (N,3) [m]
    tets: np.ndarray              # (Ne,4) corner nodes, for display
    disp: np.ndarray              # (N,3) [m]
    node_theta: np.ndarray        # (N,) rotation about the axis [rad]
    s_mid: np.ndarray             # axial positions of the lumped inertias [m]
    theta_mid: np.ndarray         # rotation there [rad]
    s_cut: np.ndarray             # section boundaries (start, cuts, end) [m]
    theta_cut: np.ndarray         # rotation there [rad]
    k_spring: np.ndarray          # (n-1,) stiffness between inertia i and i+1 [N m/rad]
    k_section: np.ndarray         # (n,) stiffness of each section piece [N m/rad]
    k_ground: float               # clamp -> nearest lumped inertia [N m/rad]
    k_total: float                # clamp -> loaded end [N m/rad]
    ground_index: int             # inertia index next to the clamp
    n_dof: int = 0
    solve_time: float = 0.0
    info: dict = field(default_factory=dict)


# --------------------------------------------------------------------------- #
# meshing
# --------------------------------------------------------------------------- #
def _mesh(gm, geo: ShaftGeometry, st: FemSettings, say):
    import_split_brep(gm, geo)
    bb = gm.model.getBoundingBox(-1, -1)
    diag = float(np.linalg.norm(np.array(bb[3:]) - np.array(bb[:3])))
    h = st.mesh_size or auto_mesh_size(geo)
    gm.option.setNumber("Mesh.MeshSizeMax", h)
    gm.option.setNumber("Mesh.MeshSizeMin", h / 10)
    gm.option.setNumber("Mesh.MeshSizeFromCurvature", 16)
    gm.option.setNumber("Mesh.Algorithm", 6)
    gm.option.setNumber("Mesh.Algorithm3D", 10 if st.order == 1 else 1)  # HXT / Delaunay
    gm.option.setNumber("Mesh.ElementOrder", int(st.order))
    gm.option.setNumber("Mesh.HighOrderOptimize", 0)
    gm.option.setNumber("Mesh.SecondOrderLinear", 0)
    gm.option.setNumber("Mesh.Optimize", 1)
    say(f"Meshing (element size {h * 1e3:.2f} mm, order {st.order}) ...")
    gm.model.mesh.generate(3)
    etype = gm.model.mesh.getElementType("Tetrahedron", int(st.order))
    etags, enodes = gm.model.mesh.getElementsByType(etype)
    nn = 4 if st.order == 1 else 10
    enodes = enodes.astype(np.int64).reshape(-1, nn)
    ntags, ncoord, _ = gm.model.mesh.getNodes()
    ntags = ntags.astype(np.int64)
    coords = ncoord.reshape(-1, 3)
    # keep only nodes used by tetrahedra, renumber 0..N-1
    used = np.unique(enodes)
    lookup = np.full(ntags.max() + 1, -1, dtype=np.int64)
    lookup[ntags] = np.arange(len(ntags))
    nodes = coords[lookup[used]]
    renum = np.full(ntags.max() + 1, -1, dtype=np.int64)
    renum[used] = np.arange(len(used))
    conn = renum[enodes]
    # reference-element data
    lc, w = gm.model.mesh.getIntegrationPoints(etype, "Gauss2" if st.order == 2 else "Gauss1")
    _, bf, _ = gm.model.mesh.getBasisFunctions(etype, lc, "GradLagrange")
    dN = np.asarray(bf).reshape(len(w), nn, 3)
    return nodes, conn, dN, np.asarray(w), h, diag


def auto_mesh_size(geo: ShaftGeometry):
    """Default element size: length/40, but at least ~4 elements across the
    smallest section diameter (estimated from J and mass)."""
    L = geo.length
    radii = []
    for sec in geo.sections:
        if sec.mass > 0:
            radii.append(np.sqrt(2 * sec.J / sec.mass))   # radius of an equivalent solid cylinder
    rmin = min(radii) if radii else L / 10
    return float(min(L / 40, rmin / 2))


# --------------------------------------------------------------------------- #
# stiffness matrix
# --------------------------------------------------------------------------- #
def _elasticity(E, nu):
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    D = np.zeros((6, 6))
    D[:3, :3] = lam
    D[np.arange(3), np.arange(3)] += 2 * mu
    D[3:, 3:] = np.eye(3) * mu
    return D


def assemble_stiffness(nodes, conn, dN, w, E, nu, chunk=8000, say=None):
    nn = conn.shape[1]
    ndof_e = 3 * nn
    D = _elasticity(E, nu)
    ndof = 3 * len(nodes)
    rows_l, cols_l, vals_l = [], [], []
    ne = len(conn)
    vol_total = 0.0
    bad = 0
    for c0 in range(0, ne, chunk):
        cc = conn[c0:c0 + chunk]
        X = nodes[cc]                                     # (m, nn, 3)
        Ke = np.zeros((len(cc), ndof_e, ndof_e))
        for g in range(len(w)):
            J = np.einsum("ni,mnj->mij", dN[g], X)        # (m,3,3)  dX_j/dxi_i
            det = np.linalg.det(J)
            bad += int(np.sum(det <= 0))
            Jinv = np.linalg.inv(J)
            G = np.einsum("mij,nj->mni", Jinv, dN[g])     # (m, nn, 3)  dN/dx
            B = np.zeros((len(cc), 6, ndof_e))
            gx, gy, gz = G[..., 0], G[..., 1], G[..., 2]
            B[:, 0, 0::3] = gx
            B[:, 1, 1::3] = gy
            B[:, 2, 2::3] = gz
            B[:, 3, 0::3] = gy
            B[:, 3, 1::3] = gx
            B[:, 4, 1::3] = gz
            B[:, 4, 2::3] = gy
            B[:, 5, 0::3] = gz
            B[:, 5, 2::3] = gx
            wd = w[g] * np.abs(det)
            vol_total += wd.sum()
            Ke += np.einsum("m,mki,kl,mlj->mij", wd, B, D, B, optimize=True)
        dofs = (3 * cc[:, :, None] + np.arange(3)).reshape(len(cc), ndof_e)
        rows_l.append(np.repeat(dofs, ndof_e, axis=1).ravel())
        cols_l.append(np.tile(dofs, (1, ndof_e)).ravel())
        vals_l.append(Ke.ravel())
        if say:
            say(f"Assembling stiffness {min(c0 + chunk, ne)}/{ne} elements ...")
    K = sp.coo_matrix((np.concatenate(vals_l), (np.concatenate(rows_l), np.concatenate(cols_l))),
                      shape=(ndof, ndof)).tocsr()
    K.sum_duplicates()
    return K, vol_total, bad


def _rigid_modes(xyz):
    """6 rigid-body modes (near null space) for algebraic multigrid."""
    n = len(xyz)
    B = np.zeros((3 * n, 6))
    x, y, z = (xyz - xyz.mean(0)).T
    for k in range(3):
        B[k::3, k] = 1.0
    B[0::3, 3], B[1::3, 3] = -y, x
    B[1::3, 4], B[2::3, 4] = -z, y
    B[0::3, 5], B[2::3, 5] = z, -x
    return B


def _solve(A, b, xyz=None, say=None):
    """Sparse SPD solve: PARDISO (pypardiso) if installed, else algebraic
    multigrid preconditioned CG (pyamg) for big systems, else SuperLU."""
    say = say or (lambda m: None)
    try:
        import pypardiso  # optional, fastest
        say("Solving (PARDISO) ...")
        return pypardiso.spsolve(A.tocsr(), b)
    except ImportError:
        pass
    if A.shape[0] > 30000 and xyz is not None:
        try:
            import pyamg
            say("Solving (AMG-preconditioned CG) ...")
            ml = pyamg.smoothed_aggregation_solver(A.tocsr(), B=_rigid_modes(xyz),
                                                   strength=("symmetric", {"theta": 0.0}),
                                                   max_coarse=500)
            res = []
            x = ml.solve(b, tol=1e-10, accel="cg", maxiter=500, residuals=res)
            if res and res[-1] <= 1e-8 * max(res[0], 1e-300):
                return x
            say("AMG did not converge - falling back to the direct solver ...")
        except ImportError:
            pass
    say("Solving (SuperLU direct) ...")
    return spla.spsolve(A.tocsc(), b, permc_spec="MMD_AT_PLUS_A")


# --------------------------------------------------------------------------- #
# rotations
# --------------------------------------------------------------------------- #
def rotation_at(s0, s, r, ut, band, min_nodes=12):
    """Least-squares rigid rotation of the cross-section at axial position s0.

    Fits  u_t = r * (theta + theta' (s - s0))  over all nodes with
    |s - s0| <= band (band widened automatically if it holds too few nodes).
    """
    ok = r > 1e-12 * max(r.max(), 1e-30)
    for _ in range(12):
        sel = ok & (np.abs(s - s0) <= band)
        if sel.sum() >= min_nodes:
            break
        band *= 1.5
    rr, ds, uu = r[sel], s[sel] - s0, ut[sel]
    A = np.column_stack([rr, rr * ds])
    coef, *_ = np.linalg.lstsq(A, uu, rcond=None)
    return float(coef[0])


def run_torsion(geo: ShaftGeometry, st: FemSettings, progress=None) -> TorsionResult:
    say = progress or (lambda m: None)
    t0 = time.time()
    with gmsh_session() as gm:
        nodes, conn, dN, w, h, diag = _mesh(gm, geo, st, say)
    say(f"Mesh: {len(nodes)} nodes, {len(conn)} tetrahedra")
    axis = geo.axis
    s = axis.s_of(nodes)
    rad = axis.radial(nodes)
    r = np.linalg.norm(rad, axis=1)
    tang = np.cross(axis.direction, rad)
    with np.errstate(invalid="ignore", divide="ignore"):
        tang = np.where(r[:, None] > 0, tang / np.maximum(r, 1e-300)[:, None], 0.0)

    band = st.end_band or max(1e-6 * geo.length, 0.002 * geo.length)
    smin, smax = s.min(), s.max()
    start_nodes = np.where(s <= smin + band)[0]
    end_nodes = np.where(s >= smax - band)[0]
    if st.fixed_end == "start":
        fix_nodes, load_nodes = start_nodes, end_nodes
    else:
        fix_nodes, load_nodes = end_nodes, start_nodes
    if len(fix_nodes) < 3 or len(load_nodes) < 3:
        raise RuntimeError("Too few nodes on the shaft end faces - increase the end band.")

    K, vol, bad = assemble_stiffness(nodes, conn, dN, w, st.E, st.nu, say=say)
    say(f"Stiffness matrix: {K.shape[0]} DOF, {K.nnz} non-zeros")

    # prescribed displacements: clamp = 0, loaded end = rigid rotation by phi0
    phi0 = 1e-3
    ndof = 3 * len(nodes)
    up = np.zeros(ndof)
    presc = np.zeros(ndof, dtype=bool)
    for n_ in fix_nodes:
        presc[3 * n_:3 * n_ + 3] = True
    rot_disp = phi0 * r[load_nodes, None] * tang[load_nodes]        # small-angle rigid rotation
    for k, n_ in enumerate(load_nodes):
        presc[3 * n_:3 * n_ + 3] = True
        up[3 * n_:3 * n_ + 3] = rot_disp[k]
    free = ~presc
    rhs = -K[free][:, presc] @ up[presc]
    Kff = K[free][:, free]
    uf = _solve(Kff, rhs, nodes[free[0::3]], say)   # whole nodes are prescribed
    u = up.copy()
    u[free] = uf
    # reaction torque on the loaded face
    f = K @ u
    fl = f.reshape(-1, 3)[load_nodes]
    T_ref = float(np.sum(np.cross(rad[load_nodes], fl) @ axis.direction))
    if abs(T_ref) < 1e-300:
        raise RuntimeError("Zero reaction torque - check mesh / boundary conditions.")
    scale = st.torque / T_ref
    u *= scale
    disp = u.reshape(-1, 3)
    ut = np.einsum("ij,ij->i", disp, tang)

    # ---- rotations at the lumped-inertia positions and at the cut planes
    slice_band = max(0.6 * h, 1e-4 * geo.length)
    s_mid = geo.mid_positions(st.mid_from_com)
    theta_mid = np.array([rotation_at(x, s, r, ut, slice_band) for x in s_mid])
    s_cut = np.array([geo.s_min] + list(geo.cuts) + [geo.s_max])
    theta_cut = np.empty(len(s_cut))
    for i, x in enumerate(s_cut):
        if i == 0 or i == len(s_cut) - 1:
            if (i == 0) == (st.fixed_end == "start"):
                theta_cut[i] = 0.0
            else:
                theta_cut[i] = phi0 * scale
            continue
        theta_cut[i] = rotation_at(x, s, r, ut, slice_band)

    # nodal rotation field (display); nodes close to the axis carry no rotation
    # information (u_t / r is noise there) -> take the axial interpolation
    order_ = np.argsort(s_cut)
    theta_interp = np.interp(s, s_cut[order_], theta_cut[order_])
    with np.errstate(invalid="ignore", divide="ignore"):
        node_theta = np.where(r > 0.05 * r.max(), ut / np.maximum(r, 1e-300), theta_interp)

    T = abs(st.torque)
    with np.errstate(divide="ignore"):
        k_spring = T / np.abs(np.diff(theta_mid))
        k_section = T / np.abs(np.diff(theta_cut))
    if st.fixed_end == "start":
        gi = 0
        k_ground = T / abs(theta_mid[0])
    else:
        gi = len(theta_mid) - 1
        k_ground = T / abs(theta_mid[-1])
    k_total = T / abs(phi0 * scale)

    res = TorsionResult(settings=st, nodes=nodes, tets=conn[:, :4].copy(), disp=disp,
                        node_theta=node_theta, s_mid=s_mid, theta_mid=theta_mid,
                        s_cut=s_cut, theta_cut=theta_cut, k_spring=k_spring,
                        k_section=k_section, k_ground=float(k_ground),
                        k_total=float(k_total), ground_index=gi, n_dof=ndof,
                        solve_time=time.time() - t0,
                        info={"mesh_size": h, "elements": len(conn), "nodes": len(nodes),
                              "mesh_volume": vol, "inverted_elements": bad,
                              "slice_band": slice_band,
                              "fixed_nodes": len(fix_nodes), "loaded_nodes": len(load_nodes)})
    say(f"FEM finished in {res.solve_time:.1f} s - total stiffness {k_total:.4g} N m/rad")
    return res
