"""FEM modal analysis of the 3-D shaft - check of the lumped torsional model.

* Same tetrahedral mesh as the static torsion FEM (step 2), consistent mass
  matrix, shift-invert Lanczos (``scipy.sparse.linalg.eigsh``) with a PARDISO
  or SuperLU factorisation.
* Boundary conditions: free-free, or the static-FEM clamp face fixed.
* Each FE mode is classified by where its kinetic energy goes: rotation about
  the shaft axis (torsion), along the axis (axial) or across it (bending).
* The section rotations of each FE mode (least-squares rotation at the lumped
  inertia positions) are compared with the lumped mode shapes (MAC), and each
  elastic lumped mode is paired with the FE mode that matches it best.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .fem import FemSettings, _mesh, _rigid_modes, assemble_stiffness, rotation_at
from .geometry import ShaftGeometry
from .gmsh_session import gmsh_session


@dataclass
class ModalResult:
    freq: np.ndarray              # (m,) FE natural frequencies [Hz], rigid-body modes included
    energy: np.ndarray            # (m,4) kinetic energy share: torsion, axial, bending, section deformation
    kind: list                    # (m,) "rigid", "torsion", "bending", "axial", "local" or "coupled"
    theta_mid: np.ndarray         # (n_sections, m) rotation at the lumped inertia positions
    nodes: np.ndarray             # (N,3) [m]
    tets: np.ndarray              # (Ne,4) corner nodes, for display
    shapes: np.ndarray            # (N,3,m) float32 mode shapes, max |u| = 1 per mode
    bc: str = "free"              # "free" or "clamped"
    n_dof: int = 0
    solve_time: float = 0.0
    info: dict = field(default_factory=dict)


def _mass_quadrature(gm, order):
    """Shape functions and gradients at points that integrate N_i N_j exactly."""
    etype = gm.model.mesh.getElementType("Tetrahedron", int(order))
    lc, w = gm.model.mesh.getIntegrationPoints(etype, "Gauss4" if order == 2 else "Gauss2")
    _, f, _ = gm.model.mesh.getBasisFunctions(etype, lc, "Lagrange")
    _, df, _ = gm.model.mesh.getBasisFunctions(etype, lc, "GradLagrange")
    nn = 10 if order == 2 else 4
    return np.asarray(f).reshape(len(w), nn), np.asarray(df).reshape(len(w), nn, 3), np.asarray(w)


def assemble_mass(nodes, conn, N, dN, w, rho, chunk=20000):
    """Consistent mass matrix (3 DOF per node, x/y/z interleaved)."""
    nn = conn.shape[1]
    rows, cols, vals = [], [], []
    for c0 in range(0, len(conn), chunk):
        cc = conn[c0:c0 + chunk]
        X = nodes[cc]
        Me = np.zeros((len(cc), nn, nn))
        for g in range(len(w)):
            det = np.abs(np.linalg.det(np.einsum("ni,mnj->mij", dN[g], X)))
            Me += (rho * w[g] * det)[:, None, None] * np.outer(N[g], N[g])[None]
        rows.append(np.repeat(cc, nn, axis=1).ravel())
        cols.append(np.tile(cc, (1, nn)).ravel())
        vals.append(Me.ravel())
    n = len(nodes)
    Ms = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                       shape=(n, n)).tocsr()
    Ms.sum_duplicates()
    return sp.kron(Ms, sp.eye(3), format="csr"), Ms


def _inverse(A, say):
    """LinearOperator x -> A^-1 x (A symmetric positive definite, sparse)."""
    try:
        from pypardiso import PyPardisoSolver
        say("Factorising (PARDISO) ...")
        A = A.tocsr()
        solver = PyPardisoSolver()
        solver.factorize(A)
        solve = lambda b: solver.solve(A, b)                     # noqa: E731
    except ImportError:
        say("Factorising (SuperLU) ...")
        lu = spla.splu(A.tocsc(), permc_spec="MMD_AT_PLUS_A")
        solve = lu.solve
    return spla.LinearOperator(A.shape, matvec=solve, dtype=float)


def _energy_split(nodes, s, axis, U, Ms, h):
    """Kinetic energy share of every mode: torsion, axial, bending, cross-section deformation.

    The shaft is cut into thin axial slices; in each slice the displacement is
    fitted by a rigid-body motion (3 translations, 3 rotations).  Rotation about
    the shaft axis is torsion, translation along it axial motion, the other four
    are bending (lateral translation and tilting); the rest deforms the section.
    """
    a = axis.direction
    e1 = axis.perpendicular()
    e2 = np.cross(a, e1)
    R = np.array([e1, e2, a])                                  # rows: local frame
    nb = int(np.clip((s.max() - s.min()) / max(h, 1e-12), 4, 200))
    edges = np.linspace(s.min(), s.max(), nb + 1)
    b = np.clip(np.searchsorted(edges, s, side="right") - 1, 0, nb - 1)
    nm = U.shape[2]
    parts = np.zeros((4,) + U.shape)                           # torsion, axial, bending, rest
    for k in range(nb):
        idx = np.where(b == k)[0]
        if len(idx) < 3:
            parts[3][idx] = U[idx]
            continue
        x = (nodes[idx] - nodes[idx].mean(0)) @ R.T             # local coordinates
        n = len(idx)
        B = np.zeros((n, 3, 6))                                 # local displacement per rigid DOF
        B[:, 0, 0] = B[:, 1, 1] = B[:, 2, 2] = 1.0              # translations e1, e2, a
        for c, w in enumerate(np.eye(3)):                       # rotations about e1, e2, a
            B[:, :, 3 + c] = np.cross(w, x)
        uloc = np.einsum("ij,njm->nim", R, U[idx])              # (n,3,m) in the local frame
        coef, *_ = np.linalg.lstsq(B.reshape(3 * n, 6), uloc.reshape(3 * n, nm), rcond=None)
        groups = {0: [5], 1: [2], 2: [0, 1, 3, 4]}
        fit = np.zeros_like(uloc)
        for gi, cols in groups.items():
            comp = np.einsum("nic,cm->nim", B[:, :, cols], coef[cols])
            parts[gi][idx] = np.einsum("ji,njm->nim", R, comp)  # back to global
            fit += comp
        parts[3][idx] = np.einsum("ji,njm->nim", R, uloc - fit)

    def q(F):                                                   # (N,3,m) -> (m,) kinetic energy
        return np.einsum("nim,nim->m", F, (Ms @ F.reshape(len(nodes), -1)).reshape(F.shape))
    e = np.column_stack([q(p) for p in parts])
    return np.maximum(e, 0) / np.maximum(e.sum(1, keepdims=True), 1e-300)


def run_modal(geo: ShaftGeometry, st: FemSettings, n_modes=20, bc="free", progress=None) -> ModalResult:
    """Lowest ``n_modes`` elastic FE modes (plus the rigid-body modes if bc == "free")."""
    say = progress or (lambda m: None)
    t0 = time.time()
    with gmsh_session() as gm:
        nodes, conn, dN, w, h, _ = _mesh(gm, geo, st, say)
        Nm, dNm, wm = _mass_quadrature(gm, st.order)
    say(f"Mesh: {len(nodes)} nodes, {len(conn)} tetrahedra")
    K, vol, _ = assemble_stiffness(nodes, conn, dN, w, st.E, st.nu, say=say)
    say("Assembling mass matrix ...")
    M, Ms = assemble_mass(nodes, conn, Nm, dNm, wm, geo.density)

    axis = geo.axis
    s = axis.s_of(nodes)
    rad = axis.radial(nodes)
    r = np.linalg.norm(rad, axis=1)
    tang = np.cross(axis.direction, rad) / np.maximum(r, 1e-300)[:, None]
    ndof = 3 * len(nodes)
    free = np.ones(ndof, dtype=bool)
    n_rigid = 6
    if bc == "clamped":
        band = st.end_band or max(1e-6 * geo.length, 0.002 * geo.length)
        fix = np.where(s <= s.min() + band)[0] if st.fixed_end == "start" else np.where(s >= s.max() - band)[0]
        for d in range(3):
            free[3 * fix + d] = False
        n_rigid = 0
    Kf, Mf = K[free][:, free], M[free][:, free]

    # shift below zero: K - sigma M is positive definite even with rigid-body modes
    sigma = -(2 * np.pi * 1.0) ** 2
    Ainv = _inverse((Kf - sigma * Mf).tocsr(), say)
    if n_rigid:
        # the 6 rigid-body modes are known exactly; keep them out of the Lanczos
        # space (they would dominate the shifted operator by ~1e6 and spoil it)
        Rb = _rigid_modes(nodes)
        Rb = Rb @ np.linalg.inv(np.linalg.cholesky(Rb.T @ (Mf @ Rb)).T)     # M-orthonormal
        MR = Mf @ Rb
        project = lambda y: y - Rb @ (MR.T @ y)                             # noqa: E731
        OPinv = spla.LinearOperator(Kf.shape, matvec=lambda x: project(Ainv @ x), dtype=float)
        v0 = project(np.random.default_rng(0).standard_normal(Kf.shape[0]))
    else:
        OPinv, v0 = Ainv, None
    nev = min(n_modes, Kf.shape[0] - n_rigid - 2)
    say(f"Lanczos: {nev} modes of {Kf.shape[0]} DOF ...")
    lam, vec = spla.eigsh(Kf, k=nev, M=Mf, sigma=sigma, which="LM", OPinv=OPinv, v0=v0, tol=1e-10)
    order = np.argsort(lam)
    lam, vec = lam[order], vec[:, order]
    if n_rigid:
        lam, vec = np.concatenate([np.zeros(n_rigid), lam]), np.hstack([Rb, vec])
        nev += n_rigid
    freq = np.sqrt(np.maximum(lam, 0.0)) / (2 * np.pi)

    u_all = np.zeros((ndof, nev))
    u_all[free] = vec
    U = u_all.reshape(len(nodes), 3, nev)

    say("Classifying modes ...")
    e = _energy_split(nodes, s, axis, U, Ms, h)
    fmax = max(freq.max(), 1e-30)
    kind = []
    for j in range(nev):
        if bc == "free" and freq[j] < 1e-3 * fmax:
            kind.append("rigid")
        elif e[j].max() < 0.6:
            kind.append("coupled")
        else:                                   # local = the cross-sections deform (disc / web modes)
            kind.append(["torsion", "axial", "bending", "local"][int(np.argmax(e[j]))])
    ut = np.einsum("nim,ni->nm", U, tang)

    slice_band = max(0.6 * h, 1e-4 * geo.length)
    s_mid = geo.mid_positions(st.mid_from_com)
    theta = np.array([[rotation_at(x, s, r, ut[:, j], slice_band) for j in range(nev)] for x in s_mid])
    shapes = (U / np.maximum(np.linalg.norm(U, axis=1).max(0), 1e-300)).astype(np.float32)
    res = ModalResult(freq=freq, energy=e, kind=kind, theta_mid=theta, nodes=nodes, tets=conn[:, :4].copy(),
                      shapes=shapes, bc=bc, n_dof=int(free.sum()), solve_time=time.time() - t0,
                      info={"mesh_size": h, "elements": len(conn), "nodes": len(nodes), "mesh_volume": vol})
    say(f"FEM modal analysis finished in {res.solve_time:.1f} s")
    return res


def mac(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float((a @ b) ** 2 / max((a @ a) * (b @ b), 1e-300))


def match_modes(modal: ModalResult, f_lumped, phi_lumped):
    """Pair the elastic lumped modes one-to-one with FE modes.

    Score of a pair: MAC of the section rotations x torsion share of the FE
    mode; the pairing maximises the total score (each FE mode used once).
    Returns rows (lumped #, f lumped, FE #, f FE, deviation [%], MAC, FE torsion share).
    """
    from scipy.optimize import linear_sum_assignment
    fmax = max(np.max(f_lumped), 1e-30)
    lum = [i for i in range(len(f_lumped)) if f_lumped[i] >= 1e-3 * fmax]     # skip rigid-body modes
    fe = [j for j in range(len(modal.freq)) if modal.kind[j] != "rigid"]
    if not lum or not fe:
        return []
    score = np.array([[mac(phi_lumped[:, i], modal.theta_mid[:, j]) * modal.energy[j, 0] for j in fe]
                      for i in lum])
    ri, ci = linear_sum_assignment(-score)
    rows = []
    for a, b in sorted(zip(ri, ci)):
        i, j = lum[a], fe[b]
        rows.append((i + 1, float(f_lumped[i]), j + 1, float(modal.freq[j]),
                     100.0 * (f_lumped[i] / modal.freq[j] - 1.0),
                     mac(phi_lumped[:, i], modal.theta_mid[:, j]), float(modal.energy[j, 0])))
    return rows


KIND_COLORS = {"torsion": "#4C78A8", "bending": "#BBBBBB", "axial": "#F58518", "local": "#54A24B",
               "coupled": "#B279A2"}


def plot_comparison(fig, geo, modal: ModalResult, f_lumped, phi_lumped, rows, n_shapes=5):
    """Frequency chart (all FE modes by type, lumped modes, pairs) and paired mode shapes."""
    fig.clear()
    ax = fig.add_subplot(211)
    for kind, c in KIND_COLORS.items():
        fx = [f for f, k in zip(modal.freq, modal.kind) if k == kind]
        if fx:
            ax.vlines(fx, 0, 1, color=c, lw=2, label=f"FE {kind}")
    fe = [f for f in f_lumped if f >= 1e-3 * max(f_lumped)]
    ax.vlines(fe, 1.1, 2.1, color="#E45756", lw=2, label="lumped model")
    for i, fl, j, ff, d, m, _ in rows:
        ax.plot([ff, fl], [1.0, 1.1], color="black", lw=0.6)
        ax.text(fl, 2.15, f"{d:+.1f}%", ha="center", fontsize=8)
    ax.set_yticks([0.5, 1.6], ["FE", "lumped"])
    ax.set_ylim(0, 3.0)
    ax.set_xlabel("frequency [Hz]")
    ax.set_title(f"Natural frequencies: 3-D FEM ({modal.bc}, all mode types) and lumped torsional model")
    ax2 = fig.add_subplot(212)
    x = (geo.mid_positions() - geo.s_min) * 1e3
    import matplotlib.pyplot as plt
    for n, (i, fl, j, ff, d, m, _) in enumerate(rows[:n_shapes]):
        c = plt.cm.tab10(n)
        pl, pf = phi_lumped[:, i - 1], modal.theta_mid[:, j - 1]
        pl = pl / pl[np.argmax(np.abs(pl))]
        pf = pf * (pf @ pl) / max(pf @ pf, 1e-300)                 # best fit scale
        ax2.plot(x, pl, "o-", color=c, label=f"lumped {i}: {fl:.0f} Hz")
        ax2.plot(x, pf, "s--", color=c, mfc="none", label=f"FE {j}: {ff:.0f} Hz (MAC {m:.2f})")
    ax2.axhline(0, color="black", lw=0.8)
    ax2.set_xlabel("position along shaft [mm]")
    ax2.set_ylabel("section rotation [-]")
    ax2.set_title("Torsional mode shapes at the lumped inertia positions (solid: lumped, dashed: FE)")
    ax2.grid(alpha=0.3)
    return ax, ax2
