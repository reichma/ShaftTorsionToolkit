"""Checks against closed-form results:  python -m pytest tests -q"""
import os
import sys

import gmsh
import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from shaft_torsion import geometry as g                       # noqa: E402
from shaft_torsion.dynamics import Load, TorsionalModel       # noqa: E402
from shaft_torsion.fem import FemSettings, run_torsion       # noqa: E402

D, L, RHO, E, NU = 0.04, 0.4, 7850.0, 210e9, 0.3
G = E / (2 * (1 + NU))
JP = np.pi * D ** 4 / 32


@pytest.fixture(scope="module")
def tilted_cylinder(tmp_path_factory):
    """Cylinder d40 x 400 mm along an oblique axis (tests the axis handling)."""
    path = str(tmp_path_factory.mktemp("cad") / "tilted.step")
    if not gmsh.isInitialized():
        gmsh.initialize(interruptible=False)
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("tilted")
    d = np.array([1.0, 2.0, 2.0]) / 3.0
    gmsh.model.occ.addCylinder(10, -5, 3, *(400 * d), 20)
    gmsh.model.occ.synchronize()
    gmsh.write(path)
    gmsh.model.remove()
    return path, d


def test_axis_and_inertia(tilted_cylinder):
    path, d = tilted_cylinder
    info = g.load_step(path)
    assert abs(abs(info.auto_axis.direction @ d) - 1) < 1e-6
    s0, s1 = g.axial_range(info.surface[0], info.auto_axis)
    assert abs((s1 - s0) - L) < 1e-6
    geo = g.split_model(path, info.auto_axis, [s0 + 0.1, s0 + 0.25], RHO, (s0, s1))
    for sec in geo.sections:
        assert sec.J == pytest.approx(RHO * JP * sec.length, rel=1e-6)
    assert geo.total_J == pytest.approx(RHO * JP * L, rel=1e-6)


def test_torsion_stiffness(tilted_cylinder):
    path, _ = tilted_cylinder
    info = g.load_step(path)
    s0, s1 = g.axial_range(info.surface[0], info.auto_axis)
    geo = g.split_model(path, info.auto_axis, [s0 + 0.1, s0 + 0.25], RHO, (s0, s1))
    res = run_torsion(geo, FemSettings(E=E, nu=NU, torque=500.0, mesh_size=0.008, order=2))
    mids = np.array([0.05, 0.175, 0.325])
    k_exact = G * JP / np.diff(mids)
    assert np.allclose(res.k_spring, k_exact, rtol=2e-3)
    assert res.k_total == pytest.approx(G * JP / L, rel=2e-3)
    assert res.k_ground == pytest.approx(G * JP / 0.05, rel=5e-3)
    assert np.allclose(res.k_section, G * JP / np.array([0.1, 0.15, 0.15]), rtol=2e-3)


def test_two_inertia_frequency_and_step():
    m = TorsionalModel([0.1, 0.3], [1e4])
    f, phi, _ = m.modes()
    assert f[1] == pytest.approx(np.sqrt(1e4 * 0.4 / 0.03) / (2 * np.pi), rel=1e-9)
    single = TorsionalModel([0.1], [], ground_node=0, k_ground=1e4)
    out = single.simulate([Load(0, "step", 10.0)], 0.05, 1e-5)
    exact = 10 / 1e4 * (1 - np.cos(np.sqrt(1e5) * out["t"]))
    assert np.abs(out["theta"][0] - exact).max() < 1e-7
    single.c_ground = 20.0
    _, zeta = single.complex_modes()
    assert zeta[0] == pytest.approx(20 / (2 * np.sqrt(1e4 * 0.1)), rel=1e-6)


def test_loads():
    t = np.array([0.0, 0.25, 0.5, 1.0])
    assert np.allclose(Load(kind="ramp", amplitude=4, t_start=0, t_end=0.5)(t), [0, 2, 4, 4])
    assert np.allclose(Load(kind="pulse", amplitude=1, t_start=0.2, t_end=0.6)(t), [0, 1, 1, 0])
    assert np.allclose(Load(kind="expression", expression="2*t")(t), 2 * t)


def test_spring_to_ground_per_inertia():
    from shaft_torsion.project import model_from_dict
    # two inertias, each held by its own mount: closed-form 2-DOF frequencies
    J1, J2, k, k1, k2 = 0.1, 0.3, 1e4, 2e3, 5e3
    m = TorsionalModel([J1, J2], [k], k_abs=[k1, k2])
    K = np.array([[k + k1, -k], [-k, k + k2]])
    w2 = np.linalg.eigvals(np.linalg.solve(np.diag([J1, J2]), K)).real
    w2.sort()
    f, _, _ = m.modes()
    assert np.allclose(f, np.sqrt(w2) / (2 * np.pi), rtol=1e-9)
    # a mount on the middle inertia: no rigid-body mode any more
    mid = TorsionalModel([0.1, 0.2, 0.1], [1e4, 1e4], k_abs=[0, 3e3, 0])
    assert mid.modes()[0][0] > 1.0
    # support torque balances a static load
    mid.c_abs = np.array([0, 50.0, 0])
    out = mid.simulate([Load(0, "constant", 10.0)], 3.0, 1e-3)
    assert out["support_torque"][1, -1] == pytest.approx(10.0, rel=1e-3)
    # saved and loaded with the project
    assert np.allclose(model_from_dict(m.to_dict()).k_abs, [k1, k2])
    assert np.allclose(model_from_dict({"J": [1.0], "k": []}).k_abs, [0.0])   # old project files


def test_load_step_after_second_order_fem(tilted_cylinder):
    """A FEM run leaves Mesh.ElementOrder = 2 in gmsh; the next import must still
    give a display surface (regression: surface was None)."""
    from shaft_torsion.gmsh_session import gmsh_session
    with gmsh_session() as gm:
        gm.option.setNumber("Mesh.ElementOrder", 2)
    info = g.load_step(tilted_cylinder[0])
    assert info.surface is not None and len(info.surface[1]) > 0


def test_fem_modal_torsion_frequencies(tilted_cylinder):
    """3-D FE modes of the cylinder: torsion frequencies n c / 2L (free) and c / 4L (clamped),
    bending and axial modes recognised, lumped model paired with the right FE mode."""
    from shaft_torsion.modal import match_modes, run_modal
    path, _ = tilted_cylinder
    info = g.load_step(path)
    s0, s1 = g.axial_range(info.surface[0], info.auto_axis)
    geo = g.split_model(path, info.auto_axis, [s0 + 0.1 * i for i in range(1, 4)], RHO, (s0, s1))
    st = FemSettings(E=E, nu=NU, mesh_size=0.008, order=2)
    c = np.sqrt(G / RHO)
    free = run_modal(geo, st, n_modes=10, bc="free")
    assert free.kind[:6] == ["rigid"] * 6
    tors = free.freq[[k == "torsion" for k in free.kind]]
    assert tors[0] == pytest.approx(c / (2 * L), rel=2e-3)
    assert "bending" in free.kind and "axial" in free.kind
    clamped = run_modal(geo, st, n_modes=6, bc="clamped")
    tors = clamped.freq[[k == "torsion" for k in clamped.kind]]
    assert tors[0] == pytest.approx(c / (4 * L), rel=2e-3)
    # lumped model of 4 equal sections, free-free: its 1st elastic mode is the FE torsion mode
    res = run_torsion(geo, st)
    lumped = TorsionalModel([s.J for s in geo.sections], res.k_spring)
    f, phi, _ = lumped.modes()
    rows = match_modes(free, f, phi)
    i, fl, j, ff, dev, mac_, et = rows[0]
    assert free.kind[j - 1] == "torsion" and mac_ > 0.99
    assert abs(dev) < 5.0                     # lumped (4 rigid discs) vs. continuous shaft
