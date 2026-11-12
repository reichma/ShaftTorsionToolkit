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
