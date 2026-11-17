"""Check the lumped torsional model against a 3-D FEM modal analysis.

    python examples/fem_modal_check.py examples/crankshaft_4cyl_project.json [output_folder] [--modes 40]

Splits the shaft as in the project, takes J, k from the project's vibration
model (or from a fresh static FEM run if the project has none), and computes
the lowest FE modes of the full 3-D solid with the same boundary condition
(free-free, or clamped at the FEM clamp face if the model is grounded there).

Prints every FE mode with its type (torsion / bending / axial / coupled) and
pairs each elastic lumped mode with the FE mode whose section rotations match
best (MAC).  Writes fem_modal_check.png and fem_modes.csv into the output
folder (default: <project>_results).
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np                                                  # noqa: E402

from shaft_torsion import geometry as g                             # noqa: E402
from shaft_torsion.dynamics import TorsionalModel                   # noqa: E402
from shaft_torsion.fem import run_torsion                           # noqa: E402
from shaft_torsion.modal import match_modes, plot_comparison, run_modal             # noqa: E402
from shaft_torsion.project import fem_settings, load_project, make_axis, model_from_results  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("project")
    ap.add_argument("out", nargs="?")
    ap.add_argument("--modes", type=int, default=40, help="number of elastic FE modes")
    a = ap.parse_args(argv)
    proj = load_project(a.project)
    out = os.path.abspath(a.out or os.path.splitext(a.project)[0] + "_results")
    os.makedirs(out, exist_ok=True)

    info = g.load_step(proj["step_file"])
    ax = proj["axis"]
    axis = make_axis(info, ax.get("mode"), ax.get("origin_mm"), ax.get("direction"))
    s_range = g.axial_range(info.surface[0], axis)
    geo = g.split_model(proj["step_file"], axis, [s_range[0] + c / 1000.0 for c in proj["cuts_mm"]],
                        proj["density"], s_range)
    st = fem_settings(proj["fem"])
    md = proj.get("model")
    if md:
        lumped = TorsionalModel(md["J"], md["k"], ground_node=md.get("ground_node", -1),
                                k_ground=md.get("k_ground", 0.0))
        if any(md.get("k_abs") or []):
            print("note: springs to ground (k to ground) are not in the FE model - left out of the comparison")
    else:
        lumped = model_from_results(geo, run_torsion(geo, st), ground=False)
    bc = "clamped" if lumped.ground_node >= 0 else "free"
    f_l, phi_l, _ = lumped.modes()

    modal = run_modal(geo, st, n_modes=a.modes, bc=bc,
                      progress=lambda t: None if t.startswith("Assembling stiff") else print(t))
    print(f"\nFE modes ({bc}, {modal.n_dof} DOF, {modal.info['elements']} elements):")
    print(f"{'#':>3} {'f [Hz]':>9}  {'type':8} {'torsion':>8} {'axial':>6} {'bending':>8} {'section':>8}")
    for j, (f, k, e) in enumerate(zip(modal.freq, modal.kind, modal.energy)):
        print(f"{j + 1:3d} {f:9.1f}  {k:8} {e[0]:8.3f} {e[1]:6.3f} {e[2]:8.3f} {e[3]:8.3f}")
    rows = match_modes(modal, f_l, phi_l)
    print("\nLumped model vs. FE (pairs by MAC of the section rotations):")
    print(f"{'lumped':>6} {'f [Hz]':>9} {'FE #':>5} {'f FE [Hz]':>10} {'diff':>7} {'MAC':>6} {'torsion':>8}")
    for i, fl, j, ff, d, m, et in rows:
        print(f"{i:6d} {fl:9.1f} {j:5d} {ff:10.1f} {d:+6.1f}% {m:6.3f} {et:8.3f}")
    if modal.freq.max() < max(f_l):
        print(f"\nnote: the {a.modes} FE modes reach {modal.freq.max():.0f} Hz only - "
              f"higher lumped modes are paired with lower FE modes; use --modes to get more")

    with open(os.path.join(out, "fem_modes.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["FE mode", "f [Hz]", "type", "torsion", "axial", "bending", "section"])
        for j, (f, k, e) in enumerate(zip(modal.freq, modal.kind, modal.energy)):
            w.writerow([j + 1, f"{f:.6g}", k] + [f"{x:.4f}" for x in e])
    plot(os.path.join(out, "fem_modal_check.png"), geo, modal, f_l, phi_l, rows)
    print("written to", out)
    return 0


def plot(path, geo, modal, f_l, phi_l, rows):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.figure import Figure
    fig = Figure(figsize=(12, 8))
    ax, ax2 = plot_comparison(fig, geo, modal, f_l, phi_l, rows)
    ax.legend(fontsize=8, ncol=6, loc="upper center")
    ax2.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=130)


if __name__ == "__main__":
    sys.exit(main())
