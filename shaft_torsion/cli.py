"""Batch run without GUI:

    python -m shaft_torsion.cli project.json [output_folder]

Runs split -> FEM -> vibration model -> time simulation (if loads are
defined) and writes CSV tables and PNG plots into the output folder.
"""
from __future__ import annotations

import os
import sys

import numpy as np


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        print(__doc__)
        return 1
    from .project import load_project, run_pipeline
    proj = load_project(argv[0])
    out = argv[1] if len(argv) > 1 else os.path.splitext(argv[0])[0] + "_results"
    os.makedirs(out, exist_ok=True)
    geo, res, model, sim, files = run_pipeline(proj, out)
    print("\nSections:")
    print(f"{'#':>3} {'start mm':>10} {'end mm':>10} {'mass kg':>10} {'J kg m2':>12} {'theta deg':>10} {'k->next':>12}")
    for i, s in enumerate(geo.sections):
        k = f"{res.k_spring[i]:12.5g}" if i < len(res.k_spring) else ""
        print(f"{i + 1:3d} {(s.s_start - geo.s_min) * 1e3:10.2f} {(s.s_end - geo.s_min) * 1e3:10.2f} "
              f"{s.mass:10.4g} {s.J:12.5g} {np.degrees(res.theta_mid[i]):10.5f} {k}")
    print(f"total stiffness {res.k_total:.5g} N m/rad, clamp->inertia {res.ground_index + 1}: {res.k_ground:.5g}")
    f, _, zeta = model.modes()
    print("natural frequencies [Hz]:", ", ".join(f"{x:.4g}" for x in f))
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        if sim is not None:
            fig, ax = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
            ax[0].plot(sim["t"], np.degrees(sim["theta"].T))
            ax[0].set_ylabel("rotation [deg]")
            ax[1].plot(sim["t"], sim["spring_torque"].T)
            ax[1].set_ylabel("spring torque [N m]")
            ax[1].set_xlabel("time [s]")
            for a in ax:
                a.grid(alpha=0.3)
            fig.savefig(os.path.join(out, "time_response.png"), dpi=130)
            files.append(os.path.join(out, "time_response.png"))
    except ImportError:
        pass
    print("written:", *files, sep="\n  ")
    return 0


if __name__ == "__main__":
    sys.exit(main())
