"""Run a project through the GUI and save the 3-D results.

The GUI window is placed outside the visible screen area while it renders
(Qt's "offscreen" platform gives VTK no render window on Windows).

    python examples/render_3d.py examples/gear_shaft_mounted_sweep.json [output_folder] [--modes 4]

Does what the buttons in the GUI do (load, split, FEM, modes, FRF, time
simulation, animation) and writes into the output folder
(default: <project>_results next to the project file):

  split_model.png              sections in colour
  fem_rotation.png             static FEM rotation (step 2)
  mode_<j>.png / mode_<j>.gif  3-D mode shapes (rigid sections, coloured by rotation)
  sim_peak.png                 3-D view at the largest rotation of the time simulation
  sim_peak.gif                 3-D animation around that instant
  modes.png, frf.png, time_response.png   the 2-D plots of tab 3
"""
from __future__ import annotations

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np                                            # noqa: E402
from PySide6.QtWidgets import QApplication, QMessageBox       # noqa: E402

ERRORS = []


def _dialog(kind):
    def show(_parent, title, text, *a, **k):
        print(f"[{kind}] {title}: {text}")
        if kind != "info":
            ERRORS.append(f"{title}: {text}")
        return QMessageBox.Ok
    return show


# dialogs would block without a window - print them instead
QMessageBox.information = staticmethod(_dialog("info"))
QMessageBox.warning = staticmethod(_dialog("warning"))
QMessageBox.critical = staticmethod(_dialog("error"))


def wait(app, win, what):
    t0 = time.time()
    app.processEvents()
    while win._threads:
        app.processEvents()
        time.sleep(0.05)
    app.processEvents()
    if ERRORS:
        raise SystemExit(f"{what} failed: {ERRORS[-1]}")
    print(f"{what}: {time.time() - t0:.1f} s")


def write_gif(win, path, frames):
    win.plotter.open_gif(path, fps=20)
    for i in frames:
        win._goto_frame(i)
        win.plotter.write_frame()
    win.plotter.mwriter.close()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("project")
    ap.add_argument("out", nargs="?")
    ap.add_argument("--modes", type=int, default=4, help="number of mode shapes to render")
    a = ap.parse_args(argv)
    project = os.path.abspath(a.project)
    out = os.path.abspath(a.out or os.path.splitext(project)[0] + "_results")
    os.makedirs(out, exist_ok=True)

    app = QApplication.instance() or QApplication([])
    from shaft_torsion.gui import MainWindow
    win = MainWindow()
    win.resize(1600, 1000)
    win.show()
    win.move(-4000, -4000)             # keep the window out of sight

    win.load_project(project)
    wait(app, win, "load STEP")
    win.run_split()
    wait(app, win, "split")
    win.view_tabs.setCurrentIndex(0)
    win.plotter.screenshot(os.path.join(out, "split_model.png"))
    win.run_fem()                       # also builds the vibration model from the project values
    wait(app, win, "FEM")
    win.view_tabs.setCurrentIndex(0)
    win.plotter.screenshot(os.path.join(out, "fem_rotation.png"))

    model = win.current_model()
    f, _, zeta = model.modes()
    print("natural frequencies [Hz]:", ", ".join(f"{x:.5g}" for x in f))
    print("damping ratios       [-]:", ", ".join(f"{x:.3g}" for x in zeta))
    win.run_modes()
    win.fig.savefig(os.path.join(out, "modes.png"), dpi=130)
    win.run_frf()
    win.fig.savefig(os.path.join(out, "frf.png"), dpi=130)

    win.cmb_disp.setCurrentIndex(0)    # rigid sections
    win.cmb_color.setCurrentIndex(0)   # coloured by rotation angle
    win.chk_auto.setChecked(True)
    win.cmb_src.setCurrentIndex(1)     # mode shapes
    for j in range(1, min(a.modes, len(f)) + 1):
        win.sp_mode.setValue(j)
        win.prepare_animation()
        win.view_tabs.setCurrentIndex(0)
        win._goto_frame(15)            # phase 90 deg = largest deflection
        win.plotter.screenshot(os.path.join(out, f"mode_{j}.png"))
        write_gif(win, os.path.join(out, f"mode_{j}.gif"), range(0, 60, 2))
        print(f"mode {j}: {f[j - 1]:.5g} Hz rendered")

    if win.loads():
        win.run_sim()
        wait(app, win, "time simulation")
        win.plot_sim()
        win.fig.savefig(os.path.join(out, "time_response.png"), dpi=130)
        win.cmb_src.setCurrentIndex(0)  # time simulation
        win.prepare_animation()
        win.view_tabs.setCurrentIndex(0)
        th = win.anim_frames
        ipk = int(np.argmax(np.abs(th).max(1)))
        win._goto_frame(ipk)
        win.plotter.screenshot(os.path.join(out, "sim_peak.png"))
        write_gif(win, os.path.join(out, "sim_peak.gif"), range(max(ipk - 40, 0), min(ipk + 40, len(th))))
        print(f"time simulation: largest rotation at {win.anim_labels[ipk]}")

    win.close()
    print("written to", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
