"""PySide6 / PyVista GUI for the shaft torsion toolkit.

Start with  ``python run_gui.py``  (or ``python -m shaft_torsion.gui``).
"""
from __future__ import annotations

import json
import os
import sys
import traceback

import numpy as np
from PySide6.QtCore import QObject, QThread, Qt, QTimer, Signal
from PySide6.QtGui import QAction, QDoubleValidator
from PySide6.QtWidgets import (QAbstractItemView, QApplication, QCheckBox, QComboBox, QDoubleSpinBox,
                               QFileDialog, QFormLayout, QGridLayout, QGroupBox, QHBoxLayout,
                               QHeaderView, QLabel, QLineEdit, QMainWindow, QMessageBox,
                               QPlainTextEdit, QPushButton, QScrollArea, QSlider, QSpinBox,
                               QSplitter, QTableWidget, QTableWidgetItem, QTabWidget, QVBoxLayout,
                               QWidget)

os.environ.setdefault("QT_API", "pyside6")
import pyvista as pv                                                    # noqa: E402
from pyvistaqt import QtInteractor                                      # noqa: E402
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT  # noqa: E402
from matplotlib.figure import Figure                                    # noqa: E402

from . import geometry as g                                             # noqa: E402
from . import project as prj                                            # noqa: E402
from .dynamics import LOAD_KINDS, Load, TorsionalModel                  # noqa: E402
from .fem import FemSettings, auto_mesh_size, run_torsion               # noqa: E402

SECTION_COLORS = ["#4C78A8", "#F58518", "#54A24B", "#E45756", "#72B7B2", "#EECA3B",
                  "#B279A2", "#FF9DA6", "#9D755D", "#BAB0AC"]


# --------------------------------------------------------------------------- #
# small widgets / helpers
# --------------------------------------------------------------------------- #
class SciEdit(QLineEdit):
    """Line edit for floats in any notation (1e6, 0.25 ...)."""

    def __init__(self, value=0.0, width=90):
        super().__init__()
        v = QDoubleValidator()
        v.setNotation(QDoubleValidator.ScientificNotation)
        self.setValidator(v)
        self.setMaximumWidth(width)
        self.setValue(value)

    def value(self):
        try:
            return float(self.text().replace(",", "."))
        except ValueError:
            return 0.0

    def setValue(self, v):
        self.setText(f"{v:.6g}")


def dspin(val, lo, hi, dec=3, step=1.0, suffix=""):
    w = QDoubleSpinBox()
    w.setRange(lo, hi)
    w.setDecimals(dec)
    w.setSingleStep(step)
    w.setValue(val)
    if suffix:
        w.setSuffix(" " + suffix)
    return w


def make_table(headers, editable=True):
    t = QTableWidget(0, len(headers))
    t.setHorizontalHeaderLabels(headers)
    t.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
    t.verticalHeader().setVisible(False)
    if not editable:
        t.setEditTriggers(QAbstractItemView.NoEditTriggers)
    t.setMinimumHeight(120)
    return t


def set_row(table, r, values, editable_cols=()):
    for c, v in enumerate(values):
        if isinstance(v, float):
            txt = f"{v:.6g}"
        else:
            txt = str(v)
        it = QTableWidgetItem(txt)
        if c not in editable_cols:
            it.setFlags(it.flags() & ~Qt.ItemIsEditable)
        table.setItem(r, c, it)


def cell_float(table, r, c, default=0.0):
    it = table.item(r, c)
    if it is None:
        return default
    try:
        return float(it.text().replace(",", "."))
    except ValueError:
        return default


class Worker(QObject):
    finished = Signal(object)
    failed = Signal(str)
    progress = Signal(str)

    def __init__(self, fn, *args, **kw):
        super().__init__()
        self.fn, self.args, self.kw = fn, args, kw

    def run(self):
        try:
            res = self.fn(*self.args, progress=self.progress.emit, **self.kw)
            self.finished.emit(res)
        except Exception as exc:                    # noqa: BLE001
            self.failed.emit(f"{exc}\n\n{traceback.format_exc()}")


# --------------------------------------------------------------------------- #
# main window
# --------------------------------------------------------------------------- #
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Shaft Torsion Toolkit - inertia, FEM stiffness, torsional vibration")
        self.resize(1500, 950)
        self.info = None          # StepInfo
        self.axis = None
        self.s_range = None
        self.geo = None           # ShaftGeometry
        self.res = None           # TorsionResult
        self.model = None         # TorsionalModel
        self.sim = None
        self._threads = []
        self.anim_frames = None   # (nframes, nsec) rotation [rad]
        self.anim_times = None
        self.anim_index = 0
        self.timer = QTimer(self)
        self.timer.timeout.connect(self._advance_frame)

        # ---- layout: left controls | right 3-D + plots
        split = QSplitter(Qt.Horizontal)
        self.setCentralWidget(split)
        left = QWidget()
        lv = QVBoxLayout(left)
        self.tabs = QTabWidget()
        lv.addWidget(self.tabs, 1)
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumHeight(130)
        lv.addWidget(self.log)
        split.addWidget(left)

        self.view_tabs = QTabWidget()
        self.plotter = QtInteractor(self)
        self.plotter.set_background("white")
        self.view_tabs.addTab(self.plotter.interactor, "3D model")
        pw = QWidget()
        pl = QVBoxLayout(pw)
        self.fig = Figure(figsize=(8, 6), tight_layout=True)
        self.canvas = FigureCanvasQTAgg(self.fig)
        pl.addWidget(NavigationToolbar2QT(self.canvas, pw))
        pl.addWidget(self.canvas)
        self.view_tabs.addTab(pw, "Plots")
        split.addWidget(self.view_tabs)
        split.setSizes([560, 940])

        self._build_tab_geometry()
        self._build_tab_fem()
        self._build_tab_dynamics()
        self._build_tab_animation()
        self._build_menu()
        self.statusBar().showMessage("Open a STEP file to start  (File > Open STEP)")
        self.plotter.add_text("Open a STEP file  (File > Open STEP)", font_size=12, color="gray")

    # ------------------------------------------------------------------ #
    def _scroll(self, w, title):
        sa = QScrollArea()
        sa.setWidgetResizable(True)
        sa.setWidget(w)
        self.tabs.addTab(sa, title)

    def _build_menu(self):
        m = self.menuBar().addMenu("&File")
        for text, fn, sc in [("Open STEP ...", self.open_step, "Ctrl+O"),
                             ("Load project ...", self.load_project, "Ctrl+L"),
                             ("Save project ...", self.save_project, "Ctrl+S"),
                             ("Export results (CSV) ...", self.export_results, "Ctrl+E"),
                             ("Save screenshot ...", self.save_screenshot, None),
                             ("Quit", self.close, "Ctrl+Q")]:
            a = QAction(text, self)
            if sc:
                a.setShortcut(sc)
            a.triggered.connect(fn)
            m.addAction(a)
        v = self.menuBar().addMenu("&View")
        for text, fn in [("Reset camera", lambda: self.plotter.reset_camera()),
                         ("View along axis", self.view_along_axis),
                         ("Isometric", lambda: self.plotter.view_isometric())]:
            a = QAction(text, self)
            a.triggered.connect(fn)
            v.addAction(a)

    # ------------------------------------------------------------------ #
    # tab 1 - geometry
    # ------------------------------------------------------------------ #
    def _build_tab_geometry(self):
        w = QWidget()
        v = QVBoxLayout(w)
        b = QPushButton("Open STEP file ...")
        b.clicked.connect(self.open_step)
        v.addWidget(b)
        self.lbl_file = QLabel("no file")
        self.lbl_file.setWordWrap(True)
        v.addWidget(self.lbl_file)

        gb = QGroupBox("Shaft axis")
        f = QGridLayout(gb)
        self.cmb_axis = QComboBox()
        self.cmb_axis.addItems(["auto (principal axis)", "X through COM", "Y through COM",
                                "Z through COM", "custom"])
        f.addWidget(QLabel("Axis"), 0, 0)
        f.addWidget(self.cmb_axis, 0, 1, 1, 3)
        self.ed_org = [SciEdit(0, 80) for _ in range(3)]
        self.ed_dir = [SciEdit(x, 80) for x in (1, 0, 0)]
        f.addWidget(QLabel("Origin [mm]"), 1, 0)
        f.addWidget(QLabel("Direction"), 2, 0)
        for i in range(3):
            f.addWidget(self.ed_org[i], 1, i + 1)
            f.addWidget(self.ed_dir[i], 2, i + 1)
        b = QPushButton("Apply axis")
        b.clicked.connect(self.apply_axis)
        f.addWidget(b, 3, 0, 1, 2)
        self.lbl_len = QLabel("")
        f.addWidget(self.lbl_len, 3, 2, 1, 2)
        self.cmb_axis.currentIndexChanged.connect(self.apply_axis)
        v.addWidget(gb)

        gb = QGroupBox("Material")
        f = QFormLayout(gb)
        self.sp_rho = dspin(7850, 1, 50000, 1, 10, "kg/m³")
        f.addRow("Density", self.sp_rho)
        v.addWidget(gb)

        gb = QGroupBox("Cut planes (distance from shaft start)")
        gv = QVBoxLayout(gb)
        self.tbl_cuts = make_table(["position [mm]"])
        self.tbl_cuts.itemChanged.connect(lambda *_: self.show_preview())
        gv.addWidget(self.tbl_cuts)
        h = QHBoxLayout()
        for text, fn in [("Add", self.add_cut), ("Remove", self.remove_cut), ("Clear", self.clear_cuts),
                         ("Sort", self.sort_cuts)]:
            bb = QPushButton(text)
            bb.clicked.connect(fn)
            h.addWidget(bb)
        gv.addLayout(h)
        h = QHBoxLayout()
        self.sp_nsec = QSpinBox()
        self.sp_nsec.setRange(2, 200)
        self.sp_nsec.setValue(6)
        bb = QPushButton("Equal sections")
        bb.clicked.connect(self.equal_cuts)
        h.addWidget(QLabel("No. of sections"))
        h.addWidget(self.sp_nsec)
        h.addWidget(bb)
        gv.addLayout(h)
        self.chk_pick = QCheckBox("Pick cut positions with the mouse (right click on the model)")
        self.chk_pick.toggled.connect(self.toggle_pick)
        gv.addWidget(self.chk_pick)
        v.addWidget(gb)

        b = QPushButton("Split model && compute inertias")
        b.setStyleSheet("font-weight: bold; padding: 6px")
        b.clicked.connect(self.run_split)
        v.addWidget(b)
        self.tbl_sec = make_table(["#", "start [mm]", "end [mm]", "mass [kg]", "J [kg·m²]"], False)
        v.addWidget(self.tbl_sec)
        self.lbl_sec = QLabel("")
        v.addWidget(self.lbl_sec)
        v.addStretch()
        self._scroll(w, "1 Geometry")

    # ------------------------------------------------------------------ #
    # tab 2 - FEM
    # ------------------------------------------------------------------ #
    def _build_tab_fem(self):
        w = QWidget()
        v = QVBoxLayout(w)
        gb = QGroupBox("Static torsion analysis")
        f = QFormLayout(gb)
        self.sp_E = dspin(210, 0.001, 2000, 3, 1, "GPa")
        self.sp_nu = dspin(0.3, 0.0, 0.499, 3, 0.01)
        self.ed_T = SciEdit(1000.0)
        self.cmb_fix = QComboBox()
        self.cmb_fix.addItems(["shaft start (load at end)", "shaft end (load at start)"])
        self.sp_h = dspin(0, 0, 1e4, 2, 0.5, "mm")
        self.sp_h.setSpecialValueText("auto")
        self.cmb_order = QComboBox()
        self.cmb_order.addItems(["quadratic tetrahedra (recommended)", "linear tetrahedra (fast, too stiff)"])
        self.sp_band = dspin(0, 0, 1e4, 3, 0.1, "mm")
        self.sp_band.setSpecialValueText("auto")
        self.cmb_lump = QComboBox()
        self.cmb_lump.addItems(["middle of section", "centre of mass of section"])
        f.addRow("Young's modulus E", self.sp_E)
        f.addRow("Poisson ratio ν", self.sp_nu)
        f.addRow("Torque [N·m]", self.ed_T)
        f.addRow("Fixed (clamped) end", self.cmb_fix)
        f.addRow("Element size", self.sp_h)
        f.addRow("Element type", self.cmb_order)
        f.addRow("End-face band", self.sp_band)
        f.addRow("Lumped inertia at", self.cmb_lump)
        v.addWidget(gb)
        b = QPushButton("Mesh && run torsion FEM")
        b.setStyleSheet("font-weight: bold; padding: 6px")
        b.clicked.connect(self.run_fem)
        v.addWidget(b)
        self.tbl_fem = make_table(["#", "pos [mm]", "θ [deg]", "k section [N·m/rad]",
                                   "k → next [N·m/rad]"], False)
        v.addWidget(self.tbl_fem)
        self.lbl_fem = QLabel("")
        self.lbl_fem.setWordWrap(True)
        v.addWidget(self.lbl_fem)
        h = QHBoxLayout()
        b = QPushButton("Show FEM result")
        b.clicked.connect(self.show_fem)
        h.addWidget(b)
        h.addWidget(QLabel("deformation scale"))
        self.ed_warp = SciEdit(0.0, 80)
        self.ed_warp.setToolTip("0 = automatic")
        h.addWidget(self.ed_warp)
        v.addLayout(h)
        v.addStretch()
        self._scroll(w, "2 FEM stiffness")

    # ------------------------------------------------------------------ #
    # tab 3 - torsional vibration model
    # ------------------------------------------------------------------ #
    def _build_tab_dynamics(self):
        w = QWidget()
        v = QVBoxLayout(w)
        b = QPushButton("Build model from steps 1 + 2")
        b.clicked.connect(self.build_model)
        v.addWidget(b)
        v.addWidget(QLabel("Inertias (editable - e.g. add attached parts):"))
        self.tbl_J = make_table(["#", "J [kg·m²]", "damping to ground [N·m·s/rad]"])
        v.addWidget(self.tbl_J)
        v.addWidget(QLabel("Springs / dampers between inertias (editable):"))
        self.tbl_k = make_table(["between", "k [N·m/rad]", "c [N·m·s/rad]"])
        v.addWidget(self.tbl_k)

        gb = QGroupBox("Boundary condition")
        f = QFormLayout(gb)
        self.cmb_bc = QComboBox()
        self.cmb_bc.addItems(["free - free (shaft floating)", "grounded at the FEM clamp"])
        self.ed_kg = SciEdit(0.0)
        self.ed_cg = SciEdit(0.0)
        f.addRow("Type", self.cmb_bc)
        f.addRow("k ground [N·m/rad]", self.ed_kg)
        f.addRow("c ground [N·m·s/rad]", self.ed_cg)
        v.addWidget(gb)

        gb = QGroupBox("Damping helpers")
        f = QGridLayout(gb)
        self.ed_beta = SciEdit(1e-5, 80)
        b = QPushButton("c = β·k")
        b.clicked.connect(self.damping_beta)
        f.addWidget(QLabel("β [s]"), 0, 0)
        f.addWidget(self.ed_beta, 0, 1)
        f.addWidget(b, 0, 2)
        self.ed_zeta = SciEdit(0.02, 80)
        b = QPushButton("c for ζ at 1st elastic mode")
        b.clicked.connect(self.damping_zeta)
        f.addWidget(QLabel("ζ [-]"), 1, 0)
        f.addWidget(self.ed_zeta, 1, 1)
        f.addWidget(b, 1, 2)
        v.addWidget(gb)

        b = QPushButton("Natural frequencies && mode shapes")
        b.clicked.connect(self.run_modes)
        v.addWidget(b)
        self.tbl_modes = make_table(["mode", "f [Hz]", "ζ [-]"], False)
        v.addWidget(self.tbl_modes)

        gb = QGroupBox("Loads (torques on inertias)")
        gv = QVBoxLayout(gb)
        self.tbl_loads = make_table(["inertia #", "type", "amplitude [N·m]", "f [Hz]", "f end [Hz]",
                                     "phase [°]", "offset [N·m]", "t start [s]", "t end [s]",
                                     "expression f(t)"])
        self.tbl_loads.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        self.tbl_loads.setMinimumHeight(140)
        gv.addWidget(self.tbl_loads)
        h = QHBoxLayout()
        for text, fn in [("Add load", self.add_load), ("Remove load", self.remove_load)]:
            bb = QPushButton(text)
            bb.clicked.connect(fn)
            h.addWidget(bb)
        gv.addLayout(h)
        v.addWidget(gb)

        gb = QGroupBox("Time simulation (Newmark)")
        f = QGridLayout(gb)
        self.ed_tend = SciEdit(0.5, 80)
        self.ed_dt = SciEdit(0.0, 80)
        self.ed_dt.setToolTip("0 = automatic (20 steps per period of the highest mode)")
        self.cmb_plot = QComboBox()
        self.cmb_plot.addItems(["rotation angles", "twist between inertias", "spring torques",
                                "angular velocities", "applied loads"])
        f.addWidget(QLabel("t end [s]"), 0, 0)
        f.addWidget(self.ed_tend, 0, 1)
        f.addWidget(QLabel("Δt [s] (0 = auto)"), 0, 2)
        f.addWidget(self.ed_dt, 0, 3)
        b = QPushButton("Run simulation")
        b.setStyleSheet("font-weight: bold")
        b.clicked.connect(self.run_sim)
        f.addWidget(b, 1, 0, 1, 2)
        f.addWidget(self.cmb_plot, 1, 2, 1, 2)
        self.cmb_plot.currentIndexChanged.connect(lambda *_: self.plot_sim())
        v.addWidget(gb)

        gb = QGroupBox("Frequency response")
        f = QGridLayout(gb)
        self.sp_frf_node = QSpinBox()
        self.sp_frf_node.setRange(1, 999)
        self.ed_fmin = SciEdit(1.0, 70)
        self.ed_fmax = SciEdit(2000.0, 70)
        f.addWidget(QLabel("excite inertia #"), 0, 0)
        f.addWidget(self.sp_frf_node, 0, 1)
        f.addWidget(QLabel("f min/max [Hz]"), 1, 0)
        f.addWidget(self.ed_fmin, 1, 1)
        f.addWidget(self.ed_fmax, 1, 2)
        b = QPushButton("Plot FRF")
        b.clicked.connect(self.run_frf)
        f.addWidget(b, 0, 2)
        v.addWidget(gb)
        v.addStretch()
        self._scroll(w, "3 Vibration model")

    # ------------------------------------------------------------------ #
    # tab 4 - animation
    # ------------------------------------------------------------------ #
    def _build_tab_animation(self):
        w = QWidget()
        v = QVBoxLayout(w)
        gb = QGroupBox("What to animate")
        f = QFormLayout(gb)
        self.cmb_src = QComboBox()
        self.cmb_src.addItems(["time simulation", "mode shape", "static FEM rotation (step 2)"])
        self.sp_mode = QSpinBox()
        self.sp_mode.setRange(1, 999)
        self.cmb_disp = QComboBox()
        self.cmb_disp.addItems(["rigid sections (lumped model)", "smooth twist (interpolated)"])
        self.cmb_color = QComboBox()
        self.cmb_color.addItems(["rotation angle", "section colours"])
        self.chk_auto = QCheckBox("automatic (max. 20° shown)")
        self.chk_auto.setChecked(True)
        self.ed_amp = SciEdit(100.0, 90)
        self.chk_ptr = QCheckBox("show rotation pointers")
        self.chk_ptr.setChecked(True)
        f.addRow("Source", self.cmb_src)
        f.addRow("Mode #", self.sp_mode)
        f.addRow("Display", self.cmb_disp)
        f.addRow("Colour", self.cmb_color)
        f.addRow("Amplification", self.chk_auto)
        f.addRow("  factor", self.ed_amp)
        f.addRow("", self.chk_ptr)
        v.addWidget(gb)
        b = QPushButton("Prepare animation")
        b.setStyleSheet("font-weight: bold; padding: 6px")
        b.clicked.connect(self.prepare_animation)
        v.addWidget(b)
        h = QHBoxLayout()
        self.btn_play = QPushButton("▶ Play")
        self.btn_play.clicked.connect(self.toggle_play)
        h.addWidget(self.btn_play)
        h.addWidget(QLabel("speed"))
        self.sp_speed = QSpinBox()
        self.sp_speed.setRange(1, 100)
        self.sp_speed.setValue(1)
        h.addWidget(self.sp_speed)
        self.chk_loop = QCheckBox("loop")
        self.chk_loop.setChecked(True)
        h.addWidget(self.chk_loop)
        v.addLayout(h)
        self.slider = QSlider(Qt.Horizontal)
        self.slider.valueChanged.connect(self._goto_frame)
        v.addWidget(self.slider)
        self.lbl_frame = QLabel("")
        v.addWidget(self.lbl_frame)
        b = QPushButton("Export animation (GIF / MP4) ...")
        b.clicked.connect(self.export_animation)
        v.addWidget(b)
        v.addStretch()
        self._scroll(w, "4 Animation")

    # ================================================================== #
    # background jobs
    # ================================================================== #
    def msg(self, text):
        self.log.appendPlainText(text)
        self.statusBar().showMessage(text.splitlines()[0] if text else "")

    def run_bg(self, fn, on_done, *args, **kw):
        """Run ``fn`` in a worker thread; ``on_done(result)`` runs in the GUI thread."""
        if self._threads:
            QMessageBox.information(self, "Busy", "Another calculation is still running.")
            return
        th = QThread()
        wk = Worker(fn, *args, **kw)
        wk.moveToThread(th)
        th.started.connect(wk.run)
        # bound methods of this (GUI-thread) object -> queued connections
        wk.progress.connect(self.msg)
        wk.finished.connect(self._job_done)
        wk.failed.connect(self._job_failed)
        self._job_cb = on_done
        self._threads.append((th, wk))
        QApplication.setOverrideCursor(Qt.WaitCursor)
        th.start()

    def _job_end(self):
        th, wk = self._threads.pop()
        th.quit()
        th.wait()
        QApplication.restoreOverrideCursor()

    def _job_done(self, res):
        self._job_end()
        try:
            self._job_cb(res)
        except Exception as exc:                    # noqa: BLE001
            self._job_failed(f"{exc}\n\n{traceback.format_exc()}", ended=True)

    def _job_failed(self, text, ended=False):
        if not ended:
            self._job_end()
        self.msg("ERROR: " + text.splitlines()[0])
        QMessageBox.critical(self, "Error", text[:3000])

    # ================================================================== #
    # step 1
    # ================================================================== #
    def open_step(self, path=None, after=None):
        if not path:
            path, _ = QFileDialog.getOpenFileName(self, "Open STEP", "",
                                                  "CAD (*.step *.stp *.STEP *.STP *.brep);;All (*)")
        if not path:
            return
        self.msg(f"Loading {path} ...")

        def done(info):
            self.info = info
            self.geo = self.res = self.model = self.sim = None
            self.lbl_file.setText(f"<b>{os.path.basename(path)}</b><br>{info.n_volumes} solid(s), "
                                  f"volume {info.volume * 1e9:.4g} mm³, "
                                  f"mass {info.volume * self.sp_rho.value():.4g} kg")
            self.apply_axis()
            self.msg(f"Loaded {os.path.basename(path)}")
            if after:
                after()

        self.run_bg(lambda progress: g.load_step(path), done)

    def apply_axis(self, *_):
        if self.info is None:
            return
        i = self.cmb_axis.currentIndex()
        mode = ["auto", "x", "y", "z", "custom"][i]
        try:
            self.axis = prj.make_axis(self.info, mode, [e.value() for e in self.ed_org],
                                      [e.value() for e in self.ed_dir])
        except Exception as exc:                    # noqa: BLE001
            QMessageBox.warning(self, "Axis", str(exc))
            return
        if mode != "custom":
            for k in range(3):
                self.ed_org[k].setValue(round(self.axis.origin[k] * 1e3, 6))
                self.ed_dir[k].setValue(round(self.axis.direction[k], 9))
        self.s_range = g.axial_range(self.info.surface[0], self.axis)
        self.lbl_len.setText(f"shaft length {(self.s_range[1] - self.s_range[0]) * 1e3:.2f} mm")
        self.geo = None
        self.show_preview()

    def cuts_mm(self):
        out = []
        for r in range(self.tbl_cuts.rowCount()):
            it = self.tbl_cuts.item(r, 0)
            if it and it.text().strip():
                try:
                    out.append(float(it.text().replace(",", ".")))
                except ValueError:
                    pass
        return out

    def set_cuts(self, values):
        self.tbl_cuts.blockSignals(True)
        self.tbl_cuts.setRowCount(0)
        for x in values:
            r = self.tbl_cuts.rowCount()
            self.tbl_cuts.insertRow(r)
            self.tbl_cuts.setItem(r, 0, QTableWidgetItem(f"{x:.3f}"))
        self.tbl_cuts.blockSignals(False)
        self.show_preview()

    def add_cut(self):
        L = (self.s_range[1] - self.s_range[0]) * 1e3 if self.s_range else 100
        cuts = self.cuts_mm()
        self.set_cuts(cuts + [L / 2 if not cuts else min(max(cuts) + L / 10, L * 0.99)])

    def remove_cut(self):
        rows = sorted({i.row() for i in self.tbl_cuts.selectedIndexes()}, reverse=True)
        if not rows and self.tbl_cuts.rowCount():
            rows = [self.tbl_cuts.rowCount() - 1]
        for r in rows:
            self.tbl_cuts.removeRow(r)
        self.show_preview()

    def clear_cuts(self):
        self.set_cuts([])

    def sort_cuts(self):
        self.set_cuts(sorted(self.cuts_mm()))

    def equal_cuts(self):
        if not self.s_range:
            return
        L = (self.s_range[1] - self.s_range[0]) * 1e3
        n = self.sp_nsec.value()
        self.set_cuts([L * i / n for i in range(1, n)])

    def toggle_pick(self, on):
        try:
            if on:
                self.plotter.enable_point_picking(callback=self._picked, show_message=False,
                                                  left_clicking=False, show_point=False)
            else:
                self.plotter.disable_picking()
        except Exception as exc:                    # noqa: BLE001
            self.msg(f"Picking not available: {exc}")

    def _picked(self, point, *_):
        if point is None or self.axis is None:
            return
        s = float(self.axis.s_of(np.asarray(point)[None] / 1e3)[0])
        self.set_cuts(sorted(self.cuts_mm() + [round((s - self.s_range[0]) * 1e3, 2)]))

    def run_split(self):
        if self.info is None:
            QMessageBox.information(self, "Split", "Open a STEP file first.")
            return
        cuts = sorted(self.cuts_mm())
        cuts_abs = [self.s_range[0] + c / 1e3 for c in cuts]
        rho = self.sp_rho.value()

        def done(geo):
            self.geo = geo
            self.res = self.model = self.sim = None
            self.tbl_sec.setRowCount(len(geo.sections))
            for r, s in enumerate(geo.sections):
                set_row(self.tbl_sec, r, [s.index + 1, (s.s_start - geo.s_min) * 1e3,
                                          (s.s_end - geo.s_min) * 1e3, s.mass, s.J])
            self.lbl_sec.setText(f"total mass {geo.total_mass:.5g} kg,  total J {geo.total_J:.5g} kg·m²")
            self._build_display()
            self.show_split()
            self.plot_sections()
            self.sp_h.setValue(0)
            self.msg(f"Split into {len(geo.sections)} sections "
                     f"(automatic element size would be {auto_mesh_size(geo) * 1e3:.2f} mm)")

        self.run_bg(g.split_model, done, self.info.path, self.axis, cuts_abs, rho, self.s_range)

    # ================================================================== #
    # step 2
    # ================================================================== #
    def fem_settings(self):
        return FemSettings(E=self.sp_E.value() * 1e9, nu=self.sp_nu.value(), torque=self.ed_T.value() or 1.0,
                           mesh_size=self.sp_h.value() / 1e3, order=2 if self.cmb_order.currentIndex() == 0 else 1,
                           fixed_end="start" if self.cmb_fix.currentIndex() == 0 else "end",
                           end_band=self.sp_band.value() / 1e3, mid_from_com=self.cmb_lump.currentIndex() == 1)

    def run_fem(self):
        if self.geo is None:
            QMessageBox.information(self, "FEM", "Split the model first (tab 1).")
            return

        def done(res):
            self.res = res
            geo = self.geo
            n = len(geo.sections)
            self.tbl_fem.setRowCount(n)
            for i in range(n):
                set_row(self.tbl_fem, i, [i + 1, (res.s_mid[i] - geo.s_min) * 1e3,
                                          np.degrees(res.theta_mid[i]), res.k_section[i],
                                          res.k_spring[i] if i < n - 1 else "-"])
            inf = res.info
            self.lbl_fem.setText(
                f"Mesh: {inf['nodes']} nodes, {inf['elements']} elements, size {inf['mesh_size'] * 1e3:.2f} mm, "
                f"{res.n_dof} DOF ({res.solve_time:.1f} s).  Mesh volume / CAD volume = "
                f"{inf['mesh_volume'] / sum(s.volume for s in geo.sections):.4f}<br>"
                f"<b>Total stiffness clamp → loaded end: {res.k_total:.5g} N·m/rad</b><br>"
                f"Clamp → inertia {res.ground_index + 1}: {res.k_ground:.5g} N·m/rad")
            self.show_fem()
            self.plot_fem()
            self.build_model()

        self.run_bg(run_torsion, done, self.geo, self.fem_settings())

    # ================================================================== #
    # step 3
    # ================================================================== #
    def build_model(self):
        if self.geo is None:
            QMessageBox.information(self, "Model", "Split the model first (tab 1).")
            return
        n = len(self.geo.sections)
        self.tbl_J.setRowCount(n)
        for i, s in enumerate(self.geo.sections):
            set_row(self.tbl_J, i, [i + 1, s.J, 0.0], editable_cols=(1, 2))
        self.tbl_k.setRowCount(n - 1)
        for i in range(n - 1):
            k = self.res.k_spring[i] if self.res is not None else 0.0
            set_row(self.tbl_k, i, [f"{i + 1} – {i + 2}", k, 0.0], editable_cols=(1, 2))
        if self.res is not None:
            self.ed_kg.setValue(self.res.k_ground)
            self.cmb_bc.setItemText(1, f"grounded at the FEM clamp (inertia {self.res.ground_index + 1})")
        self.sp_frf_node.setMaximum(n)
        pm = getattr(self, "_pending_model", None)
        if pm and len(pm.get("J", [])) == n:          # values saved in the project file
            for i in range(n):
                self.tbl_J.item(i, 1).setText(f"{pm['J'][i]:.6g}")
                self.tbl_J.item(i, 2).setText(f"{pm['c_abs'][i]:.6g}")
            for i in range(n - 1):
                if pm["k"][i] > 0:
                    self.tbl_k.item(i, 1).setText(f"{pm['k'][i]:.6g}")
                self.tbl_k.item(i, 2).setText(f"{pm['c'][i]:.6g}")
            self.cmb_bc.setCurrentIndex(1 if pm.get("ground_node", -1) >= 0 else 0)
            if pm.get("k_ground"):
                self.ed_kg.setValue(pm["k_ground"])
            self.ed_cg.setValue(pm.get("c_ground", 0.0))
            if self.res is not None:
                self._pending_model = None
            self.msg("Inertia / stiffness / damping values taken from the project file")
            return
        self.msg("Vibration model built from sections and FEM stiffness" if self.res is not None
                 else "Vibration model built - run the FEM to fill in the stiffness values")

    def current_model(self):
        n = self.tbl_J.rowCount()
        if n == 0:
            raise ValueError("Build the vibration model first.")
        J = [cell_float(self.tbl_J, i, 1) for i in range(n)]
        ca = [cell_float(self.tbl_J, i, 2) for i in range(n)]
        k = [cell_float(self.tbl_k, i, 1) for i in range(n - 1)]
        c = [cell_float(self.tbl_k, i, 2) for i in range(n - 1)]
        if any(x <= 0 for x in k):
            raise ValueError("All spring stiffnesses must be > 0 (run the FEM, tab 2).")
        gnode = -1
        if self.cmb_bc.currentIndex() == 1:
            gnode = self.res.ground_index if self.res is not None else 0
        self.model = TorsionalModel(J, k, c, ca, ground_node=gnode, k_ground=self.ed_kg.value(),
                                    c_ground=self.ed_cg.value())
        return self.model

    def damping_beta(self):
        beta = self.ed_beta.value()
        for i in range(self.tbl_k.rowCount()):
            self.tbl_k.item(i, 2).setText(f"{beta * cell_float(self.tbl_k, i, 1):.6g}")
        if self.cmb_bc.currentIndex() == 1:
            self.ed_cg.setValue(beta * self.ed_kg.value())

    def damping_zeta(self):
        try:
            m = self.current_model()
        except ValueError as exc:
            QMessageBox.warning(self, "Damping", str(exc))
            return
        f, _, _ = m.modes()
        fe = f[f > 1e-3 * f.max()]
        if not len(fe):
            return
        self.ed_beta.setValue(2 * self.ed_zeta.value() / (2 * np.pi * fe[0]))
        self.damping_beta()
        self.msg(f"Stiffness-proportional damping β = {self.ed_beta.value():.4g} s "
                 f"→ ζ = {self.ed_zeta.value()} at {fe[0]:.4g} Hz")

    def run_modes(self):
        try:
            m = self.current_model()
        except ValueError as exc:
            QMessageBox.warning(self, "Modes", str(exc))
            return
        f, phi, zeta = m.modes()
        self.tbl_modes.setRowCount(len(f))
        for j in range(len(f)):
            set_row(self.tbl_modes, j, [j + 1, f[j], zeta[j]])
        self.sp_mode.setMaximum(len(f))
        self.plot_modes(f, phi)

    def add_load(self, ld=None):
        ld = ld if isinstance(ld, Load) else Load(node=max(self.tbl_J.rowCount() - 1, 0), kind="step",
                                                  amplitude=100.0)
        r = self.tbl_loads.rowCount()
        self.tbl_loads.insertRow(r)
        cmb = QComboBox()
        cmb.addItems(LOAD_KINDS)
        cmb.setCurrentText(ld.kind)
        self.tbl_loads.setCellWidget(r, 1, cmb)
        vals = [ld.node + 1, None, ld.amplitude, ld.frequency, ld.f_end, ld.phase_deg, ld.offset,
                ld.t_start, ld.t_end, ld.expression]
        for c, x in enumerate(vals):
            if c == 1:
                continue
            self.tbl_loads.setItem(r, c, QTableWidgetItem(f"{x:.6g}" if isinstance(x, float) else str(x)))

    def remove_load(self):
        rows = sorted({i.row() for i in self.tbl_loads.selectedIndexes()}, reverse=True)
        if not rows and self.tbl_loads.rowCount():
            rows = [self.tbl_loads.rowCount() - 1]
        for r in rows:
            self.tbl_loads.removeRow(r)

    def loads(self):
        out = []
        for r in range(self.tbl_loads.rowCount()):
            it = self.tbl_loads.item(r, 9)
            out.append(Load(node=int(cell_float(self.tbl_loads, r, 0, 1)) - 1,
                            kind=self.tbl_loads.cellWidget(r, 1).currentText(),
                            amplitude=cell_float(self.tbl_loads, r, 2), frequency=cell_float(self.tbl_loads, r, 3),
                            f_end=cell_float(self.tbl_loads, r, 4), phase_deg=cell_float(self.tbl_loads, r, 5),
                            offset=cell_float(self.tbl_loads, r, 6), t_start=cell_float(self.tbl_loads, r, 7),
                            t_end=cell_float(self.tbl_loads, r, 8), expression=it.text() if it else ""))
        return out

    def run_sim(self):
        try:
            m = self.current_model()
            loads = self.loads()
            if not loads:
                raise ValueError("Define at least one load.")
            dt = self.ed_dt.value() or m.recommended_dt()
            t_end = self.ed_tend.value()
            if t_end / dt > 2e6:
                raise ValueError(f"{t_end / dt:.3g} time steps - increase Δt or shorten the simulation.")
            loads[0](np.array([0.0]))       # validate expressions early
        except Exception as exc:                    # noqa: BLE001
            QMessageBox.warning(self, "Simulation", str(exc))
            return
        self.msg(f"Simulating {t_end} s with Δt = {dt:.3g} s ...")

        def done(sim):
            self.sim = sim
            self.msg(f"Simulation finished ({len(sim['t'])} steps). Max spring torque "
                     f"{np.abs(sim['spring_torque']).max():.4g} N·m")
            self.plot_sim()

        self.run_bg(lambda progress: m.simulate(loads, t_end, dt), done)

    def run_frf(self):
        try:
            m = self.current_model()
        except ValueError as exc:
            QMessageBox.warning(self, "FRF", str(exc))
            return
        fmin, fmax = max(self.ed_fmin.value(), 1e-3), self.ed_fmax.value()
        freqs = np.geomspace(fmin, fmax, 1500)
        H = m.frf(self.sp_frf_node.value() - 1, freqs)
        self.fig.clear()
        ax = self.fig.add_subplot(211)
        for i in range(m.n):
            ax.loglog(freqs, np.abs(H[i]), label=f"θ{i + 1}", color=SECTION_COLORS[i % 10])
        ax.set_ylabel("|θ / M| [rad/(N·m)]")
        ax.set_title(f"Receptance for a unit torque on inertia {self.sp_frf_node.value()}")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=7, ncol=4)
        ax2 = self.fig.add_subplot(212, sharex=ax)
        dH = np.abs(np.diff(H, axis=0) * m.k[:, None])
        for i in range(m.n - 1):
            ax2.loglog(freqs, dH[i], label=f"k{i + 1}-{i + 2}", color=SECTION_COLORS[i % 10])
        ax2.set_xlabel("frequency [Hz]")
        ax2.set_ylabel("|spring torque / M| [-]")
        ax2.grid(True, which="both", alpha=0.3)
        ax2.legend(fontsize=7, ncol=4)
        self.canvas.draw()
        self.view_tabs.setCurrentIndex(1)

    # ================================================================== #
    # 3-D display
    # ================================================================== #
    def _axis_mm(self):
        return g.Axis(self.axis.origin * 1e3, self.axis.direction)

    def show_preview(self):
        """Raw (unsplit) model, axis and the planned cut planes."""
        if self.info is None or self.axis is None:
            return
        self.timer.stop()
        p = self.plotter
        p.clear()
        pts, tris = self.info.surface
        mesh = pv.PolyData(pts * 1e3, np.hstack([np.full((len(tris), 1), 3), tris]).ravel())
        p.add_mesh(mesh, color="#B0B8C4", smooth_shading=True, opacity=0.6 if self.cuts_mm() else 1.0,
                   pickable=True)
        self._draw_axis(p, pts)
        ax = self._axis_mm()
        rmax = np.linalg.norm(g.Axis(self.axis.origin, self.axis.direction).radial(pts), axis=1).max() * 1e3
        for c in self.cuts_mm():
            s = self.s_range[0] * 1e3 + c
            ctr = ax.origin + s * ax.direction
            p.add_mesh(pv.Disc(center=ctr, inner=0, outer=rmax * 1.15, normal=ax.direction, c_res=48),
                       color="#E45756", opacity=0.45, pickable=False)
        p.add_text(f"{len(self.cuts_mm())} cut plane(s) - press 'Split model'", font_size=10, color="black")
        p.reset_camera()
        self.view_tabs.setCurrentIndex(0)

    def _draw_axis(self, p, pts_m):
        ax = self._axis_mm()
        s0, s1 = self.s_range[0] * 1e3, self.s_range[1] * 1e3
        L = s1 - s0
        a = ax.origin + (s0 - 0.08 * L) * ax.direction
        b = ax.origin + (s1 + 0.08 * L) * ax.direction
        p.add_mesh(pv.Line(a, b), color="black", line_width=2, pickable=False)
        p.add_mesh(pv.Arrow(start=b - 0.02 * L * ax.direction, direction=ax.direction, scale=0.06 * L), color="black",
                   pickable=False)

    def _build_display(self):
        """Merged surface mesh of all sections with section index per point."""
        geo = self.geo
        P, F, S, off = [], [], [], 0
        for i in range(len(geo.sections)):
            if i not in geo.surface:
                continue
            pts, tris = geo.surface[i]
            P.append(pts * 1e3)
            F.append(np.hstack([np.full((len(tris), 1), 3), tris + off]))
            S.append(np.full(len(pts), i))
            off += len(pts)
        self.d_p0 = np.vstack(P)
        self.d_sec = np.concatenate(S)
        self.d_s = self._axis_mm().s_of(self.d_p0)
        self.d_mesh = pv.PolyData(self.d_p0.copy(), np.vstack(F).ravel())
        self.d_mesh.point_data["section"] = self.d_sec + 1
        self.d_mesh.point_data["rotation [deg]"] = np.zeros(len(self.d_p0))
        # pointers: radial line per section at the lumped position
        ax = self._axis_mm()
        perp = ax.perpendicular()
        rad = np.linalg.norm(ax.radial(self.d_p0), axis=1)
        lines, pp = [], []
        r = rad.max()
        for i, sec in enumerate(geo.sections):
            sm = sec.s_mid * 1e3
            c = ax.origin + sm * ax.direction
            pp += [c, c + perp * r * 1.35]
            lines.append([2, 2 * i, 2 * i + 1])
        self.ptr_p0 = np.array(pp)
        self.ptr_sec = np.repeat(np.arange(len(geo.sections)), 2)
        self.ptr_mesh = pv.PolyData(self.ptr_p0.copy(), lines=np.array(lines).ravel())

    def show_split(self, color_by="section", clim=None):
        if self.geo is None:
            return
        p = self.plotter
        p.clear()
        if color_by == "section":
            cols = np.array([[int(SECTION_COLORS[i % 10][k:k + 2], 16) for k in (1, 3, 5)]
                             for i in self.d_sec], dtype=np.uint8)
            self.d_mesh.point_data["rgb"] = cols
            p.add_mesh(self.d_mesh, scalars="rgb", rgb=True, smooth_shading=True)
        else:
            p.add_mesh(self.d_mesh, scalars="rotation [deg]", cmap="coolwarm", clim=clim,
                       smooth_shading=True, scalar_bar_args={"title": "rotation [deg] (true)",
                                                             "color": "black"})
        if self.chk_ptr.isChecked():
            p.add_mesh(self.ptr_mesh, color="#222222", line_width=5, render_lines_as_tubes=True)
            ax = self._axis_mm()
            tips = self.ptr_p0[1::2] + ax.perpendicular() * 0.06 * np.linalg.norm(
                self.ptr_p0[1] - self.ptr_p0[0])
            p.add_point_labels(tips, [str(i + 1) for i in range(len(self.geo.sections))],
                               font_size=14, point_size=1, shape_opacity=0.5, always_visible=True,
                               name="labels")
        self._draw_axis(p, None)
        p.reset_camera()
        self.view_tabs.setCurrentIndex(0)

    def view_along_axis(self):
        if self.axis is None:
            return
        ax = self._axis_mm()
        c = np.array(self.plotter.center)
        L = (self.s_range[1] - self.s_range[0]) * 1e3
        self.plotter.camera_position = [tuple(c - ax.direction * 3 * L), tuple(c), tuple(ax.perpendicular())]

    def show_fem(self):
        r = self.res
        if r is None:
            return
        self.timer.stop()
        cells = np.hstack([np.full((len(r.tets), 1), 4), r.tets]).ravel()
        grid = pv.UnstructuredGrid(cells, np.full(len(r.tets), pv.CellType.TETRA), r.nodes * 1e3)
        grid.point_data["rotation [deg]"] = np.degrees(r.node_theta)
        grid.point_data["disp"] = r.disp * 1e3
        try:
            surf = grid.extract_surface(algorithm="dataset_surface")
        except TypeError:                           # older PyVista
            surf = grid.extract_surface()
        scale = self.ed_warp.value()
        if scale <= 0:
            L = self.geo.length * 1e3
            dmax = np.abs(r.disp).max() * 1e3
            scale = 0.05 * L / max(dmax, 1e-30)
        warped = surf.warp_by_vector("disp", factor=scale)
        p = self.plotter
        p.clear()
        p.add_mesh(warped, scalars="rotation [deg]", cmap="jet", show_edges=True, edge_color="#555555",
                   line_width=0.3, scalar_bar_args={"title": "rotation [deg]", "color": "black"})
        p.add_mesh(surf, style="wireframe", color="gray", opacity=0.08)
        p.add_text(f"Static torsion, T = {r.settings.torque:g} N·m, deformation ×{scale:.3g}",
                   font_size=10, color="black")
        p.reset_camera()
        self.view_tabs.setCurrentIndex(0)

    # ================================================================== #
    # plots
    # ================================================================== #
    def plot_sections(self):
        geo = self.geo
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        x0 = geo.s_min
        for s in geo.sections:
            ax.bar((s.s_mid - x0) * 1e3, s.J, width=s.length * 1e3 * 0.95,
                   color=SECTION_COLORS[s.index % 10], edgecolor="black")
            ax.text((s.s_mid - x0) * 1e3, s.J, f"{s.J:.3g}", ha="center", va="bottom", fontsize=8)
        ax.set_xlabel("position along shaft [mm]")
        ax.set_ylabel("J about shaft axis [kg·m²]")
        ax.set_title("Mass moment of inertia per section")
        ax.grid(True, alpha=0.3)
        self.canvas.draw()

    def plot_fem(self):
        r, geo = self.res, self.geo
        self.fig.clear()
        ax = self.fig.add_subplot(211)
        sn = (geo.axis.s_of(r.nodes) - geo.s_min) * 1e3
        sel = np.random.default_rng(0).choice(len(sn), min(len(sn), 20000), replace=False)
        ax.plot(sn[sel], np.degrees(r.node_theta[sel]), ".", ms=1, color="#BBBBBB", label="FE nodes")
        ax.plot((r.s_cut - geo.s_min) * 1e3, np.degrees(r.theta_cut), "s--", color="#E45756", label="cut planes")
        ax.plot((r.s_mid - geo.s_min) * 1e3, np.degrees(r.theta_mid), "o-", color="#4C78A8",
                label="lumped inertia positions")
        ax.set_ylabel("rotation [deg]")
        ax.set_title(f"Static torsion, T = {r.settings.torque:g} N·m")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
        ax2 = self.fig.add_subplot(212)
        xm = (r.s_mid - geo.s_min) * 1e3
        centers = 0.5 * (xm[1:] + xm[:-1])
        ax2.bar(centers, r.k_spring, width=np.diff(xm) * 0.9, color="#54A24B", edgecolor="black")
        for x, k in zip(centers, r.k_spring):
            ax2.text(x, k, f"{k:.3g}", ha="center", va="bottom", fontsize=8)
        ax2.set_yscale("log")
        ax2.set_xlabel("position along shaft [mm]")
        ax2.set_ylabel("k between inertias [N·m/rad]")
        ax2.grid(True, alpha=0.3)
        self.canvas.draw()

    def plot_modes(self, f, phi):
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        n = phi.shape[0]
        x = np.arange(1, n + 1)
        if self.geo is not None and len(self.geo.sections) == n:
            x = (self.geo.mid_positions() - self.geo.s_min) * 1e3
        for j in range(min(len(f), 8)):
            ax.plot(x, phi[:, j], "o-", label=f"mode {j + 1}: {f[j]:.4g} Hz", color=SECTION_COLORS[j % 10])
        ax.axhline(0, color="black", lw=0.8)
        ax.set_xlabel("position along shaft [mm]" if self.geo is not None else "inertia #")
        ax.set_ylabel("relative rotation [-]")
        ax.set_title("Mode shapes")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
        self.canvas.draw()
        self.view_tabs.setCurrentIndex(1)

    def plot_sim(self):
        sim = self.sim
        if sim is None:
            return
        t = sim["t"]
        k = self.cmb_plot.currentIndex()
        data, ylabel, lab = [
            (np.degrees(sim["theta"]), "rotation [deg]", "θ{}"),
            (np.degrees(sim["twist"]), "twist [deg]", "Δθ {}-{}"),
            (sim["spring_torque"], "spring torque [N·m]", "M {}-{}"),
            (sim["omega"] * 30 / np.pi, "speed [rpm]", "n{}"),
            (sim["load"], "applied torque [N·m]", "M{}"),
        ][k]
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        for i in range(data.shape[0]):
            if k == 4 and not np.any(data[i]):
                continue
            ax.plot(t, data[i], lw=1, color=SECTION_COLORS[i % 10], label=lab.format(i + 1, i + 2))
        ax.set_xlabel("time [s]")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7, ncol=4)
        self.canvas.draw()
        self.view_tabs.setCurrentIndex(1)

    # ================================================================== #
    # animation
    # ================================================================== #
    def prepare_animation(self):
        if self.geo is None:
            QMessageBox.information(self, "Animation", "Split the model first (tab 1).")
            return
        src = self.cmb_src.currentIndex()
        n = len(self.geo.sections)
        if src == 0:
            if self.sim is None:
                QMessageBox.information(self, "Animation", "Run a time simulation first (tab 3).")
                return
            th = self.sim["theta"].T
            t = self.sim["t"]
            step = max(1, len(t) // 1500)
            frames, times = th[::step], t[::step]
            labels = [f"t = {x:.4f} s" for x in times]
        elif src == 1:
            try:
                m = self.current_model()
            except ValueError as exc:
                QMessageBox.warning(self, "Animation", str(exc))
                return
            f, phi, _ = m.modes()
            j = min(self.sp_mode.value(), len(f)) - 1
            ph = np.linspace(0, 2 * np.pi, 60, endpoint=False)
            frames = np.outer(np.sin(ph), phi[:, j]) * np.radians(1.0)   # 1 deg reference amplitude
            labels = [f"mode {j + 1}: {f[j]:.4g} Hz, phase {np.degrees(p):.0f}°" for p in ph]
        else:
            if self.res is None:
                QMessageBox.information(self, "Animation", "Run the FEM first (tab 2).")
                return
            ramp = np.concatenate([np.linspace(0, 1, 30), np.linspace(1, 0, 30)])
            frames = np.outer(ramp, self.res.theta_mid)
            labels = [f"static torque × {x:.2f}" for x in ramp]
        if frames.shape[1] != n:
            QMessageBox.warning(self, "Animation", "Model size does not match the number of sections.")
            return
        # remove rigid-body drift (free-free model) only for display of mode shapes? keep true motion.
        vmax = np.abs(frames).max()
        if self.chk_auto.isChecked():
            amp = np.radians(20.0) / max(vmax, 1e-300)
            self.ed_amp.setValue(amp)
        self.anim_amp = self.ed_amp.value()
        self.anim_frames = frames
        self.anim_labels = labels
        self.anim_clim = [-np.degrees(vmax), np.degrees(vmax)]
        self.show_split("section" if self.cmb_color.currentIndex() == 1 else "rotation",
                        clim=self.anim_clim)
        self.anim_text = self.plotter.add_text("", position="upper_left", font_size=10, color="black")
        self.slider.blockSignals(True)
        self.slider.setRange(0, len(frames) - 1)
        self.slider.setValue(0)
        self.slider.blockSignals(False)
        self._goto_frame(0)
        self.msg(f"Animation ready: {len(frames)} frames, amplification ×{self.anim_amp:.4g}")

    def _angles(self, th_sec):
        """Per display point rotation for the section rotations th_sec."""
        if self.cmb_disp.currentIndex() == 0:
            return th_sec[self.d_sec]
        sm = self.geo.mid_positions() * 1e3
        order = np.argsort(sm)
        return np.interp(self.d_s, sm[order], th_sec[order])

    def _goto_frame(self, i):
        if self.anim_frames is None:
            return
        i = int(np.clip(i, 0, len(self.anim_frames) - 1))
        self.anim_index = i
        th = self.anim_frames[i]
        ang = self._angles(th)
        ax = self._axis_mm()
        self.d_mesh.points = ax.rotate(self.d_p0, ang * self.anim_amp)
        self.d_mesh.point_data["rotation [deg]"] = np.degrees(ang)
        self.ptr_mesh.points = ax.rotate(self.ptr_p0, th[self.ptr_sec] * self.anim_amp)
        try:
            self.anim_text.SetText(2, self.anim_labels[i] + f"   (amplification ×{self.anim_amp:.3g})")
        except Exception:                           # noqa: BLE001
            pass
        self.lbl_frame.setText(self.anim_labels[i])
        self.plotter.render()

    def _advance_frame(self):
        if self.anim_frames is None:
            self.timer.stop()
            return
        i = self.anim_index + self.sp_speed.value()
        if i >= len(self.anim_frames):
            if not self.chk_loop.isChecked():
                self.toggle_play()
                return
            i = 0
        self.slider.setValue(i)

    def toggle_play(self):
        if self.timer.isActive():
            self.timer.stop()
            self.btn_play.setText("▶ Play")
        else:
            if self.anim_frames is None:
                self.prepare_animation()
                if self.anim_frames is None:
                    return
            self.timer.start(40)
            self.btn_play.setText("⏸ Pause")

    def export_animation(self):
        if self.anim_frames is None:
            QMessageBox.information(self, "Export", "Prepare an animation first.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export animation", "torsion.gif",
                                              "GIF (*.gif);;MP4 (*.mp4)")
        if not path:
            return
        self.timer.stop()
        try:
            if path.lower().endswith(".mp4"):
                self.plotter.open_movie(path, framerate=25)
            else:
                self.plotter.open_gif(path, fps=25)
            nf = len(self.anim_frames)
            for i in range(0, nf, max(1, nf // 300)):
                self._goto_frame(i)
                self.plotter.write_frame()
            self.plotter.mwriter.close()
            self.msg(f"Animation written to {path}")
        except Exception as exc:                    # noqa: BLE001
            QMessageBox.warning(self, "Export", f"Export failed ({exc}).\n"
                                "GIF needs 'imageio', MP4 needs 'imageio-ffmpeg'.")

    # ================================================================== #
    # project files
    # ================================================================== #
    def project_dict(self):
        d = json.loads(json.dumps(prj.DEFAULT_PROJECT))
        d["step_file"] = self.info.path if self.info else ""
        d["axis"] = {"mode": ["auto", "x", "y", "z", "custom"][self.cmb_axis.currentIndex()],
                     "origin_mm": [e.value() for e in self.ed_org],
                     "direction": [e.value() for e in self.ed_dir]}
        d["density"] = self.sp_rho.value()
        d["cuts_mm"] = sorted(self.cuts_mm())
        d["fem"] = {"E_GPa": self.sp_E.value(), "nu": self.sp_nu.value(), "torque": self.ed_T.value(),
                    "fixed_end": "start" if self.cmb_fix.currentIndex() == 0 else "end",
                    "mesh_size_mm": self.sp_h.value(), "order": 2 if self.cmb_order.currentIndex() == 0 else 1,
                    "end_band_mm": self.sp_band.value(),
                    "lumped_at": "com" if self.cmb_lump.currentIndex() == 1 else "middle"}
        try:
            d["model"] = self.current_model().to_dict()
        except Exception:                           # noqa: BLE001
            d["model"] = None
        d["loads"] = [ld.to_dict() for ld in self.loads()]
        d["simulation"] = {"t_end": self.ed_tend.value(), "dt": self.ed_dt.value()}
        d["frf"] = {"node": self.sp_frf_node.value(), "f_min": self.ed_fmin.value(),
                    "f_max": self.ed_fmax.value(), "n": 1500}
        return d

    def save_project(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save project", "shaft_project.json", "Project (*.json)")
        if path:
            prj.save_project(path, self.project_dict())
            self.msg(f"Project saved: {path}")

    def load_project(self, path=None):
        if not path:
            path, _ = QFileDialog.getOpenFileName(self, "Load project", "", "Project (*.json)")
        if not path:
            return
        d = prj.load_project(path)
        a = d["axis"]
        self.cmb_axis.blockSignals(True)
        self.cmb_axis.setCurrentIndex(["auto", "x", "y", "z", "custom"].index(a.get("mode", "auto")))
        self.cmb_axis.blockSignals(False)
        for k in range(3):
            self.ed_org[k].setValue(a["origin_mm"][k])
            self.ed_dir[k].setValue(a["direction"][k])
        self.sp_rho.setValue(d["density"])
        fe = d["fem"]
        self.sp_E.setValue(fe["E_GPa"])
        self.sp_nu.setValue(fe["nu"])
        self.ed_T.setValue(fe["torque"])
        self.cmb_fix.setCurrentIndex(0 if fe["fixed_end"] == "start" else 1)
        self.sp_h.setValue(fe["mesh_size_mm"])
        self.cmb_order.setCurrentIndex(0 if int(fe["order"]) == 2 else 1)
        self.sp_band.setValue(fe["end_band_mm"])
        self.cmb_lump.setCurrentIndex(1 if fe.get("lumped_at") == "com" else 0)
        self.tbl_loads.setRowCount(0)
        for ld in d["loads"]:
            self.add_load(Load(**ld))
        self.ed_tend.setValue(d["simulation"]["t_end"])
        self.ed_dt.setValue(d["simulation"]["dt"])
        self.ed_fmin.setValue(d["frf"]["f_min"])
        self.ed_fmax.setValue(d["frf"]["f_max"])
        self.sp_frf_node.setValue(int(d["frf"].get("node", 1)))
        self._pending_model = d.get("model")

        def after():
            self.set_cuts(d["cuts_mm"])
            self.msg(f"Project loaded: {path}  - press 'Split model' and 'Run FEM' to recompute.")

        if d["step_file"]:
            self.open_step(d["step_file"], after=after)
        else:
            after()

    def export_results(self):
        folder = QFileDialog.getExistingDirectory(self, "Export results to folder")
        if not folder:
            return
        model = None
        try:
            model = self.current_model()
        except Exception:                           # noqa: BLE001
            pass
        files = prj.export_results(folder, self.geo, self.res, model, self.sim)
        self.fig.savefig(os.path.join(folder, "current_plot.png"), dpi=150)
        self.msg("Exported: " + ", ".join(os.path.basename(f) for f in files) + ", current_plot.png")

    def save_screenshot(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save screenshot", "shaft.png", "PNG (*.png)")
        if path:
            self.plotter.screenshot(path)

    def closeEvent(self, ev):
        self.timer.stop()
        try:
            self.plotter.close()
        except Exception:                           # noqa: BLE001
            pass
        super().closeEvent(ev)


def main(argv=None):
    argv = sys.argv if argv is None else argv
    app = QApplication.instance() or QApplication(argv)
    win = MainWindow()
    win.show()
    if len(argv) > 1:
        arg = argv[1]
        if arg.lower().endswith(".json"):
            win.load_project(arg)
        else:
            win.open_step(arg)
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
