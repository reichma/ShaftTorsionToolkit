# Shaft Torsion Toolkit

A Python tool with a GUI that takes a shaft from a STEP file and does three things:

1. **Split and inertia.** You set cut planes along the shaft axis. The tool splits the CAD model at those planes and computes each section's mass moment of inertia about the axis.
2. **Torsional stiffness (FEM).** It meshes the split model with quadratic tetrahedra, clamps one end and applies a torque to the other end through a rigid end plate. It then reads the rotation in the middle of every section and computes the stiffness between neighbouring sections.
3. **Torsional vibration model.** It builds a model from the inertias of step 1 and the springs of step 2. You add dampers between the inertias, damping to ground and loads. It gives natural frequencies, mode shapes, a time simulation and frequency response.

It then **animates** the time response, a mode shape or the static FEM twist on the split CAD model.

Full manual (HTML, works offline): open [`docs/index.html`](docs/index.html) in a browser.

```
STEP ──► 1 split + inertia J_i ──► 2 FEM: θ at the middle of each section ──► k_i = T / Δθ
                                                                        │
         4 animation on split CAD ◄── 3 J θ'' + C θ' + K θ = M(t) ◄─────┘
```

## Installation

Windows: run `setup_venv.bat` once, then start with `run_gui.bat`. You can pass a `.step` file or a project `.json` as an argument.

Any OS:

```bash
python -m venv .venv && . .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python run_gui.py examples/stepped_shaft_project.json
```

On Linux, gmsh needs the system library `libGLU` (`apt install libglu1-mesa`).

**Windows: `ImportError: DLL load failed while importing QtCore`.** PySide6 cannot load its Qt DLLs. `run_gui.bat` then prints a diagnosis. Usual fixes, most likely first:

1. **Venv made from Anaconda Python** (error text: *procedure not found* / *Prozedur wurde nicht gefunden*): Anaconda's own `msvcp140.dll` is loaded first and is too old for PySide6 6.8+. Run `.venv\Scripts\python -m pip install "PySide6<6.8"`, or delete `.venv` and recreate it from a python.org Python.
1. Install the latest [Microsoft Visual C++ Redistributable (x64)](https://aka.ms/vs/17/release/vc_redist.x64.exe) and restart.
1. Use a 64-bit Python from python.org (not Anaconda or the Microsoft Store), delete `.venv` and run `setup_venv.bat` again.
1. Remove folders with their own Qt or `msvcp140.dll` (Anaconda, other Qt apps) from `PATH`.
1. Reinstall: `.venv\Scripts\python -m pip install --force-reinstall PySide6`.

## Using the GUI

### Tab 1: Geometry

* **Open STEP file.** Assemblies with several solids work. Touching solids are fused at their contact faces (bonded).
* **Shaft axis.** By default the tool finds the axis automatically: it is the principal axis of inertia that stands out from the other two. If the axis is within 1° of X, Y or Z, it snaps to that axis. You can also choose X, Y or Z through the centre of mass, or enter a custom origin and direction in mm.
* **Cut planes.** Positions are in mm from the shaft start. You can type them, use *Equal sections*, or right-click on the model with *Pick* enabled. The planes are shown in red before splitting.
* **Split model & compute inertias.** This gives a table with mass and J [kg·m²] for each section, and colours each section differently in the 3D view.

### Tab 2: FEM stiffness

| Setting | Meaning |
|---|---|
| E, ν, torque | Material and the torque applied to the loaded end. The model is linear, so the torque value only scales the angles. |
| Fixed end | Which end is clamped (all DOFs fixed). The other end is twisted through a rigid coupling (RBE2-like). |
| Element size | `auto` = length/40, but at least ~4 elements across the thinnest section. |
| Element type | Quadratic tetrahedra are recommended. Linear tetrahedra are 10–15 % too stiff in torsion. |
| End-face band | Nodes within this distance of an end belong to the clamp or load face. `auto` = 0.2 % of the length. |
| Lumped inertia at | Where the section rotation is measured: middle of the section, or its centre of mass. |

How the stiffness is computed:

* The rotation at an axial position is the least-squares rigid rotation of all nodes in a thin slice there: `u_t = r·(θ + θ'·Δs)`.
* The spring between inertias i and i+1 is `k_i = T / |θ_mid,i+1 − θ_mid,i|`.
* The table also shows each section's own stiffness (between its cut planes), the stiffness from the clamp to the nearest inertia, and the total stiffness.

### Tab 3: Vibration model

* **J, k and c are editable.** You can add attached parts, a coupling stiffness, and so on.
* **Damping:**
  * a damper `c` per spring (between inertias)
  * damping to ground per inertia (for example, bearings)
  * *c = β·k* (stiffness-proportional)
  * *c for ζ at the 1st elastic mode*
* **Boundary condition:** free–free, or grounded at the FEM clamp. Grounded uses the spring from the clamp to the nearest inertia (from the FEM) plus an optional damper.
* **Loads** (torque on inertia #):
  * constant, step, ramp, pulse
  * sine (with offset and phase)
  * linear sweep
  * a formula in `t`, for example `500*sin(2*pi*50*t)*exp(-5*t)`
* **Simulation** uses Newmark average acceleration, which is unconditionally stable. With Δt = 0 the step size is automatic: 20 steps per period of the highest mode.
* **Plots:** angles, twist, spring torques, speeds, loads, and a receptance FRF.

### Tab 4: Animation

* **Source:** time simulation, a mode shape, or the static FEM rotation.
* **Display:**
  * *rigid sections*: each CAD section turns by the rotation of its inertia.
  * *smooth twist*: the rotation is interpolated along the axis.
* **Colouring:** by true rotation angle or by section. Pointers show the rotation of each section.
* **Amplification** is automatic (20° maximum shown) or set by hand.
* **Export** to GIF (needs `imageio`) or MP4 (needs `imageio-ffmpeg`).

### File menu

* Save or load a project (`.json`). The STEP path, axis, cuts, material, FEM settings, model values and loads are stored. Results are recomputed after loading.
* Export results as CSV:
  * `sections.csv`
  * `torsion_fem.csv`
  * `modes.csv`
  * `time_response.csv`

  The CSV files use `;` as the separator.
* Save a screenshot.

## Batch mode (no GUI)

```bash
python -m shaft_torsion.cli examples/stepped_shaft_project.json results_folder
```

This runs the full chain and writes the CSV files plus `time_response.png`.

Project file notes:

* Lengths are in mm. Everything else is SI.
* Load `node` is **0-based** in the JSON. The GUI shows it 1-based.
* `"model": null` means the model is taken from the inertias and the FEM results.

## Verification

Run `python -m pytest tests -q`. The tests check against closed-form results:

* **Section inertias** of a Ø40 mm cylinder along an oblique axis: `ρ·π·d⁴/32·L`, accurate to 1e-6.
* **FEM stiffness** (between inertias, per section, clamp to inertia, total) against `G·Ip/L`: within 0.2 %. With a 6 mm mesh the error is 0.001 %.
* **Two-inertia natural frequency**, the **step response** of a single grounded inertia, and the **damping ratio**.

## Notes / limits

* The model is linear-elastic and small-strain with one isotropic material. Use the editable tables for parts made of other materials.
* Positions are measured along the axis. Cut planes are always normal to the shaft axis.
* The rigid coupling at the loaded end also blocks warping of that face, which matters only for non-circular end sections. The springs are measured between section middles, away from the ends.
* Large meshes: install `pypardiso` (fastest). Without it, AMG-preconditioned CG (`pyamg`) is used, then SuperLU.

## Files

| File | Content |
|---|---|
| `shaft_torsion/geometry.py` | STEP import, axis, splitting, inertias, display tessellation |
| `shaft_torsion/fem.py` | Mesh, tet4/tet10 stiffness assembly, boundary conditions, solver, section rotations, stiffness |
| `shaft_torsion/dynamics.py` | Lumped model, loads, modes, Newmark simulation, FRF |
| `shaft_torsion/gui.py` | PySide6 + PyVista GUI |
| `shaft_torsion/project.py`, `cli.py` | Project JSON, CSV export, batch run |
| `examples/` | `plain_shaft.step`, `stepped_shaft.step` (generator script included), example project |
