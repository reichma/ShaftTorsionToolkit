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

* **With Anaconda / Miniconda installed** (found via an active conda shell, `PATH`, or the usual install folders), `setup_venv.bat` creates a conda environment in `.conda_env` with Python 3.12 from **conda-forge** and installs `requirements.txt` into it with pip. conda-forge's Python brings a current C++ runtime next to its `python.exe`, so PySide6 loads; a venv made from Anaconda's own Python cannot load current PySide6 (see below). Only the free conda-forge channel is used.
* **Without conda** it creates a plain venv in `.venv`. `setup_venv.bat venv` forces this.
* `run_gui.bat` uses `.venv` if it exists, otherwise `.conda_env`. Delete either folder and run `setup_venv.bat` again to rebuild it.

Any OS:

```bash
python -m venv .venv && . .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python run_gui.py examples/stepped_shaft_project.json
```

On Linux, gmsh needs the system library `libGLU` (`apt install libglu1-mesa`).

**Windows: `ImportError: DLL load failed while importing QtCore`.** PySide6 cannot load its Qt DLLs. `setup_venv.bat` tests this at the end (`python run_gui.py --check`) and prints a diagnosis with the versions of every `msvcp140.dll` involved. On Python ≤ 3.12 it then retries with `PySide6<6.8`, which works with older C++ runtimes; for Python 3.13+ no such build exists. Usual fixes, most likely first:

1. Install the latest [Microsoft Visual C++ Redistributable (x64)](https://aka.ms/vs/17/release/vc_redist.x64.exe) and restart.
1. With Anaconda: delete `.venv` and run `setup_venv.bat` again. It builds `.conda_env` from conda-forge instead of a venv from Anaconda's Python.
1. Without Anaconda: use a 64-bit **Python 3.12** from python.org (not the Microsoft Store), delete `.venv` and run `setup_venv.bat` again. The script prefers Python 3.12, then 3.11, then 3.13 when the `py` launcher has them.
1. Remove folders with their own Qt or `msvcp140.dll` (Anaconda, MiKTeX, other Qt apps) from `PATH`.
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
  * damping to ground per inertia (for example, bearing drag)
  * *c = β·k* (stiffness-proportional)
  * *c for ζ at the 1st elastic mode*
* **Spring to ground per inertia** (*k to ground* column), for example an elastic mount or a coupling to the housing at any position along the shaft. *c = β·k* also sets the damper of every inertia that has such a spring.
* **Boundary condition:** free–free, or grounded at the FEM clamp. Grounded uses the spring from the clamp to the nearest inertia (from the FEM) plus an optional damper.
* **Loads** (torque on inertia #):
  * constant, step, ramp, pulse
  * sine (with offset and phase)
  * linear sweep
  * a formula in `t`, for example `500*sin(2*pi*50*t)*exp(-5*t)`
* **Simulation** uses Newmark average acceleration, which is unconditionally stable. With Δt = 0 the step size is automatic: 20 steps per period of the highest mode.
* **Plots:** angles, twist, spring torques, speeds, loads, and a receptance FRF.
* **Legend:** click an entry to hide or show that channel; double-click to show only that channel, and double-click it again to show all. The y axis rescales to the visible channels.

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

Saving the 3-D views and animations without opening a window:

```bash
python examples/render_3d.py examples/gear_shaft_mounted_sweep.json [output_folder] [--modes 4]
```

This runs the project through the GUI, with the window placed outside the visible screen, and writes PNG screenshots and GIF animations of the split model, the FEM rotation, the mode shapes and the time simulation at its largest rotation, plus the 2-D plots.

### Example: mounted gear shaft under a frequency sweep

`examples/gear_shaft_mounted_sweep.json` uses the gear shaft without the FEM clamp (free–free). The shaft hangs between two elastic couplings, entered as *k to ground*: 5000 N·m/rad on the drive end (inertia 1) and 10 000 N·m/rad on the output end (inertia 7). The modes are 168 Hz (the shaft rocks on the couplings), 3.99 kHz (the output end against the gear), and 8.57 and 9.18 kHz (the drive end and output end against the gear).

A 20 N·m sweep on the output end excites them in two parts. A slow sweep from 20 to 400 Hz builds up the 168 Hz resonance to about 15 times the static rotation. A fast sweep from 1 to 10 kHz then excites the higher modes. A single linear sweep from 20 Hz to 10 kHz passes the first mode in about 1 ms, which is too short for the resonance to build up. The gear is a node of modes 2–4, so a load on the gear would excite only mode 1.

The J and k values in the file come from a FEM run with the automatic mesh. Load the file in the GUI, then split and run the FEM; the model is then built from the values in the file.

### Example: crankshafts with 3, 4 and 5 cylinders

`examples/make_crankshaft.py <3|4|5>` builds an inline crankshaft as one STEP solid and writes a project file for it. The crankshaft has a damper pulley at the front, a main bearing between every pair of cylinders, webs with counterweights and a flywheel. The crankshafts are generic (stroke 86 mm, cylinder pitch 88 mm, no fillets or oil bores), not copies of real engines. The crank pins follow the usual firing orders: 3 cylinders 0/120/240°, 4 cylinders 0/180/180/0°, 5 cylinders 0/144/216/288/72°. The cut planes lie in the middle of the main journals, which gives one inertia for the pulley end, one per crank throw and one for the flywheel end.

The project files `crankshaft_<n>cyl_project.json` contain a free–free model with stiffness-proportional damping (ζ = 2 % at the first elastic mode). A 100 N·m sweep acts on the first crank throw, from 0.6 to 1.4 times the first elastic frequency over 1 s. It starts as a cosine (`phase_deg` 90): a sine start gives the free–free shaft a net angular impulse, and the whole shaft then turns away as a rigid body.

| | 3 cylinders | 4 cylinders | 5 cylinders |
|---|---|---|---|
| total J (of which flywheel) [kg·m²] | 0.228 (0.188) | 0.238 (0.188) | 0.248 (0.188) |
| elastic modes [Hz] | 427, 1051, 1616, 2269 | 360, 837, 1397, 2145, 2610 | 311, 742, 1256, 1613, 2131, 2292 |
| largest rotation in the sweep | 0.72° at 435 Hz | 0.87° at 368 Hz | 1.00° at 319 Hz |
| largest torque in the shaft | 1852 N·m | 2040 N·m | 2095 N·m |

In the first elastic mode, the pulley end swings against the flywheel, which barely moves. The longer the crankshaft, the softer it is, so this mode drops from 427 Hz to 311 Hz. `python examples/render_3d.py examples/crankshaft_4cyl_project.json` renders the 3-D mode shapes and the sweep response; it skips the 0 Hz rigid-body mode.

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
| `examples/` | `plain_shaft.step`, `stepped_shaft.step` (generator script included), example project; `gear_shaft.step` + `gear_shaft_project.json`: 420 mm gearbox output shaft with a keyed gear (3-solid assembly, generator `make_gear_shaft.py`); `gear_shaft_mounted_sweep.json`: the same shaft on two elastic couplings under a frequency sweep; `render_3d.py`: saves the 3-D results of a project as PNG and GIF; `make_crankshaft.py` + `crankshaft_<3|4|5>cyl.step` / `_project.json`: inline crankshafts with pulley and flywheel |
