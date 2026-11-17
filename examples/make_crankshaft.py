"""Inline crankshaft with pulley and flywheel as a STEP solid (mm), for 3, 4 or 5 cylinders.

    python make_crankshaft.py 4 [output.step] [project.json]

Writes the STEP file (default: crankshaft_<n>cyl.step next to this script) and
a project file with the cut planes (default: crankshaft_<n>cyl_project.json).

Layout along +X (generic passenger-car size, not a specific engine):
  0 -  50  nose d32 with damper pulley d150 x 25 (x = 5..30)
  then per cylinder: main journal d58 x 26, web 19, crank pin d48 x 24, web 19
  last main journal, flange d100 x 12, flywheel d300 x 30
Stroke 86 mm (crank radius 43), cylinder pitch 88 mm. Every web carries a
counterweight (R75 half disc) opposite its crank pin.

Crank pin angles (firing order):
  3 cyl  0, 120, 240           (1-3-2,     240 deg firing interval)
  4 cyl  0, 180, 180, 0        (1-3-4-2,   180 deg)
  5 cyl  0, 144, 216, 288, 72  (1-2-4-5-3, 144 deg)

The cut planes lie in the middle of the main journals, so the lumped model has
one inertia for the pulley end, one per crank throw and one for the flywheel end.
"""
import json
import math
import os
import sys

import gmsh

PIN_ANGLES = {3: [0, 120, 240], 4: [0, 180, 180, 0], 5: [0, 144, 216, 288, 72]}
R_CRANK, R_MAIN, R_PIN = 43.0, 29.0, 24.0
W_MAIN, W_WEB, W_PIN = 26.0, 19.0, 24.0
NOSE = 50.0


def build(n, out_step):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(f"crankshaft_{n}cyl")
    occ = gmsh.model.occ

    def cyl(x0, x1, r, y=0.0, z=0.0):
        return (3, occ.addCylinder(x0, y, z, x1 - x0, 0, 0, r))

    def web(x0, angle):
        """Web with crank-pin boss and counterweight, pin direction rotated by angle about X."""
        x1 = x0 + W_WEB
        parts = [cyl(x0, x1, 36.0), cyl(x0, x1, 30.0, y=R_CRANK),
                 (3, occ.addBox(x0, 0, -30, W_WEB, R_CRANK, 60))]
        disc = cyl(x0, x1, 75.0)
        half = (3, occ.addBox(x0 - 1, -80, -80, W_WEB + 2, 72, 160))     # side opposite the pin
        cw, _ = occ.intersect([disc], [half])
        v, _ = occ.fuse(parts[:1], parts[1:] + cw)
        occ.rotate(v, 0, 0, 0, 1, 0, 0, math.radians(angle))
        return v

    solids = [cyl(0, NOSE, 16.0), cyl(5, 30, 75.0)]                     # nose + damper pulley
    x = NOSE
    for a in PIN_ANGLES[n]:
        solids.append(cyl(x, x + W_MAIN, R_MAIN))                        # main journal
        x += W_MAIN
        solids += web(x, a)
        x += W_WEB
        pin = cyl(x, x + W_PIN, R_PIN, y=R_CRANK)
        occ.rotate([pin], 0, 0, 0, 1, 0, 0, math.radians(a))
        solids.append(pin)
        x += W_PIN
        solids += web(x, a)
        x += W_WEB
    solids.append(cyl(x, x + W_MAIN, R_MAIN))                            # last main journal
    x += W_MAIN
    solids.append(cyl(x, x + 12, 50.0))                                  # flange
    solids.append(cyl(x + 12, x + 42, 150.0))                            # flywheel
    length = x + 42
    v, _ = occ.fuse(solids[:1], solids[1:])
    occ.synchronize()
    if len(v) != 1:
        raise RuntimeError(f"crankshaft is not one solid ({len(v)} volumes)")
    gmsh.write(out_step)
    gmsh.finalize()
    pitch = W_MAIN + 2 * W_WEB + W_PIN
    cuts = [NOSE + W_MAIN / 2 + i * pitch for i in range(n + 1)]       # middle of every main journal
    return cuts, length


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    n = int(argv[0]) if argv else 4
    if n not in PIN_ANGLES:
        raise SystemExit("number of cylinders must be 3, 4 or 5")
    here = os.path.dirname(os.path.abspath(__file__))
    step = argv[1] if len(argv) > 1 else os.path.join(here, f"crankshaft_{n}cyl.step")
    proj = argv[2] if len(argv) > 2 else os.path.join(here, f"crankshaft_{n}cyl_project.json")
    cuts, length = build(n, step)
    project = {
        "step_file": os.path.relpath(step, os.path.dirname(os.path.abspath(proj))).replace("\\", "/"),
        "axis": {"mode": "custom", "origin_mm": [0, 0, 0], "direction": [1, 0, 0]},
        "density": 7850.0,
        "cuts_mm": cuts,
        "fem": {"E_GPa": 210.0, "nu": 0.3, "torque": 1000.0, "fixed_end": "start",
                "mesh_size_mm": 0.0, "order": 2, "end_band_mm": 0.0, "lumped_at": "middle"},
        "model": None,
        "loads": [],
        "simulation": {"t_end": 0.1, "dt": 0.0},
        "frf": {"node": 1, "f_min": 10.0, "f_max": 5000.0, "n": 1500},
    }
    with open(proj, "w", newline="\n") as fh:
        json.dump(project, fh, indent=2)
    print(f"{n} cylinders: length {length:.0f} mm, cuts {cuts}\n  {step}\n  {proj}")


if __name__ == "__main__":
    main()
