"""Create example STEP files (units: mm).

  plain_shaft.step   - solid cylinder d = 40 mm, L = 400 mm (analytical check)
  stepped_shaft.step - stepped shaft with a flywheel disc, a gear-like hub,
                       a centre bore and a keyway, axis along X
"""
import os
import sys

import gmsh

here = os.path.dirname(os.path.abspath(__file__))


def plain(path):
    gmsh.model.add("plain")
    gmsh.model.occ.addCylinder(0, 0, 0, 400, 0, 0, 20)
    gmsh.model.occ.synchronize()
    gmsh.write(path)
    gmsh.model.remove()


def stepped(path):
    occ = gmsh.model.occ
    gmsh.model.add("stepped")
    parts = [
        occ.addCylinder(0, 0, 0, 120, 0, 0, 20),        # journal d40
        occ.addCylinder(120, 0, 0, 40, 0, 0, 110),      # flywheel disc d220 x 40
        occ.addCylinder(160, 0, 0, 260, 0, 0, 25),      # main shaft d50
        occ.addCylinder(420, 0, 0, 60, 0, 0, 45),       # hub d90
        occ.addCylinder(480, 0, 0, 120, 0, 0, 17.5),    # output stub d35
    ]
    body, _ = occ.fuse([(3, parts[0])], [(3, p) for p in parts[1:]])
    bore = occ.addCylinder(-1, 0, 0, 100, 0, 0, 8)                  # centre bore d16
    key = occ.addBox(500, -5, 13.5, 101, 10, 10)                    # keyway 10 x 4 mm
    body, _ = occ.cut(body, [(3, bore), (3, key)])
    occ.synchronize()
    gmsh.write(path)
    gmsh.model.remove()


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else here
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    plain(os.path.join(out, "plain_shaft.step"))
    stepped(os.path.join(out, "stepped_shaft.step"))
    gmsh.finalize()
    print("written to", out)
