"""Gearbox output shaft (realistic test case) with keyed spur-gear blank, as a 3-solid STEP assembly (mm).

Shaft along +X (0..420 mm):
  0- 60  d30  drive end, keyway 8x4 (DIN 6885) x=8..52, 1 mm chamfer
 60- 85  d35  bearing seat, circlip groove d33 x 1.6 at x=80
 85- 95  d45  shoulder (fillets R1.5 / R1)
 95-165  d42  gear seat, keyway 12x5, x=110..150
165-310  d40  body
310-320  d45  shoulder
320-345  d35  bearing seat, circlip groove at x=339
345-420  d25  output end, centre bore d10 x 50 from the right end, 1 mm chamfer
Gear blank x=105..155: hub d70 (50 wide), web 15 thick, rim d130..d150 (40 wide),
bore d42 with keyway 12 x 3.  Feather key 12x8x40 (form A) bonded in both keyways.
"""
import os
import sys

import gmsh

# usage: python make_gear_shaft.py [output.step]   (default: gear_shaft.step next to this script)
out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "gear_shaft.step")
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("gear_shaft")
occ = gmsh.model.occ

def cyl(x0, x1, r):
    return (3, occ.addCylinder(x0, 0, 0, x1 - x0, 0, 0, r))

def keyshape(x0, x1, w, y0, y1):
    """Rounded-end key / keyway (form A) between x0..x1, width w, height y0..y1."""
    r = w / 2
    b = (3, occ.addBox(x0 + r, y0, -r, x1 - x0 - 2 * r, y1 - y0, w))
    c1 = (3, occ.addCylinder(x0 + r, y0, 0, 0, y1 - y0, 0, r))
    c2 = (3, occ.addCylinder(x1 - r, y0, 0, 0, y1 - y0, 0, r))
    v, _ = occ.fuse([b], [c1, c2])
    return v

steps = [(0, 60, 15), (60, 85, 17.5), (85, 95, 22.5), (95, 165, 21), (165, 310, 20),
         (310, 320, 22.5), (320, 345, 17.5), (345, 420, 12.5)]
parts = [cyl(*s) for s in steps]
shaft, _ = occ.fuse([parts[0]], parts[1:])
occ.synchronize()
# fillets R1.5 at the concave shoulder corners (x=85, 95, 310, 320)
edges, radii = [], []
for x, r, rf in [(85, 17.5, 1.5), (95, 21, 1.0), (310, 20, 1.5), (320, 17.5, 1.5)]:
    for d, t in gmsh.model.getEntitiesInBoundingBox(x - 0.01, -r - 0.01, -r - 0.01, x + 0.01, r + 0.01, r + 0.01, 1):
        edges.append(t)
        radii.append(rf)
shaft = occ.fillet([shaft[0][1]], edges, radii)
occ.synchronize()
# chamfers at both free ends
ends = []
for x, r in [(0, 15), (420, 12.5)]:
    ends += [t for _, t in gmsh.model.getEntitiesInBoundingBox(x - 0.01, -r - 0.01, -r - 0.01, x + 0.01, r + 0.01, r + 0.01, 1)]
faces = []
for x, r in [(0, 15), (420, 12.5)]:
    faces += [t for _, t in gmsh.model.getEntitiesInBoundingBox(x - 0.01, -r - 0.01, -r - 0.01, x + 0.01, r + 0.01, r + 0.01, 2)]
shaft = occ.chamfer([shaft[0][1]], ends, faces, [1.0] * len(ends))
occ.synchronize()
# grooves, keyways, centre bore
tools = [(3, occ.addCylinder(80, 0, 0, 1.6, 0, 0, 17.5)), (3, occ.addCylinder(339, 0, 0, 1.6, 0, 0, 17.5))]
tools[0] = occ.cut([tools[0]], [cyl(80, 81.6, 16.5)])[0][0]
tools[1] = occ.cut([tools[1]], [cyl(339, 340.6, 16.5)])[0][0]
tools += keyshape(8, 52, 8, 11, 16)          # d30: t1 = 4
tools += keyshape(110, 150, 12, 16, 22)      # d42: t1 = 5
tools.append(cyl(370, 421, 5))               # centre bore d10 x 50
shaft, _ = occ.cut(shaft, tools)

# gear blank
hub = cyl(105, 155, 35)
web = cyl(122.5, 137.5, 66)
rim = cyl(110, 150, 75)
gear, _ = occ.fuse([hub], [web, rim])
gear, _ = occ.cut(gear, [cyl(110.01, 121, 65), cyl(139, 149.99, 65)])
gear, _ = occ.cut(gear, [cyl(104, 156, 21)] + keyshape(104, 156, 12, 20, 24))
# feather key 12x8x40 in shaft keyway (y 16..21) + gear keyway (y 21..24)
key = keyshape(110, 150, 12, 16, 24)
occ.synchronize()
print("volumes:", occ.getEntities(3))
for d, t in occ.getEntities(3):
    print(t, round(occ.getMass(3, t), 1), "mm3")
gmsh.write(out)
gmsh.finalize()
