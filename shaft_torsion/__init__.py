"""Shaft torsion toolkit.

Workflow
--------
1. ``geometry``  - load a STEP file, define the shaft axis and cut positions,
   split the CAD model and compute the mass moment of inertia of every section.
2. ``fem``       - mesh the split model, clamp one end, apply a torque at the
   other end and derive the torsional stiffness between the sections.
3. ``dynamics``  - lumped torsional vibration model (inertias, springs,
   dampers, loads): natural frequencies, mode shapes, time response, FRF.
4. ``gui``       - PySide6/PyVista front end that animates the result on the
   split CAD model.

All internal quantities are SI (m, kg, s, N, rad).
"""

__version__ = "1.0.0"
