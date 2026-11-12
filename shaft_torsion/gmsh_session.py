"""Small helper that gives every step a clean gmsh session."""
from __future__ import annotations

import contextlib
import threading

import gmsh

_LOCK = threading.Lock()


@contextlib.contextmanager
def gmsh_session(verbose: bool = False):
    """Initialise gmsh (thread safe, no signal handler) and start an empty model.

    Geometry is always imported in metres (``Geometry.OCCTargetUnit = M``),
    whatever unit the STEP file uses.
    """
    with _LOCK:
        if not gmsh.isInitialized():
            gmsh.initialize(interruptible=False)
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.clear()
        gmsh.option.setNumber("General.Terminal", 1 if verbose else 0)
        gmsh.option.setNumber("General.Verbosity", 2 if verbose else 0)
        gmsh.option.setString("Geometry.OCCTargetUnit", "M")
        gmsh.model.add("shaft")
        try:
            yield gmsh
        finally:
            gmsh.clear()
