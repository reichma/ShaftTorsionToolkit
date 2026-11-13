"""Start the Shaft Torsion Toolkit GUI:  python run_gui.py [file.step | project.json]"""
import os
import platform
import struct
import sys


def _qt_diagnosis(err):
    """Explain the usual causes of 'DLL load failed while importing QtCore' on Windows."""
    lines = ["", "=" * 72,
             "The Qt library (PySide6) could not be loaded:",
             f"    {type(err).__name__}: {err}", "",
             f"Python     : {sys.version.split()[0]} ({struct.calcsize('P') * 8}-bit) {sys.executable}",
             f"OS         : {platform.platform()}"]
    try:
        from importlib.metadata import version
        lines.append(f"PySide6    : {version('PySide6')}")
    except Exception:
        lines.append("PySide6    : not installed")
    if "conda" in sys.version.lower() or os.path.exists(os.path.join(sys.prefix, "conda-meta")) \
            or os.environ.get("CONDA_PREFIX"):
        lines.append("NOTE       : this looks like an Anaconda/conda Python (conflicting Qt DLLs).")
    if "WindowsApps" in sys.executable:
        lines.append("NOTE       : this is the Microsoft Store Python.")
    hits = []
    for d in os.environ.get("PATH", "").split(os.pathsep):
        try:
            names = os.listdir(d)
        except OSError:
            continue
        found = [n for n in names if n.lower() in ("qt6core.dll", "qt5core.dll", "msvcp140.dll")]
        if found and "system32" not in d.lower():
            hits.append(f"    {d}  ({', '.join(found)})")
    if hits:
        lines += ["", "Folders on PATH with their own Qt or C++ runtime DLLs (can load the wrong DLL):"] + hits
    lines += ["",
              "Fixes, most likely first:",
              "  1. Install the latest Microsoft Visual C++ Redistributable (x64):",
              "     https://aka.ms/vs/17/release/vc_redist.x64.exe  - then restart the PC.",
              "  2. Use a 64-bit Python from python.org (not Anaconda / not the Store),",
              "     delete the .venv folder and run setup_venv.bat again.",
              "  3. Remove the folders listed above from PATH (or start from a clean",
              "     command prompt) so Python does not pick up foreign Qt / runtime DLLs.",
              "  4. Reinstall Qt:  .venv\\Scripts\\python -m pip install --force-reinstall PySide6",
              "  5. Windows 7/8 are not supported by current PySide6 (needs Windows 10+).",
              "=" * 72, ""]
    return "\n".join(lines)


if __name__ == "__main__":
    try:
        import PySide6.QtCore  # noqa: F401  - load Qt first to give a clear message
    except ImportError as exc:
        print(_qt_diagnosis(exc), file=sys.stderr)
        sys.exit(1)
    from shaft_torsion.gui import main
    sys.exit(main())
