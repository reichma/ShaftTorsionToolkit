"""Start the Shaft Torsion Toolkit GUI:  python run_gui.py [file.step | project.json]"""
import os
import platform
import struct
import sys


def _dll_version(path):
    """File version of a Windows DLL as 'a.b.c.d' (empty string elsewhere)."""
    try:
        import ctypes
        from ctypes import wintypes
        ver = ctypes.windll.version
        size = ver.GetFileVersionInfoSizeW(path, None)
        buf = ctypes.create_string_buffer(size)
        ver.GetFileVersionInfoW(path, 0, size, buf)
        ptr, n = ctypes.c_void_p(), wintypes.UINT()
        ver.VerQueryValueW(buf, "\\", ctypes.byref(ptr), ctypes.byref(n))
        ms, ls = ctypes.cast(ptr, ctypes.POINTER(wintypes.DWORD * 4)).contents[2:4]
        return f"{ms >> 16}.{ms & 0xFFFF}.{ls >> 16}.{ls & 0xFFFF}"
    except Exception:
        return "?"


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
    conda = ("conda" in sys.version.lower() or "anaconda" in sys.base_prefix.lower()
             or os.path.exists(os.path.join(sys.base_prefix, "conda-meta")))
    if conda:
        lines.append("NOTE       : the venv is based on an Anaconda/conda Python.")
    if "WindowsApps" in sys.executable:
        lines.append("NOTE       : this is the Microsoft Store Python.")
    # Windows searches the folder of the real python.exe first, so its C++ runtime wins.
    exe_dir = os.path.dirname(os.path.realpath(getattr(sys, "_base_executable", sys.executable)))
    for d in dict.fromkeys([exe_dir, sys.base_prefix]):
        dll = os.path.join(d, "msvcp140.dll")
        if os.path.isfile(dll):
            lines += [f"C++ runtime: {dll}  version {_dll_version(dll)}",
                      "             Windows loads this copy first; PySide6 >= 6.8 needs 14.40 or newer."]
    hits, seen = [], set()
    for d in os.environ.get("PATH", "").split(os.pathsep):
        key = os.path.normcase(os.path.normpath(d)) if d else ""
        if (not key or key in seen or "system32" in key
                or key.endswith(os.path.normcase(os.path.join("site-packages", "PySide6")))):
            continue
        seen.add(key)
        try:
            names = os.listdir(d)
        except OSError:
            continue
        found = [n for n in names if n.lower() in ("qt6core.dll", "qt5core.dll", "msvcp140.dll")]
        if found:
            hits.append(f"    {d}  ({', '.join(found)})")
    if hits:
        lines += ["", "Other folders on PATH with their own Qt or C++ runtime DLLs:"] + hits
    lines += ["", "Fixes, most likely first:"]
    if conda:
        lines += ["  * Anaconda's own msvcp140.dll is too old for this PySide6. Either",
                  "      - install an older Qt into the venv:",
                  "          .venv\\Scripts\\python -m pip install \"PySide6<6.8\"",
                  "      - or delete .venv and recreate it from a python.org Python",
                  "        (install it, then run setup_venv.bat, which prefers the py launcher)."]
    lines += ["  * Install the latest Microsoft Visual C++ Redistributable (x64) and restart:",
              "      https://aka.ms/vs/17/release/vc_redist.x64.exe",
              "  * Use a 64-bit python.org Python (not Anaconda / not the Store),",
              "    delete the .venv folder and run setup_venv.bat again.",
              "  * Remove the folders listed above from PATH, or start from a clean prompt.",
              "  * Reinstall Qt:  .venv\\Scripts\\python -m pip install --force-reinstall PySide6",
              "  * Current PySide6 needs Windows 10 or newer.",
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
