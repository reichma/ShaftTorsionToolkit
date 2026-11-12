"""Project file (JSON) and the complete batch pipeline used by the CLI.

The GUI writes the same JSON with *File > Save project*.  Lengths in the
project file are millimetres, everything else SI.
"""
from __future__ import annotations

import csv
import json
import os

import numpy as np

from . import geometry as g
from .dynamics import Load, TorsionalModel
from .fem import FemSettings, run_torsion

DEFAULT_PROJECT = {
    "step_file": "",
    "axis": {"mode": "auto", "origin_mm": [0, 0, 0], "direction": [1, 0, 0]},
    "density": 7850.0,
    "cuts_mm": [],
    "fem": {"E_GPa": 210.0, "nu": 0.3, "torque": 1000.0, "fixed_end": "start",
            "mesh_size_mm": 0.0, "order": 2, "end_band_mm": 0.0, "lumped_at": "middle"},
    "model": None,
    "loads": [],
    "simulation": {"t_end": 0.5, "dt": 0.0},
    "frf": {"node": 1, "f_min": 1.0, "f_max": 2000.0, "n": 800},
}


def make_axis(info: g.StepInfo, mode="auto", origin_mm=None, direction=None) -> g.Axis:
    mode = (mode or "auto").lower()
    if mode == "auto":
        return info.auto_axis
    if mode in ("x", "y", "z"):
        return g.Axis(info.com, np.eye(3)["xyz".index(mode)])
    return g.Axis(np.asarray(origin_mm, float) / 1000.0, direction)


def fem_settings(d) -> FemSettings:
    return FemSettings(E=d["E_GPa"] * 1e9, nu=d["nu"], torque=d["torque"],
                       mesh_size=d.get("mesh_size_mm", 0) / 1000.0, order=int(d.get("order", 2)),
                       fixed_end=d.get("fixed_end", "start"), end_band=d.get("end_band_mm", 0) / 1000.0,
                       mid_from_com=d.get("lumped_at", "middle") == "com")


def model_from_results(geo, res, ground=True, c=None, c_abs=None) -> TorsionalModel:
    n = len(geo.sections)
    return TorsionalModel(J=[s.J for s in geo.sections], k=res.k_spring,
                          c=np.zeros(n - 1) if c is None else c,
                          c_abs=np.zeros(n) if c_abs is None else c_abs,
                          ground_node=res.ground_index if ground else -1,
                          k_ground=res.k_ground)


def model_from_dict(d) -> TorsionalModel:
    return TorsionalModel(J=d["J"], k=d["k"], c=d.get("c"), c_abs=d.get("c_abs"),
                          ground_node=d.get("ground_node", -1), k_ground=d.get("k_ground", 0.0),
                          c_ground=d.get("c_ground", 0.0))


def load_project(path):
    with open(path, encoding="utf-8") as fh:
        d = json.load(fh)
    out = json.loads(json.dumps(DEFAULT_PROJECT))
    for k, v in d.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k].update(v)
        else:
            out[k] = v
    if out["step_file"] and not os.path.isabs(out["step_file"]):
        out["step_file"] = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(path)),
                                                         out["step_file"]))
    return out


def save_project(path, d):
    d = dict(d)
    try:
        d["step_file"] = os.path.relpath(d["step_file"], os.path.dirname(os.path.abspath(path)))
    except ValueError:      # different drive on Windows
        pass
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(d, fh, indent=2)


def _write_csv(path, header, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter=";")
        w.writerow(header)
        w.writerows(rows)


def export_results(folder, geo=None, res=None, model=None, sim=None):
    """Write CSV tables of every available result into ``folder``."""
    os.makedirs(folder, exist_ok=True)
    written = []
    if geo is not None:
        p = os.path.join(folder, "sections.csv")
        _write_csv(p, ["section", "start_mm", "end_mm", "length_mm", "com_mm", "mass_kg", "J_kgm2"],
                   [[s.index + 1, (s.s_start - geo.s_min) * 1e3, (s.s_end - geo.s_min) * 1e3,
                     s.length * 1e3, (s.s_com - geo.s_min) * 1e3, s.mass, s.J] for s in geo.sections])
        written.append(p)
    if res is not None and geo is not None:
        p = os.path.join(folder, "torsion_fem.csv")
        rows = []
        for i, s in enumerate(geo.sections):
            rows.append([i + 1, (res.s_mid[i] - geo.s_min) * 1e3, np.degrees(res.theta_mid[i]),
                         res.k_section[i], res.k_spring[i] if i < len(res.k_spring) else ""])
        _write_csv(p, ["section", "lumped_position_mm", "theta_deg", "k_section_Nm_per_rad",
                       "k_to_next_Nm_per_rad"], rows)
        with open(p, "a", encoding="utf-8") as fh:
            fh.write(f"# torque;{res.settings.torque}\n# k_ground (clamp->inertia {res.ground_index + 1});"
                     f"{res.k_ground}\n# k_total;{res.k_total}\n")
        written.append(p)
    if model is not None:
        f, phi, zeta = model.modes()
        p = os.path.join(folder, "modes.csv")
        _write_csv(p, ["mode", "f_Hz", "zeta"] + [f"phi_{i + 1}" for i in range(model.n)],
                   [[j + 1, f[j], zeta[j]] + list(phi[:, j]) for j in range(len(f))])
        written.append(p)
    if sim is not None:
        p = os.path.join(folder, "time_response.csv")
        n = sim["theta"].shape[0]
        hdr = (["t_s"] + [f"theta_{i + 1}_rad" for i in range(n)] +
               [f"torque_{i + 1}_{i + 2}_Nm" for i in range(n - 1)] + [f"load_{i + 1}_Nm" for i in range(n)])
        data = np.vstack([sim["t"], sim["theta"], sim["spring_torque"], sim["load"]]).T
        _write_csv(p, hdr, data.tolist())
        written.append(p)
    return written


def run_pipeline(project, out_dir, progress=print):
    """Batch: split -> FEM -> modes -> time simulation -> CSV + PNG."""
    info = g.load_step(project["step_file"])
    ax = project["axis"]
    axis = make_axis(info, ax.get("mode"), ax.get("origin_mm"), ax.get("direction"))
    s_range = g.axial_range(info.surface[0], axis)
    cuts = [s_range[0] + c / 1000.0 for c in project["cuts_mm"]]
    geo = g.split_model(project["step_file"], axis, cuts, project["density"], s_range,
                        out_brep=os.path.join(out_dir, "split_model.brep"), progress=progress)
    res = run_torsion(geo, fem_settings(project["fem"]), progress=progress)
    if project.get("model"):
        model = model_from_dict(project["model"])
    else:
        model = model_from_results(geo, res)
    loads = [Load(**ld) for ld in project.get("loads", [])]
    sim = None
    if loads:
        dt = project["simulation"].get("dt") or model.recommended_dt()
        sim = model.simulate(loads, project["simulation"]["t_end"], dt)
    files = export_results(out_dir, geo, res, model, sim)
    return geo, res, model, sim, files
