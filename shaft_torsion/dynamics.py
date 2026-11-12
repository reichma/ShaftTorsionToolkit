"""Step 3 - lumped-parameter torsional vibration model.

    J * theta'' + C * theta' + K * theta = M(t)

* ``J``   inertias of the sections (step 1)
* ``k``   springs between neighbouring inertias (step 2)
* ``c``   dampers between neighbouring inertias (user)
* ``c_abs`` damping of each inertia to ground (bearings, user)
* optional spring/damper from one inertia to ground (the FEM clamp)
* loads ``M_i(t)``: constant, step, ramp, sine, sweep, pulse or a formula
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict

import numpy as np
import scipy.linalg as sla

LOAD_KINDS = ["constant", "step", "ramp", "sine", "sweep", "pulse", "expression"]


@dataclass
class Load:
    """External torque on inertia ``node`` [N m] as a function of time.

    constant   : M = A
    step       : M = A for t >= t_start (and t < t_end if t_end > t_start)
    ramp       : M rises linearly from 0 at t_start to A at t_end, then stays
    sine       : M = offset + A sin(2 pi f (t - t_start) + phase), active from t_start
    sweep      : linear chirp from f to f_end between t_start and t_end
    pulse      : M = A between t_start and t_end, 0 otherwise
    expression : python/numpy formula in t, e.g. "500*sin(2*pi*30*t)*exp(-t)"
    """
    node: int = 0
    kind: str = "step"
    amplitude: float = 100.0
    frequency: float = 10.0       # [Hz]
    f_end: float = 100.0          # [Hz] (sweep)
    phase_deg: float = 0.0
    offset: float = 0.0
    t_start: float = 0.0
    t_end: float = 0.0
    expression: str = ""

    def __call__(self, t):
        t = np.asarray(t, float)
        A, k = self.amplitude, self.kind
        on = t >= self.t_start
        if k == "constant":
            return np.full_like(t, A)
        if k == "step":
            m = on if self.t_end <= self.t_start else on & (t < self.t_end)
            return np.where(m, A, 0.0)
        if k == "ramp":
            dur = max(self.t_end - self.t_start, 1e-12)
            return A * np.clip((t - self.t_start) / dur, 0.0, 1.0)
        if k == "sine":
            ph = np.radians(self.phase_deg)
            return np.where(on, self.offset + A * np.sin(2 * np.pi * self.frequency * (t - self.t_start) + ph), 0.0)
        if k == "sweep":
            dur = max(self.t_end - self.t_start, 1e-12)
            tau = np.clip(t - self.t_start, 0, dur)
            rate = (self.f_end - self.frequency) / dur
            arg = 2 * np.pi * (self.frequency * tau + 0.5 * rate * tau ** 2) + np.radians(self.phase_deg)
            return np.where(on & (t <= self.t_end), A * np.sin(arg), 0.0)
        if k == "pulse":
            return np.where(on & (t < self.t_end), A, 0.0)
        if k == "expression":
            ns = {n: getattr(np, n) for n in ("sin", "cos", "tan", "exp", "log", "sqrt", "abs",
                                               "pi", "where", "clip", "minimum", "maximum",
                                               "sign", "heaviside", "tanh", "arctan")}
            ns["t"] = t
            val = eval(self.expression or "0", {"__builtins__": {}}, ns)
            return np.broadcast_to(np.asarray(val, float), t.shape).copy()
        raise ValueError(f"unknown load kind {k!r}")

    def to_dict(self):
        return asdict(self)


@dataclass
class TorsionalModel:
    J: np.ndarray                         # (n,) [kg m^2]
    k: np.ndarray                         # (n-1,) [N m/rad]
    c: np.ndarray = None                  # (n-1,) [N m s/rad]
    c_abs: np.ndarray = None              # (n,) damping to ground [N m s/rad]
    ground_node: int = -1                 # inertia connected to ground; -1 = free-free
    k_ground: float = 0.0
    c_ground: float = 0.0
    names: list = field(default_factory=list)

    def __post_init__(self):
        self.J = np.asarray(self.J, float)
        n = len(self.J)
        self.k = np.asarray(self.k, float)
        self.c = np.zeros(n - 1) if self.c is None else np.asarray(self.c, float)
        self.c_abs = np.zeros(n) if self.c_abs is None else np.asarray(self.c_abs, float)
        if len(self.k) != n - 1 or len(self.c) != n - 1 or len(self.c_abs) != n:
            raise ValueError("inconsistent model sizes")
        if np.any(self.J <= 0):
            raise ValueError("all inertias must be > 0")

    @property
    def n(self):
        return len(self.J)

    def matrices(self):
        n = self.n
        M = np.diag(self.J)
        K = np.zeros((n, n))
        C = np.diag(self.c_abs.copy())
        for i in range(n - 1):
            for mat, val in ((K, self.k[i]), (C, self.c[i])):
                mat[i, i] += val
                mat[i + 1, i + 1] += val
                mat[i, i + 1] -= val
                mat[i + 1, i] -= val
        if 0 <= self.ground_node < n:
            K[self.ground_node, self.ground_node] += self.k_ground
            C[self.ground_node, self.ground_node] += self.c_ground
        return M, C, K

    # ------------------------------------------------------------------ #
    def modes(self):
        """Undamped natural frequencies [Hz] and mode shapes (columns,
        normalised to max |phi| = 1), plus modal damping ratios estimated
        from the modal damping matrix (diagonal terms)."""
        M, C, K = self.matrices()
        w2, phi = sla.eigh(K, M)
        w2 = np.clip(w2, 0, None)
        w = np.sqrt(w2)
        phi = phi / np.abs(phi).max(0)
        for j in range(phi.shape[1]):
            if phi[np.argmax(abs(phi[:, j])), j] < 0:
                phi[:, j] *= -1
        mm = np.einsum("ij,ik,kj->j", phi, M, phi)
        cm = np.einsum("ij,ik,kj->j", phi, C, phi)
        with np.errstate(divide="ignore", invalid="ignore"):
            zeta = np.where(w > 1e-9, cm / (2 * mm * w), 0.0)
        return w / (2 * np.pi), phi, zeta

    def complex_modes(self):
        """Damped eigenvalues (state-space). Returns frequency [Hz] and damping ratio."""
        M, C, K = self.matrices()
        n = self.n
        Minv = np.linalg.inv(M)
        A = np.block([[np.zeros((n, n)), np.eye(n)], [-Minv @ K, -Minv @ C]])
        lam = np.linalg.eigvals(A)
        lam = lam[np.imag(lam) >= -1e-9]
        lam = lam[np.argsort(np.abs(lam))]
        wd = np.abs(lam)
        with np.errstate(divide="ignore", invalid="ignore"):
            zeta = np.where(wd > 1e-9, -np.real(lam) / wd, 0.0)
        return wd / (2 * np.pi), zeta

    # ------------------------------------------------------------------ #
    def simulate(self, loads, t_end, dt, theta0=None, omega0=None, beta=0.25, gamma=0.5):
        """Newmark time integration (average acceleration: unconditionally
        stable, no numerical damping).  Returns a dict with t, theta, omega,
        alpha (n, nt), spring torques (n-1, nt) and applied loads (n, nt)."""
        M, C, K = self.matrices()
        n = self.n
        nt = int(round(t_end / dt)) + 1
        t = np.linspace(0, dt * (nt - 1), nt)
        F = np.zeros((n, nt))
        for ld in loads:
            if not 0 <= ld.node < n:
                raise ValueError(f"load on inertia {ld.node + 1} does not exist")
            F[ld.node] += ld(t)
        th = np.zeros((n, nt))
        om = np.zeros((n, nt))
        al = np.zeros((n, nt))
        th[:, 0] = 0 if theta0 is None else theta0
        om[:, 0] = 0 if omega0 is None else omega0
        al[:, 0] = np.linalg.solve(M, F[:, 0] - C @ om[:, 0] - K @ th[:, 0])
        Keff = M / (beta * dt ** 2) + C * gamma / (beta * dt) + K
        lu = sla.lu_factor(Keff)
        a0, a1, a2 = 1 / (beta * dt ** 2), 1 / (beta * dt), 1 / (2 * beta) - 1
        b0, b1, b2 = gamma / (beta * dt), gamma / beta - 1, dt * (gamma / (2 * beta) - 1)
        for i in range(nt - 1):
            x, v, a = th[:, i], om[:, i], al[:, i]
            rhs = (F[:, i + 1] + M @ (a0 * x + a1 * v + a2 * a) + C @ (b0 * x + b1 * v + b2 * a))
            xn = sla.lu_solve(lu, rhs)
            an = a0 * (xn - x) - a1 * v - a2 * a
            vn = v + dt * ((1 - gamma) * a + gamma * an)
            th[:, i + 1], om[:, i + 1], al[:, i + 1] = xn, vn, an
        twist = np.diff(th, axis=0)
        dtw = np.diff(om, axis=0)
        torque = self.k[:, None] * twist + self.c[:, None] * dtw
        out = {"t": t, "theta": th, "omega": om, "alpha": al, "twist": twist,
               "spring_torque": torque, "load": F}
        if 0 <= self.ground_node < n:
            g = self.ground_node
            out["ground_torque"] = self.k_ground * th[g] + self.c_ground * om[g]
        return out

    def frf(self, node, freqs):
        """Rotation amplitude of every inertia [rad per N m] for a harmonic unit
        torque on inertia ``node``.  Returns complex array (n, nf)."""
        M, C, K = self.matrices()
        f = np.zeros(self.n)
        f[node] = 1.0
        out = np.empty((self.n, len(freqs)), complex)
        for j, fr in enumerate(freqs):
            w = 2 * np.pi * fr
            A = K - w ** 2 * M + 1j * w * C
            try:
                out[:, j] = np.linalg.solve(A, f)
            except np.linalg.LinAlgError:
                out[:, j] = np.nan
        return out

    def recommended_dt(self, points_per_period=20):
        f, _, _ = self.modes()
        fmax = max(f.max(), 1e-6)
        return 1.0 / (points_per_period * fmax)

    def to_dict(self):
        return {"J": self.J.tolist(), "k": self.k.tolist(), "c": self.c.tolist(),
                "c_abs": self.c_abs.tolist(), "ground_node": self.ground_node,
                "k_ground": self.k_ground, "c_ground": self.c_ground}
