# codeaster-precice-adapter test helper
# Copyright (C) 2026 CFD FEA SERVICE S.R.L.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Mock solid for testing the coupling core without code_aster.

Euler-Bernoulli cantilever (perpendicular-flap geometry: x in [-W/2, W/2],
y in [0, H], clamped at y = 0) integrated with the trapezoidal Newmark rule.
Coupling vertices sit on both faces of the flap; forces at the same height are
summed onto the beam node, and the beam deflection is returned for both faces.
"""
import numpy as np
from precice_aster.core import SolverInterface


class BeamState:
    __slots__ = ("u", "v", "a", "t")

    def __init__(self, u, v, a, t):
        self.u, self.v, self.a, self.t = u, v, a, t


class MockBeam(SolverInterface):
    mesh_name = "Solid-Mesh"
    read_data = ("Force",)
    write_data = ("Displacement",)

    def __init__(self, E=4e6, rho=3000.0, W=0.1, H=1.0, depth=1.0, nel=40):
        self.nel = nel
        I = depth * W ** 3 / 12.0
        A = depth * W
        L = H / nel
        n = nel + 1
        K = np.zeros((2 * n, 2 * n)); M = np.zeros_like(K)
        ke = E * I / L ** 3 * np.array([[12, 6 * L, -12, 6 * L], [6 * L, 4 * L * L, -6 * L, 2 * L * L],
                                        [-12, -6 * L, 12, -6 * L], [6 * L, 2 * L * L, -6 * L, 4 * L * L]])
        me = rho * A * L / 420 * np.array([[156, 22 * L, 54, -13 * L], [22 * L, 4 * L * L, 13 * L, -3 * L * L],
                                           [54, 13 * L, 156, -22 * L], [-13 * L, -3 * L * L, -22 * L, 4 * L * L]])
        for e in range(nel):
            d = slice(2 * e, 2 * e + 4)
            K[d, d] += ke; M[d, d] += me
        self.free = np.arange(2, 2 * n)  # clamp w and theta at y = 0
        self.K = K[np.ix_(self.free, self.free)]; self.M = M[np.ix_(self.free, self.free)]
        self.y = np.linspace(0, H, n)
        self.W, self.H = W, H
        self.solves = 0

    def vertices(self, dims):
        assert dims == 2
        left = np.column_stack([np.full_like(self.y, -self.W / 2), self.y])
        right = np.column_stack([np.full_like(self.y, self.W / 2), self.y])
        self.nv = 2 * len(self.y)
        return np.vstack([left, right])

    def initial_state(self):
        z = np.zeros(len(self.free))
        return BeamState(z, z.copy(), z.copy(), 0.0)

    def initial_write_values(self, state):
        return {"Displacement": np.zeros((self.nv, 2))}

    def _load(self, forces):
        n = len(self.y)
        fx = forces[:n, 0] + forces[n:, 0]  # both faces onto the beam node
        F = np.zeros(2 * n); F[0::2] = fx
        return F[self.free]

    def solve(self, s, t, dt, read):
        self.solves += 1
        F = self._load(np.asarray(read["Force"]))
        beta, gamma = 0.25, 0.5
        Keff = self.K + self.M / (beta * dt * dt)
        rhs = F + self.M @ (s.u / (beta * dt * dt) + s.v / (beta * dt) + (0.5 / beta - 1) * s.a)
        u = np.linalg.solve(Keff, rhs)
        a = (u - s.u) / (beta * dt * dt) - s.v / (beta * dt) - (0.5 / beta - 1) * s.a
        v = s.v + dt * ((1 - gamma) * s.a + gamma * a)
        w = np.concatenate([[0.0], u[0::2]])
        disp = np.zeros((self.nv, 2)); disp[: len(w), 0] = w; disp[len(w):, 0] = w
        return BeamState(u, v, a, t + dt), {"Displacement": disp}

    def tip(self, s):
        return float(s.u[-2])
