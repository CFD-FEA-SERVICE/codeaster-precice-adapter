# codeaster-precice-adapter test: fake fluid for the CHT smoke test
# Copyright (C) 2026 CFD FEA SERVICE S.R.L.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Writes a fixed sink temperature (300 K) and heat transfer coefficient (100 W/m2/K)
on the plate top y = 0.25, x in [0, 1]; prints what the solid sends back."""
import numpy as np
import precice

p = precice.Participant("Fluid", "../precice-config.xml", 0, 1)
mesh = "Fluid-Interface"
x = np.linspace(0, 1, 51)
ids = p.set_mesh_vertices(mesh, np.column_stack([x, np.full_like(x, 0.25)]))
T_inf, h_f = 300.0, 100.0
if p.requires_initial_data():
    p.write_data(mesh, "Sink-Temperature-Fluid", ids, np.full(len(x), T_inf))
    p.write_data(mesh, "Heat-Transfer-Coefficient-Fluid", ids, np.full(len(x), h_f))
p.initialize()
while p.is_coupling_ongoing():
    dt = p.get_max_time_step_size()
    Ts = p.read_data(mesh, "Sink-Temperature-Solid", ids, dt)
    hs = p.read_data(mesh, "Heat-Transfer-Coefficient-Solid", ids, dt)
    p.write_data(mesh, "Sink-Temperature-Fluid", ids, np.full(len(x), T_inf))
    p.write_data(mesh, "Heat-Transfer-Coefficient-Fluid", ids, np.full(len(x), h_f))
    p.advance(dt)
print(f"FLUID RESULT: solid wall T mean {Ts.mean():.6f} (min {Ts.min():.6f}, max {Ts.max():.6f}), "
      f"solid HTC {hs.mean():.4g}")
p.finalize()
