# codeaster-precice-adapter: couple code_aster with other solvers through preCICE
# Copyright (C) 2026 CFD FEA SERVICE S.R.L.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Conjugate heat transfer with code_aster as the solid participant.

Robin-Robin coupling, as in the preCICE OpenFOAM CHT module:
  read  : Sink-Temperature-Fluid, Heat-Transfer-Coefficient-Fluid  -> ECHANGE load
  write : Sink-Temperature-Solid  = wall temperature
          Heat-Transfer-Coefficient-Solid = lambda / delta
All data live on the coupling nodes; the ECHANGE coefficients are FORMULE
objects reading nodal interpolants, so the load is built once.
"""

import os
import numpy as np

from .core import SolverInterface, CouplingLoop, collect_precice_outputs
from .interpolation import NodalInterpolant
from . import aster_tools as at


class ThermalState:
    __slots__ = ("result", "t")

    def __init__(self, result, t):
        self.result, self.t = result, t


class AsterCHTSolid(SolverInterface):
    """code_aster thermal participant.

    Parameters
    ----------
    conductivity : solid thermal conductivity (W/m/K) used for the written HTC
    delta : distance used for HTC = conductivity / delta (small, like a wall layer)
    steady : solve THER_LINEAIRE stationary at every coupling step
    initial_temperature : uniform initial temperature (transient runs)
    nonlinear : use THER_NON_LINE instead of THER_LINEAIRE
    """

    def __init__(self, mesh, model, material, interface_group, loads=(), conductivity=1.0,
                 delta=1e-5, steady=True, initial_temperature=293.15, nonlinear=False,
                 mesh_name="Solid-Interface-Nodes",
                 read_names=("Sink-Temperature-Fluid", "Heat-Transfer-Coefficient-Fluid"),
                 write_names=("Sink-Temperature-Solid", "Heat-Transfer-Coefficient-Solid"),
                 output_dir=None, output_every=0, ther_options=None, log=print):
        self.mesh, self.model, self.material = mesh, model, material
        self.group = interface_group
        self.node_group = (interface_group + "_PRCN")[:24]
        self.loads = list(loads)
        self.k, self.delta = float(conductivity), float(delta)
        self.steady, self.T0, self.nonlinear = steady, float(initial_temperature), nonlinear
        self.mesh_name = mesh_name
        self.read_data, self.write_data = tuple(read_names), tuple(write_names)
        self.output_dir, self.output_every = output_dir, int(output_every)
        self.ther_options = dict(ther_options or {})
        self.log = log

    def vertices(self, dimensions):
        mdim = at.mesh_dimension(self.mesh)
        if mdim != dimensions:
            raise ValueError(f"preCICE mesh is {dimensions}D but the code_aster mesh is {mdim}D")
        self.dim = dimensions
        self.mesh = at.ensure_node_group(self.mesh, self.group, self.node_group)
        self.keys, coords = at.node_coordinates(self.mesh, self.node_group)
        coords = coords[:, :dimensions]
        # nearest-node interpolation also at the Gauss points of the boundary faces
        self.t_sink = NodalInterpolant(coords, values=np.full(len(coords), self.T0))
        self.h_conv = NodalInterpolant(coords, values=np.zeros(len(coords)))
        self._formulas = (at.coordinate_formula(self.h_conv, dimensions),
                          at.coordinate_formula(self.t_sink, dimensions))
        self.robin = at.cmd("AFFE_CHAR_THER_F")(MODELE=self.model, ECHANGE=dict(
            GROUP_MA=self.group, COEF_H=self._formulas[0], TEMP_EXT=self._formulas[1]))
        self.log(f"[precice-aster] CHT: {len(self.keys)} coupling nodes on '{self.group}'")
        if self.output_dir:
            os.makedirs(self.output_dir, exist_ok=True)
        return coords

    def initial_state(self):
        return ThermalState(None, 0.0)

    def _write(self, temps):
        n = len(self.keys)
        return {self.write_data[0]: np.asarray(temps, float).reshape(n),
                self.write_data[1]: np.full(n, self.k / self.delta)}

    def initial_write_values(self, state):
        return self._write(np.full(len(self.keys), self.T0))

    def solve(self, state, t, dt, read_values):
        self.t_sink.update(np.asarray(read_values[self.read_data[0]], float).reshape(-1))
        self.h_conv.update(np.asarray(read_values[self.read_data[1]], float).reshape(-1))
        t1 = t + dt
        excit = [dict(CHARGE=l) for l in self.loads] + [dict(CHARGE=self.robin)]
        command = at.cmd("THER_NON_LINE" if self.nonlinear else "THER_LINEAIRE")
        kw = dict(MODELE=self.model, CHAM_MATER=self.material, EXCIT=excit)
        if self.steady:
            # code_aster 17: steady solve = TYPE_CALCUL='STAT', INCREMENT still mandatory:
            # a one-instant list (end of the window) gives exactly one steady solve
            kw["TYPE_CALCUL"] = "STAT"
            kw["INCREMENT"] = dict(LIST_INST=at.cmd("DEFI_LIST_REEL")(VALE=(t1,)))
        else:
            kw["TYPE_CALCUL"] = "TRAN"
            kw["INCREMENT"] = dict(LIST_INST=at.time_list(t, t1, 1))
            kw["PARM_THETA"] = 1.0
            if state.result is None:
                kw["ETAT_INIT"] = dict(VALE=self.T0)
            else:
                kw["ETAT_INIT"] = dict(EVOL_THER=state.result, INST=state.t,
                                       CRITERE="ABSOLU", PRECISION=max(dt * 1e-6, 1e-12))
        kw.update(self.ther_options)
        result = command(**kw)
        temps = at.nodal_result_values(result, "TEMP", ("TEMP",), self.node_group, self.keys,
                                       inst=None)[:, 0]
        return ThermalState(result, t1), self._write(temps)

    def on_time_window_complete(self, state, t, window):
        if window == 1 or window % 10 == 0:
            self.log(f"[precice-aster] window {window}: t = {t:g}, sink T fluid "
                     f"[{self.t_sink.values.min():.3f}, {self.t_sink.values.max():.3f}]")
        if self.output_dir and self.output_every and window % self.output_every == 0:
            at.write_med(state.result, os.path.join(self.output_dir, f"thermal_{window:05d}.rmed"),
                         fields=("TEMP",))


def run_cht(mesh, model, material, interface_group, loads=(), participant=None,
            config_file=None, **options):
    """Build the CHT solid participant and run the coupled simulation."""
    participant = participant or os.environ.get("PRECICE_PARTICIPANT", "Solid")
    config_file = config_file or os.environ.get("PRECICE_CONFIG_FILE")
    if not config_file:
        raise ValueError("give config_file= or set PRECICE_CONFIG_FILE (absolute path)")
    solver = AsterCHTSolid(mesh, model, material, interface_group, loads, **options)
    try:
        return CouplingLoop(solver, participant, config_file).run()
    finally:
        collect_precice_outputs()
