# codeaster-precice-adapter: couple code_aster with other solvers through preCICE
# Copyright (C) 2026 CFD FEA SERVICE S.R.L.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Fluid-structure interaction with code_aster as the solid participant.

Per time window:
  1. preCICE ``Force`` (nodal forces on the coupling nodes) is pushed into
     interpolants read by FORMULE-based ``FORCE_NODALE`` loads (AFFE_CHAR_MECA_F);
  2. one DYNA_NON_LINE step (Newmark) is solved from the last accepted state;
  3. ``Displacement`` (total DEPL) on the coupling nodes is sent back.

Implicit coupling: every solve starts from an immutable state (result concept +
time); the core loop keeps or discards the new state, so a preCICE checkpoint is
just a reference to the state at the start of the window.
"""

import os
import numpy as np

from .core import SolverInterface, CouplingLoop, collect_precice_outputs
from .interpolation import NodalInterpolant
from . import aster_tools as at


class SolidState:
    """Immutable solid state: a DYNA_NON_LINE result and the instant to restart from."""
    __slots__ = ("result", "t")

    def __init__(self, result, t):
        self.result, self.t = result, t


class AsterFSISolid(SolverInterface):
    """code_aster structural participant (FORCE in, DISPLACEMENT out).

    Parameters
    ----------
    mesh, model, material : code_aster mesh, model and CHAM_MATER concepts
    interface_group : cell group (GROUP_MA) of the wetted surface
    loads : list of code_aster loads always applied (supports, gravity, ...)
    mesh_name, force_name, displacement_name : preCICE names
    deformation : COMPORTEMENT DEFORMATION ('PETIT' linear, 'GROT_GDEP' large rotations)
    newmark : (beta, gamma) of the Newmark scheme; (0.25, 0.5) = trapezoidal rule
    substeps : code_aster steps per coupling step
    output_dir, output_every : MED snapshots of the accepted state every N windows
    dyna_options : extra DYNA_NON_LINE keywords (SOLVEUR, CONVERGENCE, NEWTON, ...)
    """

    def __init__(self, mesh, model, material, interface_group, loads=(),
                 mesh_name="Solid-Mesh", force_name="Force", displacement_name="Displacement",
                 deformation="PETIT", relation="ELAS", newmark=(0.25, 0.5), substeps=1,
                 output_dir=None, output_every=0, dyna_options=None, log=print):
        self.mesh, self.model, self.material = mesh, model, material
        self.group = interface_group
        self.node_group = (interface_group + "_PRCN")[:24]
        self.loads = list(loads)
        self.mesh_name = mesh_name
        self.read_data, self.write_data = (force_name,), (displacement_name,)
        self.force_name, self.disp_name = force_name, displacement_name
        self.deformation, self.relation = deformation, relation
        self.beta, self.gamma = newmark
        self.substeps = int(substeps)
        self.output_dir, self.output_every = output_dir, int(output_every)
        self.dyna_options = dict(dyna_options or {})
        self.log = log
        self.dim = None
        self.keys = None
        self.comps = None

    # -- set-up -------------------------------------------------------------------------
    def vertices(self, dimensions):
        mdim = at.mesh_dimension(self.mesh)
        if mdim != dimensions:
            raise ValueError(f"preCICE mesh {self.mesh_name} is {dimensions}D but the code_aster "
                             f"mesh is {mdim}D: use a {dimensions}D model (e.g. D_PLAN for 2D)")
        self.dim = dimensions
        self.comps = ("DX", "DY", "DZ")[:dimensions]
        self.mesh = at.ensure_node_group(self.mesh, self.group, self.node_group)
        self.keys, coords = at.node_coordinates(self.mesh, self.node_group)
        coords = coords[:, :dimensions]
        # one interpolant per force component, exact at the coupling nodes
        hmin = _min_spacing(coords)
        self.force_interp = [NodalInterpolant(coords, tol=0.25 * hmin) for _ in range(dimensions)]
        cmd_ = at.cmd("AFFE_CHAR_MECA_F")
        comp_kw = ("FX", "FY", "FZ")[:dimensions]
        self._formulas = [at.coordinate_formula(f, dimensions) for f in self.force_interp]
        self.coupling_load = cmd_(MODELE=self.model, FORCE_NODALE=dict(
            GROUP_NO=self.node_group, **{k: f for k, f in zip(comp_kw, self._formulas)}))
        self.coords = coords
        self.log(f"[precice-aster] FSI: {len(self.keys)} coupling nodes on '{self.group}', "
                 f"min spacing {hmin:.3g}")
        if self.output_dir:
            os.makedirs(self.output_dir, exist_ok=True)
        return coords

    def initial_state(self):
        return SolidState(None, 0.0)

    def initial_write_values(self, state):
        return {self.disp_name: np.zeros((len(self.keys), self.dim))}

    # -- one coupling step ----------------------------------------------------------------
    def solve(self, state, t, dt, read_values):
        forces = np.asarray(read_values[self.force_name], dtype=float).reshape(-1, self.dim)
        for c in range(self.dim):
            self.force_interp[c].update(forces[:, c])

        t1 = t + dt
        kw = dict(
            MODELE=self.model, CHAM_MATER=self.material,
            EXCIT=[dict(CHARGE=l) for l in self.loads] + [dict(CHARGE=self.coupling_load)],
            COMPORTEMENT=dict(RELATION=self.relation, DEFORMATION=self.deformation, TOUT="OUI"),
            INCREMENT=dict(LIST_INST=at.time_list(t, t1, self.substeps), INST_INIT=t, INST_FIN=t1),
            SCHEMA_TEMPS=dict(SCHEMA="NEWMARK", FORMULATION="DEPLACEMENT",
                              BETA=self.beta, GAMMA=self.gamma),
        )
        if state.result is not None:
            kw["ETAT_INIT"] = dict(EVOL_NOLI=state.result, INST=state.t,
                                   CRITERE="ABSOLU", PRECISION=_eps(dt))
        kw.update(self.dyna_options)
        misses = sum(fi.misses for fi in self.force_interp)
        result = at.cmd("DYNA_NON_LINE")(**kw)
        if sum(fi.misses for fi in self.force_interp) > misses:
            raise RuntimeError(
                "code_aster evaluated the coupling force away from the coupling nodes "
                "(FORMULE called at points farther than 1/4 of the node spacing): the "
                "force would be lost. Check that the interface group is the loaded group.")

        disp = at.nodal_result_values(result, "DEPL", self.comps, self.node_group, self.keys,
                                      inst=t1, precision=_eps(dt))
        return SolidState(result, t1), {self.disp_name: disp}

    # -- monitoring -----------------------------------------------------------------------
    def on_time_window_complete(self, state, t, window):
        f = self.force_interp
        ftot = [float(np.sum(fi.values)) for fi in f]
        if window == 1 or window % 10 == 0:
            self.log(f"[precice-aster] window {window}: t = {t:g}, total coupling force "
                     + ", ".join(f"{c}={v:.4g}" for c, v in zip(("FX", "FY", "FZ"), ftot)))
        if self.output_dir and self.output_every and window % self.output_every == 0:
            at.write_med(state.result, os.path.join(self.output_dir, f"solid_{window:05d}.rmed"),
                         inst=t, fields=("DEPL",))


def _min_spacing(coords):
    n = len(coords)
    if n < 2:
        return 1.0
    sample = coords if n <= 2000 else coords[np.random.default_rng(0).choice(n, 2000, replace=False)]
    best = np.inf
    for i in range(len(sample)):
        d = np.linalg.norm(coords - sample[i], axis=1)
        d = d[d > 0]
        if d.size:
            best = min(best, d.min())
    return float(best) if np.isfinite(best) else 1.0


def _eps(dt):
    return max(abs(dt) * 1e-6, 1e-12)


def run_fsi(mesh, model, material, interface_group, loads=(), participant=None,
            config_file=None, **options):
    """Build the FSI solid participant and run the coupled simulation.

    ``participant`` defaults to $PRECICE_PARTICIPANT (or 'Solid');
    ``config_file`` defaults to $PRECICE_CONFIG_FILE (use an absolute path:
    code_aster runs in a temporary directory).
    """
    participant = participant or os.environ.get("PRECICE_PARTICIPANT", "Solid")
    config_file = config_file or os.environ.get("PRECICE_CONFIG_FILE")
    if not config_file:
        raise ValueError("give config_file= or set PRECICE_CONFIG_FILE (absolute path)")
    solver = AsterFSISolid(mesh, model, material, interface_group, loads, **options)
    try:
        return CouplingLoop(solver, participant, config_file).run()
    finally:
        collect_precice_outputs()
