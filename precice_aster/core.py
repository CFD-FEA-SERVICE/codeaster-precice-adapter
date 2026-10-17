# codeaster-precice-adapter: couple code_aster with other solvers through preCICE
# Copyright (C) 2026 CFD FEA SERVICE S.R.L.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Solver-agnostic preCICE coupling loop.

The loop owns the preCICE participant and the checkpoint logic; a solver only
has to implement :class:`SolverInterface`.  The solver never mutates its
"committed" state in place: ``solve`` returns a new state, and the loop decides
whether to keep it (time step accepted) or to go back to the checkpoint
(implicit coupling sub-iteration).  This keeps the code_aster side free of any
preCICE bookkeeping.
"""

import glob
import math
import os
import shutil
import time as _time


class SolverInterface:
    """What a coupled solver must provide.

    Subclasses implement these methods; states are opaque objects for the loop.
    """

    #: name of the preCICE mesh the solver provides
    mesh_name = None
    #: preCICE data names read by the solver (list of str)
    read_data = ()
    #: preCICE data names written by the solver (list of str)
    write_data = ()
    #: largest time step the solver accepts (None = take preCICE's)
    max_dt = None

    def vertices(self, dimensions):
        """Coordinates of the coupling vertices, array (n, dimensions)."""
        raise NotImplementedError

    def initial_state(self):
        """State at t = 0."""
        raise NotImplementedError

    def initial_write_values(self, state):
        """dict name -> array, sent if preCICE requires initial data."""
        raise NotImplementedError

    def solve(self, state, t, dt, read_values):
        """Advance from ``state`` at time ``t`` by ``dt`` with the coupling data
        ``read_values`` (dict name -> array).  Return ``(new_state, write_values)``.
        Must not modify ``state``."""
        raise NotImplementedError

    def on_time_window_complete(self, state, t, window):
        """Called once a time window is accepted (output, monitoring)."""

    def finalize(self, state, t):
        """Called at the end of the coupled run."""


class CouplingLoop:
    """Drive a :class:`SolverInterface` with a preCICE participant."""

    def __init__(self, solver, participant_name, config_file, rank=0, size=1,
                 precice_module=None, log=print):
        if precice_module is None:
            import precice as precice_module  # pyprecice
        self.precice = precice_module
        self.solver = solver
        self.participant_name = participant_name
        self.config_file = config_file
        self.rank, self.size = rank, size
        self.log = log
        self.stats = {"windows": 0, "solves": 0, "restores": 0, "wall": 0.0}

    def run(self):
        s = self.solver
        p = self.precice.Participant(self.participant_name, self.config_file,
                                     self.rank, self.size)
        dims = p.get_mesh_dimensions(s.mesh_name)
        coords = s.vertices(dims)
        if coords.ndim != 2 or coords.shape[1] != dims:
            raise ValueError(f"vertices must be (n, {dims}), got {coords.shape}")
        ids = p.set_mesh_vertices(s.mesh_name, coords)
        self.log(f"[precice-aster] participant {self.participant_name}: "
                 f"{len(ids)} vertices on mesh {s.mesh_name} ({dims}D)")

        state = s.initial_state()
        if p.requires_initial_data():
            for name, values in s.initial_write_values(state).items():
                p.write_data(s.mesh_name, name, ids, values)
        p.initialize()

        t, window, checkpoint = 0.0, 0, None
        t0_wall = _time.time()
        while p.is_coupling_ongoing():
            if p.requires_writing_checkpoint():
                checkpoint = (state, t)

            dt = p.get_max_time_step_size()
            if s.max_dt is not None:
                dt = min(dt, s.max_dt)

            read_values = {name: p.read_data(s.mesh_name, name, ids, dt)
                           for name in s.read_data}
            new_state, write_values = s.solve(state, t, dt, read_values)
            self.stats["solves"] += 1
            for name in s.write_data:
                p.write_data(s.mesh_name, name, ids, write_values[name])

            p.advance(dt)

            if p.requires_reading_checkpoint():
                if checkpoint is None:
                    raise RuntimeError("preCICE asked to read a checkpoint that was never written")
                state, t = checkpoint
                self.stats["restores"] += 1
            else:
                state, t = new_state, t + dt
                if p.is_time_window_complete():
                    window += 1
                    self.stats["windows"] = window
                    s.on_time_window_complete(state, t, window)

        self.stats["wall"] = _time.time() - t0_wall
        s.finalize(state, t)
        p.finalize()
        st = self.stats
        avg = st["solves"] / max(st["windows"], 1)
        self.log(f"[precice-aster] done: t = {t:g}, {st['windows']} windows, "
                 f"{st['solves']} solves ({avg:.2f} per window), "
                 f"{st['restores']} checkpoint restores, {st['wall']:.1f} s")
        return state, t


def nearly(a, b, rel=1e-9):
    """Float comparison used for time instants."""
    return math.isclose(a, b, rel_tol=rel, abs_tol=rel)


def collect_precice_outputs(dest=None, src="."):
    """Copy preCICE output files (watch points, iteration/convergence logs, profiling,
    exports) from ``src`` to ``dest`` (default $PRECICE_CASE_DIR).

    code_aster runs in a temporary directory that is deleted at the end of the job,
    and preCICE writes these files in the current directory.
    """
    dest = dest or os.environ.get("PRECICE_CASE_DIR")
    if not dest or os.path.abspath(dest) == os.path.abspath(src):
        return []
    os.makedirs(dest, exist_ok=True)
    copied = []
    for path in glob.glob(os.path.join(src, "precice-*")):
        target = os.path.join(dest, os.path.basename(path))
        if os.path.isdir(path):
            shutil.copytree(path, target, dirs_exist_ok=True)
        else:
            shutil.copy2(path, target)
        copied.append(target)
    return copied
