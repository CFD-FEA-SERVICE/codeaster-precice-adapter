# codeaster-precice-adapter: couple code_aster with other solvers through preCICE
# Copyright (C) 2026 CFD FEA SERVICE S.R.L.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Nodal interpolants evaluated by code_aster FORMULE objects.

code_aster loads such as ``FORCE_NODALE`` (AFFE_CHAR_MECA_F) or ``ECHANGE``
(AFFE_CHAR_THER_F) accept functions of the coordinates.  We build one FORMULE
per component that calls a :class:`NodalInterpolant`; preCICE data are pushed
into the interpolant with :meth:`NodalInterpolant.update` before every solve, so
the load objects are created once and only their values change.
"""

import numpy as np


class NodalInterpolant:
    """Value at the coupling vertex nearest to a point.

    Used at the coupling vertices themselves (nodal forces, nodal Robin data) the
    lookup is exact; elsewhere (Gauss points of boundary faces) it is a
    nearest-vertex interpolation.  ``tol`` (absolute distance) turns the lookup
    strict: a query farther than ``tol`` from every vertex returns ``default``.
    """

    def __init__(self, coords, values=None, tol=None, default=0.0):
        self.coords = np.asarray(coords, dtype=float)
        if self.coords.ndim != 2:
            raise ValueError("coords must be (n, dim)")
        self.n, self.dim = self.coords.shape
        self.values = np.zeros(self.n) if values is None else np.asarray(values, float).copy()
        self.tol = tol
        self.default = float(default)
        self.calls = 0
        self.misses = 0
        # simple uniform-grid bucket index for O(1)-ish lookups: code_aster calls the
        # interpolant once per node and per component at every load evaluation, so a
        # brute-force search would cost O(n) per call
        span = self.coords.max(0) - self.coords.min(0)
        span[span == 0] = 1.0
        self._h = float(np.max(span)) / max(int(np.sqrt(self.n)), 1)
        self._origin = self.coords.min(0)
        self._buckets = {}
        for i, key in enumerate(self._keys(self.coords)):
            self._buckets.setdefault(key, []).append(i)

    def _keys(self, pts):
        return [tuple(k) for k in np.floor((pts - self._origin) / self._h).astype(int)]

    def update(self, values):
        """Replace the vertex values (same order as ``coords``); geometry is unchanged."""
        values = np.asarray(values, dtype=float)
        if values.shape != (self.n,):
            raise ValueError(f"expected {self.n} values, got {values.shape}")
        self.values[:] = values

    def nearest(self, point):
        """(index, distance) of the coupling vertex nearest to ``point``."""
        p = np.asarray(point, dtype=float)[: self.dim]
        key = tuple(np.floor((p - self._origin) / self._h).astype(int))
        best, best_d = -1, np.inf
        # search the bucket and its neighbours; fall back to brute force
        rng = range(-1, 2)
        cand = []
        if self.dim == 2:
            for a in rng:
                for b in rng:
                    cand += self._buckets.get((key[0] + a, key[1] + b), [])
        else:
            for a in rng:
                for b in rng:
                    for c in rng:
                        cand += self._buckets.get((key[0] + a, key[1] + b, key[2] + c), [])
        idx = np.array(cand, dtype=int) if cand else np.arange(self.n)
        d = np.linalg.norm(self.coords[idx] - p, axis=1)
        j = int(np.argmin(d))
        best, best_d = int(idx[j]), float(d[j])
        if best_d > self._h and cand:  # point outside the searched neighbourhood
            d = np.linalg.norm(self.coords - p, axis=1)
            best = int(np.argmin(d)); best_d = float(d[best])
        return best, best_d

    def __call__(self, *xyz):
        self.calls += 1
        i, d = self.nearest(xyz)
        if self.tol is not None and d > self.tol:
            self.misses += 1
            return self.default
        return float(self.values[i])
