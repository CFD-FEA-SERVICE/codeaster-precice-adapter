# codeaster-precice-adapter test helper
# Copyright (C) 2026 CFD FEA SERVICE S.R.L.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Minimal stand-in for the code_aster commands used by precice_aster.

It checks keyword usage (unknown keywords raise), evaluates FORMULE objects the
way code_aster would (expression + NOM_PARA + Python context) and integrates the
same Euler-Bernoulli cantilever as tests/mock_beam.py, so a coupled run with the
stub must reproduce the mock-beam run.  It does NOT validate code_aster itself.
"""
import numpy as np
from mock_beam import MockBeam

CALLS = []


def _check(name, kw, allowed):
    bad = set(kw) - set(allowed)
    if bad:
        raise TypeError(f"{name}: unexpected keywords {sorted(bad)}")
    CALLS.append(name)


class Mesh:
    def __init__(self, coords, cell_groups, names, dim=2):
        self.coords, self.cell_groups, self.names, self.dim = coords, cell_groups, names, dim
        self.node_groups = {}

    def getDimension(self):
        return self.dim

    def hasGroupOfNodes(self, name):
        return name in self.node_groups


class Table:
    def __init__(self, d):
        self._d = d

    def EXTR_TABLE(self):
        d = self._d
        class _T:
            def values(self_inner):
                return d
        return _T()


class Formula:
    def __init__(self, VALE, NOM_PARA, **ctx):
        self.expr, self.params, self.ctx = VALE, tuple(NOM_PARA), ctx

    def __call__(self, *args):
        env = dict(self.ctx); env.update(dict(zip(self.params, args)))
        return eval(self.expr, {"__builtins__": {}}, env)


class Load:
    def __init__(self, kind, kw):
        self.kind, self.kw = kind, kw


class Result:
    def __init__(self):
        self.states = {}  # inst -> (u, v, a) beam dofs
        self.nodal = {}   # inst -> {node: (dx, dy)}


class Stub:
    def __init__(self, mesh):
        self.mesh = mesh
        self.beam = MockBeam(nel=40)
        self.beam.vertices(2)

    def __call__(self, name):
        return getattr(self, name)

    # -- commands ---------------------------------------------------------------------
    def DEFI_GROUP(self, **kw):
        _check("DEFI_GROUP", kw, ("reuse", "MAILLAGE", "CREA_GROUP_NO"))
        g = kw["CREA_GROUP_NO"]
        nodes = sorted({n for c in self.mesh.cell_groups[g["GROUP_MA"]] for n in c})
        self.mesh.node_groups[g["NOM"]] = nodes
        return self.mesh

    def CREA_CHAMP(self, **kw):
        _check("CREA_CHAMP", kw, ("TYPE_CHAM", "OPERATION", "MAILLAGE", "NOM_CHAM", "MODELE", "AFFE"))
        if kw["OPERATION"] == "AFFE":
            assert kw["TYPE_CHAM"] == "NOEU_DEPL_R" and set(kw["AFFE"]["VALE"]) == {0.0}
            return ("ZERO",)
        assert kw["TYPE_CHAM"] == "NOEU_GEOM_R" and kw["NOM_CHAM"] == "GEOMETRIE"
        return ("GEOM", kw["MAILLAGE"])

    def POST_RELEVE_T(self, **kw):
        _check("POST_RELEVE_T", kw, ("ACTION",))
        a = kw["ACTION"]
        _check("POST_RELEVE_T/ACTION", a, ("OPERATION", "INTITULE", "CHAM_GD", "RESULTAT", "NOM_CHAM",
                                            "GROUP_NO", "TOUT_CMP", "NOM_CMP", "INST", "CRITERE",
                                            "PRECISION", "TOUT_ORDRE"))
        nodes = self.mesh.node_groups[a["GROUP_NO"]]
        xyz = self.mesh.coords[nodes]
        d = {"NOEUD": [self.mesh.names[n] for n in nodes],
             "COOR_X": list(xyz[:, 0]), "COOR_Y": list(xyz[:, 1]), "COOR_Z": [0.0] * len(nodes)}
        if "CHAM_GD" in a:
            d.update(X=list(xyz[:, 0]), Y=list(xyz[:, 1]), Z=[0.0] * len(nodes))
            return Table(d)
        res = a["RESULTAT"]
        insts = sorted(res.nodal) if a.get("TOUT_ORDRE") == "OUI" else \
            [min(res.nodal, key=lambda t: abs(t - a["INST"]))]
        if "INST" in a:
            assert abs(insts[0] - a["INST"]) <= a["PRECISION"] * 1.0001, "instant not found"
        out = {k: [] for k in ["NOEUD", "INST", "COOR_X", "COOR_Y", "COOR_Z"] + list(a["NOM_CMP"])}
        for t in insts:
            for i, n in enumerate(nodes):
                out["NOEUD"].append(self.mesh.names[n]); out["INST"].append(t)
                out["COOR_X"].append(xyz[i, 0]); out["COOR_Y"].append(xyz[i, 1]); out["COOR_Z"].append(0.0)
                vals = res.nodal[t][n]
                for c, v in zip(a["NOM_CMP"], vals):
                    out[c].append(v)
        return Table(out)

    def FORMULE(self, **kw):
        CALLS.append("FORMULE")
        return Formula(**kw)

    def AFFE_CHAR_MECA_F(self, **kw):
        _check("AFFE_CHAR_MECA_F", kw, ("MODELE", "FORCE_NODALE"))
        _check("FORCE_NODALE", kw["FORCE_NODALE"], ("GROUP_NO", "FX", "FY", "FZ"))
        return Load("force", kw)

    def DEFI_LIST_REEL(self, **kw):
        _check("DEFI_LIST_REEL", kw, ("DEBUT", "INTERVALLE", "VALE"))
        if "VALE" in kw:
            return list(kw["VALE"])
        i = kw["INTERVALLE"]
        return list(np.linspace(kw["DEBUT"], i["JUSQU_A"], i["NOMBRE"] + 1))

    def DYNA_NON_LINE(self, **kw):
        _check("DYNA_NON_LINE", kw, ("MODELE", "CHAM_MATER", "EXCIT", "COMPORTEMENT", "INCREMENT",
                                     "SCHEMA_TEMPS", "ETAT_INIT", "SOLVEUR", "CONVERGENCE", "NEWTON"))
        st = kw["SCHEMA_TEMPS"]
        assert st["SCHEMA"] == "NEWMARK" and st["FORMULATION"] == "DEPLACEMENT"
        times = kw["INCREMENT"]["LIST_INST"]
        # nodal forces from the FORMULE-based load, evaluated at the group nodes
        load = [e["CHARGE"] for e in kw["EXCIT"] if isinstance(e["CHARGE"], Load)][0]
        fn = load.kw["FORCE_NODALE"]
        nodes = self.mesh.node_groups[fn["GROUP_NO"]]
        f = np.array([[fn["FX"](*self.mesh.coords[n, :2]), fn["FY"](*self.mesh.coords[n, :2])] for n in nodes])
        # map nodal forces to beam nodes (same layout as MockBeam.vertices)
        b = self.beam
        y = self.mesh.coords[nodes, 1]
        forces = np.zeros((2 * len(b.y), 2))
        for k, n in enumerate(nodes):
            j = int(np.argmin(np.abs(b.y - y[k])))
            side = 0 if self.mesh.coords[n, 0] < 0 else len(b.y)
            forces[side + j] += f[k]
        if "ETAT_INIT" in kw and "EVOL_NOLI" not in kw["ETAT_INIT"]:
            ei = kw["ETAT_INIT"]
            assert ei["DEPL"] == ei["VITE"] == ei["ACCE"] == ("ZERO",)
            s = b.initial_state()
        elif "ETAT_INIT" in kw:
            ei = kw["ETAT_INIT"]; prev = ei["EVOL_NOLI"]
            t0 = min(prev.states, key=lambda t: abs(t - ei["INST"]))
            assert abs(t0 - ei["INST"]) <= ei["PRECISION"] * 1.0001
            u, v, a = prev.states[t0]
            from mock_beam import BeamState
            s = BeamState(u, v, a, t0)
        else:
            s = b.initial_state()
        res = Result()
        res.states[times[0]] = (s.u, s.v, s.a)
        for t0, t1 in zip(times[:-1], times[1:]):
            s, _ = b.solve(s, t0, t1 - t0, {"Force": forces})
            res.states[t1] = (s.u, s.v, s.a)
            w = np.concatenate([[0.0], s.u[0::2]])
            res.nodal[t1] = {n: (float(np.interp(self.mesh.coords[n, 1], b.y, w)), 0.0) for n in nodes}
        return res


def flap_mesh():
    """Coupling nodes on both faces of the flap, at the mock-beam stations."""
    y = np.linspace(0, 1, 41)
    left = np.column_stack([np.full_like(y, -0.05), y])
    right = np.column_stack([np.full_like(y, 0.05), y])
    coords = np.vstack([left, right])
    names = [f"N{i+1}" for i in range(len(coords))]
    segs = [(i, i + 1) for i in range(40)] + [(41 + i, 42 + i) for i in range(40)]
    return Mesh(coords, {"flap": segs}, names, dim=2)
