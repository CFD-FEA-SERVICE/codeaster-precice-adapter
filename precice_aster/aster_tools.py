# codeaster-precice-adapter: couple code_aster with other solvers through preCICE
# Copyright (C) 2026 CFD FEA SERVICE S.R.L.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Thin helpers around code_aster commands.

Only long-standing, documented commands are used (CREA_CHAMP, DEFI_GROUP,
POST_RELEVE_T, FORMULE, DEFI_LIST_REEL, DEFI_FICHIER, IMPR_RESU).  Nodal values
are read through POST_RELEVE_T tables and matched to the coupling vertices by
node name (or by coordinates), so the adapter does not depend on the internal
node numbering of the Python API.
"""

import numpy as np

# Commands are looked up lazily so the module can be imported (and tested with
# a stub) without code_aster.  ``set_command_provider`` lets tests inject one.
_provider = None


def set_command_provider(provider):
    """provider(name) -> code_aster command object (tests only)."""
    global _provider
    _provider = provider


def _to_factor(value, F):
    """Plain dicts (and lists/tuples of dicts) -> code_aster _F factor keywords.

    Some commands (e.g. DYNA_NON_LINE, STAT_NON_LINE) inspect factor keywords
    as ``keywords["INCREMENT"][0]`` before the syntax is normalised, which only
    works with _F objects, not with plain dicts.
    """
    if isinstance(value, dict) and not isinstance(value, F):
        return F(**{k: _to_factor(v, F) for k, v in value.items()})
    if isinstance(value, (list, tuple)) and value and all(isinstance(v, dict) for v in value):
        return tuple(_to_factor(v, F) for v in value)
    return value


class _LocalF(dict):
    """Fallback with the behaviour of code_aster's _F: ``f[0] is f``."""

    def __getitem__(self, key):
        if isinstance(key, int):
            if key == 0:
                return self
            raise IndexError(key)
        return dict.__getitem__(self, key)


_F_CLASS = None


def _factor_class(commands_module):
    """code_aster's _F class: its location changed across versions."""
    global _F_CLASS
    if _F_CLASS is None:
        import builtins
        import importlib
        candidates = [getattr(commands_module, "_F", None), getattr(builtins, "_F", None)]
        for mod in ("code_aster.Cata.Syntax", "code_aster.Cata.Language.Syntax",
                    "code_aster.Supervis", "code_aster.Utilities", "code_aster.CA"):
            try:
                candidates.append(getattr(importlib.import_module(mod), "_F", None))
            except Exception:  # noqa: BLE001  (module absent in this version)
                pass
        # first candidate that is a dict subclass; otherwise our own equivalent, which
        # is enough for what code_aster needs (factor keyword access as kw[0])
        found = next((c for c in candidates if isinstance(c, type) and issubclass(c, dict)), None)
        _F_CLASS = found or _LocalF
    return _F_CLASS


def cmd(name):
    """code_aster command ``name``, wrapped so that plain dicts work as factor keywords.

    The adapter writes ``INCREMENT=dict(...)``, ``EXCIT=[dict(...), ...]``; the wrapper
    converts them to code_aster _F objects before the call.
    """
    if _provider is not None:
        return _provider(name)
    import code_aster.Commands as C  # noqa: N812  (code_aster >= 15)
    command, F = getattr(C, name), _factor_class(C)

    def call(**keywords):
        return command(**{k: _to_factor(v, F) for k, v in keywords.items()})

    call.__name__ = name
    return call


COORD_COLS = ("COOR_X", "COOR_Y", "COOR_Z")


def table_to_dict(table):
    """Return {column: list} for a code_aster table concept."""
    if hasattr(table, "EXTR_TABLE"):
        return table.EXTR_TABLE().values()
    if hasattr(table, "getValues"):  # newer API, kept as fallback
        return table.getValues()
    raise TypeError("cannot read code_aster table: no EXTR_TABLE/getValues")


def mesh_dimension(mesh):
    """Mesh dimension (2 for a plane mesh such as the .mail examples, 3 otherwise)."""
    try:
        return int(mesh.getDimension())
    except Exception:  # pragma: no cover - very old API
        return 3


def ensure_node_group(mesh, cell_group, node_group):
    """Create node group ``node_group`` from the nodes of cell group ``cell_group``."""
    # reuse the group if it exists already (e.g. a second participant on the same mesh)
    try:
        if hasattr(mesh, "hasGroupOfNodes") and mesh.hasGroupOfNodes(node_group):
            return mesh
    except Exception:
        pass
    DEFI_GROUP = cmd("DEFI_GROUP")
    return DEFI_GROUP(reuse=mesh, MAILLAGE=mesh,
                      CREA_GROUP_NO=dict(NOM=node_group, GROUP_MA=cell_group))


def _rows(tab_dict, value_cols):
    """List of (key, coords, values) for every row of an extraction table.

    The key identifies the node: its name (NOEUD column) when the table has it,
    otherwise its rounded coordinates.  Keys are how values are matched to the
    coupling vertices, independently of code_aster's internal node numbering.
    """
    n = len(tab_dict[value_cols[0]])
    have_node = "NOEUD" in tab_dict and tab_dict["NOEUD"] and tab_dict["NOEUD"][0] is not None
    coords = np.array([[float(tab_dict[c][i]) if c in tab_dict and tab_dict[c][i] is not None else 0.0
                        for c in COORD_COLS] for i in range(n)])
    keys = [str(tab_dict["NOEUD"][i]).strip() for i in range(n)] if have_node else \
           [tuple(np.round(coords[i], 9)) for i in range(n)]
    vals = np.array([[float(tab_dict[c][i]) for c in value_cols] for i in range(n)])
    return keys, coords, vals


def node_coordinates(mesh, node_group):
    """(keys, coords (n,3)) of the nodes of ``node_group``, in table order."""
    CREA_CHAMP, POST_RELEVE_T = cmd("CREA_CHAMP"), cmd("POST_RELEVE_T")
    # geometry as a nodal field, then extracted on the group like any result field
    geom = CREA_CHAMP(TYPE_CHAM="NOEU_GEOM_R", OPERATION="EXTR",
                      MAILLAGE=mesh, NOM_CHAM="GEOMETRIE")
    tab = POST_RELEVE_T(ACTION=dict(OPERATION="EXTRACTION", INTITULE="PRECICE_GEOM",
                                    CHAM_GD=geom, GROUP_NO=node_group, TOUT_CMP="OUI"))
    d = table_to_dict(tab)
    col = "COOR_X" if "COOR_X" in d else "X"
    keys, coords, _ = _rows(d, [col])
    if "COOR_X" not in d:  # geometry field components are X, Y, Z
        coords = np.array([[float(d.get(c, [0.0] * len(keys))[i]) for c in ("X", "Y", "Z")]
                           for i in range(len(keys))])
        keys = [tuple(np.round(coords[i], 9)) for i in range(len(keys))] if "NOEUD" not in d else keys
    # one row per node (a table can list a node more than once)
    seen, ukeys, ucoords = set(), [], []
    for k, c in zip(keys, coords):
        if k not in seen:
            seen.add(k); ukeys.append(k); ucoords.append(c)
    return ukeys, np.array(ucoords)


def nodal_result_values(result, field, components, node_group, keys, inst=None, precision=1e-9):
    """Values (len(keys), len(components)) of ``field`` at ``inst`` on ``node_group``,
    ordered like ``keys``.  ``inst=None`` takes the last stored instant."""
    POST_RELEVE_T = cmd("POST_RELEVE_T")
    action = dict(OPERATION="EXTRACTION", INTITULE="PRECICE_" + field, RESULTAT=result,
                  NOM_CHAM=field, GROUP_NO=node_group, NOM_CMP=tuple(components))
    if inst is None:
        action["TOUT_ORDRE"] = "OUI"
    else:
        action.update(INST=inst, CRITERE="ABSOLU", PRECISION=precision)
    d = table_to_dict(POST_RELEVE_T(ACTION=action))
    rkeys, rcoords, vals = _rows(d, list(components))
    # inst=None: steady thermal results, where only the last computed instant matters
    if inst is None and "INST" in d:  # keep the last instant per node
        insts = np.array([float(x) for x in d["INST"]])
        last = insts.max()
        sel = np.isclose(insts, last, rtol=0, atol=max(abs(last), 1.0) * 1e-12)
        rkeys = [k for k, s in zip(rkeys, sel) if s]
        vals = vals[sel]
    # reorder the table rows to the order of the preCICE vertices
    lookup = {k: v for k, v in zip(rkeys, vals)}
    missing = [k for k in keys if k not in lookup]
    if missing:
        raise KeyError(f"{len(missing)} coupling nodes missing from {field} extraction "
                       f"(first: {missing[0]!r})")
    return np.array([lookup[k] for k in keys])


def time_list(t0, t1, nsteps=1):
    """DEFI_LIST_REEL from t0 to t1 in ``nsteps`` equal steps."""
    DEFI_LIST_REEL = cmd("DEFI_LIST_REEL")
    return DEFI_LIST_REEL(DEBUT=t0, INTERVALLE=dict(JUSQU_A=t1, NOMBRE=int(nsteps)))


def coordinate_formula(interpolant, dim):
    """FORMULE f(X, Y[, Z]) that evaluates ``interpolant`` at the point."""
    FORMULE = cmd("FORMULE")
    # the interpolant is passed in the FORMULE context under the name used in VALE;
    # code_aster calls it with the coordinates of each point where the load is needed
    params = ("X", "Y", "Z")[:dim]
    expr = "precice_interp(" + ", ".join(params) + ")"
    return FORMULE(VALE=expr, NOM_PARA=params, precice_interp=interpolant)


def write_med(result, path, unit=81, inst=None, fields=None):
    """Write ``result`` (optionally one instant) to a MED file."""
    DEFI_FICHIER, IMPR_RESU = cmd("DEFI_FICHIER"), cmd("IMPR_RESU")
    DEFI_FICHIER(ACTION="ASSOCIER", UNITE=unit, FICHIER=path, TYPE="BINARY")
    resu = dict(RESULTAT=result)
    if inst is not None:
        resu.update(INST=inst, CRITERE="ABSOLU", PRECISION=1e-9)
    if fields:
        resu["NOM_CHAM"] = tuple(fields)
    IMPR_RESU(FORMAT="MED", UNITE=unit, RESU=resu)
    DEFI_FICHIER(ACTION="LIBERER", UNITE=unit)
