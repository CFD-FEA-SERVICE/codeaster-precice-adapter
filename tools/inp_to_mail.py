#!/usr/bin/env python3
# codeaster-precice-adapter: couple code_aster with other solvers through preCICE
# Copyright (C) 2026 CFD FEA SERVICE S.R.L.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Convert an Abaqus/CalculiX tetrahedral mesh (.inp/.msh + .nam + .sur) to code_aster .mail.

usage: inp_to_mail.py out.mail mesh.inp [sets.nam ...] [surface.sur ...]

* elements: C3D4 -> TETRA4, C3D10 -> TETRA10 (same node order in both codes;
  elements with negative volume are reoriented)
* *NSET, NSET=Nfoo         -> node group  'foo'  (leading 'N' dropped)
* *SURFACE, TYPE=ELEMENT   -> cell group  'foo'  of TRIA3/TRIA6 faces (leading 'S' dropped),
                              built from the element faces S1..S4
* all volume elements      -> cell group 'solid'
Group names are cut to 24 characters (code_aster limit).
"""
import re
import sys
from collections import OrderedDict

import numpy as np

# Abaqus tetra faces (corner indices, then mid-side indices), 0-based
FACES4 = {"S1": (0, 1, 2), "S2": (0, 3, 1), "S3": (1, 3, 2), "S4": (2, 3, 0)}
FACES10 = {"S1": (0, 1, 2, 4, 5, 6), "S2": (0, 3, 1, 7, 8, 4),
           "S3": (1, 3, 2, 8, 9, 5), "S4": (2, 3, 0, 9, 7, 6)}


def read_blocks(paths):
    """Yield (keyword line, [data lines]) for every *KEYWORD block in the files."""
    for path in paths:
        key, data = None, []
        with open(path) as f:
            for raw in f:
                line = raw.strip()
                if not line or line.startswith("**"):
                    continue
                if line.startswith("*"):
                    if key is not None:
                        yield key, data
                    key, data = line, []
                else:
                    data.append(line)
        if key is not None:
            yield key, data


def opt(key, name):
    m = re.search(name + r"\s*=\s*([^,\s]+)", key, re.I)
    return m.group(1) if m else None


def ints(lines):
    out = []
    for ln in lines:
        out += [int(v) for v in ln.replace(",", " ").split()]
    return out


def main(argv):
    out, inputs = argv[0], argv[1:]
    nodes, elems, etype = OrderedDict(), OrderedDict(), None
    elsets, nsets, surfaces = {}, OrderedDict(), OrderedDict()
    for key, data in read_blocks(inputs):
        k = key.split(",")[0].lower()
        if k == "*node":
            for ln in data:
                v = [x for x in ln.replace(",", " ").split()]
                nodes[int(v[0])] = tuple(float(x) for x in v[1:4])
        elif k == "*element":
            t = opt(key, "type").upper()
            if t not in ("C3D4", "C3D10"):
                raise SystemExit(f"unsupported element type {t}")
            etype = etype or t
            n = 4 if t == "C3D4" else 10
            flat = ints(data)
            for i in range(0, len(flat), n + 1):
                elems[flat[i]] = flat[i + 1:i + 1 + n]
        elif k == "*elset":
            elsets[opt(key, "elset")] = ints(data)
        elif k == "*nset":
            nsets[opt(key, "nset")] = ints(data)
        elif k == "*surface":
            surfaces[opt(key, "name")] = [tuple(x.strip() for x in ln.split(",")) for ln in data]

    faces_tab = FACES4 if etype == "C3D4" else FACES10
    face_cells, face_groups = [], OrderedDict()
    # faces are built from the original (Abaqus) connectivity, before any reorientation
    for sname, entries in surfaces.items():
        gname = (sname[1:] if sname[0].upper() == "S" else sname)[:24]
        members = []
        for ref, side in entries:
            ids = elsets.get(ref, None)
            ids = ids if ids is not None else [int(ref)]
            for eid in ids:
                con = elems[eid]
                members.append(len(face_cells))
                face_cells.append([con[j] for j in faces_tab[side.upper()]])
        face_groups[gname] = members

    xyz = {i: np.array(c) for i, c in nodes.items()}
    flipped = 0
    for eid, con in elems.items():
        a, b, c, d = (xyz[con[j]] for j in range(4))
        if np.dot(np.cross(b - a, c - a), d - a) < 0:      # reorient: swap nodes 2 and 3
            con[1], con[2] = con[2], con[1]
            if len(con) == 10:  # mid-sides 5<->7 (1-2 <-> 3-1), 9<->10 (2-4 <-> 3-4)
                con[4], con[6] = con[6], con[4]
                con[8], con[9] = con[9], con[8]
            flipped += 1

    nn = len(nodes)
    vt = "TETRA4" if etype == "C3D4" else "TETRA10"
    ft = "TRIA3" if etype == "C3D4" else "TRIA6"
    with open(out, "w") as f:
        f.write("TITRE\n converted by precice_aster inp_to_mail.py\nFINSF\n%\n")
        f.write("COOR_3D\n")
        for i, (x, y, z) in nodes.items():
            f.write(f" N{i} {x:.12e} {y:.12e} {z:.12e}\n")
        f.write("FINSF\n%\n")
        f.write(f"{vt}\n")
        for eid, con in elems.items():
            f.write(f" M{eid} " + " ".join(f"N{c}" for c in con) + "\n")
        f.write("FINSF\n%\n")
        if face_cells:
            f.write(f"{ft}\n")
            for i, con in enumerate(face_cells):
                f.write(f" F{i + 1} " + " ".join(f"N{c}" for c in con) + "\n")
            f.write("FINSF\n%\n")
        f.write("GROUP_MA\n solid\n")
        for k, eid in enumerate(elems):
            f.write(f" M{eid}" + ("\n" if k % 8 == 7 else ""))
        f.write("\nFINSF\n%\n")
        for g, members in face_groups.items():
            f.write(f"GROUP_MA\n {g}\n")
            for k, m in enumerate(members):
                f.write(f" F{m + 1}" + ("\n" if k % 8 == 7 else ""))
            f.write("\nFINSF\n%\n")
        for name, ids in nsets.items():
            g = (name[1:] if name[0].upper() == "N" else name)[:24]
            f.write(f"GROUP_NO\n {g}\n")
            for k, n in enumerate(ids):
                f.write(f" N{n}" + ("\n" if k % 8 == 7 else ""))
            f.write("\nFINSF\n%\n")
        f.write("FIN\n")
    print(f"{out}: {nn} nodes, {len(elems)} {vt}, {len(face_cells)} {ft}, reoriented {flipped}; "
          f"cell groups solid + {list(face_groups)}, node groups {[n[1:] for n in nsets]}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    main(sys.argv[1:])
