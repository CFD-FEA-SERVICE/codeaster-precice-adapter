#!/usr/bin/env python3
# codeaster-precice-adapter: couple code_aster with other solvers through preCICE
# Copyright (C) 2026 CFD FEA SERVICE S.R.L.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Structured 2D rectangle mesh (QUAD8 + SEG3 edges) in code_aster native format.

usage: rect_mesh.py out.mail x0 x1 y0 y1 nx ny  [edge groups]
edge groups: name=side[+side...] with side in bottom,right,top,left
example (perpendicular flap):
  rect_mesh.py flap.mail -0.05 0.05 0 1 4 40 flap=left+top+right base=bottom
Every rectangle also gets the cell group 'solid' (all QUAD8).
"""
import sys


def build(x0, x1, y0, y1, nx, ny):
    """Return node dict (i,j)->(name,x,y), quads, edges per side (quadratic grid indices)."""
    nodes, names = {}, {}
    k = 0
    for j in range(2 * ny + 1):
        for i in range(2 * nx + 1):
            if i % 2 == 1 and j % 2 == 1:
                continue  # QUAD8: no centre node
            k += 1
            x = x0 + (x1 - x0) * i / (2 * nx)
            y = y0 + (y1 - y0) * j / (2 * ny)
            nodes[(i, j)] = (f"N{k}", x, y)
    quads = []
    for j in range(ny):
        for i in range(nx):
            I, J = 2 * i, 2 * j
            ids = [(I, J), (I + 2, J), (I + 2, J + 2), (I, J + 2),       # corners, counter-clockwise
                   (I + 1, J), (I + 2, J + 1), (I + 1, J + 2), (I, J + 1)]  # mid-side nodes
            quads.append([nodes[p][0] for p in ids])
    edges = {"bottom": [], "top": [], "left": [], "right": []}
    for i in range(nx):
        I = 2 * i
        edges["bottom"].append([nodes[(I, 0)][0], nodes[(I + 2, 0)][0], nodes[(I + 1, 0)][0]])
        edges["top"].append([nodes[(I + 2, 2 * ny)][0], nodes[(I, 2 * ny)][0], nodes[(I + 1, 2 * ny)][0]])
    for j in range(ny):
        J = 2 * j
        edges["right"].append([nodes[(2 * nx, J)][0], nodes[(2 * nx, J + 2)][0], nodes[(2 * nx, J + 1)][0]])
        edges["left"].append([nodes[(0, J + 2)][0], nodes[(0, J)][0], nodes[(0, J + 1)][0]])
    return nodes, quads, edges


def write_mail(path, nodes, quads, edges, groups, title="precice_aster rectangle"):
    lines = ["TITRE", f" {title}", "FINSF", "%", "COOR_2D"]
    for name, x, y in nodes.values():
        lines.append(f" {name:<8s} {x: .12e} {y: .12e}")
    lines += ["FINSF", "%", "QUAD8"]
    qnames = []
    for q, conn in enumerate(quads, 1):
        qnames.append(f"Q{q}")
        lines.append(f" Q{q:<7d} " + " ".join(conn))
    lines += ["FINSF", "%", "SEG3"]
    side_cells = {s: [] for s in edges}
    e = 0
    for side, segs in edges.items():
        for conn in segs:
            e += 1
            side_cells[side].append(f"S{e}")
            lines.append(f" S{e:<7d} " + " ".join(conn))
    lines += ["FINSF", "%"]

    def group(name, cells):
        out = ["GROUP_MA", f" {name}"]
        for i in range(0, len(cells), 8):
            out.append(" " + " ".join(cells[i:i + 8]))
        return out + ["FINSF", "%"]

    lines += group("solid", qnames)
    for gname, sides in groups.items():
        cells = [c for s in sides for c in side_cells[s]]
        lines += group(gname, cells)
    lines.append("FIN")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def main(argv):
    if len(argv) < 8:
        print(__doc__); return 1
    out = argv[1]
    x0, x1, y0, y1 = map(float, argv[2:6])
    nx, ny = int(argv[6]), int(argv[7])
    groups = {}
    for g in argv[8:]:
        name, sides = g.split("=")
        groups[name] = sides.split("+")
    nodes, quads, edges = build(x0, x1, y0, y1, nx, ny)
    write_mail(out, nodes, quads, edges, groups)
    print(f"{out}: {len(nodes)} nodes, {len(quads)} QUAD8, groups solid, {', '.join(groups)}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
