# codeaster-precice-adapter test: precice_aster.cht with a thermal stub
# Copyright (C) 2026 CFD FEA SERVICE S.R.L.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Thermal stub: plate of thickness L = 0.25, conductivity k = 100, bottom at 310 K,
Robin condition on top (COEF_H, TEMP_EXT read from the FORMULE objects).
Exact 1D steady solution: T_top = (k/L*310 + h*T_ext) / (k/L + h)."""
import sys, os
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "..")); sys.path.insert(0, os.path.join(here, ".."))
import numpy as np
from precice_aster import aster_tools
from precice_aster.cht import AsterCHTSolid
from precice_aster.core import CouplingLoop
import aster_stub
from aster_stub import Mesh, Table, Formula, Load, _check


class ThermalStub(aster_stub.Stub):
    k, L, Tb = 100.0, 0.25, 310.0

    def AFFE_CHAR_THER_F(self, **kw):
        _check("AFFE_CHAR_THER_F", kw, ("MODELE", "ECHANGE"))
        _check("ECHANGE", kw["ECHANGE"], ("GROUP_MA", "COEF_H", "TEMP_EXT"))
        return Load("robin", kw)

    def THER_LINEAIRE(self, **kw):
        _check("THER_LINEAIRE", kw, ("MODELE", "CHAM_MATER", "EXCIT", "ETAT_INIT", "INCREMENT", "PARM_THETA"))
        assert kw["ETAT_INIT"] == {"STATIONNAIRE": "OUI"}
        robin = [e["CHARGE"] for e in kw["EXCIT"] if isinstance(e["CHARGE"], Load)][0].kw["ECHANGE"]
        nodes = sorted({n for c in self.mesh.cell_groups[robin["GROUP_MA"]] for n in c})
        res = aster_stub.Result()
        res.nodal[0.0] = {}
        for n in nodes:
            x, y = self.mesh.coords[n, :2]
            h, Te = robin["COEF_H"](x, y), robin["TEMP_EXT"](x, y)
            res.nodal[0.0][n] = ((self.k / self.L * self.Tb + h * Te) / (self.k / self.L + h),)
        return res


x = np.linspace(0, 1, 41)
coords = np.column_stack([x, np.full_like(x, 0.25)])
mesh = Mesh(coords, {"interface": [(i, i + 1) for i in range(40)]}, [f"N{i+1}" for i in range(41)], dim=2)
aster_tools.set_command_provider(ThermalStub(mesh))
solid = AsterCHTSolid(mesh, "MODEL", "MATER", "interface", loads=["T_BOTTOM"], conductivity=100.0,
                      delta=1e-5, steady=True, initial_temperature=310.0)
CouplingLoop(solid, "Solid", os.path.abspath(sys.argv[1])).run()
print("expected wall T:", (400 * 310 + 100 * 300) / 500)
