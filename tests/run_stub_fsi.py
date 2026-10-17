# codeaster-precice-adapter test: precice_aster.fsi with a code_aster stub, coupled through preCICE
# Copyright (C) 2026 CFD FEA SERVICE S.R.L.
# SPDX-License-Identifier: GPL-3.0-or-later
import sys, os, json
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..")); sys.path.insert(0, here)
from precice_aster import aster_tools
from precice_aster.fsi import AsterFSISolid
from precice_aster.core import CouplingLoop
import aster_stub

mesh = aster_stub.flap_mesh()
aster_tools.set_command_provider(aster_stub.Stub(mesh))
solid = AsterFSISolid(mesh, model="MODEL", material="MATER", interface_group="flap",
                      loads=["BC"], log=lambda *a: print(*a, flush=True))
hist = []
orig = solid.on_time_window_complete
def track(state, t, w):
    orig(state, t, w)
    r = state.result
    hist.append((t, max(v[0] for v in r.nodal[max(r.nodal)].values())))
solid.on_time_window_complete = track
loop = CouplingLoop(solid, "Solid", os.path.abspath(sys.argv[1]))
state, t = loop.run()
from collections import Counter
json.dump({"stats": loop.stats, "history": hist, "calls": Counter(aster_stub.CALLS)},
          open("stub_result.json", "w"))
print("tip displacement at end:", hist[-1][1])
