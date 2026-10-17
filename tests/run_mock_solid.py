# codeaster-precice-adapter test: mock beam coupled through the real preCICE library
# Copyright (C) 2026 CFD FEA SERVICE S.R.L.
# SPDX-License-Identifier: GPL-3.0-or-later
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from precice_aster.core import CouplingLoop
from mock_beam import MockBeam

config = sys.argv[1]
beam = MockBeam()
history = []


class Tracking(MockBeam):
    pass


orig = beam.on_time_window_complete
def on_window(state, t, window):
    history.append((t, beam.tip(state)))
beam.on_time_window_complete = on_window

loop = CouplingLoop(beam, "Solid", config)
state, t = loop.run()
json.dump({"stats": loop.stats, "history": history}, open("mock_result.json", "w"))
print("tip displacement at end:", beam.tip(state))
