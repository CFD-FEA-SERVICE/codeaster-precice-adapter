# codeaster-precice-adapter: couple code_aster with other solvers through preCICE
# Copyright (C) 2026 CFD FEA SERVICE S.R.L.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""preCICE adapter for code_aster (version 17 and later).

Modules
-------
core           solver-agnostic preCICE coupling loop with checkpointing
interpolation  nodal interpolants used to feed preCICE data into code_aster loads
fsi            fluid-structure interaction: reads Force, writes Displacement
cht            conjugate heat transfer: Robin-Robin coupling (sink temperature + HTC)

The code_aster specific modules are imported lazily, so the core can be
tested without code_aster installed.
"""

__version__ = "0.1.0"

from .core import CouplingLoop, SolverInterface  # noqa: F401
from .interpolation import NodalInterpolant  # noqa: F401
