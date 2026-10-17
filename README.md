# codeaster-precice-adapter

Couple **code_aster** (17.x) with OpenFOAM, CalculiX, SU2, FEniCS, … through **preCICE v3**.

| Module | code_aster side | Reads | Writes |
|---|---|---|---|
| `precice_aster.fsi` | `DYNA_NON_LINE` (Newmark), one call per coupling step | `Force` (nodal) | `Displacement` |
| `precice_aster.cht` | `THER_LINEAIRE` / `THER_NON_LINE`, steady or transient | `Sink-Temperature-*`, `Heat-Transfer-Coefficient-*` (Robin) | wall `T`, `k/δ` |

This is a clean rewrite. It shares no code with the original `precice/code_aster-adapter`, which covers CHT only, targets code_aster 14 and is licensed GPL-2.0.

License: **GPL-3.0-or-later** (see `LICENSE`). © 2026 CFD FEA SERVICE S.R.L.
The `fluid-openfoam` folders and the flap `precice-config.xml` in `examples/` come from the
[preCICE tutorials](https://github.com/precice/tutorials) (commit 3d47c2c, LGPL-3.0). Changes: standalone `run.sh`/`clean.sh` scripts, and for the heated plate 1 cell in z and 4000 iterations.

## How it works

```
 preCICE  ──Force──▶  NodalInterpolant ──▶ FORMULE(X,Y[,Z]) ──▶ AFFE_CHAR_MECA_F(FORCE_NODALE)
                                                                         │
 preCICE  ◀─Displacement── POST_RELEVE_T(DEPL @ t+dt) ◀── DYNA_NON_LINE(ETAT_INIT = checkpoint)
```

* The coupling mesh is made of the **nodes** of a code_aster cell group, for example `flap`. The adapter creates the node group `<group>_PRCN` itself.
* The load is built **once**. Its components are `FORMULE` objects that read Python arrays, so new preCICE data only refreshes the arrays.
* **Implicit coupling** works without code_aster save/restore. Every coupling step is a new `DYNA_NON_LINE` call restarted from the stored result (`ETAT_INIT=_F(EVOL_NOLI=…, INST=…)`). Reverting to a checkpoint means restarting from the older result.
* `core.CouplingLoop` holds all the preCICE logic (checkpoints, sub-cycling, time-window bookkeeping). It knows nothing about code_aster. Any participant that implements `SolverInterface` can use it.

## Requirements

* code_aster ≥ 17 (tested target: conda-forge `code-aster=17.3.17` openmpi build)
* `pyprecice` and `libprecice` **3.3.x in the same environment as code_aster**
* numpy

Docker / Apptainer (cloudHPC image), added after the micromamba line:

```dockerfile
COPY codeaster-precice-adapter /opt/codeaster-precice-adapter
ENV PRECICE_ASTER_ROOT=/opt/codeaster-precice-adapter
```

(or `RUN git clone <repo-url> /opt/codeaster-precice-adapter`). No build step is needed; it is pure Python.

## Using it in a .comm file

```python
import os, sys
sys.path.insert(0, os.environ.get("PRECICE_ASTER_ROOT", "/opt/codeaster-precice-adapter"))
from code_aster.Commands import *
DEBUT()
mesh  = LIRE_MAILLAGE(UNITE=20, FORMAT="MED")
model = AFFE_MODELE(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"))
...
from precice_aster.fsi import run_fsi
run_fsi(mesh, model, mater, interface_group="wet", loads=[clamp])   # reads $PRECICE_CONFIG_FILE
FIN()
```

Main options of `run_fsi` / `AsterFSISolid`:
`mesh_name="Solid-Mesh"`, `force_name="Force"`, `displacement_name="Displacement"`,
`deformation="PETIT"|"GROT_GDEP"`, `relation="ELAS"`, `newmark=(beta, gamma)`, `substeps=1`,
`output_dir`, `output_every` (MED every N windows), `dyna_options={...}`. The last one is passed to `DYNA_NON_LINE`, for example `SOLVEUR` or `CONVERGENCE`.

`run_cht` / `AsterCHTSolid`: `conductivity`, `delta` (HTC = k/δ), `steady=True`,
`initial_temperature`, `nonlinear=False`, `mesh_name="Solid-Interface-Nodes"`, `read_names`, `write_names`.

Environment variables:

| Variable | Meaning |
|---|---|
| `PRECICE_CONFIG_FILE` | absolute path of `precice-config.xml` (code_aster runs in a temp dir) |
| `PRECICE_PARTICIPANT` | participant name, default `Solid` |
| `PRECICE_CASE_DIR` | case folder (used by the examples for `results/`) |

**The `exchange-directory` of `<m2n:sockets>` must be absolute.** The example `run.sh` scripts
rewrite it on each run.

## Examples

| Folder | Coupling | Compare with |
|---|---|---|
| `examples/selftest` | code_aster only, **no preCICE** | beam theory |
| `examples/perpendicular-flap` | OpenFOAM `pimpleFoam` ↔ code_aster FSI (2D `D_PLAN`) | CalculiX / solids4foam tutorial runs |
| `examples/heated-plate` | OpenFOAM `buoyantSimpleFoam` ↔ code_aster CHT (2D `PLAN`) | `chtMultiRegionSimpleFoam` OpenFOAM-only case |

Run a coupled case by starting both participants, in any order:

```bash
cd examples/perpendicular-flap
./fluid-openfoam/run.sh &          # or: ./fluid-openfoam/run.sh -parallel &
./solid-codeaster/run.sh
```

`solid-codeaster/run.sh` calls `micromamba run -r /opt/conda -n aster run_aster`. Override it with
`ASTER_RUN="run_aster"` when code_aster is already on the PATH.

## Test plan (new installation)

1. **Self-test, code_aster only, about 1 min, 1 core.**
   `cd examples/selftest && micromamba run -r /opt/conda -n aster run_aster selftest.export`.
   Then `grep SELFTEST selftest.mess`. The run must end with `PRECICE_ASTER SELFTEST: PASS`. It checks five things:
   * 169 coupling nodes
   * static tip deflection within 5 % of Timoshenko theory
   * load update (FORMULE re-read)
   * bit-identical restart from a checkpoint
   * first period within 8 % of Euler-Bernoulli
2. **Perpendicular flap, 5 s, 500 windows.** Tip displacement is in `solid-codeaster/precice-Solid-watchpoint-Flap-Tip.log`. Compare the peak (≈0.12–0.14 m) and the period (≈1.25 s) with solids4foam and CalculiX.
3. **Heated plate, 400 windows, 4000 SIMPLE iterations.** Compare the interface temperature along x and the heat flow with the OpenFOAM-only case. Reference values: interface T 303.7–307.2 K, Q ≈ 586 W per 0.4 m depth. Divide by 0.4 to compare with this 2D case per metre.

## Tests without code_aster

`tests/` holds tests that run the coupling logic against preCICE without code_aster: a mock beam and a code_aster command stub, with the preCICE tutorial fake fluid. Results:

* Flap: 500 windows and 274 checkpoint restores; the FSI module reproduces the mock solver to 3e-11.
* CHT: the wall temperature matches the analytic value exactly.

These tests validate the preCICE side and the command keywords used; the self-test validates code_aster.

## Limitations (v0.1)

* Serial code_aster (`mpi_nbcpu 1`) for the coupling process. The linear solver can still use threads (`dyna_options={"SOLVEUR": ...}`).
* Nodal forces only (conservative mapping on the fluid side, as in the preCICE tutorials).
* With `GROT_GDEP` the force is applied at the initial node positions, as a follower-free nodal load.
* Each coupling step creates a new result concept. Memory grows with the number of steps, which is fine for 10³–10⁴ steps on small meshes.
