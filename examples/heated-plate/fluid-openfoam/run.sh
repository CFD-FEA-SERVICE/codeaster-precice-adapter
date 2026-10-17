#!/usr/bin/env bash
# Fluid participant (OpenFOAM + preCICE OpenFOAM adapter). Serial by default;
# "./run.sh -parallel" decomposes with system/decomposeParDict and runs with mpirun.
set -euo pipefail
cd "$(dirname "$0")"
here=$(pwd); case_root=$(dirname "$here")
# same absolute exchange directory as the code_aster side (it runs in a temp folder)
sed -i -E "s#exchange-directory=\"[^\"]*\"#exchange-directory=\"$case_root\"#" "$case_root/precice-config.xml"
blockMesh > log.blockMesh 2>&1
solver=$(foamDictionary -entry application -value system/controlDict)
if [ "${1:-}" = "-parallel" ]; then
    decomposePar -force > log.decomposePar 2>&1
    n=$(foamDictionary -entry numberOfSubdomains -value system/decomposeParDict)
    mpirun -np "$n" "$solver" -parallel 2>&1 | tee log.fluid
    reconstructPar > log.reconstructPar 2>&1
else
    "$solver" 2>&1 | tee log.fluid
fi
