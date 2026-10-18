#!/usr/bin/env bash
# Fluid participant (OpenFOAM pimpleFoam + preCICE adapter), mesh in constant/polyMesh.
# Serial by default; "./run.sh -parallel" decomposes and runs with mpirun.
set -euo pipefail
cd "$(dirname "$0")"
here=$(pwd); case_root=$(dirname "$here")
sed -i -E "s#exchange-directory=\"[^\"]*\"#exchange-directory=\"$case_root\"#" "$case_root/precice-config.xml"
if [ "${1:-}" = "-parallel" ]; then
    decomposePar -force > log.decomposePar 2>&1
    n=$(foamDictionary -entry numberOfSubdomains -value system/decomposeParDict)
    mpirun -np "$n" pimpleFoam -parallel 2>&1 | tee log.fluid
    reconstructPar > log.reconstructPar 2>&1
else
    pimpleFoam 2>&1 | tee log.fluid
fi
