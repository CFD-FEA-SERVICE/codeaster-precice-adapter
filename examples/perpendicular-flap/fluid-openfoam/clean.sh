#!/usr/bin/env sh
cd "$(dirname "$0")"
foamListTimes -rm > /dev/null 2>&1 || true
rm -rf processor* constant/polyMesh log.* precice-*-iterations.log precice-*-convergence.log \
       precice-*-watchpoint-*.log precice-profiling postProcessing
