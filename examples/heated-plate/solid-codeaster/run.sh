#!/usr/bin/env bash
# Solid participant of the heated-plate case: code_aster through precice_aster.
# Start the fluid participant (e.g. ../fluid-openfoam) in parallel, in any order.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
case_root=$(dirname "$here")
cfg="$case_root/precice-config.xml"

# code_aster runs in a temporary folder: preCICE needs an absolute exchange directory.
# Both participants read the same file, so it is edited in place (idempotent).
sed -i -E "s#exchange-directory=\"[^\"]*\"#exchange-directory=\"$case_root\"#" "$cfg"

export PRECICE_CONFIG_FILE="$cfg"
export PRECICE_PARTICIPANT="${PRECICE_PARTICIPANT:-Solid}"
export PRECICE_CASE_DIR="$here"
export PRECICE_ASTER_ROOT="${PRECICE_ASTER_ROOT:-/opt/codeaster-precice-adapter}"

cd "$here"
[ -f plate.mail ] || python3 "$PRECICE_ASTER_ROOT/tools/rect_mesh.py" plate.mail 0 1 0 0.25 80 20 interface=top bottom=bottom left=left right=right

# run_aster from the conda environment (override with ASTER_RUN=...)
${ASTER_RUN:-/usr/local/bin/micromamba run -r /opt/conda -n aster run_aster} plate.export 2>&1 | tee log.solid
