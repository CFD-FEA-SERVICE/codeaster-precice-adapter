#!/usr/bin/env bash
# Run a fake fluid and a solid test participant side by side (local tests only).
# usage: run_pair.sh <case dir> <fluid script> <solid script> [timeout s]
set -u
case_dir=$1; fluid=$2; solid=$3; tmo=${4:-900}
cd "$case_dir" || exit 1
( cd fluid && timeout "$tmo" python -I "$fluid" > ../fluid.log 2>&1 ) &
fpid=$!
( cd solid && timeout "$tmo" python "$solid" ../precice-config.xml > ../solid.log 2>&1 )
srv=$?
wait $fpid; frv=$?
echo "fluid exit $frv, solid exit $srv"
