#!/usr/bin/env bash
#
# Run gate 1 (params contract + preprocess seam) across every scenario.
#
#   ./run_matrix.sh [scenario ...]      # default: all scenarios in scenarios.json
#
# The app must already be up (`make app.start`). Each scenario drives a fresh
# browser session, so budget ~5-15 min per scenario depending on fixture size.
# Failures do not stop the batch -- the point is a full report, not a fast exit.
#
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OPG_MAIN="$(cd "$HERE/../.." && pwd)"
# see Makefile: R_LIBS_USER replaces the personal library, so append it back
R_PERSONAL="$HOME/R/x86_64-pc-linux-gnu-library/4.3"
PB_EDGY="$HERE/libs/pb-edgy:$R_PERSONAL"

if [ $# -gt 0 ]; then
  SCENARIOS=("$@")
else
  mapfile -t SCENARIOS < <(Rscript -e '
    cfg <- jsonlite::fromJSON("'"$HERE"'/scenarios.json", simplifyDataFrame = FALSE)
    cat(vapply(cfg$scenarios, function(s) s$id, ""), sep = "\n")' 2>/dev/null)
fi

echo "== matrix: ${#SCENARIOS[@]} scenario(s)"
for s in "${SCENARIOS[@]}"; do
  echo ""
  echo "================ $s ================"
  # A stale raw_ dir from a crashed run would be picked up by the poller.
  rm -rf "$OPG_MAIN"/data/USER_INPUT/raw_*

  if ! node "$HERE/drive_app.mjs" "$s" > "$HERE/runs/driver-$s.log" 2>&1; then
    echo "  CAPTURE FAILED -- see runs/driver-$s.log"
    grep -E "^\[drive\] (FAILED|page state)" "$HERE/runs/driver-$s.log" | tail -2
    continue
  fi

  R_LIBS_USER="$PB_EDGY" Rscript "$HERE/run_script.R" "$s" \
    > "$HERE/runs/script-$s.log" 2>&1 || { echo "  PATH B FAILED"; continue; }

  R_LIBS_USER="$PB_EDGY" Rscript "$HERE/compare.R" "$s" 1 2>&1 \
    | grep -E "^\s+\[|^== gate"
done

echo ""
echo "== done -> $HERE/runs/report.csv"
