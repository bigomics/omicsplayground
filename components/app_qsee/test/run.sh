#!/usr/bin/env bash
## Launch the standalone QSEE (batch-effects/QC) test app.
##
## Usage: ./run.sh [port]
## A/B comparison toggles (see test/app.R): QSEE_LAZY=false QSEE_PURGE=false ./run.sh

set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

PORT="${1:-8793}"

Rscript -e "shiny::runApp('.', port = ${PORT}, launch.browser = TRUE, host = '127.0.0.1')"
