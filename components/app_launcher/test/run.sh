#!/usr/bin/env bash
## Launch the standalone Smart Tools app-launcher test app.
##
## Usage: ./run.sh [port]

set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

PORT="${1:-8792}"

Rscript -e "shiny::runApp('.', port = ${PORT}, launch.browser = TRUE, host = '127.0.0.1')"
