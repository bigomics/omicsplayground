#!/usr/bin/env bash
## Launch the standalone Methylome Profiler test app.
##
## Usage: ./run.sh [port]
##        METHYLOME_PGX=/path/to/dataset.pgx ./run.sh

set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

PORT="${1:-8793}"

Rscript -e "shiny::runApp('.', port = ${PORT}, launch.browser = TRUE, host = '127.0.0.1')"
