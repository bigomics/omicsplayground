#!/usr/bin/env bash
## Launch the standalone Methylome Profiler test app.
## Usage: ./run.sh [port]
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
PORT="${1:-8794}"
Rscript -e "shiny::runApp('.', port = ${PORT}, launch.browser = FALSE, host = '127.0.0.1')"
