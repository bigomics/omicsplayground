#!/usr/bin/env bash
# Build bigomics/omicsplayground:trixie-edgy on shibuya, from scratch.
#
# Not run by CI. This is the host-specific orchestration script used on the
# shibuya VM (tokyo's Proxmox VM 101) to hand-build the trixie image locally,
# because bigomics/omicsplayground:trixie-edgy is not published to Docker Hub.
# See shiny-deploy/12-shinyproxy-docker-service-swarm/tokyo/AGENTS.md gotcha 6
# on that deploy for the full context.
#
# Expects $HOME/trixie-build/{playbase,omicsplayground} to be checkouts of
# those two repos (this script does not clone them) and $HOME/trixie-build/.pat
# to contain a GITHUB_PAT.
#
# Built here, not on the Proxmox host: enabling docker there sets iptables
# FORWARD DROP and makes every VM unreachable, which would take down this very
# deploy. Building here also means the image lands directly in the daemon
# ShinyProxy uses, with no transfer.
#
# REQUIREMENT: host CPU must expose AVX/AVX2. PPM ships tiledb as a binary built
# for a modern ISA, and TileDBArray dyn.loads it during install -- without AVX
# that is an instant SIGILL, which fails TileDBArray -> playbase (an Import) ->
# the whole build.
#   grep -w avx2 /proc/cpuinfo   # must print something
set -euo pipefail

ROOT=$HOME/trixie-build
BRANCH=edgy                     # opg code branch cloned INSIDE the image
PB_BRANCH=edgy                  # playbase branch. Explicit: Dockerfile.{,update}
                                 # install exactly this (fail-fast if missing)
                                 # instead of the same-name check with silent
                                 # fallback to playbase@main. Decouples the two
                                 # repos -- no same-named playbase alias needed.
SQUASH=$HOME/.local/bin/docker-squash
export DOCKER_BUILDKIT=1

GITHUB_PAT=$(cat "$ROOT/.pat"); export GITHUB_PAT
test -n "$GITHUB_PAT"

step(){ echo; echo "=================== $* ($(date -Is)) ==================="; echo; }

#---------------------------------------------------------------------
# playbase base.
# NOT `make docker`: its docker.squash guard is always-true and always exits 1,
# and docker.pkg pipes to tee without pipefail so a failed build returns 0.
#---------------------------------------------------------------------
cd "$ROOT/playbase"
cp dev/rspm.R dev/Rprofile          # gitignored; without it everything builds from source
grep -q '2026-07-15' dev/Rprofile   # PPM pin must be present or lattice drifts to 0.23-1

step "1/5  playbase OS layer"
if docker image inspect playbase-os >/dev/null 2>&1; then echo "reusing playbase-os"; else
  docker build --no-cache -f dev/Dockerfile.os -t playbase-os . 2>&1 | tee "$ROOT/os.log"
fi

step "2/5  playbase R base layer"
if docker image inspect playbase-rbase >/dev/null 2>&1; then echo "reusing playbase-rbase"; else
  docker build -f dev/Dockerfile.rbase -t playbase-rbase . 2>&1 | tee "$ROOT/rbase.log"
fi

step "3/5  playbase package layer + squash"
if docker image inspect bigomics/playbase:trixie >/dev/null 2>&1; then
  echo "reusing bigomics/playbase:trixie"
else
  umask 077; printf 'GITHUB_PAT=%s\n' "$GITHUB_PAT" > dev/Renviron
  docker build -f dev/Dockerfile -t playbase-pkg . 2>&1 | tee "$ROOT/pkg.log"
  rm -f dev/Renviron
  docker image inspect playbase-pkg >/dev/null
  # Squash for PAT hygiene, not size: the token is baked into an intermediate
  # layer and this docker group has two members (xavier, kwee).
  $SQUASH playbase-pkg -t bigomics/playbase:trixie
  docker image inspect bigomics/playbase:trixie >/dev/null
  docker rmi playbase-pkg
fi

#---------------------------------------------------------------------
# Plain `docker build` for the opg stages. On an earlier build host (docker
# 26.1.5) `quarto install tinytex` died under AppArmor with a deno SIGILL-style
# panic, which forced a privileged buildx builder plus a local registry to feed
# it the base image. Docker 29.7.1 does NOT reproduce that -- verified by
# running the exact quarto+tinytex install under default confinement -- so the
# registry and the container driver are both dropped. That also avoids pushing
# the squashed 20.5GB single-layer image, whose one ~10GB blob timed out.
#---------------------------------------------------------------------
step "4/5  omicsplayground stage 1 (docker/Dockerfile)"
cd "$ROOT/omicsplayground"
docker build --no-cache --progress plain \
  --build-arg BRANCH=$BRANCH \
  --build-arg PLAYBASE_BRANCH=$PB_BRANCH \
  --build-arg BASE_IMAGE=bigomics/playbase:trixie \
  -f docker/Dockerfile \
  -t bigomics/omicsplayground:trixie-edgy-base . 2>&1 | tee "$ROOT/stage1.log"
docker image inspect bigomics/omicsplayground:trixie-edgy-base >/dev/null

step "5/5  omicsplayground stage 2 (docker/Dockerfile.update)"
docker build --progress plain \
  --build-arg BRANCH=$BRANCH \
  --build-arg PLAYBASE_BRANCH=$PB_BRANCH \
  --build-arg BASE_IMAGE=bigomics/omicsplayground:trixie-edgy-base \
  --build-arg CACHEBUST_CODE="$(date +%s)" \
  --secret id=GITHUB_PAT,env=GITHUB_PAT \
  -f docker/Dockerfile.update \
  -t bigomics/omicsplayground:trixie-edgy . 2>&1 | tee "$ROOT/stage2.log"
docker image inspect bigomics/omicsplayground:trixie-edgy >/dev/null

#---------------------------------------------------------------------
step "SMOKE TEST"
#---------------------------------------------------------------------
docker run --rm --entrypoint sh bigomics/omicsplayground:trixie-edgy \
  -c '. /etc/os-release; echo "$PRETTY_NAME"; quarto --version'
# Asserts the packages LOAD, not just that they are installed: an earlier run
# passed an installed-only check while methylumi/lumi/wateRmelon were unloadable
# because lattice had been upgraded out from under them. tiledb/TileDBArray are
# the AVX canaries.
docker run --rm --entrypoint Rscript bigomics/omicsplayground:trixie-edgy -e '
  cat(R.version.string, "| BioC", as.character(BiocManager::version()), "\n")
  cat("lattice", as.character(packageVersion("lattice")),
      "| exports parallel:", "parallel" %in% getNamespaceExports("lattice"), "\n")
  need <- c("collapse","isotree","WGCNAplus","omicsai","playbase","tiledb",
            "TileDBArray","PCSF","visNetwork","wateRmelon","methylumi","lumi")
  bad <- character(0)
  for (x in need) {
    r <- tryCatch({ suppressMessages(loadNamespace(x)); "LOADS" },
                  error = function(e) paste("FAILS:", conditionMessage(e)))
    cat(sprintf("%-12s %s\n", x, r)); if (r != "LOADS") bad <- c(bad, x)
  }
  cat("playbase", as.character(packageVersion("playbase")), "\n")
  if (length(bad)) { cat("FAILED TO LOAD:", paste(bad, collapse=" "), "\n"); quit(status=1) }'

rm -f "$ROOT/.pat"
step "BUILD COMPLETE"
