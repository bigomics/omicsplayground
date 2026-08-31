##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Run the full application with the Apps launcher enabled, for testing the
## Methylome Profiler mini app.
##
##   make methylome
##   make methylome pgx=~/Downloads/GSE43976-methyl-mini.pgx
##   make methylome port=3939
##
## The Apps tab is gated on opt$DEVMODE, which is only readable from
## etc/OPTIONS. Rather than editing that tracked file, we override opt in
## memory after global.R has parsed it.

port <- as.integer(Sys.getenv("METHYLOME_PORT", "3838"))
pgx_file <- Sys.getenv("METHYLOME_PGX", "")

## The dataset the panel guide and the capability deck are both built from.
## GSE87571: 149 whole-blood donors aged 15-83, real continuous age, both sexes.
DEMO <- "GSE87571-methyl-sub150.pgx"

## Where the demo dataset might already be sitting. Checked in order, first
## hit wins; data/ first so a dataset already in the library is left alone.
demo_candidates <- function() {
  c(
    file.path("data", DEMO),
    file.path(path.expand("~/Downloads/methylome-profiler"), DEMO),
    file.path(path.expand("~/Downloads"), DEMO),
    Sys.glob(file.path("/mnt/c/Users/*/Downloads/methylome-profiler", DEMO)),
    Sys.glob(file.path("/mnt/c/Users/*/Downloads", DEMO))
  )
}

## .pgx files are read with load(), not readRDS(). A pgx written with
## saveRDS() is silently skipped by the folder scan ("magic number 'X'"), so
## convert it rather than letting it vanish from the library.
ensure_rdata <- function(f) {
  ok <- tryCatch({
    e <- new.env(); load(f, envir = e); TRUE
  }, error = function(e) FALSE)
  if (!ok) {
    message("[methylome] ", basename(f), " is RDS; converting to save() format")
    pgx <- readRDS(f)
    save(pgx, file = f, compress = TRUE)
  }
  invisible(f)
}

stage <- function(src) {
  dest <- file.path("data", basename(src))
  if (normalizePath(src, mustWork = FALSE) != normalizePath(dest, mustWork = FALSE)) {
    file.copy(src, dest, overwrite = TRUE)
    ## Drop the whole index, not just the csv: a stale datasets-sigdb.h5
    ## makes the rescan die in H5Dcreate and no dataset list is written.
    unlink(Sys.glob("data/datasets-*"))
    message("[methylome] staged ", dest)
  }
  ensure_rdata(dest)
}

if (nzchar(pgx_file)) {
  if (!file.exists(pgx_file)) stop("no such pgx file: ", pgx_file)
  stage(pgx_file)
} else {
  ## Default: mount the methylation demo dataset without being asked.
  hit <- Filter(file.exists, demo_candidates())
  if (length(hit)) {
    stage(hit[1])
  } else {
    message(
      "[methylome] ", DEMO, " not found in data/ or Downloads.\n",
      "[methylome] Build it with:  Rscript components/app_methylome/test/build-test-pgx.Rscript\n",
      "[methylome] or point at one: make methylome pgx=/path/to/dataset.pgx"
    )
  }
}

setwd("components/app/R")
source("global.R")

opt$DEVMODE <- TRUE ## show the Apps launcher, which carries the Methylome tile
message("[methylome] DEVMODE enabled for this session only; etc/OPTIONS untouched")

source("ui.R")
source("server.R")

message("[methylome] Library -> load a methylomics dataset; it opens the Methylome app directly")

shinyApp(
  ui = app_ui,
  server = app_server,
  uiPattern = ".*",
  options = list(
    launch.browser = FALSE,
    host = "0.0.0.0",
    port = port
  )
)
