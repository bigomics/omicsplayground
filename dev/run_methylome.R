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

## Stage a dataset into data/ if one was named, and drop the folder index so
## the library rescans it on first open.
if (nzchar(pgx_file)) {
  if (!file.exists(pgx_file)) stop("no such pgx file: ", pgx_file)
  dest <- file.path("data", basename(pgx_file))
  file.copy(pgx_file, dest, overwrite = TRUE)

  ## .pgx files are read with load(), not readRDS(). A pgx written with
  ## saveRDS() is silently skipped by the folder scan ("magic number 'X'"),
  ## so convert it here rather than letting it disappear from the library.
  is_rdata <- tryCatch({
    e <- new.env(); load(dest, envir = e); TRUE
  }, error = function(e) FALSE)
  if (!is_rdata) {
    message("[methylome] ", basename(dest), " is RDS; converting to save() format")
    pgx <- readRDS(dest)
    save(pgx, file = dest, compress = TRUE)
    rm(pgx)
  }
  unlink(Sys.glob("data/datasets-*"))
  message("[methylome] staged ", dest)
}

setwd("components/app/R")
source("global.R")

opt$DEVMODE <- TRUE ## show the Apps launcher, which carries the Methylome tile
message("[methylome] DEVMODE enabled for this session only; etc/OPTIONS untouched")

source("ui.R")
source("server.R")

message("[methylome] Library -> load a methylomics dataset -> Apps -> Methylome Profiler")

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
