##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Minimal standalone Shiny app to run the Methylome Profiler (app_methylome)
## in isolation, outside the full Playground application.
##
## Run with:
##   ./run.sh              # picks a methylomics pgx out of data/
##   ./run.sh 8793         # ...on another port
##   METHYLOME_PGX=/path/to/x.pgx ./run.sh
## or open this file in RStudio and click "Run App".
##
## Unlike app_idat's harness, this one cannot get away with plain shiny: the
## panels are PlotModule/TableModule, so the platform's module code has to be
## present. 00SourceAll.R is the same loader global.R uses, so what runs here
## is what runs in the app.
##

library(shiny)
library(bslib)
library(bigdash)
library(DT)

OPG <- normalizePath("../../..")
options(shiny.maxRequestSize = 999 * 1024^2)

## Globals the panels read that global.R would otherwise define. Only the two
## the methylome modules actually touch; everything else comes from the source
## tree below.
SCROLLY_MODAL <<- "55vh"
TABLE_HEIGHT_MODAL <<- "75vh"
WATERMARK <<- FALSE
opt <<- list(DEVMODE = TRUE, WATERMARK = FALSE, ENABLE_AI = FALSE)

## Resource paths, so the stylesheet and icons resolve as they do in the app.
shiny::addResourcePath("custom", file.path(OPG, "components/assets"))
shiny::addResourcePath("static", file.path(OPG, "components/assets"))

## The platform's modules (PlotModule, TableModule, withTooltip, bs_alert...)
## plus every app_methylome file.
local({
  wd <- getwd()
  setwd(file.path(OPG, "components"))
  on.exit(setwd(wd))
  source("00SourceAll.R", chdir = TRUE)
})

## A methylomics dataset to open. Env var wins; otherwise take the first
## methylomics pgx in data/, so a fresh clone with the demo staged just works.
pgx_file <- Sys.getenv("METHYLOME_PGX", "")
if (!nzchar(pgx_file)) {
  cand <- Sys.glob(file.path(OPG, "data", "*.pgx"))
  for (f in cand) {
    e <- new.env()
    ok <- tryCatch({ load(f, envir = e); TRUE }, error = function(e) FALSE)
    if (ok && identical(tolower(as.character(e$pgx$datatype)[1]), "methylomics")) {
      pgx_file <- f
      break
    }
  }
}
if (!nzchar(pgx_file) || !file.exists(pgx_file)) {
  stop("No methylomics .pgx found. Build one with test/build-test-pgx.Rscript, ",
       "or point METHYLOME_PGX at a file.", call. = FALSE)
}
message("[test] loading ", basename(pgx_file))
.e <- new.env(); load(pgx_file, envir = .e)
PGX_STATIC <- .e$pgx
message("[test] ", PGX_STATIC$name, ": ", nrow(PGX_STATIC$X), " probes x ",
        ncol(PGX_STATIC$X), " samples",
        if (!is.null(PGX_STATIC$meth$cells)) {
          paste0(", ", length(PGX_STATIC$meth$cells$counts), " stored cell panels")
        } else ", no stored cell composition")

ui <- bslib::page_fluid(
  padding = 0,
  theme = bigdash::big_theme(),
  bigdash::dependencies(),
  shiny::tags$head(shiny::tags$link(rel = "stylesheet", href = "custom/styles.min.css")),
  shinyjs::useShinyjs(),
  div(class = "omicspanel fullheight-page", style = "height:100vh;",
      methylome_ui("methylome"))
)

server <- function(input, output, session) {
  ## The app hands the module a reactiveValues PGX; a plain list is accepted
  ## too (methylome_server wraps it), which is all a static harness needs.
  methylome_server("methylome", pgx = PGX_STATIC)
}

shinyApp(ui, server)
