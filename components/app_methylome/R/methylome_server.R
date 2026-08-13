##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Methylome Profiler server. Takes the loaded PGX like every other board;
## the app owns no upload path of its own.

methylome_server <- function(id = "methylome", pgx) {
  shiny::moduleServer(id, function(input, output, session) {
    ## Accept either a reactive or a plain pgx so the standalone harness in
    ## ../test/app.R can pass a readRDS() result directly.
    PGX <- if (shiny::is.reactive(pgx)) pgx else shiny::reactive(pgx)

    methylome_ledger_server("ledger", PGX)
    methylome_age_server("age", PGX)
    methylome_character_server("character", PGX)
    methylome_ewas_server("ewas", PGX)
  })
}
