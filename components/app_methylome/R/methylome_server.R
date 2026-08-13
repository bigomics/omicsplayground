##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Methylome Profiler server. Takes the loaded PGX like every other board; the
## app owns no upload path of its own.

methylome_server <- function(id = "methylome", pgx, watermark = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    ## Accept a reactive, a reactiveValues PGX, or a plain pgx list so the
    ## standalone harness can pass a readRDS() result directly.
    PGX <- if (shiny::is.reactive(pgx)) pgx else shiny::reactive(pgx)

    ## Only contrast-driven screen; everything else is per-sample.
    r_contrast <- shiny::reactive({
      p <- PGX()
      if (is.null(p$gx.meta)) NULL else names(p$gx.meta$meta)[1]
    })

    methylome_table_ledger_server("ledger_tbl", PGX)
    methylome_plot_betadist_server("ledger_dens", PGX, watermark = watermark)

    methylome_plot_clocks_server("age_clocks", PGX, watermark = watermark)
    methylome_plot_agegroup_server("age_group", PGX, watermark = watermark)
    methylome_table_coverage_server("age_cov", PGX)

    methylome_plot_context_server("char_island", PGX, what = "island", watermark = watermark)
    methylome_plot_context_server("char_gene", PGX, what = "gene", watermark = watermark)

    methylome_plot_manhattan_server("ewas_manhattan", PGX, r_contrast, watermark = watermark)
    methylome_plot_qq_server("ewas_qq", PGX, r_contrast, watermark = watermark)
    methylome_plot_enrichment_server("ewas_enrich", PGX, r_contrast, watermark = watermark)
  })
}
