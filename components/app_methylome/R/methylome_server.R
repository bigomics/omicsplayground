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

    ## Clock selection: the family checkboxes drive the individual list, and
    ## the individual list is what is actually computed - so a user can tick a
    ## family and then drop one clock out of it.
    shiny::observeEvent(input$clock_families, {
      sel <- unlist(MP_CLOCK_FAMILIES[input$clock_families], use.names = FALSE)
      shiny::updateCheckboxGroupInput(session, "clocks", selected = sel)
    }, ignoreNULL = FALSE)

    r_clocks <- shiny::reactive({
      cl <- input$clocks
      if (is.null(cl) || !length(cl)) MP_CLOCK_ALL else cl
    })
    r_mincov <- shiny::reactive({
      v <- input$min_cov
      if (is.null(v) || is.na(v)) 0.8 else v
    })

    ## methylclock takes ~20s on a full cohort, so the clock set is computed
    ## once and shared across the four panels, and only on demand: reacting to
    ## every checkbox tick would refit ten clocks each time. Loading a new
    ## dataset also triggers it, otherwise the panels would keep showing the
    ## previous dataset's ages.
    applied <- shiny::reactiveVal(NULL)
    r_clockset <- shiny::eventReactive(
      list(input$recompute_clocks, PGX()),
      {
        shiny::req(PGX())
        a <- list(clocks = r_clocks(), min_cov = r_mincov())
        applied(a)
        mp_clock_set(PGX(), a$clocks, a$min_cov)
      },
      ignoreNULL = FALSE
    )

    ## Tell the user when the panels are showing something older than the
    ## current settings, rather than letting them wonder if it worked.
    output$clock_stale <- shiny::renderUI({
      a <- applied()
      if (is.null(a)) return(NULL)
      stale <- !identical(sort(a$clocks), sort(r_clocks())) ||
        !isTRUE(all.equal(a$min_cov, r_mincov()))
      if (!stale) return(NULL)
      shiny::div(
        style = "margin-top:6px; font-size:11.5px; color:#8a5a06;",
        "Settings changed - press Recompute to apply."
      )
    })

    methylome_plot_agecor_server("age_cor", r_clockset, watermark = watermark)
    methylome_plot_clocks_server("age_clocks", r_clockset, watermark = watermark)
    methylome_plot_agegroup_server("age_group", PGX, r_clockset, watermark = watermark)
    methylome_table_coverage_server("age_cov", r_clockset)

    methylome_plot_context_server("char_island", PGX, what = "island", watermark = watermark)
    methylome_plot_context_server("char_gene", PGX, what = "gene", watermark = watermark)

    ## Threshold shared by the Manhattan line, the context enrichment and the
    ## hit table. Defaults hold until the settings panel has initialised.
    r_thresh <- shiny::reactive({
      v <- input$thresh_value
      if (is.null(v) || is.na(v) || v <= 0 || v > 1) v <- 0.05
      db <- input$min_dbeta
      if (is.null(db) || is.na(db) || db < 0) db <- 0
      list(
        type = if (is.null(input$thresh_type)) "q" else input$thresh_type,
        value = v,
        min_dbeta = db
      )
    })

    r_topn <- shiny::reactive({
      n <- input$top_n
      if (is.null(n) || is.na(n) || n < 1) 6 else min(as.integer(n), 12)
    })

    methylome_plot_manhattan_server("ewas_manhattan", PGX, r_contrast, r_thresh,
                                    watermark = watermark)
    methylome_plot_qq_server("ewas_qq", PGX, r_contrast, watermark = watermark)
    methylome_plot_enrichment_server("ewas_enrich", PGX, r_contrast, r_thresh,
                                     watermark = watermark)
    methylome_table_hits_server("ewas_hits", PGX, r_contrast, r_thresh)
    methylome_plot_stripcharts_server("ewas_strips", PGX, r_contrast, r_thresh,
                                      r_topn, watermark = watermark)
  })
}
