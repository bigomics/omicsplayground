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

    ## The app hands us a reactiveValues, which is not itself reactive: PGX()
    ## returns the same object every time and never invalidates. Reading its
    ## fields inside a reactive does create a dependency, so derive an explicit
    ## change signal and trigger off that instead - otherwise the pickers are
    ## populated once, before any dataset exists, and never again.
    r_dataset <- shiny::reactive({
      p <- PGX()
      shiny::req(p)
      paste(p$name, ncol(p$X), paste(colnames(p$contrasts), collapse = ","))
    })

    ## ---------------------------------------------------------- choices --
    ## Populate the model pickers from whatever the loaded dataset offers.
    shiny::observeEvent(r_dataset(), {
      p <- PGX()
      shiny::req(p)
      cmps <- colnames(p$contrasts)
      if (is.null(cmps) && !is.null(p$gx.meta)) cmps <- names(p$gx.meta$meta)
      shiny::updateSelectInput(session, "ewas_contrast", choices = cmps,
                               selected = cmps[1])
      shiny::updateSelectizeInput(session, "ewas_covars",
                                  choices = mp_model_vars(p), selected = character(0),
                                  server = TRUE)
      phe <- mp_model_vars(p)
      shiny::updateSelectInput(session, "comp_pheno", choices = phe, selected = phe[1])
    })

    ## ----------------------------------------------------- cell fractions --
    ## Deconvolution is a projection over ~500 CpGs but still not instant, and
    ## it is a prerequisite for the adjusted model, so it is explicit.
    applied_ref <- shiny::reactiveVal(NULL)
    r_cells <- shiny::eventReactive(
      list(input$run_deconv, r_dataset()),
      {
        shiny::req(PGX())
        applied_ref(input$deconv_ref)
        mp_cell_counts(PGX(), if (is.null(input$deconv_ref)) MP_DECONV_REFS[1] else input$deconv_ref)
      },
      ignoreNULL = FALSE
    )
    output$deconv_stale <- shiny::renderUI({
      a <- applied_ref()
      if (is.null(a) || identical(a, input$deconv_ref)) return(NULL)
      shiny::div(style = "margin-top:6px; font-size:11.5px; color:#8a5a06;",
                 "Reference changed - press Estimate composition to apply.")
    })

    ## ------------------------------------------------------- EWAS model --
    r_ewas <- shiny::eventReactive(
      list(input$run_ewas, r_dataset()),
      {
        p <- PGX()
        shiny::req(p)
        cmp <- input$ewas_contrast
        if (is.null(cmp)) {
          cmp <- colnames(p$contrasts)[1]
          if (is.null(cmp) && !is.null(p$gx.meta)) cmp <- names(p$gx.meta$meta)[1]
        }
        cf <- if (isTRUE(input$ewas_adjust_cells)) r_cells() else NULL
        fit <- mp_fit_ewas(p, cmp, covars = input$ewas_covars,
                           cellfracs = cf, mask = input$ewas_mask)
        d <- mp_ewas_annotate(p, fit$table)
        list(data = d, contrast = cmp, meta = fit,
             bacon = mp_bacon(stats::qnorm(fit$table$p / 2) *
                              sign(fit$table$logFC_M)))
      },
      ignoreNULL = FALSE
    )

    ## The fitted model, printed so the user can see what was actually tested.
    output$ewas_model <- shiny::renderUI({
      m <- r_ewas()$meta
      shiny::div(
        style = "margin-top:8px; font-size:11.5px; color:#697586; line-height:1.5;",
        shiny::tags$b("Model: "), m$formula, shiny::tags$br(),
        sprintf("n = %d  (%s)", m$n,
                paste(names(m$groups), m$groups, collapse = " / ")),
        shiny::tags$br(),
        sprintf("%s probes masked", format(m$masked, big.mark = ",")),
        if (length(m$dropped_covars)) {
          shiny::tags$div(style = "color:#8a5a06;",
            paste("dropped:", paste(m$dropped_covars, collapse = ", ")))
        }
      )
    })

    r_contrast <- shiny::reactive(r_ewas()$contrast)

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
      list(input$recompute_clocks, r_dataset()),
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

    methylome_plot_composition_server("comp_bars", PGX, r_cells, watermark = watermark)
    methylome_table_composition_server("comp_tbl", r_cells)
    methylome_plot_compgroup_server("comp_group", PGX, r_cells,
                                    shiny::reactive(input$comp_pheno),
                                    watermark = watermark)

    methylome_plot_manhattan_server("ewas_manhattan", r_ewas, r_thresh,
                                    watermark = watermark)
    methylome_plot_qq_server("ewas_qq", r_ewas, watermark = watermark)
    methylome_plot_enrichment_server("ewas_enrich", r_ewas, r_thresh,
                                     watermark = watermark)
    methylome_table_hits_server("ewas_hits", r_ewas, r_thresh)
    methylome_plot_stripcharts_server("ewas_strips", PGX, r_ewas, r_thresh,
                                      r_topn, watermark = watermark)
  })
}
