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
      ## One picker, two groups: a two-group contrast or a continuous exposure
      ## are the same model to limma, so they are the same choice to the user.
      cont <- setdiff(mp_continuous_vars(p), cmps)
      ch <- list()
      if (length(cmps)) ch[["Contrasts"]] <- cmps
      if (length(cont)) ch[["Continuous variables"]] <- cont
      shiny::updateSelectInput(session, "ewas_contrast", choices = ch,
                               selected = if (length(cmps)) cmps[1] else cont[1])
      shiny::updateSelectizeInput(session, "ewas_covars",
                                  choices = mp_model_vars(p), selected = character(0),
                                  server = TRUE)
      phe <- mp_model_vars(p)
      shiny::updateSelectInput(session, "comp_pheno", choices = phe, selected = phe[1])
      ## Acceleration is split by group, so only the categorical columns.
      cat_vars <- Filter(function(k) {
        u <- unique(as.character(p$samples[[k]]))
        u <- u[!is.na(u) & u != ""]
        length(u) >= 2 && length(u) <= 8
      }, phe)
      shiny::updateSelectInput(session, "acc_pheno", choices = cat_vars,
                               selected = cat_vars[1])
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
    ## The contrast picker is populated by an observer, so it is still NULL on
    ## the first pass; include it in the trigger or the initial fit runs with
    ## no contrast and the failure sticks. Refitting on a contrast change is
    ## also the right behaviour - it is a different question, not a tweak.
    r_ewas <- shiny::eventReactive(
      list(input$run_ewas, r_dataset(), input$ewas_contrast),
      {
        p <- PGX()
        shiny::req(p)
        cmp <- input$ewas_contrast
        if (is.null(cmp) || !nzchar(cmp)) {
          cmp <- colnames(p$contrasts)[1]
          if (is.null(cmp) && !is.null(p$gx.meta)) cmp <- names(p$gx.meta$meta)[1]
        }
        ## Wait rather than fail if nothing is resolvable yet.
        shiny::req(!is.null(cmp) && nzchar(cmp))
        ## Refuse rather than fail open. Asking for cell adjustment and
        ## silently getting an unadjusted model is the one wrong answer this
        ## board exists to prevent - and it is invisible, because the formula
        ## line simply does not mention the fractions.
        cf <- NULL
        if (isTRUE(input$ewas_adjust_cells)) {
          cf <- r_cells()
          shiny::validate(shiny::need(!is.null(cf) && nrow(cf) > 0,
            "Cell-composition adjustment was requested but no proportions are available. Estimate them on the Cell composition screen, or untick the adjustment."))
        }
        fit <- mp_fit_ewas(p, cmp, covars = input$ewas_covars,
                           cellfracs = cf, mask = input$ewas_mask)
        d <- mp_ewas_annotate(p, fit$table)
        list(data = d, contrast = cmp, meta = fit,
             ## The moderated t is already a signed test statistic. Rebuilding
             ## one as qnorm(p/2)*sign(logFC) inverted every z - qnorm(p/2) is
             ## always negative - which flipped the sign of the reported bias.
             bacon = mp_bacon(fit$t))
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
                if (is.null(m$groups)) m$desc
                else paste(names(m$groups), m$groups, collapse = " / ")),
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

    ## Only clocks that survived the coverage floor can be the acceleration
    ## clock, and the list changes with every recompute - keep the current
    ## choice when it is still usable rather than snapping back to the first.
    shiny::observeEvent(r_clockset(), {
      cl <- r_clockset()
      if (is.null(cl) || !ncol(cl$age)) return()
      cur <- shiny::isolate(input$acc_clock)
      shiny::updateSelectInput(
        session, "acc_clock", choices = colnames(cl$age),
        selected = if (!is.null(cur) && cur %in% colnames(cl$age)) cur else colnames(cl$age)[1]
      )
    })

    methylome_plot_agecor_server("age_cor", r_clockset, watermark = watermark)
    methylome_plot_clocks_server("age_clocks", r_clockset, watermark = watermark)
    methylome_plot_agegroup_server("age_group", PGX, r_clockset,
                                   shiny::reactive(input$acc_pheno),
                                   shiny::reactive(input$acc_clock),
                                   r_cells,
                                   shiny::reactive(isTRUE(input$acc_intrinsic)),
                                   watermark = watermark)
    methylome_table_coverage_server("age_cov", r_clockset)

    methylome_plot_context_server("char_island", PGX, what = "island", watermark = watermark)
    methylome_plot_context_server("char_gene", PGX, what = "gene", watermark = watermark)

    ## ------------------------------------------------- methylome landscape --
    ## Moved from the Epigenomics board. Its plot modules read pgx$X directly
    ## rather than taking a reactive, so they get the object this board was
    ## handed, not PGX() - passing the evaluated value would freeze them on
    ## whichever dataset happened to be loaded when the board was wired.
    shiny::observeEvent(r_dataset(), {
      p <- PGX()
      shiny::req(p$X, p$samples, p$genes)
      kk <- grep("chr|chromosome|chromosomes|chrom|chroms", tolower(colnames(p$genes)))
      chroms <- unique(stats::na.omit(sub("(p|q|cen).*", "", as.character(p$genes[, kk[1]]))))
      chroms <- unique(sub("^chr", "", chroms))
      sex.chr <- intersect(c("X", "Y"), chroms)
      chroms <- paste0("chr", c(sort(as.numeric(setdiff(chroms, sex.chr))), sex.chr))
      shiny::updateSelectizeInput(session, "select_chromosome", choices = chroms,
                                  selected = utils::head(chroms, 4), server = TRUE)
      shiny::updateSelectInput(session, "data_samplefilter", choices = colnames(p$X))
      Y <- p$samples
      pheno <- colnames(Y)[sapply(Y, function(v) any(!is.na(v) & v != ""))]
      pheno <- pheno[!grepl("cell_cycle", pheno, ignore.case = TRUE)]
      shiny::updateSelectInput(session, "select_pheno",
                               choices = c("<ungrouped>", pheno), selected = "<ungrouped>")
    })

    ## The ideogram is drawn per chromosome, so it is capped; the boxplot just
    ## gets narrower columns and is not.
    ## Autosomes numerically, then X and Y. as.numeric() alone turned "X" and
    ## "Y" into NA, so selecting a sex chromosome silently drew nothing.
    mp_sort_chroms <- function(ch) {
      ch <- unique(sub("^chr", "", as.character(ch)))
      num <- suppressWarnings(as.numeric(ch))
      c(as.character(sort(num[!is.na(num)])), intersect(c("X", "Y"), ch))
    }
    r_chroms_capped <- shiny::reactive({
      ch <- mp_sort_chroms(input$select_chromosome)
      shiny::validate(
        shiny::need(length(ch) > 0, "Select at least one chromosome in the settings panel."),
        shiny::need(length(ch) < 7, "Select at most six chromosomes at a time - the ideogram becomes unreadable beyond that.")
      )
      ch
    })
    r_chroms_all <- shiny::reactive({
      ch <- mp_sort_chroms(input$select_chromosome)
      shiny::validate(shiny::need(length(ch) > 0,
        "Select at least one chromosome in the settings panel."))
      ch
    })
    r_land_samples <- shiny::reactive({
      ss <- rownames(PGX()$samples)
      drop <- input$data_samplefilter
      if (!is.null(drop)) {
        drop <- intersect(ss, drop[!is.na(drop) & drop != ""])
        if (length(drop)) ss <- setdiff(ss, drop)
      }
      shiny::validate(shiny::need(length(ss) > 0, "No samples remaining after filtering."))
      ss
    })
    r_land_pheno <- shiny::reactive(input$select_pheno)

    ## The landscape screen's "Global beta distribution" was a second beta
    ## density, of the per-CpG mean across samples. It has been dropped in
    ## favour of the per-sample one on the Sample ledger: averaging across
    ## samples pulls the two modes toward the middle, and collapsing samples
    ## makes the panel unable to answer the question a beta density is for -
    ## did any one sample fail. It also ignored the phenotype selector while
    ## documenting that it split by group.
    dataview_table_beta_server("methyltable", pgx, r.samples = r_land_samples,
                               r.pheno = r_land_pheno, scrollY = "30vh")
    epigenomics_plot_boxplot_beta_server("boxplotBeta", pgx,
                                         r.chromosome = r_chroms_all,
                                         r.samples = r_land_samples,
                                         r.pheno = r_land_pheno, watermark = watermark)
    epigenomics_plot_methylIdeogram_server("methylIdeogram", pgx,
                                           r.chromosome = r_chroms_capped,
                                           r.samples = r_land_samples,
                                           r.pheno = r_land_pheno, watermark = watermark)

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

    methylome_plot_composition_server("comp_bars", PGX, r_cells,
                                     shiny::reactive(input$comp_pheno),
                                     watermark = watermark)
    methylome_table_composition_server("comp_tbl", r_cells)
    methylome_plot_compgroup_server("comp_group", PGX, r_cells,
                                    shiny::reactive(input$comp_pheno),
                                    watermark = watermark)

    methylome_plot_manhattan_server("ewas_manhattan", r_ewas, r_thresh,
                                    watermark = watermark)
    methylome_plot_volcano_server("ewas_volcano", r_ewas, r_thresh,
                                  watermark = watermark)
    methylome_plot_hitmap_server("ewas_hitmap", PGX, r_ewas, r_thresh,
                                 watermark = watermark)
    methylome_plot_qq_server("ewas_qq", r_ewas, watermark = watermark)
    methylome_plot_enrichment_server("ewas_enrich", r_ewas, r_thresh,
                                     watermark = watermark)
    methylome_table_hits_server("ewas_hits", r_ewas, r_thresh)
    methylome_plot_stripcharts_server("ewas_strips", PGX, r_ewas, r_thresh,
                                      r_topn, watermark = watermark)

    ## ------------------------------------------- regions and gene sets --
    ## Both are minutes-scale on a full array, so both are explicit.
    ## reactiveVal rather than eventReactive: an eventReactive that has never
    ## fired simply does not run, so its table renders blank instead of saying
    ## what to press. This way the table gets NULL and shows the prompt.
    regions_val <- shiny::reactiveVal(NULL)
    shiny::observeEvent(input$run_dmr, {
      p <- PGX(); res <- r_ewas()
      gap <- input$dmr_maxgap
      if (is.null(gap) || is.na(gap) || gap < 1) gap <- 500
      shiny::withProgress(message = "Calling regions...", value = 0.4, {
        dmrs <- mp_call_dmrs(p, res, maxgap = gap)
        regions_val(list(dmrs = dmrs, genes = mp_dmr_genes(p, dmrs, res$data)))
      })
    })
    ## A new model invalidates any regions already called.
    shiny::observeEvent(r_ewas(), regions_val(NULL), ignoreInit = TRUE)
    r_regions <- shiny::reactive(regions_val())

    ## Clicking a region in the table is what selects it for the detail plot -
    ## the region is already on screen there, so a second copy of the same list
    ## in the settings panel was a picker for something the user can just point
    ## at. TableModuleServer hands back rows_selected for exactly this.
    dmr_tbl <- methylome_table_dmr_server("ewas_dmr", PGX, r_regions)
    methylome_plot_dmrregion_server("ewas_dmrplot", PGX, r_ewas, r_regions,
                                    dmr_tbl$rows_selected,
                                    watermark = watermark)

    enrich_val <- shiny::reactiveVal(NULL)
    shiny::observeEvent(input$run_gometh, {
      p <- PGX(); d <- r_ewas()$data
      sig <- mp_ewas_sig(d, r_thresh())
      arr <- mp_array_type(p)
      shiny::withProgress(message = "Testing gene sets...", value = 0.4, {
        enrich_val(mp_run_gometh(
          d$probe[sig], d$probe,
          collection = if (is.null(input$gs_collection)) "GO" else input$gs_collection,
          array_type = arr))
      })
    })
    shiny::observeEvent(r_ewas(), enrich_val(NULL), ignoreInit = TRUE)
    r_enrich <- shiny::reactive(enrich_val())
    methylome_table_enrich_server("ewas_enrichgs", r_enrich)
  })
}
