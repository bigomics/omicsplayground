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
      ## One picker, three groups: a two-group contrast, a continuous exposure
      ## and a multi-level factor are all the same machinery to limma, so they
      ## are the same choice to the user. Without the third group a three-level
      ## phenotype is unreachable - playbase builds contrasts two groups at a
      ## time, so nothing in the first two groups can express it.
      cont <- setdiff(mp_continuous_vars(p), cmps)
      cat <- setdiff(mp_categorical_vars(p), cmps)
      ch <- list()
      if (length(cmps)) ch[["Contrasts"]] <- cmps
      if (length(cont)) ch[["Continuous variables"]] <- cont
      if (length(cat)) ch[["Categorical variables"]] <- cat
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

      ## The clustering scatter is copied from the Dashboard, where its two
      ## pickers are filled by the BOARD server rather than by the module. Copy
      ## that too, or both selectInputs render empty and the panel never draws.
      ##
      ## Filled from this app's own helper, not playbase::pgx.getCategoricalPhenotypes:
      ## that one apply()s over the sample sheet and dies with "dim(X) must have
      ## a positive length" on an empty one, which is exactly what this observer
      ## sees before a dataset is loaded - and an unhandled error in an observer
      ## severs the whole session, not just the panel.
      grp <- mp_categorical_vars(p, max_levels = 12, min_per_level = 1)
      grp <- grp[grep("sample|patient", grp, invert = TRUE)]
      shiny::updateSelectInput(session, "overview_pca-hmpca.colvar",
                               choices = c(grp, "<none>"),
                               selected = if (length(grp)) grp[1] else "<none>")
      shiny::updateSelectInput(session, "overview_pca-hmpca.shapevar",
                               choices = c("<none>", grp), selected = "<none>")
    })

    ## ----------------------------------------------------- cell fractions --
    ## Every reference panel is fitted at dataset creation and stored in the
    ## pgx, so the usual path is a lookup: the screen is populated on load and
    ## a panel switch is instant. The button stays for the one case that still
    ## costs something - a pgx built before that slot existed, where
    ## mp_cell_counts() falls back to fitting live. That fit still
    ## quantile-normalises the whole matrix and runs a QP per sample, tens of
    ## seconds on a real cohort, which is why it is never dataset-triggered:
    ## as a standalone app the UI is in the DOM from page load and nothing
    ## suspends it, so a dataset-driven fit would run on every load.
    applied_ref <- shiny::reactiveVal(NULL)
    cells_val <- shiny::reactiveVal(NULL)
    ## Free when the pgx carries the fit, NULL when it does not - so this is
    ## safe to run on every dataset change and on the panel picker.
    take_stored <- function(ref) {
      pgx <- PGX()
      if (is.null(pgx)) return(FALSE)
      cc <- mp_stored_cells(pgx, if (is.null(ref)) MP_DECONV_REFS[1] else ref)
      if (is.null(cc)) return(FALSE)
      cells_val(cc)
      applied_ref(if (is.null(ref)) MP_DECONV_REFS[1] else ref)
      TRUE
    }
    shiny::observeEvent(input$run_deconv, {
      shiny::req(PGX())
      applied_ref(input$deconv_ref)
      shiny::withProgress(
        message = "Estimating cell composition...", value = 0.4,
        cells_val(mp_cell_counts(
          PGX(), if (is.null(input$deconv_ref)) MP_DECONV_REFS[1] else input$deconv_ref))
      )
    })
    ## Switching panel costs nothing when it was precomputed; when it was not,
    ## leave what is showing and let deconv_stale tell the user to press the
    ## button, which is the pre-existing behaviour.
    shiny::observeEvent(input$deconv_ref, {
      take_stored(input$deconv_ref)
    }, ignoreInit = TRUE)
    shiny::observeEvent(r_dataset(), {
      cells_val(NULL)
      applied_ref(NULL)
      take_stored(input$deconv_ref)
    }, ignoreInit = TRUE)
    ## The first dataset arrives before any observer above has fired.
    shiny::observeEvent(PGX(), {
      if (is.null(cells_val())) take_stored(input$deconv_ref)
    }, once = FALSE)
    r_cells <- shiny::reactive(cells_val())
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
    ## Also button-driven. The contrast picker is repopulated on every dataset
    ## load, so having input$ewas_contrast in the trigger fired a full limma
    ## fit the moment a dataset arrived - masking, betaToM over every probe,
    ## eBayes - before the user had asked for anything.
    r_ewas_fit <- shiny::eventReactive(
      input$run_ewas,
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
            "Cell-composition adjustment was requested but no proportions are available. Press Estimate composition on the Cell composition screen, or untick the adjustment."))
        }
        ## sva is an SVD plus IRW iterations over the whole probe matrix -
        ## seconds on a subset, minutes on a full array - so the fit gets a
        ## progress bar rather than an unexplained pause.
        ## Starts low and is driven by the fit's own per-chromosome ticks, so
        ## the bar reports where it actually is rather than sitting at 40%.
        fit <- shiny::withProgress(
          message = "Fitting model...", value = 0.05, detail = "preparing",
          mp_fit_ewas(p, cmp, covars = input$ewas_covars,
                      cellfracs = cf, mask = input$ewas_mask,
                      n_sv = if (isTRUE(input$ewas_sva)) -1L else 0L,
                      collapse = if (is.null(input$ewas_collapse)) "gene_region"
                                 else input$ewas_collapse))
        d <- mp_ewas_annotate(p, fit$table, annot = fit$collapse_annot)
        res <- list(data = d, contrast = cmp, meta = fit,
             ## The moderated t is already a signed test statistic. Rebuilding
             ## one as qnorm(p/2)*sign(logFC) inverted every z - qnorm(p/2) is
             ## always negative - which flipped the sign of the reported bias.
             ## bacon models a signed z; an F is not one, so the anova branch
             ## gets no bias estimate rather than a confident wrong one.
             bacon = if (isTRUE(fit$anova)) NULL else mp_bacon(fit$t))
        res
      },
      ## Nothing until the button is actually pressed. ignoreNULL = FALSE here
      ## would fire the fit at startup, which is the behaviour being removed.
      ignoreNULL = TRUE, ignoreInit = TRUE
    )

    ## An eventReactive that has never fired is silent, so every panel on this
    ## screen rendered as an empty card with no hint that a button was waiting.
    ## A flag in front of it turns that into one prompt, shared by all of them.
    ##
    ## Deliberately a thin wrapper and not a reactiveVal holding the result:
    ## mp_fit_ewas() refuses a rank-deficient design, a missing cell-fraction
    ## estimate and too few complete samples with validate(), and those refusals
    ## only reach the panels because they propagate out of a reactive. Moving
    ## the fit into an observer would swallow every one of them.
    ewas_asked <- shiny::reactiveVal(FALSE)
    shiny::observeEvent(input$run_ewas, ewas_asked(TRUE))
    ## A new dataset invalidates the fit that was made against the old one.
    shiny::observeEvent(r_dataset(), ewas_asked(FALSE), ignoreInit = TRUE)

    r_ewas <- shiny::reactive({
      shiny::validate(shiny::need(ewas_asked(), MP_NO_EWAS))
      r_ewas_fit()
    })

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

    ## The clocks now come with the dataset, so the panels fill in on load and
    ## both settings apply instantly - they filter a fit that used neither.
    ##
    ## A pgx built before the slot existed has nothing to filter, and fitting
    ## ten clocks takes ~40s on a full cohort. That used to run on every
    ## dataset load: as a board the UI was inserted lazily so nothing read it
    ## until you opened this screen, but as an app the whole UI is in the DOM
    ## from page load, so it fired every time - and R being single-threaded,
    ## any click during those 40s queues behind it, which is what made the
    ## Library's own loading modal look like it took forever to appear. So the
    ## fallback stays behind the button, and its raw fit is kept in a
    ## reactiveVal that a new dataset clears rather than letting the panels
    ## show the previous cohort's ages.
    live_val <- shiny::reactiveVal(NULL)
    shiny::observeEvent(input$recompute_clocks, {
      shiny::req(PGX())
      ## Nothing to do when the dataset already carries its clocks - pressing
      ## the button would refit for 40s and land on the same numbers.
      if (!is.null(mp_stored_clocks(PGX()))) return()
      shiny::withProgress(
        message = "Computing epigenetic clocks...", value = 0.4,
        live_val(playbase.epigenetics::compute_clocks(mp_beta(PGX())))
      )
    })
    shiny::observeEvent(r_dataset(), live_val(NULL), ignoreInit = TRUE)

    r_clockset <- shiny::reactive({
      p <- PGX()
      shiny::req(p)
      ## The stored path never touches mp_beta(), so the datatype check has to
      ## happen here or an RNA-seq dataset would be told to press a button.
      mp_require_methylomics(p)
      cl <- mp_stored_clocks(p)
      if (!is.null(cl)) cl$stored <- TRUE else cl <- live_val()
      mp_clock_filter(cl, r_clocks(), r_mincov(), mp_chronological_age(p))
    })

    ## Only clocks that survived the coverage floor can be the acceleration
    ## clock, and the list changes with the floor - keep the current
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
    ## The "Methylation table" that used to sit here is gone: it was the
    ## chromosome boxplot with the dispersion thrown away - same matrix, same
    ## annotation, same centromere filter - and its own help conceded it was
    ## "read for what stands out rather than for its values", which is what a
    ## boxplot does and a column of means cannot. Past ~20 samples it emitted
    ## one column per sample and was unreadable.
    ## Components, loadings and phenotype directions for the scatter. One
    ## decomposition per dataset, shared by the points and by both arrow modes.
    ## Upstream builds this in the Clustering board server and passes it in; the
    ## module takes it as an argument, so the board decides the scale - M here.
    r_pca_components <- shiny::reactive({
      p <- PGX(); shiny::req(p)
      ## 1,000 most variable by default; the settings toggle drops the cut.
      rsd <- if (isTRUE(input$structure_allfeatures)) 0L else 1000L
      shiny::withProgress(message = "Computing components...", value = 0.5,
                          mp_pca_components(p, reduce.sd = rsd))
    })

    ## Raw pgx, like every other panel here: the modules read pgx$X and
    ## pgx$cluster directly, and the app hands the board a reactiveValues.
    methylome_plot_clustpca_server(
      "overview_pca", pgx,
      selected_samples = shiny::reactive({
        p <- PGX(); shiny::req(p); colnames(p$X)
      }),
      clustmethod = shiny::reactive("pca"),
      pca_components = r_pca_components,
      ## PC1 vs PC2, fixed. Upstream lets the user pick any pair of PC1-PC5 from
      ## the board sidebar; this app has no such sidebar, and the heatmap beside
      ## the scatter already reports every component against every variable.
      pca_dims = shiny::reactive(c(1L, 2L)),
      ## Feature arrows are CpGs; the gene symbol is what makes one readable.
      labeltype = shiny::reactive("symbol"),
      watermark = watermark, parent = session$ns)
    ## Same components the scatter draws - one decomposition, two panels.
    methylome_plot_pccovar_server("overview_pccov", PGX,
                                  pca_components = r_pca_components,
                                  watermark = watermark)

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

    ## EWAS Catalog rows for the current hit list, read once and shared by the
    ## hits table and the two EWAS Catalog panels. Each read materialises
    ## ~129 MB for ~0.7s, so three consumers reading it themselves would be
    ## three reads per threshold nudge.
    ## The catalog is indexed by CpG id, so it has nothing to say about a fit
    ## collapsed to genes or gene regions. That is a fact about the catalog, not
    ## a reason to withhold the hit list: this returns `ok = FALSE` with the
    ## reason and lets each consumer decide. The dedicated catalog panels refuse;
    ## the hits table simply drops its two catalog columns and still shows the
    ## hits, which is the whole point of the screen.
    MP_CATALOG_NA <- paste(
      "The EWAS Catalog is indexed by CpG, so it has nothing to say about a fit",
      "collapsed to genes or gene regions. Set 'Test at' to probe level to use it."
    )
    r_catalog <- shiny::reactive({
      res <- r_ewas()
      if (!identical(res$meta$collapse, "probe")) {
        return(list(ok = FALSE, reason = MP_CATALOG_NA))
      }
      d <- res$data
      sig <- mp_ewas_sig(d, r_thresh())
      probes <- d$probe[sig]
      if (!length(probes)) {
        return(list(ok = FALSE, reason = "No CpG passes the current threshold."))
      }
      ## playdata::get_file() is system.file(mustWork = TRUE), so a playdata
      ## built before the catalog landed does not return "" - it throws, and
      ## the throw would surface as a red error on the main hits table rather
      ## than on the catalog panels. Degrade to "unavailable" instead.
      tr <- tryCatch(playbase.epigenetics::ewas_catalog_traits(probes),
                     error = function(e) e)
      if (inherits(tr, "error")) {
        return(list(ok = FALSE, reason = paste(
          "The bundled EWAS Catalog is not available in this installation:",
          conditionMessage(tr))))
      }
      list(ok = TRUE, probes = probes, traits = tr)
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
    methylome_plot_traitfreq_server("ewas_traitfreq", r_catalog,
                                    watermark = watermark)
    methylome_table_hits_server("ewas_hits", r_ewas, r_thresh, r_catalog)
    methylome_table_traits_server("ewas_traits", r_catalog)
    methylome_plot_stripcharts_server("ewas_strips", PGX, r_ewas, r_thresh,
                                      r_topn, watermark = watermark)

    ## ------------------------------------------- regions and gene sets --
    ## Both are minutes-scale on a full array, so both are explicit.
    ## reactiveVal rather than eventReactive: an eventReactive that has never
    ## fired simply does not run, so its table renders blank instead of saying
    ## what to press. This way the table gets NULL and shows the prompt.
    regions_val <- shiny::reactiveVal(NULL)
    shiny::observeEvent(input$run_dmr, mp_on_click({
      p <- PGX(); res <- r_ewas()
      gap <- input$dmr_maxgap
      if (is.null(gap) || is.na(gap) || gap < 1) gap <- 500
      shiny::withProgress(message = "Calling regions...", value = 0.4, {
        dmrs <- mp_call_dmrs(p, res, maxgap = gap)
        regions_val(list(dmrs = dmrs, genes = mp_dmr_genes(p, dmrs, res$data)))
      })
    }))
    ## A new model invalidates any regions already called.
    shiny::observeEvent(r_ewas_fit(), regions_val(NULL), ignoreInit = TRUE)
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
    shiny::observeEvent(input$run_gometh, mp_on_click({
      p <- PGX(); res <- r_ewas(); d <- res$data
      ## gometh's whole point is reweighting by probes per gene, and it takes
      ## CpG ids to do it. A collapsed fit has already averaged those probes
      ## away, so the correction has nothing to correct and the ids would not
      ## resolve either.
      shiny::validate(shiny::need(identical(res$meta$collapse, "probe"),
        "Bias-corrected gene-set testing needs per-probe CpG ids, and this model was fitted on collapsed features. Set 'Test at' back to probe level."))
      sig <- mp_ewas_sig(d, r_thresh())
      arr <- mp_array_type(p)
      shiny::withProgress(message = "Testing gene sets...", value = 0.4, {
        enrich_val(mp_run_gometh(
          d$probe[sig], d$probe,
          collection = if (is.null(input$gs_collection)) "GO" else input$gs_collection,
          array_type = arr))
      })
    }))
    shiny::observeEvent(r_ewas_fit(), enrich_val(NULL), ignoreInit = TRUE)
    r_enrich <- shiny::reactive(enrich_val())
    methylome_table_enrich_server("ewas_enrichgs", r_enrich)

    ## Promoter vs body. Region calling and gene-set testing are both per-probe
    ## answers and refuse on a collapsed fit, which left this sub-tab dead there;
    ## this is the question only the collapsed fit can answer, so it takes their
    ## place. Both directions refuse cleanly, so nothing has to know which fit
    ## is live.
    methylome_plot_divergence_server("ewas_divplot", r_ewas, r_thresh,
                                     watermark = watermark)
    methylome_table_divergence_server("ewas_divtbl", r_ewas, r_thresh)
  })
}
