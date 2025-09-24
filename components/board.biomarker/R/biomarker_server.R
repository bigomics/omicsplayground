##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

BiomarkerBoard <- function(id, pgx) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    fullH <- 800 ## full height of panel
    rowH <- 320 ## row height of panel
    imgH <- 260

    pdx_infotext <- tspan("<center><iframe width='560' height='315' src='https://www.youtube.com/embed/IICgZVUSrpU?si=0m9gvGU7jArsZNnW&amp;start=193' title='YouTube video player' frameborder='0' allow='accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share' referrerpolicy='strict-origin-when-cross-origin' allowfullscreen></iframe></center><br><br>The <strong>Biomarker Board</strong> performs
    the biomarker selection that can be used for classification or prediction purposes.
    <br><br>To better understand which genes, mutations, or gene sets influence
    the final phenotype the most, Playground calculates a variable importance
    score for each feature using state-of-the-art machine learning algorithms,
    including LASSO, elastic nets, random forests, and extreme gradient boosting,
    and provides the top 50 features according to cumulative ranking by the algorithms.
    By combining several methods, the platform aims to select the best possible biomarkers.
    <br><br>The phenotype of interest can be multi-categorical classes or patient
    survival data. Instead of choosing a phenotype, users can also specify a particular
    contrast from the analysis and perform biomarker selection. The platform also
    provides a heatmap of samples based on identified top features.
    <br><br>In addition, it creates a classification tree using top features and
    provides expression boxplots by phenotype classes for features present in the
    tree. The platform can also provide a survival tree analysis using top features
    and provides expression boxplots by phenotype classes for features present in
    the tree.", js = FALSE)

    ## =========================================================================
    ## ======================== OBSERVERS ======================================
    ## =========================================================================

    # Observe tabPanel change to update Settings visibility
    tab_elements <- list(
      "Feature selection" = list(
        enable = NULL,
        disable = NULL
      ),
      "Feature-set ranking" = list(
        enable = NULL,
        disable = c("pdx_target", "pdx_filter")
      )
    )
    shiny::observeEvent(input$tabs1, {
      bigdash::update_tab_elements(input$tabs1, tab_elements)
    })

    shiny::observeEvent(input$pdx_info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Biomarker Board</strong>"),
        shiny::HTML(pdx_infotext),
        easyClose = TRUE, size = "l"
      ))
    })


    shiny::observe({
      shiny::req(pgx$X)
      ct <- colnames(pgx$Y)
      shiny::updateSelectInput(session, "pdx_target", choices = ct)
    })

    shiny::observe({
      shiny::req(pgx$Y)
      ## levels for sample filter
      levels <- playbase::getLevels(pgx$Y)
      shiny::updateSelectInput(session, "pdx_samplefilter", choices = levels)
    })

    shiny::observe({
      shiny::req(pgx$X)
      if (FALSE && shiny::isolate(input$pdx_level == "geneset")) {
        gset_collections <- playbase::pgx.getGeneSetCollections(gsets = rownames(pgx$gsetX))
        ft <- names(gset_collections)
        nn <- sapply(gset_collections, function(x) sum(x %in% rownames(pgx$gsetX)))
        ft <- ft[nn >= 10]
      } else {
        ## gene level
        ft <- names(pgx$families)
      }
      ft <- sort(ft)
      ft <- sort(c("<custom>", ft))
      shiny::updateSelectInput(session, "pdx_filter", choices = ft, selected = "<all>")
      shiny::updateTextAreaInput(
        session = session,
        inputId = "pdx_select",
        placeholder = tspan("Paste your gene list", js = FALSE)
      )
    })

    # Enable or disable the run button in the UI
    # if the pdx_target overlaps with the pdx_samplefilter variable
    shiny::observeEvent(
      list(
        pgx$Y,
        input$pdx_samplefilter,
        input$pdx_target
      ),
      {
        shiny::req(pgx$Y, input$pdx_target)
        # check how many levels pgx_predicted has if it has more than
        # 1 level, then enable the run button if it has 1 level, then
        # disable the run button
        kk <- selected_samples()
        levels_filtered <- unique(pgx$Y[kk, input$pdx_target])
        if (length(levels_filtered) > 1) {
          shinyjs::enable("pdx_runbutton")
        } else {
          shinyjs::disable("pdx_runbutton")
        }
      }
    )

    is_computed <- reactiveVal(FALSE)
    observeEvent(
      {
        list(
          input$pdx_target,
          input$pdx_samplefilter,
          input$pdx_filter,
          pgx$X
        )
      },
      {
        is_computed(FALSE)
      }
    )

    ## =========================================================================
    ## ============================= REACTIVES =================================
    ## =========================================================================

    input_pdx_select <- shiny::reactive({
      gg <- input$pdx_select
      if (is.null(gg)) {
        return(NULL)
      }

      gg <- strsplit(as.character(gg), split = "[, \n\t]")[[1]]
      if (length(gg) == 0) {
        return(NULL)
      }
      if (length(gg) == 1 && gg[1] != "") gg <- c(gg, gg) ## hack to allow single gene....
      return(gg)
    }) %>% shiny::debounce(1000)

    ## get selected samples after sample filtering
    selected_samples <- shiny::eventReactive(
      list(pgx$Y, input$pdx_samplefilter),
      {
        shiny::req(pgx$Y)
        samples <- rownames(pgx$Y)
        sel <- input$pdx_samplefilter
        if (!is.null(sel) && all(sel != "")) {
          samples <- playbase::selectSamplesFromSelectedLevels(pgx$Y, sel)
        }
        samples
      }
    )

    ## calculate variable importance upon compute button
    calcVariableImportance <- shiny::eventReactive(input$pdx_runbutton, {
      ## This code also features a progress indicator.
      shiny::req(pgx$X, input$pdx_target)

      shiny::isolate(ph <- input$pdx_target)
      do.survival <- grepl("survival", ph, ignore.case = TRUE)
      if (is.null(ph)) {
        return(NULL)
      }

      ## Create a Progress object
      progress <- shiny::Progress$new()
      ## Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())
      progress$set(message = "", value = 0)

      if (!(ph %in% colnames(pgx$Y))) {
        return(NULL)
      }

      ft <- shiny::isolate(input$pdx_filter)
      if (is.null(ft)) {
        return(NULL)
      }
      shiny::isolate(sel <- input_pdx_select())

      progress$inc(0.33, detail = "Calculating variable importance. Please wait...")
      if (pgx$datatype == "multi-omics") {
        res <- playbase::pgx.compute_importance(
          pgx,
          pheno = ph,
          level = "genes",
          multiomics = 2, ## 2-pass
          filter_features = ft,
          select_features = sel,
          select_samples = selected_samples(),
          nfeatures = 50
        )
      } else {
        res <- playbase::pgx.compute_importance(
          pgx,
          pheno = ph,
          level = "genes",
          filter_features = ft,
          select_features = sel,
          select_samples = selected_samples(),
          nfeatures = 50
        )
      }

      is_computed(TRUE)
      return(res)
    })

    ## ===========================================================================
    ## =============================== PLOTS =====================================
    ## ===========================================================================

    biomarker_plot_importance_server(
      "pdx_importance",
      pgx,
      calcVariableImportance,
      is_computed,
      watermark = WATERMARK
    )

    biomarker_plot_heatmap_server(
      "pdx_heatmap",
      calcVariableImportance,
      pgx,
      reactive(input$pdx_target),
      is_computed,
      watermark = WATERMARK
    )

    biomarker_plot_decisiontree_server(
      "pdx_decisiontree",
      calcVariableImportance,
      pgx,
      is_computed,
      watermark = WATERMARK
    )

    biomarker_plot_boxplots_server(
      "pdx_boxplots",
      pgx,
      calcVariableImportance,
      is_computed,
      watermark = WATERMARK
    )

    biomarker_plot_featurerank_server(
      id = "featurerank",
      pgx = pgx,
      ft_level = shiny::reactive("gene"),
      samplefilter = shiny::reactive(input$pdx_samplefilter),
      watermark = WATERMARK
    )
  })
} ## end-of-Board
