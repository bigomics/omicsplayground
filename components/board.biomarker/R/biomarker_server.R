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

    pdx_infotext <- tspan("The <strong>Biomarker Board</strong> performs
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

    ## ================================================================================
    ## ======================= REACTIVE/OBSERVE FUNCTIONS =============================
    ## ================================================================================

    # Observe tabPanel change to update Settings visibility
    tab_elements <- list(
      "Feature selection" = list(
        enable = NULL,
        disable = NULL
      ),
      "Feature-set ranking" = list(
        enable = NULL,
        disable = c("pdx_predicted", "pdx_filter")
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

    shiny::observe({
      shiny::req(pgx$X)
      ct <- colnames(pgx$Y)
      shiny::updateSelectInput(session, "pdx_predicted", choices = ct)
    })

    shiny::observe({
      shiny::req(pgx$Y)
      ## levels for sample filter
      levels <- playbase::getLevels(pgx$Y)
      shiny::updateSelectInput(session, "pdx_samplefilter", choices = levels)
    })

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
    # if the pdx_predicted overlaps with the pdx_samplefilter variable
    shiny::observeEvent(
      list(
        pgx$Y,
        input$pdx_samplefilter,
        input$pdx_predicted
      ),
      {
        shiny::req(pgx$Y, input$pdx_predicted)
        # check how many levels pgx_predicted has if it has more than
        # 1 level, then enable the run button if it has 1 level, then
        # disable the run button
        kk <- selected_samples()
        levels_filtered <- unique(pgx$Y[kk, input$pdx_predicted])
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
          input$pdx_predicted,
          input$pdx_samplefilter,
          input$pdx_filter,
          pgx$X
        )
      },
      {
        is_computed(FALSE)
      }
    )

    calcVariableImportance <- shiny::eventReactive(input$pdx_runbutton, {
      ## This code also features a progress indicator.
      shiny::req(pgx$X, input$pdx_predicted)

      ct <- 2
      ct <- 12
      colnames(pgx$Y)
      shiny::isolate(ct <- input$pdx_predicted)
      do.survival <- grepl("survival", ct, ignore.case = TRUE)
      if (is.null(ct)) {
        return(NULL)
      }

      NFEATURES <- 50
      NFEATURES <- 60

      ## Create a Progress object
      progress <- shiny::Progress$new()
      ## Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())

      progress$set(message = "Variable importance", value = 0)

      if (!(ct %in% colnames(pgx$Y))) {
        return(NULL)
      }
      y0 <- as.character(pgx$Y[, ct])
      names(y0) <- rownames(pgx$Y)
      y0 <- y0[names(y0) %in% selected_samples()]
      y <- y0[!is.na(y0)]

      ## augment to at least 100 samples per level :)
      ii <- tapply(1:length(y), y, function(ii) sample(c(ii, ii), size = 100, replace = TRUE))
      y <- y[unlist(ii)]

      ## -------------------------------------------
      ## select features
      ## -------------------------------------------
      if (FALSE && shiny::isolate(input$pdx_level) == "geneset") {
        X <- pgx$gsetX[, names(y)]
        if (any(is.na(X))) {
          X <- X[complete.cases(X), , drop = FALSE]
        }
      } else {
        X <- pgx$X[, names(y)]
        if (any(is.na(X))) {
          X <- pgx$impX[, names(y)]
        }
      }
      X0 <- X

      ## ----------- filter with selected features
      progress$inc(1 / 10, detail = "Filtering features")

      ft <- "<all>"
      ft <- shiny::isolate(input$pdx_filter)
      if (is.null(ft)) {
        return(NULL)
      }
      shiny::isolate(sel <- input_pdx_select())

      is.family <- (ft %in% c(names(pgx$families), names(playdata::iGSETS)))

      if (ft == "<custom>" && !is.null(sel) && length(sel) > 0) {
        ## ------------- filter with user selection
        if (sel[1] != "") {
          pp <- rownames(X)[which(toupper(rownames(X)) %in% toupper(sel))]
          X <- X[pp, , drop = FALSE]
        }
      } else if (is.family) {
        pp <- rownames(X)
        if (ft %in% names(pgx$families)) {
          gg <- pgx$families[[ft]]
          pp <- playbase::filterProbes(pgx$genes, gg)
        } else if (ft %in% names(playdata::iGSETS)) {
          gg <- unlist(playdata::getGSETS(ft))
          pp <- playbase::filterProbes(pgx$genes, gg)
        }
        pp <- intersect(pp, rownames(X))
        X <- X[pp, , drop = FALSE]
      }

      ## ----------- restrict to top SD -----------
      X <- head(X[order(-apply(X, 1, sd, na.rm = TRUE)), , drop = FALSE], 10 * NFEATURES) ## top 100
      sdx0 <- matrixStats::rowSds(X, na.rm = TRUE)
      sdx1 <- 0.5 * sdx0 + 0.5 * mean(sdx0, na.rm = TRUE)
      X <- X + 0.25 * sdx1 * matrix(rnorm(length(X)), nrow(X), ncol(X)) ## add some noise

      progress$inc(4 / 10, detail = "computing scores")

      ## -------------------------------------------
      ## compute importance values
      ## -------------------------------------------
      if (do.survival) {
        time <- abs(y)
        status <- (y > 0) ## dead is positive time
        methods <- c("glmnet", "randomforest", "xgboost", "pls")
        P <- playbase::pgx.survivalVariableImportance(
          X,
          time = time, status = status, methods = methods
        )
      } else {
        methods <- c("glmnet", "randomforest", "xgboost", "splsda", "correlation", "ftest")
        X1 <- X
        y1 <- y
        names(y1) <- colnames(X1) <- paste0("x", 1:ncol(X))
        res <- playbase::pgx.variableImportance(
          X1, y1,
          methods = methods, reduce = 1000, resample = 0,
          scale = FALSE, add.noise = 0
        )
        P <- res$importance
      }
      P <- abs(P) ## sometimes negative according to sign

      P[is.na(P)] <- 0
      P[is.nan(P)] <- 0
      P <- t(t(P) / (1e-3 + apply(P, 2, max, na.rm = TRUE)))
      P <- P[order(-rowSums(P, na.rm = TRUE)), , drop = FALSE]

      R <- P
      if (nrow(R) > 1) {
        R <- (apply(P, 2, rank) / nrow(P))**4
        R <- R[order(-rowSums(R, na.rm = TRUE)), , drop = FALSE]
      }

      progress$inc(3 / 10, detail = "drawing tree")

      ## ------------------------------
      ## create partition tree
      ## ------------------------------

      R <- R[order(-rowSums(R, na.rm = TRUE)), , drop = FALSE]
      sel <- head(rownames(R), 100)
      sel <- intersect(sel, rownames(X))
      sel <- head(rownames(R), NFEATURES) ## top50 features
      tx <- t(X[sel, , drop = FALSE])


      ## formula wants clean names, so save original names
      colnames(tx) <- gsub("[: +-.,]", "_", colnames(tx))
      colnames(tx) <- gsub("[')(]", "", colnames(tx))
      colnames(tx) <- gsub("\\[|\\]", "", colnames(tx))
      orig.names <- sel
      names(orig.names) <- colnames(tx)
      jj <- names(y)

      if (do.survival) {
        time <- abs(y)
        status <- (y > 0) ## dead if positive time
        df <- data.frame(time = time + 0.001, status = status, tx)
        rf <- rpart::rpart(survival::Surv(time, status) ~ ., data = df)
      } else {
        df <- data.frame(y = y, tx)
        rf <- rpart::rpart(y ~ ., data = df)
      }
      table(rf$where)
      rf$cptable
      rf$orig.names <- orig.names

      rf.nsplit <- rf$cptable[, "nsplit"]
      if (grepl("survival", ct)) {
        MAXSPLIT <- 4 ## maximum five groups....
      } else {
        MAXSPLIT <- 1.5 * length(unique(y)) ## maximum N+1 groups
      }
      if (max(rf.nsplit) > MAXSPLIT) {
        cp.idx <- max(which(rf.nsplit <= MAXSPLIT))
        cp0 <- rf$cptable[cp.idx, "CP"]
        rf <- rpart::prune(rf, cp = cp0)
      }

      progress$inc(2 / 10, detail = "done")

      y <- y[rownames(tx)]
      colnames(tx) <- orig.names[colnames(tx)]
      res <- list(R = R, y = y, X = t(tx), rf = rf)

      is_computed(TRUE)
      return(res)
    })

    ## ================================================================================
    ## ==================================== PLOTS =====================================
    ## ================================================================================

    biomarker_plot_importance_server(
      "pdx_importance",
      calcVariableImportance,
      is_computed,
      watermark = WATERMARK
    )

    biomarker_plot_heatmap_server(
      "pdx_heatmap",
      calcVariableImportance,
      pgx,
      reactive(input$pdx_predicted),
      is_computed,
      watermark = WATERMARK
    )

    biomarker_plot_decisiontree_server(
      "pdx_decisiontree",
      calcVariableImportance,
      is_computed,
      watermark = WATERMARK
    )

    biomarker_plot_auc_server(
      "pdx_auc",
      calcVariableImportance,
      is_computed,
      watermark = WATERMARK
    )

    biomarker_plot_boxplots_server(
      "pdx_boxplots",
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
