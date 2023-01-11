##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

BiomarkerBoard <- function(id, inputData) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    fullH <- 800 ## full height of panel
    rowH <- 320 ## row height of panel
    imgH <- 260


    pdx_infotext <- strwrap("The <strong>Biomarker Board</strong> performs
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
    <br><br>In addition, it generates a classification tree using top features and
    provides expression boxplots by phenotype classes for features present in the
    tree. The platform can also provide a survival tree analysis using top features
    and provides expression boxplots by phenotype classes for features present in
    the tree.")

    ## ================================================================================
    ## ======================= REACTIVE/OBSERVE FUNCTIONS =============================
    ## ================================================================================

    shiny::observeEvent(input$pdx_info, {
      dbg("<module-biomarker::observe pdxinfo> reacted")
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Biomarker Board</strong>"),
        shiny::HTML(pdx_infotext),
        easyClose = TRUE, size = "l"
      ))
    })

    input_pdx_select <- shiny::reactive({
      dbg("[BiomarkerBoard:<input_pdx_select>]  reacted")
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
      ngs <- inputData()
      ## if(is.null(ngs)) return(NULL)
      shiny::req(ngs)
      dbg("[BiomarkerBoard::observe1] reacted")
      dbg("[BiomarkerBoard::observe1] dim(ngs$Y) = ", dim(ngs$Y))
      ct <- colnames(ngs$Y)
      ## ct <- grep("group|sample|patient|donor",ct,value=TRUE,invert=TRUE)
      ## ct <- grep("sample|patient|donor",ct,value=TRUE,invert=TRUE)
      shiny::updateSelectInput(session, "pdx_predicted", choices = ct)
    })

    shiny::observe({
      ngs <- inputData()
      shiny::req(ngs)
      ## input$pdx_runbutton
      dbg("[BiomarkerBoard::observe2] reacted")

      if (FALSE && shiny::isolate(input$pdx_level == "geneset")) {
        ft <- names(COLLECTIONS)
        nn <- sapply(COLLECTIONS, function(x) sum(x %in% rownames(ngs$gsetX)))
        ft <- ft[nn >= 10]
      } else {
        ## gene level
        ## ft <- pgx.getFamilies(ngs,nmin=10,extended=FALSE)
        ft <- names(ngs$families)
      }
      ft <- sort(ft)
      ## if(input$pdx_level == "gene") ft = sort(c("<custom>",ft))
      ft <- sort(c("<custom>", ft))
      shiny::updateSelectInput(session, "pdx_filter", choices = ft, selected = "<all>")
    })

    calcVariableImportance <- shiny::eventReactive(input$pdx_runbutton, {
      ##
      ## This code also features a progress indicator.
      ##

      ## input$pdx_runbutton
      dbg("[calcVariableImportance] reacted on runbutton")

      ngs <- inputData()
      if (is.null(ngs)) {
        return(NULL)
      }
      shiny::req(ngs, input$pdx_predicted)

      dbg("[calcVariableImportance] 0: ")

      ct <- 2
      ct <- 12
      colnames(ngs$Y)
      shiny::isolate(ct <- input$pdx_predicted)
      do.survival <- grepl("survival", ct, ignore.case = TRUE)
      if (is.null(ct)) {
        return(NULL)
      }

      dbg("[calcVariableImportance] 1: called!")

      NFEATURES <- 50
      NFEATURES <- 60

      ## Create a Progress object
      progress <- shiny::Progress$new()
      ## Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())

      progress$set(message = "Variable importance", value = 0)

      if (!(ct %in% colnames(ngs$Y))) {
        return(NULL)
      }
      y0 <- as.character(ngs$Y[, ct])
      names(y0) <- rownames(ngs$Y)
      y <- y0[!is.na(y0)]


      ## augment to 100 samples
      table(y)
      ## if(length(y)<40) y <- head(rep(y,10),100)
      if (length(y) < 100) y <- head(rep(y, 100), 100)
      table(y)

      ## -------------------------------------------
      ## select features
      ## -------------------------------------------
      ## group prediction
      if (FALSE && shiny::isolate(input$pdx_level) == "geneset") {
        X <- ngs$gsetX[, names(y)]
      } else {
        X <- ngs$X[, names(y)]
      }
      dim(X)
      X0 <- X
      length(y)

      ## ----------- filter with selected features
      progress$inc(1 / 10, detail = "Filtering features")

      ft <- "<all>"
      ft <- shiny::isolate(input$pdx_filter)
      if (is.null(ft)) {
        return(NULL)
      }
      shiny::isolate(sel <- input_pdx_select())

      is.family <- (ft %in% c(names(ngs$families), names(iGSETS)))

      if (ft == "<custom>" && !is.null(sel) && length(sel) > 0) {
        ## ------------- filter with user selection
        if (sel[1] != "") {
          dbg("[calcVariableImportance] 2: using custom list of variable ")
          ## pp <- intersect(rownames(X),sel)
          pp <- rownames(X)[which(toupper(rownames(X)) %in% toupper(sel))]
          X <- X[pp, , drop = FALSE]
        }
      } else if (is.family) {
        pp <- rownames(X)
        if (ft %in% names(ngs$families)) {
          dbg("[calcVariableImportance] 2: using ngs$families")
          gg <- ngs$families[[ft]]
          pp <- filterProbes(ngs$genes, gg)
        } else if (ft %in% names(iGSETS)) {
          dbg("[calcVariableImportance] 2: using genesets")
          gg <- unlist(getGSETS(ft))
          pp <- filterProbes(ngs$genes, gg)
        }
        pp <- intersect(pp, rownames(X))
        X <- X[pp, , drop = FALSE]
      } else {
        dbg("[calcVariableImportance] 2: using all features")
      }

      dbg("[calcVariableImportance] 3: dim.X = ", dim(X))

      ## ----------- restrict to top 100
      dim(X)
      X <- head(X[order(-apply(X, 1, sd)), , drop = FALSE], 10 * NFEATURES) ## top 100
      sdx <- mean(apply(X, 1, sd))
      X <- X + 0.25 * sdx * matrix(rnorm(length(X)), nrow(X), ncol(X)) ## add some noise
      dim(X)

      dbg("[calcVariableImportance] 4: dim.X = ", dim(X))

      progress$inc(4 / 10, detail = "computing scores")

      ## -------------------------------------------
      ## compute importance values
      ## -------------------------------------------
      if (do.survival) {
        time <- abs(y)
        status <- (y > 0) ## dead is positive time
        methods <- c("glmnet", "randomforest", "boruta", "xgboost", "pls")
        methods <- c("glmnet", "randomforest", "xgboost", "pls")
        P <- pgx.survivalVariableImportance(
          X,
          time = time, status = status, methods = methods
        )
      } else {
        methods <- c("glmnet", "randomforest", "boruta", "xgboost", "pls")
        methods <- c("glmnet", "randomforest", "xgboost", "pls")
        X1 <- X
        y1 <- y
        names(y1) <- colnames(X1) <- paste0("x", 1:ncol(X))
        P <- pgx.multiclassVariableImportance(X1, y1, methods = methods)
        ## P <- pgx.variableImportance(X1, y1, methods=methods)
      }
      P <- abs(P)
      head(P)

      P[is.na(P)] <- 0
      P[is.nan(P)] <- 0
      P <- t(t(P) / (1e-3 + apply(P, 2, max, na.rm = TRUE)))
      ## P <- pmax(P,0.1)
      P <- P[order(-rowSums(P, na.rm = TRUE)), , drop = FALSE]
      head(P)

      R <- P
      if (nrow(R) > 1) {
        R <- (apply(P, 2, rank) / nrow(P))**4
        R <- R[order(-rowSums(R)), , drop = FALSE]
        head(R)
      }

      if (FALSE && DEV) {
        is.multiomics <- any(grepl("\\[gx\\]|\\[mrna\\]", rownames(R)))
        is.multiomics
        do.multiomics <- (is.multiomics && shiny::isolate(input$pdx_multiomics))
        if (do.multiomics) {
          dbg("calcVariableImportance:: 5: EXPERIMENTAL: multi-omics weighting")
          ## EXPERIMENTAL: multi-omics weighting
          rr <- rowMeans(R)
          dtype <- ngs$genes[rownames(R), "data_type"]
          gene <- ngs$genes[rownames(R), "gene_name"]
          gg <- sort(unique(gene))
          dtypes <- unique(dtype)
          dt <- "gx"
          rho <- c()
          for (dt in dtypes) {
            r1 <- rr[which(dtype == dt)]
            g1 <- gene[which(dtype == dt)]
            r1 <- r1[match(gg, g1)]
            rho <- cbind(rho, r1)
          }
          rownames(rho) <- gg
          colnames(rho) <- dtypes
          rho[is.na(rho)] <- 0
          rho <- rho[order(-rowSums(rho, na.rm = TRUE)), ]
          dim(rho)
          barplot(t(head(rho, 50)), las = 3)
          jj <- match(gene, rownames(rho))
          rname <- rownames(R)
          R <- rho[jj, ]
          rownames(R) <- rname
        }
      }

      dbg("calcVariableImportance:: 5: drawing tree\n")
      progress$inc(3 / 10, detail = "drawing tree")

      ## ------------------------------
      ## create partition tree
      ## ------------------------------

      R <- R[order(-rowSums(R)), , drop = FALSE]
      sel <- head(rownames(R), 100)
      sel <- intersect(sel, rownames(X))
      sel <- head(rownames(R), NFEATURES) ## top50 features
      tx <- t(X[sel, , drop = FALSE])
      dim(tx)

      ## formula wants clean names, so save original names
      colnames(tx) <- gsub("[: +-.,]", "_", colnames(tx))
      colnames(tx) <- gsub("[')(]", "", colnames(tx))
      colnames(tx) <- gsub("\\[|\\]", "", colnames(tx))
      orig.names <- sel
      names(orig.names) <- colnames(tx)
      jj <- names(y)
      ## ny <- length(unique(y))
      ## if(length(jj) < ny*20) jj <- c(jj,jj,jj)

      if (do.survival) {
        time <- abs(y)
        status <- (y > 0) ## dead if positive time
        df <- data.frame(time = time + 0.001, status = status, tx)
        ## df <- df[jj,]
        rf <- rpart::rpart(survival::Surv(time, status) ~ ., data = df)
      } else {
        df <- data.frame(y = y, tx)
        ## df <- df[jj,]
        rf <- rpart::rpart(y ~ ., data = df)
      }
      table(rf$where)
      rf$cptable
      rf$orig.names <- orig.names

      rf.nsplit <- rf$cptable[, "nsplit"]
      max(rf.nsplit)
      length(unique(y))
      if (grepl("survival", ct)) {
        MAXSPLIT <- 4 ## maximum five groups....
      } else {
        MAXSPLIT <- 1.5 * length(unique(y)) ## maximum N+1 groups
      }
      if (max(rf.nsplit) > MAXSPLIT) {
        cp.idx <- max(which(rf.nsplit <= MAXSPLIT))
        cp0 <- rf$cptable[cp.idx, "CP"]
        ## rf <- rpart::prune(rf, cp=0.05)
        rf <- rpart::prune(rf, cp = cp0)
      }
      table(rf$where)

      dbg("[calcVariableImportance] done!!!\n")
      progress$inc(2 / 10, detail = "done")

      ## y <- y[rownames(ngs$samples)]
      ## tx <- tx[names(y),]
      y <- y[rownames(tx)]
      colnames(tx) <- orig.names[colnames(tx)]
      ## res <- list(P=P, R=R, y=y, X=t(tx), rf=rf)
      res <- list(R = R, y = y, X = t(tx), rf = rf)

      return(res)
    })

    ## ================================================================================
    ## ==================================== PLOTS =====================================
    ## ================================================================================

    biomarker_plot_importance_server(
      'pdx_importance',
      calcVariableImportance,
      watermark = WATERMARK
    )

    biomarker_plot_heatmap_server(
      'pdx_heatmap',
      calcVariableImportance,
      inputData,
      reactive(input$pdx_predicted),
      watermark = WATERMARK
    )

    biomarker_plot_decisiontree_server(
      'pdx_decisiontree',
      calcVariableImportance,
      watermark = WATERMARK
    )

    biomarker_plot_boxplots_server(
      'pdx_boxplots',
      calcVariableImportance,
      watermark = WATERMARK
    )

  })
} ## end-of-Board
