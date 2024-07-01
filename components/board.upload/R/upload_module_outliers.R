##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


## =============================================================================
## ==================== OUTLIERS UI/SERVER =================================
## =============================================================================


upload_module_outliers_ui <- function(id, height = "100%") {
  ns <- shiny::NS(id)

  score.infotext <-
    "The outlier z-score is calculated as the average of z-score from correlation, euclidean distance and avarage feature z-score."

  missing.infotext <-
    "Analysis of variables by plotting their significance in correlation with the phenotype against their significance in correlation with a principal component (PC) vector. Strong model variables are situate 'top right'. Batch effect variables with high PC correlation but low phenotype correlation are on the 'top left'. A well-designed experiment shows strong model variables in PC1, else it may be a sign of significant batch-effects."

  missing.options <- tagList(
    shiny::radioButtons(ns("missing_plottype"), "Plot type:", c("heatmap", "ratio plot"),
      selected = "heatmap", inline = TRUE
    ),
  )

  norm.options <- tagList(
    shiny::radioButtons(ns("norm_plottype"), "Plot type:", c("boxplot", "histogram", "density"),
      selected = "boxplot", inline = TRUE
    ),
    shiny::checkboxInput(ns("norm_zero_na"), "zero as NA", FALSE)
  )

  outlier.options <- tagList(

    shiny::checkboxInput(ns("outlier_shownames"), "show sample names", FALSE)
  )

  bec.options <- tagList(
    shiny::radioButtons(ns("bec_plottype"), "Plot type:", c("pca", "tsne", "heatmap"),
      inline = TRUE
    )
  )

  bslib::layout_columns(
    col_widths = c(2, 10),
    height = "calc(100vh - 200px)",
    heights_equal = "row",
    bslib::card(bslib::card_body(
      style = "padding: 0px;",
      bslib::accordion(
        multiple = FALSE,
        style = "background-color: #F7FAFD99;",
        bslib::accordion_panel(
          title = "1. Missing values",
          shiny::p("Replace missing values using an imputation method:\n"),
          shiny::selectInput(ns("impute_method"), NULL,
          ##  choices = c("bpca","LLS","MinDet","MinProb","NMF","RF","SVD2","zero"),
            choices = c(
                "SVDimpute (default)" = "SVD2",
                "Zero" = "zero",
                "MinDet",
                "MinProb",
                "NMF",
                "Skip imputation" = "skip_imputation"
            ),
            selected = "SVD2"
            ),
          shiny::checkboxInput(ns("zero_as_na"), label = "Treat zero as NA", value = FALSE),
          br()
        ),
        bslib::accordion_panel(
          title = "2. Normalization",
          div("Normalize data values:\n"),
          shiny::selectInput(ns("scaling_method"), NULL,
            choices = c(
                "LogCPM (default)" = "CPM",
                ## "Mean centering" = "mean.center",
                ## "Median centering" = "median.center",
                ## "TMM" = "TMM",
                ## "RLE" = "RLE",
                ## "RLE2" = "RLE2",
                "LogMaxMedian" = "logMaxMedian",
                "LogMaxSum" = "logMaxSum",
                "Skip normalization" = "Skip_normalization"
            ),
            selected = "CPM"
          ),
          shiny::checkboxInput(ns("quantile_norm"), "Quantile normalization", value = TRUE),
          br()
        ),
        bslib::accordion_panel(
          title = "3. Detect and remove outliers",
          shiny::p("Identify and remove outliers (i.e., bad samples) from your dataset.\n"),
          shiny::checkboxInput(ns("remove_outliers"), "Check and remove outliers", value = TRUE),
          ## shiny::selectInput(ns("remove_outliers_samples"), NULL,
          shiny::sliderInput(ns("outlier_threshold"), "Threshold:", 1, 12, 6, 1),
          br()
        ),
        ## bslib::accordion_panel(
        ##  title = "4. Remove unwanted technical and biogical variation",
        ##  shiny::p("Correct for technical and biogical factors:\n"),
        ##  shiny::selectInput(ns("correct_factor"), NULL,
        ##    choices = c(
        ##        "library size" = "lib",
        ##        "ribo", "cellcycle", "gender",
        ##        "Skip correction (Default)" = "Skip_correction"
        ##    ),
        ##    selected = "Skip_correction"
        ##  ),
        ##  br()
        ## ),
        bslib::accordion_panel(
          title = "4. Batch-effect correction",
          shiny::p("Remove batch effects from your data:\n"),
          shiny::selectInput(ns("bec_method"), NULL,
              choices = c("uncorrected (default)" = "uncorrected",
                          "ComBat" = "ComBat",
                          "SVA" = "SVA",
                          "RUV3" = "RUV3",
                          "NPM" = "NPmatch"),
            selected = "uncorrected"
          ),
          shiny::conditionalPanel(
            "input.bec_method == 'ComBat' || input.bec_method == 'limma'",
            ns = ns,
            shiny::textOutput(ns("bec_param_text")),
            shiny::br(),
          ),
          shiny::checkboxInput(ns("bec_preview_all"), "Preview all methods", value = TRUE),
          br()
        )
      ),
      br()
    )),

    ## ---------------------------- canvas ----------------------------------
    bslib::layout_columns(
      width = 12,
      bslib::layout_columns(
        col_widths = 6,
        row_heights = c(3, 3),
        height = "calc(100vh - 200px)",
        heights_equal = "row",
        PlotModuleUI(
          ns("plot2"),
          title = "Missing values",
          info.text = missing.infotext,
          caption = missing.infotext,
          options = missing.options,
          height = c("100%", "70vh"),
          show.maximize = FALSE
        ),
        PlotModuleUI(
          ns("plot1"),
          title = "Normalization",
          options = norm.options,
          height = c("100%", "70vh"),
          show.maximize = FALSE
        ),
        PlotModuleUI(
          ns("plot3"),
          title = "Outliers detection",
          info.text = score.infotext,
          caption = score.infotext,
          options = outlier.options,
          height = c("100%", "70vh"),
          show.maximize = FALSE
        ),
        PlotModuleUI(
          ns("plot4"),
          title = "Batch-effects correction",
          options = NULL,
          height = c("100%", "70vh"),
          show.maximize = FALSE
        )
      )
    )
  )
}


upload_module_outliers_server <- function(id, r_X, r_samples, r_contrasts,
                                          is.count = FALSE, height = 720) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      output$bec_param_text <- renderText({
        pars <- get_model_parameters()
        shiny::req(pars)
        batch.pars <- pars$batch.pars
        paste("Batch parameters:", paste(batch.pars, collapse = "+"), "\n")
      })

      get_model_parameters <- eventReactive(
        {
          list(r_X(), r_samples(), r_contrasts())
        },
        {
          shiny::req(r_X(), r_samples(), r_contrasts())

          X <- r_X()
          samples <- r_samples()
          contrasts <- r_contrasts()

          shiny::req(ncol(X) == nrow(samples))
          shiny::req(nrow(contrasts) == ncol(X))

          pars <- playbase::get_model_parameters(X, samples, pheno = NULL, contrasts = contrasts)

          list(
            batch.pars = pars$batch.pars,
            pheno.pars = pars$pheno.pars,
            batch = pars$batch,
            pheno = pars$pheno
          )
        }
      )

      ## ------------------------------------------------------------------
      ## Object reactive chain
      ## ------------------------------------------------------------------
      
      ## Impute; Remove duplicate features
      imputedX <- reactive({
          shiny::req(r_X())
          counts <- r_X()                   

          counts[which(is.nan(counts))] <- NA
          counts[which(is.infinite(counts))] <- NA

          zeros <- sum(counts==0, na.rm=TRUE)
          negs <- sum(counts<0, na.rm=TRUE) ## what about for Olink NPX?
          nmissing <- sum(is.na(counts))
          infin <- sum(is.infinite(counts))
          dbg("[outliers_server] Counts data have ", zeros, " zero values.")
          dbg("[outliers_server] Counts data have ", negs, " negative values.")
          dbg("[outliers_server] Counts data have ", nmissing, " missing values.")
          dbg("[outliers_server] Counts data have ", infin, " infinite values.")

          if(input$zero_as_na) {
              dbg("[outliers_server] Setting 0 values to NA")
              counts[which(counts==0)] <- NA
          } else {
              dbg("[outliers_server] Adding 1e-20 to counts data")
              counts <- counts + 1e-20
          }

          X <- log2(counts)
          X[playbase::is.xxl(X, z = 10)] <- NA
          nmissing <- sum(is.na(X))

          if(nmissing>0) {
              dbg("[outliers_server] X has ", nmissing, " missing values.")
              if (input$impute_method != "skip_imputation") {
                  dbg("[outliers_server] Imputing data using ", input$impute_method)
                  X <- playbase::imputeMissing(X, method = input$impute_method)
                  dbg("[outliers_server] dim.imputedX = ", dim(X))
              } else {
                  dbg("[outliers_server] Skipping imputation")
                  dbg("[outliers_server] Assigning 1e-20 to ", nmissing, "missing values")
                  counts <- 2**X
                  counts[which(is.na(counts))] <- 1e-20
                  X <- log2(counts)
              }
          } else {
              dbg("[outliers_server] No missing values detected in the data. Not imputing.")
          }
          dbg("[outliers_server] Checking for duplicated features")
          X <- playbase::counts.mergeDuplicateFeatures(X, is.counts = FALSE)
          dbg("[outliers_server] dim.X = ", dim(X))
                    
          counts <- 2 ** X
          LL <- list(counts = counts, X = X)
          LL
      })

      ## Normalize
      normalizedX <- reactive({
          shiny::req(dim(imputedX()$X), dim(imputedX()$counts))
          counts <- imputedX()$counts
          X <- imputedX()$X ## can be imputed or not (see above). log2. Can have negatives.
          dbg("[outliers_server] Normalization step: dim.X = ", dim(X))
          dbg("[outliers_server] Normalization step: dim.counts = ", dim(counts))

          ## zeros <- sum(counts==0, na.rm=TRUE)
          ## negs <- sum(counts<0, na.rm=TRUE)
          ## nmissing <- sum(is.na(counts))
          ## infin <- sum(is.infinite(counts))
          ## dbg("[outliers_server] counts data have ", zeros, " zero values.")
          ## dbg("[outliers_server] countsdata have ", negs, " negative values.")
          ## dbg("[outliers_server] counts data have ", nmissing, " missing values.")
          ## dbg("[outliers_server] counts data have ", infin, " infinite values.")

          if (input$scaling_method == "Skip_normalization") {
              dbg("[outliers_server] Skipped normalization: dim.X = ", dim(X))
              return(X)
          } else {
              shiny::withProgress(message = "Normalizing the data...", value = 0, {
                  shiny::incProgress(amount = 0.25, "Normalization...")
                  dbg("[outliers_server] Normalization: dim.X = ", dim(X))
                  dbg("[outliers_server] Normalization: Normalizing data using ", input$scaling_method)
                  normCounts <- playbase::pgx.countNormalization(counts, method = input$scaling_method)
                  if(input$scaling_method == "CPM") {
                      normCounts <- normCounts + 1
                  } else {
                      normCounts <- normCounts + 1e-20
                  }
                  X <- log2(normCounts)
                  if (input$quantile_norm) {
                      dbg("[outliers_server] Applying quantile normalization")
                      shiny::incProgress(amount = 0.25, "Quantile normalization...")
                      ## Ignore nearly 0 values
                      ## jj <- which(X < 0.01)
                      ## X[jj] <- NA
                      X <- limma::normalizeQuantiles(X)
                      ## X[jj] <- 0
                  }
              })
          }
          dbg("[outliers_server:normalizedX] dim.counts = ", dim(counts))
          dbg("[outliers_server:normalizedX] dim.normalizedX = ", dim(X))

          ## zeros <- sum(X==0, na.rm=TRUE)
          ## negs <- sum(X<0, na.rm=TRUE)
          ## nmissing <- sum(is.na(X))
          ## infin <- sum(is.infinite(X))
          ## dbg("[outliers_server] Normalized data have ", zeros, " zero values.")
          ## dbg("[outliers_server] Normalized data have ", negs, " negative values.")
          ## dbg("[outliers_server] Normalized data have ", nmissing, " missing values.")
          ## dbg("[outliers_server] Normalized data have ", infin, " infinite values.")

          X
      })

      ## Remove outliers
      cleanX <- reactive({
          shiny::req(dim(normalizedX()))
          X <- normalizedX()
          dbg("[outliers_server] Removing outliers step: dim.X = ", dim(X))
          if(input$remove_outliers) {
              dbg("[outliers_server] Removing outliers: Threshold = ",
                  input$outlier_threshold)
              res <- playbase::detectOutlierSamples(X, plot = FALSE)
              is.outlier <- (res$z.outlier > input$outlier_threshold)
              if (any(is.outlier) && !all(is.outlier)) {
                  X <- X[, which(!is.outlier), drop = FALSE]
                  ## also filter counts?
              }
          }   
          pos <- NULL
          if (NCOL(X) > 1) {
               pos <- irlba::irlba(X, nv = 2)$v
               rownames(pos) <- colnames(X)
          }
          dbg("[outliers_server] dim.cleanX = ", dim(X))
          ## zeros <- sum(X==0, na.rm=TRUE)
          ## negs <- sum(X<0, na.rm=TRUE)
          ## nmissing <- sum(is.na(X))
          ## infin <- sum(is.infinite(X))
          ## dbg("[outliers_server] cleanX data have ", zeros, " zero values.")
          ## dbg("[outliers_server] cleanX data have ", negs, " negative values.")
          ## dbg("[outliers_server] cleanX data have ", nmissing, " missing values.")
          ## dbg("[outliers_server] cleanX data have ", infin, " infinite values.")

          list(X = X, pos = pos)
      })

      ## Technical and biological effects correction
      correctedX <- shiny::reactive({
          shiny::req(dim(cleanX()$X), dim(imputedX()$counts),
                     dim(r_contrasts()), dim(r_samples()))
          
          ## Technical & biological effects correction ## ps: also includes biological
          ## shiny::incProgress(amount = 0.25, "Correcting for technical effects...")
          ## pheno <- playbase::contrasts2pheno(contrasts, samples)
          ## X <- playbase::removeTechnicalEffects(
          ##  X, samples, pheno,
          ##  p.pheno = 0.05, p.pca = 0.5, force = FALSE,
          ##  params = c("lib", "mito", "ribo", "cellcycle", "gender"),
          ##  nv = 2, k.pca = 10, xrank = NULL)        

          ## recompute chosed correction method with full
          ## matrix. previous was done on shortened matrix.
          X1 <- cleanX()$X
          counts <- imputedX()$counts
          samples <- r_samples()
          contrasts <- r_contrasts()
          kk <- intersect(colnames(X1), rownames(samples))
          kk <- intersect(kk, rownames(contrasts))
          counts <- counts[, kk, drop = FALSE]
          X1 <- X1[, kk, drop = FALSE]
          contrasts <- contrasts[kk, , drop = FALSE]
          samples <- samples[kk, , drop = FALSE]

          m <- input$bec_method
          dbg("[outliers_server]: Batch correction method = ", m)
          if(m == "uncorrected (default)") {
              m <- "uncorrected"
              dbg("[outliers_server]: Data not corrected for (potential) batch effects")
              cx <- X1
          } else {
              mm <- unique(c("uncorrected", input$bec_method))
              pars <- get_model_parameters()
              xlist <- playbase::runBatchCorrectionMethods(
                                     X = X1,
                                     batch = pars$batch,
                                     y = pars$pheno,
                                     controls = NULL,
                                     methods = mm,
                                     ntop = Inf,
                                     sc = FALSE,
                                     remove.failed = TRUE
                                 )
              shiny::removeModal()
              dbg("[outliers_server] names.xlist = ", names(xlist))
              cx <- xlist[[m]]
              dbg("[outliers_server] dim.correctedX = ", dim(cx))
          }
          ## pmax(2**cx - 1, 0)
          ## cx <- 2 ** cx - 1
          return(cx)
      })

      ## return object
      correctedCounts <- reactive({
          X <- correctedX()
          counts <- 2 ** X
          dbg("[outliers_server] dim.correctedCounts = ", dim(counts))
          counts
      })

      ## ------------------------------------------------------------------
      ## Compute reactive
      ## ------------------------------------------------------------------

      results_correction_methods <- reactive({
        ## shiny::req(cleanX()$X, r_contrasts(), r_samples())
        shiny::req(dim(cleanX()$X), dim(r_contrasts()), dim(r_samples()))

        X0 <- imputedX()$X
        X1 <- cleanX()$X
        samples <- r_samples()
        contrasts <- r_contrasts()
        kk <- intersect(colnames(X1), colnames(X0))
        kk <- intersect(kk, rownames(samples))
        kk <- intersect(kk, rownames(contrasts))
        X1 <- X1[, kk, drop = FALSE]
        X0 <- X0[, kk, drop = FALSE]
        contrasts <- contrasts[kk, , drop = FALSE]
        samples <- samples[kk, , drop = FALSE]
        xlist.init <- list("raw" = X0, "normalized" = X1)

        shiny::withProgress(message = "Comparing batch-correction methods...", value = 0.3, {
          res <- playbase::compare_batchcorrection_methods(
            X1, samples,
            pheno = NULL, contrasts = contrasts,
            methods = c("ComBat", "RUV", "SVA", "NPM"),
            ntop = 4000, xlist.init = xlist.init
          )
        })

        selected <- res$best.method
        shiny::updateSelectInput(session, "bec_method", selected = selected)

        return(res)
      })

      results_outlier_methods <- eventReactive(
        {
          list(normalizedX())
        },
        {
          X <- normalizedX()
          shiny::validate(shiny::need(!is.null(X), "no data. please upload."))
          shiny::validate(shiny::need(!is.null(nrow(X)), "no data. please upload."))

          X <- head(X[order(-matrixStats::rowSds(X)), ], 1000)
          out <- playbase::detectOutlierSamples(X, plot = FALSE)

          nb <- min(30, dim(X) / 5)
          scaledX <- t(scale(t(scale(t(X), scale = FALSE))))
          corX <- cor(t(scaledX))

          ## standard dim reduction methods
          pos <- list()
          ## pos[['tsne']] <- Rtsne::Rtsne(scaledX, check_duplicates=FALSE, perplexity=nb)$Y
          pos[["pca"]] <- irlba::irlba(scaledX, nu = 2, nv = 0)$u
          ## pos[['umap']] <- uwot::umap(scaledX, n_neighbors = ceiling(nb/2))
          for (i in 1:length(pos)) {
            rownames(pos[[i]]) <- rownames(scaledX)
            colnames(pos[[i]]) <- paste0(names(pos)[i], "_", 1:2)
          }

          out$pos <- pos
          out$corX <- corX
          out
        }
      )

      ## ------------------------------------------------------------------
      ## Plot functions
      ## ------------------------------------------------------------------

      plot_normalization <- function() {
        rX <- r_X()
        X0 <- imputedX()$X
        X1 <- normalizedX()

        if (input$norm_zero_na) {
          which.zero <- which(rX == 0)
          X0[X0 == 0] <- NA
          X1[X1 == 0] <- NA
        }

        if (input$norm_plottype == "boxplot") {
          if (ncol(X0) > 40) {
            jj <- sample(ncol(X0), 40)
            ii <- rownames(X0) ## names!
            ## just downsampling for boxplots
            if (length(ii) > 2000) ii <- sample(ii, 2000)
            X0 <- X0[ii, jj]
            X1 <- X1[ii, jj]
            rX <- rX[ii, jj]
          }

          ymax <- max(max(X0, na.rm = TRUE), max(X1, na.rm = TRUE))
          ymin <- quantile(X0[which(rX > 0)], probs = 0.001, na.rm = TRUE)
          if (ymin > 0) ymin <- 0
          dy <- 0.1 * (ymax - ymin)
          ylim <- c(ymin - dy, ymax + dy)

          par(mfrow = c(1, 2), mar = c(3.2, 6, 2, 0.5),
              mgp = c(1.8, 0.2, 0), tcl = -0.2)
          boxplot(X0,
            main = "raw", ylim = ylim, las = 2,
            ylab = "expression (log2)", xlab = "", cex.axis=0.8
          )
          boxplot(X1,
            main = "normalized", ylim = ylim, las = 2, 
            ylab = "", xlab = "", cex.axis=0.8
          )
        }

        if (input$norm_plottype == "histogram") {
          xmax0 <- quantile(X0[which(rX > 0)], probs = 0.999, na.rm = TRUE)
          xmax1 <- quantile(X1[which(rX > 0)], probs = 0.999, na.rm = TRUE)
          xmin0 <- quantile(X1[which(rX > 0)], probs = 0.001, na.rm = TRUE)
          xmin1 <- quantile(X1[which(rX > 0)], probs = 0.001, na.rm = TRUE)
          xmin0 <- min(xmin0, 0)
          xmin1 <- min(xmin1, 0)
          xlim0 <- c(xmin0, xmax0)
          xlim1 <- c(xmin1, xmax1)
          par(mfrow = c(1, 2), mar = c(3.2, 3, 2, 0.5), mgp = c(2.1, 0.8, 0))
          hist(X0,
            breaks = 70, main = "raw", xlim = xlim0,
            las = 1, xlab = "expression (log2)"
          )
          hist(X1,
            breaks = 60, main = "normalized", xlim = xlim1,
            las = 1, xlab = "expression (log2)", ylab = ""
          )
        }

        if (input$norm_plottype == "density") {
          xmax0 <- quantile(X0[which(rX > 0)], probs = 0.999, na.rm = TRUE)
          xmax1 <- quantile(X1[which(rX > 0)], probs = 0.999, na.rm = TRUE)
          xmin0 <- quantile(X1[which(rX > 0)], probs = 0.001, na.rm = TRUE)
          xmin1 <- quantile(X1[which(rX > 0)], probs = 0.001, na.rm = TRUE)
          xmin0 <- min(xmin0, 0)
          xmin1 <- min(xmin1, 0)
          xlim0 <- c(xmin0, xmax0)
          xlim1 <- c(xmin1, xmax1)

          par(mfrow = c(1, 2), mar = c(3.2, 3, 2, 0.5), mgp = c(2.1, 0.8, 0))
          playbase::gx.hist(X0,
            breaks = 70, main = "raw", xlim = xlim0,
            las = 1, xlab = "expression (log2)"
          )
          playbase::gx.hist(X1,
            breaks = 60, main = "normalized", xlim = xlim1,
            las = 1, xlab = "expression (log2)", ylab = ""
          )
        }
      }

      ## missing values
      plot_missingvalues <- function() {
          X0 <- r_X()
          X1 <- imputedX()$X
          X0 <- X0[rownames(X1), ] ## remove duplicates
          
          if(!any(is.na(X0)) & !any(is.na(X1))) {
              plot.new()
              text(0.5, 0.5, "No missing values", cex=1.2)
          } else if (any(is.na(X0)) & !any(is.na(X1))) {
              X0[!is.na(X0)] <- 2
              X0[is.na(X0)] <- 1
              par(mfrow = c(1, 2), mar = c(3.2, 3.2, 1.5, 0.5), mgp = c(2.2, 0.85, 0))
              playbase::gx.imagemap(X0, cex = -1)
              title("Missing values patterns in raw data", cex.main = 0.8)
              X1[!is.na(X1)] <- 2
              X1[is.na(X1)] <- 1
              playbase::gx.imagemap(X1, cex = -1)
              title("No missing values in imputed data", cex.main = 0.8)
          } else {
              
              ii <- which(is.na(X0))
              if (isolate(input$zero_as_na)) {
                  ii <- which(is.na(X0) | X0 == 0)
              }
              q999 <- quantile(X1, probs = 0.999, na.rm = TRUE)[1]
              X1[X1 > q999] <- NA
              h <- hist(X1, breaks = 80, plot = FALSE, las = 1)
              hh <- h$breaks

              ## set zero value to 1, NA values to 2
              X2 <- 1 * is.na(X0)
              ## X2 <- X2 + 1 * (!is.na(X0) & X0 == 0)
              if (input$zero_as_na) X2[X0 == 0] <- 1
              ## jj <- head( order(-rowMeans(is.na(X0))), 200)
              jj <- head(order(-apply(X2, 1, sd)), 200)
              X2 <- X2[jj, ]

              par(mfrow = c(1, 2), mar = c(3.2, 3.2, 0.8, 0.5), mgp = c(2.2, 0.85, 0))

              if (length(ii) > 0) {
                  hist(X1[-ii], breaks = hh, main = "", xlab = "expression (log2CPM)", las = 1)
                  hist(X1[ii], breaks = hh, add = TRUE, col = "red", las = 1)
              } else {
                  hist(X1, breaks = hh, main = "", xlab = "expression (log2CPM)", las = 1)
              }

              if (input$missing_plottype == "heatmap") {
                  if (any(X2 > 0)) {
                      ## NA heatmap
                      par(mar = c(3, 3, 2, 2), mgp = c(2.5, 0.85, 0))
                      playbase::gx.imagemap(X2, cex = -1)
                      title("missing values patterns", cex.main = 1.2)
                  } else {
                      plot.new()
                      text(0.5, 0.5, "no missing values")
                  }
              }

              if (input$missing_plottype == "ratio plot") {
                  if (any(X2 > 0)) {
                      ## NA ratio plot
                      par(mar = c(3, 3, 2, 2), mgp = c(2.0, 0.75, 0))
                      x.avg <- rowMeans(X1, na.rm = TRUE)
                      x.nar <- rowMeans(is.na(X0))
                      x.avg2 <- cut(x.avg, breaks = 20)
                      x.nar2 <- tapply(1:nrow(X0), x.avg2, function(i) mean(is.na(X0[i, , drop = FALSE])))
                      aa <- sort(unique(as.numeric(gsub(".*,|\\]", "", as.character(x.avg2)))))
                      barplot(rbind(x.nar2, 1 - x.nar2),
                              beside = FALSE, names.arg = aa, las = 1,
                              xlab = "average intensity (log2)", ylab = "missing value ratio"
                              )
                      title("missingness vs. average intensity")
                  } else {
                      plot.new()
                      text(0.5, 0.5, "no missing values")
                  }
              }
          }
      }

      ## sample outlier PCA plot
      plot.outlierPCA <- function(pos, z, z0, shownames) {
        is.outlier <- (z > z0)
        col1 <- "grey70"
        ## col1 <- res.outliers$dbscan$cluster + 1
        cex1 <- cut(nrow(pos),
          breaks = c(0, 40, 100, 250, 1000, 999999),
          c(1, 0.85, 0.7, 0.55, 0.4)
        )
        cex1 <- 3 * as.numeric(as.character(cex1))
        pos <- playbase::uscale(pos)

        ## How about plotly??
        plot(pos,
          col = col1, cex = 0.8 * cex1, pch = 20, las = 1,
          xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.1),
          xlab = "PC1", ylab = "PC2", main = "outliers"
        )

        if (shownames) {
          pos1 <- pos
          j <- which(is.outlier)
          if (length(j)) pos1 <- pos[-j, , drop = FALSE]
          text(pos1, rownames(pos1), cex = 0.85, offset = 0.8, pos = 1:4)
        }

        if (any(is.outlier)) {
          j <- which(is.outlier)
          points(pos[j, , drop = FALSE], col = "red", cex = 0.8 * cex1, lwd = 3, pch = 1)
          outlier.name <- rownames(pos)[j]
          text(pos[j, 1], pos[j, 2], outlier.name, cex = 1.0, offset = 0.8, pos = 1:4)
        }
      }

      ## sample outlier scores
      plot_outliers <- function() {
        res <- results_outlier_methods()
        z0 <- as.numeric(input$outlier_threshold)
        zscore <- res$z.outlier
        Z <- res$Z
        pos <- res$pos[["pca"]]
        ##        plottype <- input$outlier_plottype
        plottype <- "pca"
        if (plottype == "pca") {
          par(mfrow = c(1, 2), mar = c(3.2, 3, 2, 0.5), mgp = c(2.1, 0.8, 0))
          barplot(zscore,
            main = "outlier score",  ylab = "z-score",
            las = 1, ylim = c(0, max(7, 1.2 * max(Z))),
          )
          abline(h = z0, lty = 3, lwd = 1.5, col = "red")
          plot.outlierPCA(pos, zscore, z0, input$outlier_shownames)
        }

        if (plottype == "heatmap") {
          par(mfrow = c(1, 2), mar = c(0, 3, 0, 1), mgp = c(2.1, 0.8, 0))
          playbase::gx.heatmap(res$corX,
            sym = TRUE, mar = c(1, 12), keysize = 0.4,
            cexCol = 0.0001, scale = "none", key = FALSE
          )
        }
      }

      plot_correction <- function() {
        if (input$bec_preview_all == FALSE) {
          plot_before_after()
        }
        if (input$bec_preview_all == TRUE) {
          plot_all_methods()
        }
      }

      plot_all_methods <- function() {
        res <- results_correction_methods()
        pos.list <- res$pos[["tsne"]]
        pheno <- res$pheno
        xdim <- length(res$pheno)
        col1 <- factor(pheno)
        cex1 <- cut(xdim,
          breaks = c(0, 40, 100, 250, 1000, 999999),
          c(1, 0.85, 0.7, 0.55, 0.4)
        )
        cex1 <- 2.5 * as.numeric(as.character(cex1))

        par(mfrow = c(2, 3), mar = c(2, 2, 2, 1))
        for (i in 1:length(pos.list)) {
          plot(pos.list[[i]], col = col1, cex = cex1, pch = 20)
          title(names(pos.list)[i], cex.main = 1.5)
        }
      }

      plot_before_after <- function() {
        ##        out.res <- results_outlier_methods()
        res <- results_correction_methods()

        method <- input$bec_method
        if (method == "uncorrected") method <- "normalized"
        pos0 <- res$pos[["tsne"]][["normalized"]]
        pos1 <- res$pos[["tsne"]][[method]]

        kk <- intersect(rownames(pos0), rownames(pos1))
        pos0 <- pos0[kk, ]
        pos1 <- pos1[kk, ]

        ## pheno <- r_contrasts()[,1]
        pheno <- playbase::contrasts2pheno(r_contrasts(), r_samples())
        pheno <- pheno[rownames(pos0)]
        col1 <- factor(pheno)
        cex1 <- cut(nrow(pos1),
          breaks = c(0, 40, 100, 250, 1000, 999999),
          c(1, 0.85, 0.7, 0.55, 0.4)
        )
        cex1 <- 2.7 * as.numeric(as.character(cex1))

        par(mfrow = c(1, 2), mar = c(3.2, 3, 2, 0.5), mgp = c(2.1, 0.8, 0))
        plot(pos0,
          col = col1, pch = 20, cex = 1.0 * cex1, main = "before",
          xlab = "PC1", ylab = "PC2"
        )
        plot(pos1,
          col = col1, pch = 20, cex = 1.0 * cex1, main = "after",
          xlab = "PC1", ylab = "PC2"
        )
      }

      ## ------------------------------------------------------------------
      ## Plot modules
      ## ------------------------------------------------------------------

      PlotModuleServer(
        "plot1",
        plotlib = "base",
        func = plot_normalization,
        ##      func2 = plot.RENDER,
        ##      csvFunc = plot_data,
        res = c(75, 120),
        pdf.width = 12,
        pdf.height = 6,
        add.watermark = FALSE
      )

      PlotModuleServer(
        "plot2",
        plotlib = "base",
        func = plot_missingvalues,
        res = c(75, 120),
        pdf.width = 12,
        pdf.height = 6,
        add.watermark = FALSE
      )

      PlotModuleServer(
        "plot3",
        plotlib = "base",
        func = plot_outliers,
        ##      func2 = plot.RENDER,
        ##      csvFunc = plot_data,
        res = c(75, 120),
        pdf.width = 12,
        pdf.height = 6,
        add.watermark = FALSE
      )

      PlotModuleServer(
        "plot4",
        plotlib = "base",
        func = plot_correction,
        ## func2 = plot.RENDER,
        ## csvFunc = plot_data,
        res = c(75, 120),
        pdf.width = 12,
        pdf.height = 6,
        add.watermark = FALSE
      )

      correctedCounts <- shiny::reactive({
          shiny::req(imputedX()$counts)
          dbg("[outliers_server] dim.counts = ", dim(imputedX()$counts))
          return(imputedX()$counts)
      })
      
      return(
        list(
            counts = correctedCounts,
            X = correctedX,
            results = results_correction_methods
        )
      ) ## pointing to reactive
    } ## end-of-server
  )
}
