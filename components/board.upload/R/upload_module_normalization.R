##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


## =========================================================================
## ==================== NORMALIZATION UI/SERVER =================================
## =========================================================================


upload_module_normalization_ui <- function(id, height = "100%") {
  ns <- shiny::NS(id)
  uiOutput(ns("normalization"), fill = TRUE)
}

upload_module_normalization_server <- function(
    id,
    r_counts,
    r_samples,
    r_contrasts,
    upload_datatype,
    is.count = FALSE,
    height = 720) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns


      observeEvent(input$normalization_method, {
        shiny::req(input$normalization_method == "reference")
        gg <- sort(rownames(r_counts()))
        shiny::updateSelectizeInput(
          session, "ref_gene",
          choices = gg,
          selected = character(0), server = TRUE
        )
      })


      ## ------------------------------------------------------------------
      ## Object reactive chain
      ## ------------------------------------------------------------------

      ## Impute and remove duplicated features
      imputedX <- reactive({
        shiny::req(r_counts())
        shiny::req(!is.null(input$zero_as_na))
        counts <- r_counts()

        counts[which(is.nan(counts))] <- NA
        counts[which(is.infinite(counts))] <- NA

        negs <- sum(counts < 0, na.rm = TRUE)
        if (negs > 0) {
          counts <- pmax(counts, 0) ## NEED RETHINK (eg: what about Olink NPX)
        }

        if (input$zero_as_na) {
          dbg("[normalization_server:imputedX] Setting 0 values to NA")
          counts[which(counts == 0)] <- NA
        }

        m <- input$normalization_method
        prior <- ifelse(m %in% c("CPM", "CPM+quantile"), 1, 1e-4) ## NEW
        X <- log2(counts + prior) ## NEED RETHINK

        ## if (input$remove_xxl) {
        ##     dbg("[normalization_server:imputedX]: Assign NA to outlier features")
        ##     X[playbase::is.xxl(X, z = 10)] <- NA
        ## }

        nmissing <- sum(is.na(X))
        dbg("[normalization_server:imputedX] X has ", nmissing, " missing values (NAs).")
        if (nmissing > 0 && input$impute) {
          m <- input$impute_method
          dbg("[normalization_server:imputedX] Imputing data using ", m)
          X <- playbase::imputeMissing(X, method = m)
          dbg("[normalization_server:imputedX] dim.imputedX = ", dim(X))
        } else {
          dbg("[normalization_server:imputedX] No imputation.")
        }

        dbg("[normalization_server:imputedX] Checking for duplicated features")
        X <- playbase::counts.mergeDuplicateFeatures(X, is.counts = FALSE)
        X
      })

      ## Normalize
      normalizedX <- reactive({
        shiny::req(dim(imputedX()))
        X <- imputedX() ## can be imputed or not (see above). log2. Can have negatives.
        if (input$normalize) {
          m <- input$normalization_method
          dbg("[normalization_server:normalizedX] Normalizing data using ", m)
          ## NEED RETHINK: would be better to rewrite Normalization in log2-space (IK)
          ref <- NULL
          if (m == "reference") {
            ref <- input$ref_gene
            shiny::validate(shiny::need(
              isTruthy(ref),
              tspan("Please select reference gene", js = FALSE)
            ))
            shiny::req(ref)
          }
          prior <- ifelse(m %in% c("CPM", "CPM+quantile"), 1, 1e-4)
          m0 <- m
          if (m == "CPM+quantile") m0 <- "CPM"
          normCounts <- playbase::pgx.countNormalization(
            pmax(2**X - prior, 0),
            method = m0, ref = ref
          )
          X <- log2(normCounts + prior)
          ## if (FALSE && input$quantile_norm) {
          if (m == "CPM+quantile") {
            dbg("[normalization_server:normalizedX] Applying quantile normalization")
            X <- limma::normalizeQuantiles(X)
          }
        } else {
          dbg("[normalization_server:normalizedX] Skipping normalization")
        }
        return(X)
      })

      ## Remove outliers
      cleanX <- reactive({
        shiny::req(dim(normalizedX()))
        X <- normalizedX()

        if (input$remove_outliers) {
          threshold <- input$outlier_threshold
          dbg("[normalization_server:cleanX] Removing outliers: Threshold = ", threshold)
          res <- playbase::detectOutlierSamples(X, plot = FALSE)
          is.outlier <- (res$z.outlier > threshold)
          if (any(is.outlier) && !all(is.outlier)) {
            X <- X[, which(!is.outlier), drop = FALSE] ## also filter counts?
          }
        }
        pos <- NULL
        if (NCOL(X) > 1) {
          pos <- irlba::irlba(X, nv = 2)$v
          rownames(pos) <- colnames(X)
        }
        dbg("[normalization_server:cleanX] dim.cleanX = ", dim(X))
        list(X = X, pos = pos)
      })

      ## Technical and biological effects correction
      correctedX <- shiny::reactive({
        shiny::req(dim(cleanX()$X), dim(r_contrasts()), dim(r_samples()))
        X1 <- cleanX()$X
        samples <- r_samples()
        contrasts <- r_contrasts()

        ## recompute chosed correction method with full
        ## matrix. previous was done on shortened matrix.
        kk <- intersect(colnames(X1), rownames(samples))
        kk <- intersect(kk, rownames(contrasts))
        X1 <- X1[, kk, drop = FALSE]
        contrasts <- contrasts[kk, , drop = FALSE]
        samples <- samples[kk, , drop = FALSE]

        nmissing <- sum(is.na(X1))
        if (!input$batchcorrect) {
          dbg("[normalization_server:correctedX] Data not corrected for batch effects")
          if (nmissing == 0) {
            cx <- list(X = X1)
          } else {
            dbg("[normalization_server:correctedX] ", nmissing, " missing values in X1.")
            dbg("[normalization_server:correctedX] Generating an internal, SVD2-imputed matrix")
            impX1 <- playbase::imputeMissing(X1, method = "SVD2")
            cx <- list(X = X1, impX1 = impX1)
          }
        } else {
          m <- input$bec_method
          dbg("[normalization_server:correctedX] Batch correction method = ", m)
          mm <- unique(c("uncorrected", m))

          ##  pars <- get_model_parameters()
          pars <- playbase::get_model_parameters(X1, samples, pheno = NULL, contrasts)
          batch.pars <- input$bec_param
          if (any(grepl("<autodectect>", batch.pars))) {
            batch.pars <- pars$batch.pars
          }
          if (any(grepl("<none>", batch.pars))) {
            batch.pars <- ""
          }
          batch.pars <- intersect(batch.pars, colnames(samples))

          if (length(batch.pars)) {
            batch <- samples[, batch.pars, drop = FALSE] ## matrix
          } else {
            batch <- NULL
          }

          pheno <- pars$pheno
          if (nmissing == 0) {
            xlist <- playbase::runBatchCorrectionMethods(
              X = X1,
              batch = batch,
              y = pheno,
              methods = mm,
              ntop = Inf
            )
            cx <- list(X = xlist[[m]])
          } else {
            dbg("[normalization_server:correctedX] missing values in X1")
            impX1 <- playbase::imputeMissing(X1, method = "SVD2")
            xlist <- playbase::runBatchCorrectionMethods(
              X = impX1,
              batch = batch,
              y = pheno,
              methods = mm,
              ntop = Inf
            )
            bc_impX1 <- xlist[[m]] ## Batch corrected, imputed
            jj <- which(is.na(X1), arr.ind = TRUE)
            xlist[[m]][jj] <- NA ## Batch corrected, with original NAs restored
            cx <- list(X = xlist[[m]], impX1 = bc_impX1)
          }
        }
        shiny::removeModal()

        return(cx)
      })

      ## return object
      correctedCounts <- reactive({
        shiny::req(dim(correctedX()$X))
        X <- correctedX()$X
        prior <- ifelse(input$normalization_method %in% c("CPM", "CPM+quantile"), 1, 1e-4)
        dbg("[normalization_server:correctedCounts] Generating counts. Prior=", prior)
        counts <- pmax(2**X - prior, 0)
        counts
      })

      ## ------------------------------------------------------------------
      ## Compute reactive
      ## ------------------------------------------------------------------
      results_correction_methods <- reactive({
        shiny::req(dim(cleanX()$X), dim(r_contrasts()), dim(r_samples()))

        X0 <- imputedX()
        X1 <- cleanX()$X ## normalized+cleaned
        samples <- r_samples()
        contrasts <- r_contrasts()
        batch.pars <- input$bec_param

        nmissing0 <- sum(is.na(X0))
        if (nmissing0 > 0) {
          X0 <- playbase::imputeMissing(X0, method = "SVD2")
        }

        nmissing1 <- sum(is.na(X1))
        if (nmissing1 > 0) {
          X1 <- playbase::imputeMissing(X1, method = "SVD2")
        }

        kk <- intersect(colnames(X1), colnames(X0))
        kk <- intersect(kk, rownames(samples))
        kk <- intersect(kk, rownames(contrasts))
        X1 <- X1[, kk, drop = FALSE]
        X0 <- X0[, kk, drop = FALSE]
        contrasts <- contrasts[kk, , drop = FALSE]
        samples <- samples[kk, , drop = FALSE]

        if (any(grepl("<autodetect>", batch.pars))) batch.pars <- "<autodetect>"
        if (any(grepl("<none>", batch.pars))) batch.pars <- NULL

        methods <- c("ComBat", "limma", "RUV", "SVA", "NPM")
        xlist.init <- list("uncorrected" = X0, "normalized" = X1)

        shiny::withProgress(
          message = "Comparing batch-correction methods...",
          value = 0.3,
          {
            res <- playbase::compare_batchcorrection_methods(
              X1, samples,
              pheno = NULL,
              contrasts = contrasts,
              batch.pars = batch.pars,
              clust.method = "tsne",
              methods = methods,
              evaluate = FALSE, ## no score computation
              xlist.init = xlist.init
            )
          }
        )

        ## ## take out failed methods
        ## xlist.ok <- sapply(res$xlist, function(x) !any(class(x)=="try-error"))
        ## pos.ok <- sapply(res$pos, function(x) !any(class(x)=="try-error"))
        ## res$xlist <- res$xlist[which(xlist.ok && pos.ok)]
        ## res$pos <- res$pos[which(xlist.ok && pos.ok)]

        return(res)
      })

      ## Remove?
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
          pos[["pca"]] <- irlba::irlba(scaledX, nu = 2, nv = 0)$u
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
        rX <- r_counts()
        X0 <- imputedX()
        ## X1 <- normalizedX()
        X1 <- cleanX()$X
        main.tt <- ifelse(input$normalize, norm_method(), "no normalization")

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

          par(mfrow = c(1, 2), mar = c(6, 3, 2, 0.5), mgp = c(2.1, 0.8, 0))
          boxplot(
            X0,
            main = "raw",
            ylim = range(X0, na.rm = TRUE) + 0.2 * c(-1, 1) * diff(range(X0, na.rm = TRUE)),
            las = 2,
            ylab = tspan("counts (log2)", js = FALSE),
            xlab = "",
            cex.axis = 0.8,
            cex = 0.5
          )

          boxplot(
            X1,
            main = main.tt,
            ylim = range(X1, na.rm = TRUE) + 0.2 * c(-1, 1) * diff(range(X1, na.rm = TRUE)),
            las = 2,
            ylab = "",
            xlab = "",
            cex.axis = 0.8,
            cex = 0.5
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
            las = 1, xlab = tspan("counts (log2)", js = FALSE)
          )
          hist(X1,
            breaks = 60, main = main.tt, xlim = xlim1,
            las = 1, xlab = tspan("counts (log2)", js = FALSE), ylab = ""
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
            las = 1, xlab = tspan("counts (log2)", js = FALSE)
          )

          playbase::gx.hist(X1,
            breaks = 60, main = main.tt, xlim = xlim1,
            las = 1, xlab = tspan("counts (log2)", js = FALSE), ylab = ""
          )
        }
      }

      ## missing values
      plot_missingvalues <- function() {
        X0 <- r_counts()
        X1 <- imputedX()
        X0 <- X0[rownames(X1), ] ## remove duplicates

        has.zeros <- any(X0 == 0, na.rm = TRUE)
        if (!any(is.na(X0)) && !(input$zero_as_na && has.zeros)) {
          plot.new()
          text(0.5, 0.5, "No missing values", cex = 1.2)
        } else if (FALSE && any(is.na(X0)) && !any(is.na(X1))) {
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
          if (input$zero_as_na) X2[X0 == 0] <- 1
          jj <- head(order(-apply(X2, 1, sd)), 200)
          X2 <- X2[jj, ]

          par(mfrow = c(1, 2), mar = c(3.2, 3.2, 0.8, 0.5), tcl = -0.15, mgp = c(2.2, 0.2, 0))

          if (length(ii) > 0) {
            hist(X1[-ii], breaks = hh, main = "", xlab = "expression (log2)", las = 1)
            hist(X1[ii], breaks = hh, add = TRUE, col = "red", las = 1)
          } else {
            hist(X1, breaks = hh, main = "", xlab = "expression (log2)", las = 1)
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
        ## plottype <- input$outlier_plottype
        plottype <- "pca"
        if (plottype == "pca") {
          par(mfrow = c(1, 2), mar = c(3.2, 3, 2, 0.5), mgp = c(2.1, 0.8, 0))
          barplot(zscore,
            main = "outlier score", ylab = "z-score",
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
        if (input$batchcorrect) {
          plot_before_after()
        } else {
          plot_all_methods()
        }
      }

      plot_all_methods <- function() {
        res <- results_correction_methods()
        out.res <- results_outlier_methods()
        shiny::req(res)
        shiny::req(out.res)

        methods <- c("uncorrected", sort(c("ComBat", "limma", "RUV", "SVA", "NPM")))
        ## methods <- intersect(methods, names(res$pos))
        ## pos.list <- res$pos[methods]
        pos.list <- res$pos
        ## get same positions as after outlier detection
        pos0 <- out.res$pos[["pca"]]
        pos.list <- c(list("uncorrected" = pos0), pos.list)
        #        names(pos.list) <- sub("ComBat", "auto-ComBat", names(pos.list))
        #        names(pos.list) <- sub("limma", "auto-limma", names(pos.list))

        pheno <- res$pheno
        xdim <- length(pheno)
        col1 <- factor(pheno)
        cex1 <- cut(xdim,
          breaks = c(0, 40, 100, 250, 1000, 999999),
          c(1, 0.85, 0.7, 0.55, 0.4)
        )
        cex1 <- 2.5 * as.numeric(as.character(cex1))
        par(mfrow = c(2, 3), mar = c(2, 2, 2, 1))
        for (m in methods) {
          if (m %in% names(pos.list)) {
            plot(
              pos.list[[m]],
              col = col1,
              cex = cex1,
              pch = 20
            )
          } else {
            plot.new()
            text(0.45, 0.5, "method failed")
          }
          title(m, cex.main = 1.5)
        }
      }

      plot_before_after <- function() {
        out.res <- results_outlier_methods()
        res <- results_correction_methods()
        ## get same positions as after outlier detection
        ## pos0 <- res$pos[["normalized"]]
        pos0 <- out.res$pos[["pca"]]
        method <- input$bec_method

        if (!method %in% names(res$pos)) {
          plot.new()
          text(0.45, 0.5, "method failed")
          return(NULL)
        }

        if (!input$batchcorrect) {
          pos1 <- pos0
        } else {
          pos1 <- res$pos[[method]]
        }

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
        #        method <- sub("ComBat", "auto-ComBat", method)
        #        method <- sub("limma", "auto-limma", method)
        par(mfrow = c(1, 2), mar = c(3.2, 3, 2, 0.5), mgp = c(2.1, 0.8, 0))
        plot(pos0,
          col = col1, pch = 20, cex = 1.0 * cex1, las = 1,
          main = "uncorrected", xlab = "PC1", ylab = "PC2"
        )
        plot(pos1,
          col = col1, pch = 20, cex = 1.0 * cex1, las = 1,
          main = method, xlab = "PC1", ylab = "PC2"
        )
      }

      ## ------------------------------------------------------------------
      ## Plot UI
      ## ------------------------------------------------------------------


      getBatchParams <- eventReactive(
        {
          list(r_counts(), r_samples(), r_contrasts())
        },
        {
          shiny::req(dim(r_counts()), dim(r_samples()), dim(r_contrasts()))
          X <- r_counts()
          samples <- r_samples()
          contrasts <- r_contrasts()
          pars <- playbase::get_model_parameters(
            X,
            samples,
            pheno = NULL,
            contrasts = contrasts
          )
          all.pars <- setdiff(colnames(samples), pars$pheno.pars)
          all.pars <- union(all.pars, pars$batch.pars)
          names(all.pars) <- ifelse(all.pars %in% pars$batch.pars,
            paste(all.pars, "*"), all.pars
          )
          ##        all.pars <- c("<autodetect>","<none>",all.pars)
          all.pars <- c("<autodetect>", all.pars)
          return(all.pars)
        }
      )





      output$normalization <- shiny::renderUI({
        ## reactive
        batch_params <- getBatchParams()

        score.infotext <-
          "Outliers markedly deviate from the vast majority of samples. Outliers could be caused by technical factors and negatively affect data analysis. Here, outliers are identified and marked for removal should you wish so."

        missing.infotext <-
          "Missing values (MVs) reduce the completeness of biological data and hinder preprocessing steps. MVs (i.e., NA), more often populate proteomics and metabolomics data. Here, MVs are identified and their patterns in your data is shown."

        normalization.infotext <-
          "Normalization enables to standardize the data and improve their consistency, comparability and reproducibility. Boxplots of raw (unnormalized) and normalized data are shown. Normalization method can be selected on the left, under “Normalization”."

        batcheff.infotext <-
          "Batch effects (BEs) are due to technical, experimental factors that introduce unwanted variation into the measurements. Here, BEs are detected and BEs correction is shown. BE correction methods can be selected on the left, under “Batch-effects correction”."

        missing.options <- tagList(
          shiny::radioButtons(ns("missing_plottype"), "Plot type:", c("heatmap", "ratio plot"),
            selected = "heatmap", inline = TRUE
          ),
        )

        norm.options <- tagList(
          shiny::radioButtons(
            ns("norm_plottype"),
            label = "Plot type:",
            choices = c("boxplot", "histogram"),
            selected = "boxplot", inline = FALSE
          )
        )

        outlier.options <- tagList(
          shiny::checkboxInput(ns("outlier_shownames"), "show sample names", FALSE)
        )

        bec.options <- tagList(
          shiny::radioButtons(ns("bec_plottype"), "Plot type:", c("pca", "tsne", "heatmap"),
            inline = TRUE
          )
        )

        navmenu <- tagList(
          bslib::card(bslib::card_body(
            style = "padding: 0px;",
            bslib::accordion(
              multiple = FALSE,
              style = "background-color: #F7FAFD99;",
              bslib::accordion_panel(
                title = "1. Missing values",
                shiny::p("Replace missing values using an imputation method:\n"),
                shiny::checkboxInput(ns("zero_as_na"), label = "Treat zero as NA", value = FALSE),
                shiny::checkboxInput(ns("impute"), label = "Impute missing values", value = TRUE),
                shiny::conditionalPanel(
                  "input.impute == true",
                  ns = ns,
                  shiny::selectInput(ns("impute_method"), NULL,
                    choices = c(
                      "SVDimpute" = "SVD2"
                      # "Zero" = "zero",
                      # "MinDet",
                      # "MinProb"
                      # "NMF"
                    ),
                    selected = "SVD2"
                  )
                ),
                ## shiny::checkboxInput(ns("remove_xxl"), label = "Treat XXL as NA", value = FALSE),
                br()
              ),
              bslib::accordion_panel(
                title = "2. Normalization",
                shiny::checkboxInput(ns("normalize"), label = "Normalize data", value = TRUE),
                shiny::conditionalPanel(
                  "input.normalize == true",
                  ns = ns,
                  shiny::selectInput(
                    ns("normalization_method"), NULL,
                    choices = if (grepl("proteomics", upload_datatype(), ignore.case = TRUE)) {
                      c(
                        "maxMedian", "maxSum", ## "TMM",
                        "reference"
                      )
                    } else {
                      c(
                        "CPM", "CPM+quantile", ## "quantile",
                        "maxMedian", "maxSum", ## "TMM",
                        "reference"
                      )
                    },
                    selected = ifelse(grepl("proteomics", upload_datatype(), ignore.case = TRUE),
                      "maxMedian", "CPM+quantile"
                    )
                  ),
                  shiny::conditionalPanel(
                    "input.normalization_method == 'reference'",
                    ns = ns,
                    shiny::selectizeInput(
                      ns("ref_gene"), NULL,
                      choices = NULL,
                      multiple = FALSE,
                      options = list(
                        placeholder = tspan("Choose gene...", js = FALSE)
                      )
                    )
                  )
                  ## shiny::checkboxInput(ns("quantile_norm"), "Add quantile normalization", value = TRUE)
                ),
                br()
              ),
              bslib::accordion_panel(
                title = "3. Remove outliers",
                shiny::p("Automatically detect and remove outlier samples."),
                shiny::checkboxInput(ns("remove_outliers"), "remove outliers", value = FALSE),
                shiny::conditionalPanel(
                  "input.remove_outliers == true",
                  ns = ns,
                  shiny::sliderInput(
                    ns("outlier_threshold"), "Select threshold:", 1, 12, 6, 1
                  )
                ),
                br()
              ),
              bslib::accordion_panel(
                title = "4. Batch-effect correction",
                shiny::p("Automatically remove unwanted variation from your data."),
                shiny::checkboxInput(ns("batchcorrect"),
                  label = "Remove batch effects",
                  value = FALSE
                ),
                shiny::conditionalPanel(
                  "input.batchcorrect == true",
                  ns = ns,
                  shiny::selectInput(
                    ns("bec_method"),
                    label = "Select method:",
                    choices = c(
                      "ComBat",
                      "limma",
                      "NPM" = "NPM",
                      "RUV" = "RUV",
                      "SVA" = "SVA"
                    ),
                    selected = "SVA"
                  ),
                  shiny::conditionalPanel(
                    "input.bec_method == 'ComBat' || input.bec_method == 'limma'",
                    ns = ns,
                    shiny::selectizeInput(
                      ns("bec_param"),
                      label = "Batch parameter:",
                      choices = batch_params, ## reactive
                      selected = batch_params[1],
                      multiple = TRUE,
                      options = list(
                        placeholder = "Select..."
                      )
                    ),
                    shiny::br()
                  )
                ),
                br()
              )
            ),
            br()
          ))
        )

        ## ---------------------------- UI ----------------------------------
        ui <- div(
          bslib::as_fill_carrier(),
          style = "width: 100%; display: flex; ",
          bslib::layout_columns(
            col_widths = c(2, 10),
            style = "margin-bottom: 20px;",
            heights_equal = "row",
            ## ----------- menu ------------
            navmenu,
            ## ----------- canvas ------------
            bslib::layout_columns(
              col_widths = c(6, 6),
              row_heights = c(3, 3),
              heights_equal = "row",
              PlotModuleUI(
                ns("plot2"),
                title = "Missing values",
                info.text = missing.infotext,
                caption = missing.infotext,
                options = missing.options,
                height = c("auto", "100%"),
                show.maximize = FALSE
              ),
              PlotModuleUI(
                ns("plot1"),
                title = "Normalization",
                options = norm.options,
                info.text = normalization.infotext,
                height = c("auto", "100%"),
                show.maximize = FALSE
              ),
              PlotModuleUI(
                ns("plot3"),
                title = "Outliers detection",
                info.text = score.infotext,
                caption = score.infotext,
                options = outlier.options,
                height = c("auto", "100%"),
                show.maximize = FALSE
              ),
              PlotModuleUI(
                ns("plot4"),
                title = "Batch-effects correction",
                options = NULL,
                info.text = batcheff.infotext,
                height = c("auto", "100%"),
                show.maximize = FALSE
              )
            )
          ),
          ## div(shiny::selectInput(ns("normalizationUI"),NULL,choices=TRUE),style='display:none')
          div(shiny::checkboxInput(ns("normalizationUI"), NULL, TRUE), style = "visibility:hidden")
        )

        return(ui)
      })

      ## ------------------------------------------------------------------
      ## Plot modules
      ## ------------------------------------------------------------------

      PlotModuleServer(
        "plot1",
        plotlib = "base",
        func = plot_normalization,
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
        res = c(75, 120),
        pdf.width = 12,
        pdf.height = 6,
        add.watermark = FALSE
      )

      PlotModuleServer(
        "plot4",
        plotlib = "base",
        func = plot_correction,
        res = c(75, 120),
        pdf.width = 12,
        pdf.height = 6,
        add.watermark = FALSE
      )

      cX <- reactive({
        shiny::req(dim(correctedX()$X))
        cX <- correctedX()$X
        cX
      })

      impX <- reactive({
        shiny::req(dim(correctedX()$X))
        if (length(correctedX()) == 2) {
          impX <- correctedX()$impX1
        } else {
          impX <- NULL
        }
        impX
      })

      norm_method <- reactive({
        m <- input$normalization_method
        if (!input$normalize) m <- "skip_normalization"
        m
      })

      return(
        list(
          counts = correctedCounts,
          X = cX,
          impX = impX,
          norm_method = norm_method
          ## results = results_correction_methods  ## IK reallz needed??
        )
      ) ## pointing to reactive
    } ## end-of-server
  )
}
