##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

upload_module_normalization_ui <- function(id, height = "100%") {
  ns <- shiny::NS(id)
  uiOutput(ns("normalization"), fill = TRUE)
}


upload_module_normalization_server <- function(
  id,
  r_counts,
  r_samples,
  r_contrasts,
  r_annot,
  upload_datatype,
  is.olink,
  meth_type,
  is.count = FALSE,
  height = 720,
  recompute_pgx = NULL
) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns

      zero_as_na <- function() isTRUE(input$zero_as_na)

      observeEvent(input$normalization_method, {
        shiny::req(input$normalization_method == "reference")
        gg <- sort(rownames(r_counts()))
        shiny::updateSelectizeInput(session, "ref_gene",
          choices = gg,
          selected = character(0), server = TRUE
        )
      })

      ## ImputedX
      imputedX <- reactive({
        shiny::req(dim(r_counts()), !is.null(input$normalize))
        counts <- r_counts()
        samples <- r_samples()
        contrasts <- r_contrasts()
        annot <- r_annot()
        shiny::req(dim(contrasts))
        ## shiny::req(dim(r_counts()))
        ## shiny::req(!is.null(input$zero_as_na))
        ## shiny::req(!is.null(input$normalize)) ## new

        counts[which(is.nan(counts))] <- NA
        counts[which(is.infinite(counts))] <- NA

        ## Olink NPX are passed on up to here unaltered.
        if (is.olink()) {
          dbg("[normalization_server:imputedX] Olink NPX Proteomics")
          counts <- 2**counts
          shiny::updateCheckboxInput(session, "normalize", value = FALSE)
        }

        if (any(counts < 0, na.rm = TRUE)) counts <- pmax(counts, 0)

        #if (input$zero_as_na) {
        if (zero_as_na()) {
          dbg("[normalization_server:imputedX] Setting 0 values to NA")
          counts[which(counts == 0)] <- NA
        }

        is.mox <- playbase::is.multiomics(rownames(counts))
        if (is.mox) {
          X <- counts
          dtypes <- unique(sub(":.*", "", rownames(X)))
          for (i in 1:length(dtypes)) {
            ii <- grep(paste0("^", dtypes[i], ":"), rownames(counts))
            prior <- 1
            if (dtypes[i] != "gx") prior <- playbase::getPrior(counts[ii, ])
            X[ii, ] <- log2(counts[ii, ] + prior)
          }
        } else {
          is.meth.beta <- FALSE
          vv <- range(counts, na.rm = TRUE)
          is.meth.beta <- (all(vv>=0 & vv<=1) & upload_datatype() == "methylomics")
          if (is.meth.beta) {
            X <- counts
            prior <- 0
          } else {
            dbg("-------------------------------M1")
            prior0 <- playbase::getPrior(counts)
            m <- input$normalization_method
            prior <- ifelse(grepl("CPM|TMM", m), 1, prior0)
            X <- log2(counts + prior)
          }
        }
        dbg("[normalization_server:imputedX] X has ", sum(is.na(X)), " missing values (NAs).")
        dbg("[normalization_server:imputedX] X has ", sum(rowSums(is.na(X)) > 0), " rows with NAs.")

        ## Filter probes for maximum missingness as required
        if (sum(is.na(X)) > 0 && input$filtermissing) {
          f <- input$filterthreshold
          dbg(paste0("[normalization_server:imputedX] Threshold NA filter: ", f))
          sample.contrasts <- playbase::contrasts.convertToLabelMatrix(contrasts, samples)
          grp <- apply(sample.contrasts, 1, paste, collapse = "_")
          if (f >= 1) {
            grp.sum <- tapply(1:ncol(counts), grp, function(i) {
              rx <- counts[, i, drop = FALSE]
              rowSums(!is.na(rx))
            })
            maxsum <- apply(do.call(cbind, grp.sum), 1, max, na.rm = TRUE)
            sel <- (maxsum >= 3)
          } else if (f < 0) {
            grp.avg <- tapply(1:ncol(counts), grp, function(i) {
              rx <- counts[, i, drop = FALSE]
              rowMeans(!is.na(rx))
            })
            maxavg <- apply(do.call(cbind, grp.avg), 1, max, na.rm = TRUE)
            sel <- (maxavg >= 0.5) # maxavg >= abs(f)
          } else {
            sel <- (rowMeans(is.na(X)) <= f)
          }
          dbg("[normalization_server:imputedX] nrows excluded due to NA: n=", sum(!sel))
          X <- X[which(sel), , drop = FALSE]
          counts <- counts[which(sel), , drop = FALSE]
          annot <- annot[which(sel), , drop = FALSE]
        }

        ## Impute if required
        if (any(is.na(X)) & input$impute) {
          if (is.mox) {
            X <- playbase::imputeMissing.mox(X, method = input$impute_method)
          } else {
            X <- playbase::imputeMissing(X, method = input$impute_method)
          }
        }

        return(list(counts = counts, X = X, prior = prior, annot = annot))

      })

      ## Normalize
      normalizedX <- reactive({

        shiny::req(dim(imputedX()$X))
        X <- imputedX()$X ## can be imputed or not. log2. Can have negatives.
        prior <- imputedX()$prior

        if (input$normalize) {
          m <- input$normalization_method
          ref <- NULL
          if (m == "reference") {
            ref <- input$ref_gene
            shiny::validate(shiny::need(isTruthy(ref), tspan("Please select reference gene", js = FALSE)))
            shiny::req(ref)
          }
          if (upload_datatype() == "multi-omics") {
            X <- playbase::normalizeMultiOmics(X)
          } else if (upload_datatype() == "methylomics") {
            write.csv(X, "~/Desktop/XX.csv")
            nX <- try(playbase::normalizeMethylation(X, m), silent = TRUE)
            if (is.null(nX)) dbg("--------------------nX is NULL")
            if (!is.null(nX)) X=nX; rm(nX)
          } else {
            dbg("[normalization_server:normalizedX] normalizing data using", m)
            X <- playbase::normalizeExpression(X, method = m, ref = ref, prior = prior)
          }
        } else {
          dbg("[normalization_server:normalizedX] Skipping normalization")
        }

        return(X)

      })

      ## Remove outliers
      cleanX <- reactive({
        shiny::req(dim(normalizedX()), dim(imputedX()$counts))
        X <- normalizedX()
        counts <- imputedX()$counts
        is.mox <- playbase::is.multiomics(rownames(counts))
        if (input$remove_outliers) {
          threshold <- input$outlier_threshold
          dbg("[normalization_server:cleanX] Removing outliers: Threshold = ", threshold)
          if (sum(is.na(X)) > 0) {
            if (is.mox) {
              X <- playbase::imputeMissing.mox(X, method = "SVD2")
            } else {
              X <- playbase::imputeMissing(X, method = "SVD2")
            }
          }
          res <- playbase::detectOutlierSamples(X, plot = FALSE)
          is.outlier <- (res$z.outlier > threshold)
          if (any(is.outlier) && !all(is.outlier)) {
            X <- X[, which(!is.outlier), drop = FALSE]
            counts <- counts[, colnames(X), drop = FALSE]
          }
        }
        return(list(counts = counts, X = X))
      })

      correctedX <- shiny::reactive({
        shiny::req(dim(cleanX()$X))
        X <- cleanX()$X
        cx <- list(X = X)
        return(cx)
      })

      annot <- shiny::reactive({
        annot <- imputedX()$annot
        return(annot)
      })

      ## ------------------------------------------------------------------
      ## Compute reactive
      ## ------------------------------------------------------------------
      results_correction_methods <- reactive({
        shiny::req(dim(cleanX()$X), dim(r_contrasts()), dim(r_samples()))
        X0 <- imputedX()$X
        X1 <- cleanX()$X ## normalized+cleaned
        samples <- r_samples()
        contrasts <- r_contrasts()
        batch.pars <- input$bec_param

        ## Average (if any dups) for BC overview
        dups <- sum(duplicated(rownames(X0)))
        if (dups > 0) X0 <- playbase::counts.mergeDuplicateFeatures(X0, is.counts = FALSE)
        dups <- sum(duplicated(rownames(X1)))
        if (dups > 0) X1 <- playbase::counts.mergeDuplicateFeatures(X1, is.counts = FALSE)

        is.mox <- playbase::is.multiomics(rownames(X0))

        if (sum(is.na(X0)) > 0) {
          if (is.mox) {
            X0 <- playbase::imputeMissing.mox(X0, method = "SVD2")
          } else {
            X0 <- playbase::imputeMissing(X0, method = "SVD2")
          }
        }

        if (sum(is.na(X1)) > 0) {
          if (is.mox) {
            X1 <- playbase::imputeMissing.mox(X1, method = "SVD2")
          } else {
            X1 <- playbase::imputeMissing(X1, method = "SVD2")
          }
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
        if (ncol(X0) > 100) methods <- methods[methods != "NPM"]
        shiny::updateSelectInput(
          session,
          "bec_method",
          choices = methods
        )
        xlist.init <- list("uncorrected" = X0, "normalized" = X1)

        shiny::withProgress(
          message = "Comparing batch-correction methods...",
          value = 0.3,
          {
            res <- playbase::compare_batchcorrection_methods(
              X1,
              samples,
              pheno = NULL,
              contrasts = contrasts,
              batch.pars = batch.pars,
              clust.method = "pca",
              methods = methods,
              evaluate = FALSE, ## no score computation
              xlist.init = xlist.init,
              ntop = 1000
            )
          }
        )

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

          is.mox <- playbase::is.multiomics(rownames(X))

          if (sum(is.na(X)) > 0) {
            if (is.mox) {
              X <- playbase::imputeMissing.mox(X, method = "SVD2")
            } else {
              X <- playbase::imputeMissing(X, method = "SVD2")
            }
          }

          out <- playbase::detectOutlierSamples(X, plot = FALSE)

          scaledX <- playbase::double_center_scale_fast(X)
          corX <- HiClimR::fastCor(t(scaledX), optBLAS = TRUE)

          ## standard dim reduction methods
          pos <- list()
          set.seed(1234)
          pca <- irlba::irlba(scaledX, nu = 2, nv = 0)
          pos[["pca"]] <- pca$u
          for (i in 1:length(pos)) {
            rownames(pos[[i]]) <- rownames(scaledX)
            colnames(pos[[i]]) <- paste0(names(pos)[i], "_", 1:2)
          }
          pos[["pca.varexp"]] <- (pca$d^2 / sum(pca$d^2)) * 100
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
        X0 <- imputedX()$X
        X1 <- cleanX()$X
        main.tt <- ifelse(input$normalize, norm_method(), "no normalization")

        if (input$norm_plottype == "boxplot") {
          if (ncol(X1) > 40) {
            jj <- sample(ncol(X1), 40)
            ii <- rownames(X1)
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
        # X0 <- r_counts()
        X0 <- cleanX()$counts
        X1 <- imputedX()$X

        dups <- sum(duplicated(rownames(X0)))
        if (dups > 0) X0 <- playbase::counts.mergeDuplicateFeatures(X0, is.counts = TRUE)
        dups <- sum(duplicated(rownames(X1)))
        if (dups > 0) X1 <- playbase::counts.mergeDuplicateFeatures(X1, is.counts = FALSE)

        X0 <- X0[rownames(X1), , drop = FALSE]

        has.zeros <- any(X0 == 0, na.rm = TRUE)
        ## if (!any(is.na(X0)) && !(input$zero_as_na && has.zeros)) {
        if (!any(is.na(X0)) && !(zero_as_na() && has.zeros)) {
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
          ## if (isolate(input$zero_as_na)) {
          if (isolate(zero_as_na())) {
            ii <- which(is.na(X0) | X0 == 0)
          }
          q999 <- quantile(X1, probs = 0.999, na.rm = TRUE)[1]
          X1[X1 > q999] <- NA
          h <- hist(X1, breaks = 80, plot = FALSE, las = 1)
          hh <- h$breaks

          ## set zero value to 1, NA values to 2
          X2 <- 1 * is.na(X0)
          ## if (input$zero_as_na) X2[X0 == 0] <- 1
          if (zero_as_na()) X2[X0 == 0] <- 1
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
              par(mar = c(3, 3, 2, 2), mgp = c(2.5, 0.85, 0))
              playbase::gx.imagemap(X2, cex = -1, col = rev(heat.colors(64)))
              title("missing values patterns", cex.main = 1.2)
            } else {
              plot.new()
              text(0.5, 0.5, "no missing values")
            }
          }

          if (input$missing_plottype == "ratio plot") {
            if (any(X2 > 0)) {
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

          if (input$missing_plottype == "missingness per sample") {
            if (any(X2 > 0)) {
              par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), mgp = c(2.5, 0.75, 0))
              X3 <- imputedX()$X
              pct.na <- colMeans(is.na(X3)) * 100
              bp <- barplot(pct.na,
                col = "grey", xaxt = "n",
                ylab = "Missing %", ylim = c(0, max(pct.na) + 10),
                cex.lab = 1.5, las = 2
              )
              text(
                x = bp, y = par("usr")[3] - 0.02 * diff(par("usr")[3:4]),
                labels = names(pct.na), srt = 45, adj = 1, xpd = TRUE, cex = 1
              )
              title("missingness per sample")
              grid()
              rm(X3)
            } else {
              plot.new()
              text(0.5, 0.5, "no missing values")
            }
          }

          if (input$missing_plottype == "missingness across features") {
            if (any(X2 > 0)) {
              par(mfrow = c(1, 1), mar = c(2, 3.5, 2, 2), mgp = c(2.5, 0.75, 0))
              X3 <- imputedX()$X
              pct.na <- round(rowMeans(is.na(X3)) * 100)
              hh <- hist(pct.na,
                xlim = c(0, 100), col = "grey", main = "",
                las = 1, tcl = -0.1, mgp = c(2.5, 0.5, 0), yaxs = "i",
                xlab = "Missingness across features (%)", ylab = "Number of features"
              )
              abline(v = mean(pct.na), col = "red")
              abline(v = median(pct.na), col = "blue")
              xpos <- 90
              ypos <- max(hh$counts) * 0.95
              lab1 <- paste0("Mean: ", round(mean(pct.na)), "%")
              lab2 <- paste0("Median: ", round(median(pct.na)), "%")
              text(xpos, ypos, labels = lab1, col = "red")
              text(xpos, ypos - (ypos * 8 / 100), labels = lab2, col = "blue")
              title("Distribution of missing values across features")
              grid()
              rm(X3)
            } else {
              plot.new()
              text(0.5, 0.5, "no missing values")
            }
          }

          if (input$missing_plottype == "PCA of imputed data") {
            if (any(X2 > 0)) {
              X3 <- imputedX()$X
              if (input$impute) X3 <- log2(imputedX()$counts + imputedX()$prior)
              mm <- c("SVD2", "QRILC", "MinProb", "Perseus")
              imp <- list()
              is.mox <- playbase::is.multiomics(rownames(X3))
              for (i in 1:length(mm)) {
                if (is.mox) {
                  imp[[mm[i]]] <- playbase::imputeMissing.mox(X3, mm[i])
                } else {
                  imp[[mm[i]]] <- playbase::imputeMissing(X3, mm[i])
                }
              }
              scaled.imp <- lapply(imp, function(x) playbase::double_center_scale_fast(x))
              par(mfrow = c(2, 2), mar = c(4, 3, 2, 0.5), las = 1, mgp = c(2, 0.4, 0), tcl = -0.1)
              cex1 <- cut(ncol(X3),
                breaks = c(0, 40, 100, 250, 1000, 999999),
                c(1, 0.85, 0.7, 0.55, 0.4)
              )
              cex1 <- 2.7 * as.numeric(as.character(cex1))
              for (i in 1:length(scaled.imp)) {
                set.seed(1234)
                pca <- irlba::irlba(scaled.imp[[i]], nu = 2, nv = 0)
                pca.pos <- pca$u
                pca.var <- (pca$d^2 / sum(pca$d^2)) * 100
                plot(pca.pos[, 1], pca.pos[, 2],
                  col = "black", pch = 20,
                  cex = cex1, cex.lab = 1, main = names(scaled.imp)[i],
                  xlab = paste0("PC1 (", round(pca.var[1], 2), "%)"),
                  ylab = paste0("PC2 (", round(pca.var[2], 2), "%)"),
                  asp = 1
                )
                grid()
                rm(pca, pca.pos, pca.var)
                gc()
              }
              rm(X3, imp, scaled.imp)
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
        shiny::validate(shiny::need(nrow(r_samples()) > 2, "Outlier detection requires at least 3 samples."))
        res <- results_outlier_methods()
        z0 <- as.numeric(input$outlier_threshold)
        zscore <- res$z.outlier
        Z <- res$Z
        pos <- res$pos[["pca"]]
        plottype <- "pca"
        if (plottype == "pca") {
          par(mfrow = c(1, 2), mar = c(3.2, 3, 2, 0.5), mgp = c(2.1, 0.8, 0))
          Z[which(is.infinite(Z) | is.nan(Z))] <- NA
          barplot(zscore,
            main = "outlier score", ylab = "z-score",
            las = 1, ylim = c(0, max(7, 1.2 * max(Z, na.rm = TRUE))),
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
        shiny::validate(shiny::need(nrow(r_samples()) > 2, "Batch-effect correction requires at least 3 samples."))
        if (input$batchcorrect) {
          plot_before_after()
        } else {
          plot_all_methods()
        }
      }

      plot_all_methods <- function() {
        out.res <- results_outlier_methods()
        res <- results_correction_methods()
        shiny::req(res)
        shiny::req(out.res)
        samples <- r_samples()

        methods <- c("uncorrected", sort(c("ComBat", "limma", "RUV", "SVA", "NPM")))

        pos.list <- res$pos
        pos0 <- out.res$pos[["pca"]]

        pos.list <- c(list("uncorrected" = pos0), pos.list)

        colorby_var <- input$colorby_var
        colorby_var <- intersect(colorby_var, colnames(samples))
        col1 <- factor(samples[, colorby_var])

        pheno <- res$pheno
        xdim <- length(pheno)
        breaks <- c(0, 40, 100, 250, 1000, 999999)
        labs <- c(1, 0.85, 0.7, 0.55, 0.4)
        cex1 <- cut(xdim, breaks, labs)
        cex1 <- 2.5 * as.numeric(as.character(cex1))
        if (is.na(cex1)) cex1 <- 1
        
        cols <- NULL
        ncol <- length(col1)
        col1a <- as.character(unname(col1))
        c1 <- all(!is.na(as.numeric(col1a)))
        c2 <- all(grepl("[0-9]", col1a))
        is.num <- (c1 & c2)
        if (is.num) {
          col1 <- as.numeric(col1a[!is.na(col1a)])
          pal <- colorRampPalette(c("gray", "black"))(ncol)
          cols <- pal[cut(col1, breaks = ncol, include.lowest = TRUE)]
        }
        color <- if (all(!is.null(cols))) cols else col1

        if (is.num) {
          par(mfrow = c(2, 4), mar = c(3, 3, 2, 1), mgp = c(2, 0.4, 0), tcl = -0.1)
        } else {
          par(mfrow = c(2, 3), mar = c(3, 3, 2, 1), mgp = c(2, 0.4, 0), tcl = -0.1)
          cex.axis <- 1
          cex.lab <- 1
          cex.main <- 1.2
        }
        
        for (m in methods) {
          if (m %in% names(pos.list)) {
            plot(pos.list[[m]],
              col = color, cex = cex1,
              pch = 20, las = 1, xlab = "PCA_1", ylab = "PCA_2"
            )
          } else {
            plot.new()
            text(0.45, 0.5, "method failed")
          }
          title(m, cex.main = 1.5)
        }

        if (is.num) {
          plot.new()
          fields::image.plot(
            legend.only = TRUE, col = pal, zlim = range(col1),
            legend.width = 4, axis.args = list(cex.axis = 1.3, las = 1)
          )
        }
      }

      plot_before_after <- function() {
        out.res <- results_outlier_methods()
        res <- results_correction_methods()
        samples <- r_samples()

        pos0 <- out.res$pos[["pca"]]
        pos0.varexp <- out.res$pos[["pca.varexp"]]
        pos1.varexp <- res[["pca.varexp"]]
        method <- input$bec_method

        if (!method %in% names(res$pos)) {
          plot.new()
          text(0.45, 0.5, "method failed")
          return(NULL)
        }

        pos1 <- res$pos[[method]]
        if (!input$batchcorrect) pos1 <- pos0

        kk <- intersect(rownames(pos0), rownames(pos1))
        pos0 <- pos0[kk, , drop = FALSE]
        pos1 <- pos1[kk, , drop = FALSE]

        pheno <- playbase::contrasts2pheno(r_contrasts(), r_samples())
        pheno <- pheno[rownames(pos0)]
        colorby_var <- input$colorby_var
        colorby_var <- intersect(colorby_var, colnames(samples))
        samples <- samples[rownames(pos0), , drop = FALSE]
        col1 <- factor(samples[, colorby_var])
        breaks <- c(0, 40, 100, 250, 1000, 999999)
        labs <- c(1, 0.85, 0.7, 0.55, 0.4)
        cex1 <- cut(nrow(pos1), breaks, labs)
        cex1 <- 2.7 * as.numeric(as.character(cex1))

        cols <- NULL
        ncol <- length(col1)
        col1a <- as.character(unname(col1))
        c1 <- all(!is.na(as.numeric(col1a)))
        c2 <- all(grepl("[0-9]", col1a))
        is.num <- (c1 & c2)
        if (is.num) {
          col1 <- as.numeric(col1a[!is.na(col1a)])
          pal <- colorRampPalette(c("gray", "black"))(ncol)
          cols <- pal[cut(col1, breaks = ncol, include.lowest = TRUE)]
        }
        color <- if (all(!is.null(cols))) cols else col1

        if (is.num) {
          layout(matrix(c(1, 2, 3), ncol = 3), widths = c(4, 4, 2))
          cex.axis <- 1.4
          cex.lab <- 1.2
          cex.main <- 1.7
        } else {
          par(mfrow = c(1, 2), mar = c(3.2, 3, 2, 0.5), mgp = c(2.1, 0.4, 0), tcl = -0.1)
          cex.axis <- 1
          cex.lab <- 1
          cex.main <- 1.2
        }

        if (is.num) {
          par(mar = c(3.4, 3.5, 2, 0.1), mgp = c(2.3, 0.4, 0), tcl = -0.1)
        }
        plot(pos0,
          col = color, pch = 20, cex = cex1, las = 1,
          cex.axis = cex.axis, cex.lab = cex.lab,
          main = "uncorrected", cex.main = cex.main,
          xlab = paste0("PC1 (", round(pos0.varexp[1], 2), "%)"),
          ylab = paste0("PC2 (", round(pos0.varexp[2], 2), "%)")
        )

        if (is.num) {
          par(mar = c(3.4, 4.5, 2, 0.1), mgp = c(2.4, 0.4, 0), tcl = -0.1)
        }
        plot(pos1,
          col = color, pch = 20, cex = cex1, las = 1,
          cex.axis = cex.axis, cex.lab = cex.lab,
          main = method, cex.main = cex.main,
          xlab = paste0("PC1 (", round(pos1.varexp[[method]][1], 2), "%)"),
          ylab = paste0("PC2 (", round(pos1.varexp[[method]][2], 2), "%)")
        )

        if (is.num) {
          plot.new()
          par(mar = c(3, 3, 3, 4.8))
          fields::image.plot(
            legend.only = TRUE, col = pal, zlim = range(col1),
            legend.width = 4, axis.args = list(cex.axis = 1.3, las = 1)
          )
        }
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
          if (nrow(samples) < 3) {
            return(NULL)
          }
          pars <- playbase::get_model_parameters(X, samples, pheno = NULL, contrasts = contrasts)
          all.pars <- setdiff(colnames(samples), pars$pheno.pars)
          all.pars <- union(all.pars, pars$batch.pars)
          names(all.pars) <- ifelse(all.pars %in% pars$batch.pars,
            paste(all.pars, "*"), all.pars
          )
          all.pars <- c("<autodetect>", all.pars)
          return(all.pars)
        }
      )

      getMetadataVars <- eventReactive(
        {
          list(r_samples())
        },
        {
          shiny::req(dim(r_samples()))
          samples <- r_samples()
          return(colnames(samples))
        }
      )

      output$normalization <- shiny::renderUI({
        batch_params <- getBatchParams()
        metadata_vars <- getMetadataVars()

        ## -----------------------------------------------------------------
        ## Get default values from recompute_pgx if available
        ## -----------------------------------------------------------------
        pgx <- recompute_pgx()
        pgx_settings <- if (!is.null(pgx)) pgx$settings else NULL

        ## Imputation defaults
        default_zero_as_na <- FALSE
        default_impute <- DEFAULTS$qc$impute
        default_impute_method <- "SVD2"
        if (!is.null(pgx_settings$imputation_method) && is.list(pgx_settings$imputation_method)) {
          imp <- pgx_settings$imputation_method
          if (!is.null(imp$zero_as_na)) default_zero_as_na <- imp$zero_as_na
          if (!is.null(imp$imputation)) {
            default_impute <- (imp$imputation != "no_imputation")
            if (default_impute) default_impute_method <- imp$imputation
          }
        }

        ## Normalization defaults
        default_normalize <- TRUE
        default_norm_method <- 1
        if (!is.null(pgx_settings$norm_method)) {
          if (pgx_settings$norm_method == "skip_normalization") {
            default_normalize <- FALSE
          } else {
            default_normalize <- TRUE
            default_norm_method <- pgx_settings$norm_method
          }
        }

        ## Outlier removal defaults
        default_remove_outliers <- FALSE
        default_outlier_threshold <- 6
        if (!is.null(pgx_settings$remove_outliers)) {
          ro <- pgx_settings$remove_outliers
          if (is.character(ro) && ro == "no_outlier_removal") {
            default_remove_outliers <- FALSE
          } else if (is.numeric(ro)) {
            default_remove_outliers <- TRUE
            default_outlier_threshold <- ro
          }
        }

        ## Batch correction defaults
        default_batchcorrect <- FALSE
        default_bec_method <- "SVA"
        default_bec_param <- batch_params[1]
        if (!is.null(pgx_settings$bc_method)) {
          bc <- pgx_settings$bc_method
          if (is.character(bc) && bc == "no_batch_correct") {
            default_batchcorrect <- FALSE
          } else if (is.list(bc)) {
            default_batchcorrect <- TRUE
            if (!is.null(bc$method)) default_bec_method <- bc$method
            if (!is.null(bc$param)) default_bec_param <- bc$param
          }
        }

        score.infotext <-
          "Outliers markedly deviate from the vast majority of samples. Outliers could be caused by technical factors and negatively affect data analysis. Here, outliers are identified and marked for removal should you wish so."

        missing.infotext <-
          "Missing values (MVs) reduce the completeness of biological data and hinder preprocessing steps. MVs (i.e., NA), more often populate proteomics and metabolomics data. Here, MVs are identified and their patterns in your data is shown. PCA is also optionally performed on data imputed with all methods to aid comparison."

        normalization.infotext <-
          "Normalization enables to standardize the data and improve their consistency, comparability and reproducibility. Boxplots of raw (unnormalized) and normalized data are shown. Normalization method can be selected on the left, under “Normalization”."

        batcheff.infotext <-
          "Batch effects (BEs) are due to technical, experimental factors that introduce unwanted variation into the measurements. Here, BEs are detected and BEs correction is shown. BE correction methods can be selected on the left, under “Batch-effects correction”."

        missing.options <- tagList(
          shiny::radioButtons(ns("missing_plottype"), "Plot type:",
            c(
              "heatmap", "ratio plot", "missingness per sample",
              "missingness across features", "PCA of imputed data"
            ),
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
          shiny::radioButtons(
            ns("colorby_var"),
            label = "Annotate by:",
            choices = metadata_vars,
            selected = metadata_vars[1],
            inline = FALSE
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
                shiny::div(
                  style = "display: flex; align-items: center; justify-content: space-between;",
                  shiny::p("Handle missing values:\n"),
                  shiny::HTML("<a href='https://bigomics.ch/blog/imputation-of-missing-values-in-proteomics' target='_blank' class='info-link' style='margin-left: 15px;'>
                      <i class='fa-solid fa-circle-info info-icon' style='color: blue; font-size: 20px;'></i>
                      </a>")
                ),
                if (upload_datatype() != "methylomics") {
                  shiny::checkboxInput(ns("zero_as_na"), label = "Treat zero as NA", value = default_zero_as_na)
                },
                shiny::checkboxInput(ns("filtermissing"), label = "Remove NA rows", value = FALSE),
                shiny::conditionalPanel("input.filtermissing == true",
                  ns = ns,
                  shiny::selectInput(ns("filterthreshold"), NULL,
                    choices = c(
                      ">10% NA" = 0.1, ">20% NA" = 0.2, ">50% NA" = 0.5,
                      "<=3 valid in any group" = 3, "<=50% valid in any group" = -0.5
                    ),
                    selected = 0.2
                  )
                ),
                shiny::checkboxInput(ns("impute"), label = "Impute NA", value = default_impute),
                shiny::conditionalPanel("input.impute == true",
                  ns = ns,
                  shiny::selectInput(ns("impute_method"), NULL,
                    choices = c("SVDimpute" = "SVD2", "QRILC", "MinProb", "Perseus-like" = "Perseus"),
                    selected = default_impute_method
                  )
                ),
                br()
              ),
              bslib::accordion_panel(
                title = "2. Normalization",
                shiny::div(
                  style = "display: flex; align-items: center; justify-content: space-between;",
                  shiny::p("Normalize the data using one of the following methods:"),
                  shiny::HTML("<a href='https://omicsplayground.readthedocs.io/en/latest/methods/#normalization' target='_blank' class='info-link' style='margin-left: 15px;'>
                      <i class='fa-solid fa-circle-info info-icon' style='color: blue; font-size: 20px;'></i>
                      </a>")
                ),
                shiny::checkboxInput(ns("normalize"), label = "Normalize data", value = default_normalize),
                shiny::conditionalPanel(
                  "input.normalize == true",
                  ns = ns,
                  shiny::selectInput(
                    ns("normalization_method"), NULL,
                    choices = if (grepl("proteomics|metabolomics", upload_datatype(),
                      ignore.case = TRUE
                    )) {
                      c("maxMedian", "maxSum", "quantile", "reference")
                    } else if (grepl("methylomics", upload_datatype(),
                      ignore.case = TRUE
                    )) {
                      c(
                        "BMIQ", "quantile"
                      )
                    } else if (grepl("multi-omics", upload_datatype(),
                      ignore.case = TRUE
                    )) {
                      c(
                        "multi-omics median" = "median"
                        ## "multi-omics combat" = "combat"
                      )
                    } else {
                      c(
                        "CPM", "CPM+quantile", "TMM", "quantile",
                        "maxMedian", "maxSum", "reference"
                      )
                    },
                    selected = default_norm_method
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
                ),
                br()
              ),
              bslib::accordion_panel(
                title = "3. Remove outliers",
                shiny::p("Detect and remove outlier samples."),
                shiny::checkboxInput(ns("remove_outliers"), "remove outliers", value = default_remove_outliers),
                shiny::conditionalPanel("input.remove_outliers == true",
                  ns = ns,
                  shiny::sliderInput(ns("outlier_threshold"), "Select threshold:", 1, 12, default_outlier_threshold, 1)
                ),
                br()
              ),
              bslib::accordion_panel(
                title = "4. Batch-effect correction",
                shiny::div(
                  style = "display: flex; align-items: center; justify-content: space-between;",
                  shiny::p("Remove unwanted variation from your data."),
                  shiny::HTML("<a href='https://omicsplayground.readthedocs.io/en/latest/methods/#batch-correction' target='_blank' class='info-link' style='margin-left: 15px;'>
                      <i class='fa-solid fa-circle-info info-icon' style='color: blue; font-size: 20px;'></i>
                      </a>")
                ),
                shiny::checkboxInput(ns("batchcorrect"),
                  label = "Remove batch effects",
                  value = default_batchcorrect
                ),
                shiny::conditionalPanel(
                  "input.batchcorrect == true",
                  ns = ns,
                  shiny::selectInput(
                    ns("bec_method"),
                    label = "Select method:",
                    choices = c("ComBat", "limma", "NPM" = "NPM", "RUV" = "RUV", "SVA" = "SVA"),
                    selected = default_bec_method
                  ),
                  shiny::conditionalPanel(
                    "input.bec_method == 'ComBat' || input.bec_method == 'limma'",
                    ns = ns,
                    shiny::selectizeInput(
                      ns("bec_param"),
                      label = "Batch parameter:",
                      choices = batch_params, ## reactive
                      selected = default_bec_param,
                      multiple = TRUE,
                      options = list(placeholder = "Select...")
                    ),
                    shiny::br()
                  )
                ),
                br(),
                shiny::HTML("<div style='margin-top: 10px;'><a href='https://academic.oup.com/bioinformatics/article/41/3/btaf084/8042340' target='_blank' style='color: #0066cc; text-decoration: none;'>Learn about NPM <i class='fa-solid fa-external-link' style='font-size: 12px;'></i></a></div>")
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
                info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#normalization",
                show.maximize = FALSE
              ),
              PlotModuleUI(
                ns("plot3"),
                title = "Outliers detection",
                info.text = score.infotext,
                caption = score.infotext,
                options = outlier.options,
                height = c("auto", "100%"),
                info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#identification-of-outlier-samples",
                show.maximize = FALSE
              ),
              PlotModuleUI(
                ns("plot4"),
                title = "Batch-effects correction",
                options = bec.options,
                info.text = batcheff.infotext,
                height = c("auto", "100%"),
                info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#batch-correction",
                show.maximize = FALSE
              )
            )
          ),
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

      counts <- reactive({
        shiny::req(dim(cleanX()$counts))
        counts <- cleanX()$counts
        return(counts)
      })

      cX <- reactive({
        shiny::req(dim(correctedX()$X))
        cX <- correctedX()$X
        return(cX)
      })

      imputation_method <- reactive({
        ## ll <- list(zero_as_na = input$zero_as_na, imputation = input$impute_method)
        ## if (!input$impute) {
        ##   ll <- list(zero_as_na = input$zero_as_na, imputation = "no_imputation")
        ## }
        ll <- list(zero_as_na = zero_as_na(), imputation = input$impute_method)
        if (!input$impute) {
          ll <- list(zero_as_na = zero_as_na(), imputation = "no_imputation")
        }
        return(ll)
      })

      norm_method <- reactive({
        m <- input$normalization_method
        if (!input$normalize) m <- "skip_normalization"
        return(m)
      })

      remove_outliers <- reactive({
        ro <- input$outlier_threshold
        if (input$remove_outliers == FALSE) ro <- "no_outlier_removal"
        return(ro)
      })

      bc_method <- reactive({
        param <- input$bec_param
        if ("<autodetect>" %in% param && length(param) > 1) {
          param <- setdiff(param, "<autodetect>")
          shiny::updateSelectizeInput(session, "bec_param", selected = param)
        }
        ll <- list(method = input$bec_method, param = param)
        if (input$batchcorrect == FALSE) ll <- "no_batch_correct"
        return(ll)
      })

      return(
        list(
          counts = counts,
          X = cX,
          norm_method = norm_method,
          imputation_method = imputation_method,
          bc_method = bc_method,
          remove_outliers = remove_outliers,
          annot = annot
        )
      ) ## pointing to reactive
    } ## end-of-server
  )
}
