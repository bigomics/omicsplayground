# Bulk-upload preprocessing controls and visual previews.
#
# This file owns UI policy and plotting for the canonical matrix pipeline.
# Numerical stages run only through playbase.preprocess functions.

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
  is.nulisa,
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

      ## Translates a norm_method a legacy object recorded into a value the
      ## dropdown still offers. "median" named a multi-omics algorithm deleted
      ## from playbase on 2025-11-11 (f3d5b3d8); without the translation it
      ## reaches selected= as a value not in choices, and selectInput() then
      ## silently takes choices[[1]] instead of saying so.
      .normalize_selected <- function(x) {
        if (identical(x, "median")) "multiomics" else x
      }

      observeEvent(input$normalization_method, {
        shiny::req(input$normalization_method == "reference")
        gg <- sort(rownames(r_counts()))
        pgx <- recompute_pgx()
        selected <- if (!is.null(pgx)) {
          pgx$settings$preprocess$options$normalize_args$ref
        } else {
          NULL
        }
        if (is.null(selected)) selected <- character(0)
        shiny::updateSelectizeInput(session, "ref_gene",
          choices = gg,
          selected = selected, server = TRUE
        )
      })

      ## Canonical staged previews
      imputedX <- reactive({
        shiny::req(dim(r_counts()), !is.null(input$normalize))
        shiny::req(dim(r_contrasts()))
        .opg_run_preprocess_preview(
          counts = r_counts(),
          samples = r_samples(),
          contrasts = r_contrasts(),
          annot = r_annot(),
          options = preprocess(),
          through = "impute"
        )
      })

      normalizedX <- reactive({
        shiny::req(dim(imputedX()$X))
        if (identical(input$normalization_method, "reference")) {
          shiny::validate(shiny::need(
            isTruthy(input$ref_gene),
            tspan("Please select reference gene", js = FALSE)
          ))
        }
        .opg_run_preprocess_preview(
          counts = r_counts(),
          samples = r_samples(),
          contrasts = r_contrasts(),
          annot = r_annot(),
          options = preprocess(),
          through = "normalize"
        )$X
      })

      cleanX <- reactive({
        shiny::req(dim(normalizedX()))
        .opg_run_preprocess_preview(
          counts = r_counts(),
          samples = r_samples(),
          contrasts = r_contrasts(),
          annot = r_annot(),
          options = preprocess(),
          through = "outliers"
        )
      })

      correctedX <- shiny::reactive({
        shiny::req(dim(cleanX()$X))
        .opg_run_preprocess_preview(
          counts = r_counts(),
          samples = r_samples(),
          contrasts = r_contrasts(),
          annot = r_annot(),
          options = preprocess(),
          through = "batch"
        )
      })

      annot <- shiny::reactive({
        r_annot()
      })

      ## ------------------------------------------------------------------
      ## Compute reactive
      ## ------------------------------------------------------------------
      results_correction_methods <- reactive({
        shiny::req(dim(cleanX()$X), dim(r_contrasts()), dim(r_samples()))
        X0 <- imputedX()$X
        X1 <- cleanX()$X
        samples <- r_samples()
        contrasts <- r_contrasts()
        batch.pars <- input$bec_param

        ## Average (if any dups) for BC overview
        dups <- sum(duplicated(rownames(X0)))
        if (dups > 0) {
          X0 <- playbase.preprocess::pp.deduplicate(
            X0,
            method = "average",
            space = "log2"
          )$X
        }
        dups <- sum(duplicated(rownames(X1)))
        if (dups > 0) {
          X1 <- playbase.preprocess::pp.deduplicate(
            X1,
            method = "average",
            space = "log2"
          )$X
        }

        if (sum(is.na(X0)) > 0) {
          X0 <- .opg_impute(X0, method = "SVD2")
        }

        if (sum(is.na(X1)) > 0) {
          X1 <- .opg_impute(X1, method = "SVD2")
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

        ## Top-1000 most variable features by default, or all if toggled.
        ntop_features <- if (isTRUE(input$bec_full_features)) Inf else 1000

        methods <- c("ComBat", "limma", "RUV", "SVA", "NPM")
        ## NPM does not scale: drop it for many samples or for the very large
        ## feature space of methylation arrays.
        if (ncol(X0) > 100 || upload_datatype() == "methylomics") {
          methods <- methods[methods != "NPM"]
        }
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
              ntop = ntop_features,
              npc = 5
            )
          }
        )

        return(res)
      })

      ## Trim PC selector choices to the number of PCs actually computed.
      shiny::observeEvent(results_correction_methods(), {
        res <- results_correction_methods()
        shiny::req(res, res$pos, length(res$pos) > 0)
        npc_avail <- ncol(res$pos[[1]])
        out.res <- results_outlier_methods()
        if (!is.null(out.res) && !is.null(out.res$pos[["pca"]])) {
          npc_avail <- min(npc_avail, ncol(out.res$pos[["pca"]]))
        }
        choices <- paste0("PC", seq_len(max(2, npc_avail)))
        sel.x <- if (input$bec_xpc %in% choices) input$bec_xpc else choices[1]
        sel.y <- if (input$bec_ypc %in% choices) input$bec_ypc else choices[min(2, length(choices))]
        shiny::updateSelectInput(session, "bec_xpc", choices = choices, selected = sel.x)
        shiny::updateSelectInput(session, "bec_ypc", choices = choices, selected = sel.y)
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

          outlier_result <- playbase.preprocess::pp.removeOutliers(
            X,
            threshold = Inf,
            methods = c("z.correlation", "z.distance", "z.features")
          )
          X <- outlier_result$X
          out <- outlier_result$scores

          scaledX <- playbase::double_center_scale_fast(X)
          corX <- HiClimR::fastCor(t(scaledX), optBLAS = TRUE)

          ## standard dim reduction methods
          pos <- list()
          set.seed(1234)
          npc <- max(2, min(5, min(dim(scaledX)) - 1))
          pca <- irlba::irlba(scaledX, nu = npc, nv = 0)
          pos[["pca"]] <- pca$u[, seq_len(npc), drop = FALSE]
          for (i in 1:length(pos)) {
            rownames(pos[[i]]) <- rownames(scaledX)
            colnames(pos[[i]]) <- paste0(names(pos)[i], "_", seq_len(ncol(pos[[i]])))
          }
          ## total-variance denominator (see compare_batchcorrection_methods):
          ## keeps uncorrected %varexp comparable to the corrected methods'
          pos[["pca.varexp"]] <- (pca$d^2 / sum(scaledX^2)) * 100
          out$pos <- pos
          out$corX <- corX
          out
        }
      )

      ## ------------------------------------------------------------------
      ## Plot functions
      ## ------------------------------------------------------------------

      plot_normalization <- function() {
        clean <- cleanX()
        rX <- .opg_preview_counts(clean)
        X0 <- .opg_align_preview_X(imputedX(), clean)
        X1 <- clean$X
        main.tt <- ifelse(input$normalize, norm_method(), "no normalization")

        if (input$norm_plottype == "boxplot") {
          if (ncol(X1) > 40) {
            jj <- withr::with_seed(1234, sample(seq_len(ncol(X1)), 40))
            ii <- seq_len(nrow(X1))
            if (length(ii) > 2000) {
              ii <- withr::with_seed(1235, sample(ii, 2000))
            }
            X0 <- X0[ii, jj, drop = FALSE]
            X1 <- X1[ii, jj, drop = FALSE]
            rX <- rX[ii, jj, drop = FALSE]
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

      plot_missingvalues <- function() {
        preview <- imputedX()
        X0 <- .opg_preview_counts(preview)
        X1 <- preview$X

        has.zeros <- any(X0 == 0, na.rm = TRUE)
        if (!any(is.na(X0)) && !(zero_as_na() && has.zeros)) {
          plot.new()
          text(0.5, 0.5, "No missing values", cex = 1.2)
        } else {
          ii <- which(is.na(X0))
          if (isolate(zero_as_na())) {
            ii <- which(is.na(X0) | X0 == 0)
          }
          q999 <- quantile(X1, probs = 0.999, na.rm = TRUE)[1]
          X1[X1 > q999] <- NA
          h <- hist(X1, breaks = 80, plot = FALSE, las = 1)
          hh <- h$breaks

          ## set zero value to 1, NA values to 2
          X2 <- 1 * is.na(X0)
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
              preview_options <- preprocess()
              preview_options$impute <- FALSE
              X3 <- .opg_run_preprocess_preview(
                counts = r_counts(),
                samples = r_samples(),
                contrasts = r_contrasts(),
                annot = r_annot(),
                options = preview_options,
                through = "impute"
              )$X
              mm <- c("SVD2", "QRILC", "MinProb", "Perseus")
              imp <- list()
              for (i in 1:length(mm)) {
                imp[[mm[i]]] <- .opg_impute(X3, mm[i])
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

      plot_outliers <- function() {
        shiny::validate(shiny::need(nrow(r_samples()) > 2, "Outlier detection requires at least 3 samples."))
        res <- results_outlier_methods()
        z0 <- as.numeric(input$outlier_threshold)
        zscore <- res$z.outlier
        Z <- res$Z
        pos <- res$pos[["pca"]][, 1:2, drop = FALSE]
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
      }

      plot_correction <- function() {
        shiny::validate(shiny::need(nrow(r_samples()) > 2, "Batch-effects correction requires at least 3 samples."))
        bec_view <- if (is.null(input$bec_view)) "pca" else input$bec_view
        switch(bec_view,
          loadings = plot_bec_biplot("loadings"),
          pheno = plot_bec_biplot("pheno"),
          scree = plot_bec_scree(),
          if (input$batchcorrect) plot_before_after() else plot_all_methods()
        )
      }

      ## Biplot grid for the BC panel: one biplot per method (sample scores as
      ## points, arrows overlaid) on the two selected PCs. Same method list /
      ## grid as plot_all_methods. Two arrow sources:
      ##  - "loadings": top feature loadings (res$loadings, same PCA as scores)
      ##  - "pheno": each annotation's correlation with the two PCs, precomputed
      ##    in playbase (res$pheno.cor; eigencorplot-style, plain cor())
      plot_bec_biplot <- function(arrows = "loadings") {
        res <- results_correction_methods()
        shiny::req(res, res$pos)
        if (arrows == "loadings") shiny::req(res$loadings)
        samples <- r_samples()

        xpc <- as.integer(sub("PC", "", if (is.null(input$bec_xpc)) "PC1" else input$bec_xpc))
        ypc <- as.integer(sub("PC", "", if (is.null(input$bec_ypc)) "PC2" else input$bec_ypc))
        shiny::req(!is.na(xpc), !is.na(ypc))

        methods <- c("uncorrected", sort(c("ComBat", "limma", "RUV", "SVA", "NPM")))
        methods <- intersect(methods, names(res$pos))
        colorby_var <- intersect(input$colorby_var, colnames(samples))

        ## arrow endpoints for a method: named 2-col matrix in loading- or
        ## correlation-units (scaled to the score cloud later). Both come
        ## precomputed from playbase::compare_batchcorrection_methods.
        arrow_ends <- function(m, scores) {
          if (arrows == "loadings") {
            L <- res$loadings[[m]]
            if (is.null(L) || max(xpc, ypc) > ncol(L)) {
              return(NULL)
            }
            A <- L[, c(xpc, ypc), drop = FALSE]
            A[head(order(-(A[, 1]^2 + A[, 2]^2)), 8), , drop = FALSE]
          } else {
            Rc <- res$pheno.cor[[m]]
            if (is.null(Rc) || max(xpc, ypc) > ncol(Rc)) {
              return(NULL)
            }
            A <- Rc[, c(xpc, ypc), drop = FALSE]
            A <- A[stats::complete.cases(A), , drop = FALSE]
            if (!nrow(A)) {
              return(NULL)
            }
            A[head(order(-(A[, 1]^2 + A[, 2]^2)), 12), , drop = FALSE]
          }
        }

        draw_biplot <- function(m) {
          P <- res$pos[[m]]
          if (is.null(P) || max(xpc, ypc) > ncol(P)) {
            plot.new()
            text(0.45, 0.5, "method failed")
            title(m, cex.main = 1.3)
            return(invisible())
          }
          scores <- P[, c(xpc, ypc), drop = FALSE]

          smp <- samples[rownames(scores), , drop = FALSE]
          color <- factor(smp[, colorby_var])
          cex1 <- cut(nrow(scores), c(0, 40, 100, 250, 1000, 999999), c(1, 0.85, 0.7, 0.55, 0.4))
          cex1 <- 2 * as.numeric(as.character(cex1))
          if (is.na(cex1)) cex1 <- 1

          lim <- c(-1, 1) * max(abs(scores)) * 1.3
          plot(scores,
            col = color, pch = 20, cex = cex1, las = 1, xlim = lim, ylim = lim,
            xlab = paste0("PC", xpc), ylab = paste0("PC", ypc)
          )
          title(m, cex.main = 1.3)
          abline(h = 0, v = 0, lty = 3, col = "grey70")

          A <- arrow_ends(m, scores)
          if (is.null(A) || !nrow(A) || max(abs(A)) == 0) {
            return(invisible())
          }
          ## scale arrows to fill the score cloud (loadings/correlations sit on a
          ## much smaller scale than the scores)
          s <- 0.85 * max(abs(scores)) / max(abs(A))
          ax <- A[, 1] * s
          ay <- A[, 2] * s
          arrows(0, 0, ax, ay, length = 0.05, col = "#B3444488", lwd = 1.3)
          ## repel labels off each other (base-graphics; places each away from
          ## its nearest neighbour) rather than a fixed left/right offset
          plotrix::thigmophobe.labels(ax, ay, rownames(A),
            cex = 0.65, col = "#7A2E2E", offset = 0.4, xpd = NA
          )
        }

        par(mfrow = c(2, 3), mar = c(3, 3, 2, 1), mgp = c(1.9, 0.4, 0), tcl = -0.2)
        for (m in methods) draw_biplot(m)
      }

      ## Scree grid for the BC panel: per method, % variance explained per PC
      ## (bars) with the cumulative % overlaid as a line. Both are in % units on
      ## the same y-axis. Uses res$pca.varexp (now normalised to total variance).
      plot_bec_scree <- function() {
        res <- results_correction_methods()
        shiny::req(res, res$pca.varexp)

        methods <- c("uncorrected", sort(c("ComBat", "limma", "RUV", "SVA", "NPM")))
        methods <- intersect(methods, names(res$pca.varexp))

        par(mfrow = c(2, 3), mar = c(3.2, 3.4, 2, 1), mgp = c(2, 0.5, 0), tcl = -0.2)
        for (m in methods) {
          ve <- res$pca.varexp[[m]]
          if (is.null(ve)) {
            plot.new()
            text(0.45, 0.5, "method failed")
            title(m, cex.main = 1.3)
            next
          }
          ve <- ve[seq_len(min(length(ve), 10))]
          cum <- cumsum(ve)
          bp <- barplot(ve,
            col = "#8FB3D9", border = NA, ylim = c(0, min(100, max(cum) * 1.1)),
            names.arg = paste0("PC", seq_along(ve)), las = 2, cex.names = 0.7,
            ylab = "% variance explained"
          )
          title(m, cex.main = 1.3)
          lines(bp, cum, type = "b", pch = 20, col = "#B34444", lwd = 1.5)
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

        xpc <- as.integer(sub("PC", "", if (is.null(input$bec_xpc)) "PC1" else input$bec_xpc))
        ypc <- as.integer(sub("PC", "", if (is.null(input$bec_ypc)) "PC2" else input$bec_ypc))
        shiny::req(!is.na(xpc), !is.na(ypc))
        pos.list <- lapply(pos.list, function(p) {
          if (is.null(p) || ncol(p) < max(xpc, ypc)) {
            return(NULL)
          }
          p[, c(xpc, ypc), drop = FALSE]
        })

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
          if (m %in% names(pos.list) && !is.null(pos.list[[m]])) {
            plot(pos.list[[m]],
              col = color, cex = cex1,
              pch = 20, las = 1,
              xlab = paste0("PC", xpc), ylab = paste0("PC", ypc)
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

        xpc <- as.integer(sub("PC", "", if (is.null(input$bec_xpc)) "PC1" else input$bec_xpc))
        ypc <- as.integer(sub("PC", "", if (is.null(input$bec_ypc)) "PC2" else input$bec_ypc))
        shiny::req(!is.na(xpc), !is.na(ypc))
        npc.avail <- min(ncol(pos0), ncol(pos1))
        if (max(xpc, ypc) > npc.avail) {
          xpc <- min(xpc, npc.avail)
          ypc <- min(ypc, npc.avail)
        }
        pos0 <- pos0[, c(xpc, ypc), drop = FALSE]
        pos1 <- pos1[, c(xpc, ypc), drop = FALSE]

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
          xlab = paste0("PC", xpc, " (", round(pos0.varexp[xpc], 2), "%)"),
          ylab = paste0("PC", ypc, " (", round(pos0.varexp[ypc], 2), "%)")
        )

        if (is.num) {
          par(mar = c(3.4, 4.5, 2, 0.1), mgp = c(2.4, 0.4, 0), tcl = -0.1)
        }

        plot(pos1,
          col = color, pch = 20, cex = cex1, las = 1,
          cex.axis = cex.axis, cex.lab = cex.lab,
          main = method, cex.main = cex.main,
          xlab = paste0("PC", xpc, " (", round(pos1.varexp[[method]][xpc], 2), "%)"),
          ylab = paste0("PC", ypc, " (", round(pos1.varexp[[method]][ypc], 2), "%)")
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

      plot_methyl <- function() {
        X <- normalizedX()
        if (input$methyl_plottype == "Density") {
          par(mfrow = c(1, 1), mar = c(3.3, 3.2, 0.8, 0.5), las = 1, mgp = c(2.1, 0.35, 0), tcl = -0.1)
          minfi::densityPlot(X, pal = "gray60", xlab = "Beta signal", main = "", cex.lab = 1.4, cex.axis = 1.3)
          grid()
        } else if (input$methyl_plottype == "Beanplot") {
          par(mfrow = c(1, 1), mar = c(4.5, 3.3, 0.8, 0.5), las = 2, tcl = -0.1, mgp = c(2.2, 0.5, 0))
          x <- reshape2::melt(X, varnames = c("cpg", "sample"))
          ww <- c(0, 1, 1, 0)
          beanplot::beanplot(value ~ sample,
            data = x, horizontal = FALSE, what = ww, log = "",
            ylim = c(0, 1), ylab = "Beta signal", method = "stack", main = "", beanlinewd = 1,
            cex.lab = 1.2, cex.axis = 0.8, border = "gray", frame.plot = FALSE
          )
          grid()
        }
        ## if (input$infer_sex) { ## placeholder
        ##   S <- playbase::infer_sex_methyl(data = X, meth_type = meth_type())
        ##   pred_sex <- S[["pred_sex"]]
        ##   x_med <- S[["x_med"]]
        ##   y_med <- S[["y_med"]]
        ##   shiny::validate(shiny::need(is.null(pred_sex), "No X or Y-linked probes found. Could not infer sex."))
        ##   par(mfrow = c(1, 1), mar = c(5, 5, 1.5, 0.5), las = 1, mgp = c(2.5, 0.5, 0), tcl = -0.1)
        ##   cols <- ifelse(pred_sex == "F", "red", "blue")
        ##   plot(1:length(y_med), y_med, col = cols, ylim = c(0, max(y_med) + mean(y_med)),
        ##     ylab = "Median Y-linked CpG beta value", pch = 15, xlab = "", xaxt = "n", cex = 1.5)
        ##   axis(side = 1, at = 1:length(y_med), labels = FALSE)
        ##   text(x = 1:length(y_med), y = par("usr")[3] - diff(par("usr")[3:4]) * 0.03,
        ##     labels = names(y_med), srt = 45, adj = 1, xpd = TRUE, cex = 0.5)
        ##   legend("topright", legend = c("F", "M"), fill = c("red", "blue"), cex = 1.2)
        ##   grid()
        ## }
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
          safe.pars <- setdiff(colnames(samples), pars$pheno.pars)
          safe.pars <- union(safe.pars, pars$batch.pars)
          confounded.pars <- setdiff(intersect(colnames(samples), pars$pheno.pars), safe.pars)
          all.pars <- c(safe.pars, confounded.pars)
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
        pgx_options <- if (!is.null(pgx)) {
          pgx$settings$preprocess$options
        } else {
          NULL
        }

        ## Imputation defaults
        default_zero_as_na <- FALSE
        default_filter_missing <- FALSE
        default_filter_threshold <- 0.2
        if (grepl("proteomics|metabolomics", upload_datatype())) default_zero_as_na <- TRUE
        default_impute <- DEFAULTS$qc$impute
        default_impute_method <- "SVD2"
        if (is.list(pgx_options)) {
          default_zero_as_na <- isTRUE(pgx_options$zero_as_na)
          default_filter_missing <- isTRUE(pgx_options$filter_missing)
          if (!is.null(pgx_options$filter_threshold)) {
            default_filter_threshold <- pgx_options$filter_threshold
          }
          default_impute <- isTRUE(pgx_options$impute)
          if (!is.null(pgx_options$impute_method)) {
            default_impute_method <- unname(pgx_options$impute_method[[1L]])
          }
        }

        ## Normalization defaults
        default_normalize <- !(is.olink() || is.nulisa())
        default_norm_method <- 1
        if (is.list(pgx_options)) {
          default_normalize <- isTRUE(pgx_options$normalize)
          if (!is.null(pgx_options$norm_method)) {
            default_norm_method <- if (length(pgx_options$norm_method) > 1L) {
              "multiomics"
            } else {
              unname(pgx_options$norm_method[[1L]])
            }
          }
        }

        ## Outlier removal defaults
        default_remove_outliers <- FALSE
        default_outlier_threshold <- 6
        if (is.list(pgx_options)) {
          default_remove_outliers <- isTRUE(pgx_options$remove_outliers)
          if (!is.null(pgx_options$outlier_threshold)) {
            default_outlier_threshold <- pgx_options$outlier_threshold
          }
        }

        ## Batch correction defaults
        default_batchcorrect <- FALSE
        default_bec_method <- "SVA"
        default_bec_param <- batch_params[1]
        if (is.list(pgx_options)) {
          default_batchcorrect <- isTRUE(pgx_options$batch_correct)
          if (!is.null(pgx_options$batch_method)) {
            default_bec_method <- unname(pgx_options$batch_method[[1L]])
          }
          if (!is.null(pgx_options$batch)) {
            default_bec_param <- colnames(pgx_options$batch)
          }
        }

        score.infotext <-
          "Outliers markedly deviate from the vast majority of samples. Outliers could be caused by technical factors and negatively affect data analysis. Here, outliers are identified and marked for removal should you wish so."

        missing.infotext <-
          "Missing values (MVs) reduce the completeness of biological data and hinder preprocessing steps. MVs (i.e., NA), more often populate proteomics and metabolomics data. Here, MVs are identified and their patterns in your data is shown. PCA is also optionally performed on data imputed with all methods to aid comparison."

        normalization.infotext <-
          "Normalization enables to standardize the data and improve their consistency, comparability and reproducibility. Boxplots of raw (unnormalized) and normalized data are shown. Normalization method can be selected on the left, under “Normalization”."

        batcheff.infotext <-
          paste(
            "Batch effects (BEs) are technical variation in the measurements.",
            "Method comparisons use the 1,000 most variable features by default",
            "and are labelled as an approximate preview; submitted correction",
            "uses the full processed matrix."
          )

        methyl.infotext <- "Density plot of beta values. Optionally, sample-specific beanplot of beta value distribution can be plotted."

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
            ns("bec_view"),
            label = "Show:",
            choices = c(
              "Samples (PCA)" = "pca",
              "Biplot (loadings)" = "loadings",
              "Biplot (phenotypes)" = "pheno",
              "Variance explained" = "scree"
            ),
            selected = "pca",
            inline = FALSE
          ),
          shiny::radioButtons(
            ns("colorby_var"),
            label = "Annotate by:",
            choices = metadata_vars,
            selected = metadata_vars[1],
            inline = FALSE
          ),
          shiny::selectInput(
            ns("bec_xpc"),
            label = "X-axis:",
            choices = paste0("PC", 1:5),
            selected = "PC1"
          ),
          shiny::selectInput(
            ns("bec_ypc"),
            label = "Y-axis:",
            choices = paste0("PC", 1:5),
            selected = "PC2"
          )
        )

        methyl.options <- tagList(
          shiny::radioButtons(
            ns("methyl_plottype"),
            label = "Plot type:",
            choices = c("Density", "Beanplot"),
            selected = "Density", inline = FALSE
          )
        )

        navmenu <- tagList(
          bslib::card(bslib::card_body(
            style = "padding: 0px;",
            do.call(bslib::accordion, c(
              list(multiple = FALSE, style = "background-color: #F7FAFD99;"),
              Filter(Negate(is.null), list(
                if (upload_datatype() != "methylomics") {
                  bslib::accordion_panel(
                    title = "Missing values",
                    shiny::div(
                      style = "display: flex; align-items: center; justify-content: space-between;",
                      shiny::p("Handle missing values:\n"),
                      shiny::HTML("<a href='https://bigomics.ch/blog/imputation-of-missing-values-in-proteomics' target='_blank' class='info-link' style='margin-left: 15px;'>
                      <i class='fa-solid fa-circle-info info-icon' style='color: blue; font-size: 20px;'></i>
                      </a>")
                    ),
                    shiny::checkboxInput(ns("zero_as_na"), label = "Treat zero as NA", value = default_zero_as_na),
                    shiny::checkboxInput(
                      ns("filtermissing"),
                      label = "Remove NA rows",
                      value = default_filter_missing
                    ),
                    shiny::conditionalPanel("input.filtermissing == true",
                      ns = ns,
                      shiny::selectInput(ns("filterthreshold"), NULL,
                        choices = c(
                          ">10% NA" = 0.1, ">20% NA" = 0.2, ">50% NA" = 0.5,
                          "<3 valid in any group" = 3, "<50% valid in any group" = -0.5
                        ),
                        selected = default_filter_threshold
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
                  )
                },
                bslib::accordion_panel(
                  title = "Normalization",
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
                        ## The canonical pipeline applies one declared method
                        ## to each row-prefix layer.
                        c(
                          "multi-omics per-block (gx: CPM, other: maxMedian)" = "multiomics"
                        )
                      } else {
                        c(
                          "CPM", "CPM+quantile", "TMM", "quantile",
                          "maxMedian", "maxSum", "reference"
                        )
                      },
                      selected = .normalize_selected(default_norm_method)
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
                  title = "Remove outliers",
                  shiny::p("Detect and remove outlier samples."),
                  shiny::checkboxInput(ns("remove_outliers"), "remove outliers", value = default_remove_outliers),
                  shiny::conditionalPanel("input.remove_outliers == true",
                    ns = ns,
                    shiny::sliderInput(ns("outlier_threshold"), "Select threshold:", 1, 12, default_outlier_threshold, 1)
                  ),
                  br()
                ),
                bslib::accordion_panel(
                  title = "Batch-effect correction",
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
                  shiny::checkboxInput(ns("bec_full_features"),
                    label = "Use all features for exact BC preview (slower)",
                    value = FALSE
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
                    ),
                    br(),
                    shiny::HTML("<div style='margin-top: 10px;'><a href='https://academic.oup.com/bioinformatics/article/41/3/btaf084/8042340' target='_blank' style='color: #0066cc; text-decoration: none;'>Learn about NPM <i class='fa-solid fa-external-link' style='font-size: 12px;'></i></a></div>")
                  )
                )
              ))
            )),
            br()
          ))
        )

        ## ---------------------------- UI ----------------------------------
        ui <- div(
          bslib::as_fill_carrier(),
          style = "width: 100%; display: flex; ",
          bslib::layout_columns(
            col_widths = c(2, 10),
            style = "margin-bottom: 0px;",
            heights_equal = "row",
            ## ----------- menu ------------
            navmenu,
            ## ----------- canvas ------------
            bslib::layout_columns(
              col_widths = c(6, 6),
              row_heights = c(3, 3),
              heights_equal = "row",
              ## --------new
              if (upload_datatype() == "methylomics") {
                PlotModuleUI(
                  ns("plot5"),
                  title = "Distribution of Beta values",
                  info.text = methyl.infotext,
                  caption = methyl.infotext,
                  options = methyl.options,
                  height = c("auto", "100%"),
                  show.maximize = FALSE
                )
                ## --------new
              } else {
                PlotModuleUI(
                  ns("plot2"),
                  title = "Missing values",
                  info.text = missing.infotext,
                  caption = missing.infotext,
                  options = missing.options,
                  height = c("auto", "100%"),
                  show.maximize = FALSE
                )
              },
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

      PlotModuleServer(
        "plot5",
        plotlib = "base",
        func = plot_methyl,
        res = c(75, 120),
        pdf.width = 12,
        pdf.height = 6,
        add.watermark = FALSE
      )

      counts <- reactive({
        shiny::req(dim(r_counts()))
        r_counts()
      })

      cX <- reactive({
        shiny::req(dim(correctedX()$X))
        return(correctedX()$X)
      })

      imputation_method <- reactive({
        ll <- list(zero_as_na = zero_as_na(), imputation = input$impute_method)
        if (!isTRUE(input$impute)) {
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

      ## Canonical options are the single preprocessing description shared by
      ## staged previews and the submitted raw matrix.
      preprocess <- reactive({
        batch_inputs <- .opg_upload_batch_inputs(
          counts = r_counts(),
          samples = r_samples(),
          contrasts = r_contrasts(),
          selection = input$bec_param,
          enabled = isTRUE(input$batchcorrect)
        )
        .opg_upload_preprocess_options(
          counts = r_counts(),
          datatype = upload_datatype(),
          is_npx = isTRUE(is.olink()) || isTRUE(is.nulisa()),
          zero_as_na = zero_as_na(),
          filter_missing = isTRUE(input$filtermissing),
          filter_threshold = input$filterthreshold,
          impute = isTRUE(input$impute),
          impute_method = input$impute_method,
          normalize = isTRUE(input$normalize),
          norm_method = input$normalization_method,
          ref_gene = if (identical(input$normalization_method, "reference")) {
            input$ref_gene
          } else {
            NULL
          },
          remove_outliers = isTRUE(input$remove_outliers),
          outlier_threshold = input$outlier_threshold,
          batch_correct = isTRUE(input$batchcorrect),
          batch_method = if (is.null(input$bec_method)) "SVA" else input$bec_method,
          batch = batch_inputs$batch,
          target = batch_inputs$target,
          meth_type = meth_type(),
          max_features = NULL
        )
      })

      return(
        list(
          counts = counts,
          X = cX,
          norm_method = norm_method,
          imputation_method = imputation_method,
          bc_method = bc_method,
          remove_outliers = remove_outliers,
          annot = annot,
          preprocess = preprocess
        )
      ) ## pointing to reactive
    } ## end-of-server
  )
}
