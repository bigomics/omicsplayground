##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


biomarker_plot_featurerank_ui <- function(
    id,
    title,
    info.text,
    info.methods,
    info.references,
    info.extra_link,
    caption,
    label = "",
    height,
    width) {
  ns <- shiny::NS(id)

  clust_featureRank.opts <- shiny::tagList(
    withTooltip(
      shiny::radioButtons(ns("clust_featureRank_method"), "Method:",
        choices = c("p-value", "correlation", "meta"),
        inline = TRUE
      ),
      "Choose ranking method: p-value based or correlation-based.",
      placement = "right", options = list(container = "body")
    )
  )



  PlotModuleUI(
    ns("pltmod"),
    label = label,
    plotlib = "plotly",
    title = title,
    caption = caption,
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    options = clust_featureRank.opts,
    download.fmt = c("png", "pdf", "csv"),
    width = width,
    height = height
  )
}

biomarker_plot_featurerank_server <- function(id,
                                              pgx,
                                              ft_level,
                                              samplefilter,
                                              watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    calcFeatureRanking <- shiny::reactive({
      pgx <- pgx
      ft_level <- ft_level()

      shiny::req(pgx$X, pgx$Y, pgx$gsetX, pgx$genes)

      features <- X <- NULL
      if (ft_level == "geneset") {
        gset_collections <- playbase::pgx.getGeneSetCollections(gsets = rownames(pgx$gsetX))
        features <- gset_collections
        X <- pgx$gsetX
        ## if(any(is.na(X))) {
        ##     X <- X[complete.cases(X), , drop = FALSE]
        ## }
      } else {
        features <- pgx$families
        X <- pgx$X
        if (any(is.na(X))) {
          X <- pgx$impX
        }
      }

      ## ------------ intersect features, set minimum set size

      genes <- toupper(pgx$genes$symbol)
      features <- lapply(features, toupper)
      features <- lapply(features, function(f) intersect(toupper(f), genes))
      features <- features[sapply(features, length) >= 10]

      ## ------------ Just to get current samples
      samples <- playbase::selectSamplesFromSelectedLevels(pgx$Y, samplefilter())
      X <- X[, samples]

      ## the code below overwrittes user input, and should be removed

      cvar <- playbase::pgx.getCategoricalPhenotypes(pgx$Y, max.ncat = 999)
      cvar <- grep("sample|patient|years|days|months|gender",
        cvar,
        invert = TRUE, value = TRUE
      ) ## no sample IDs
      Y <- pgx$Y[colnames(X), cvar, drop = FALSE]
      kk <- which(apply(Y, 2, function(y) length(unique(y)) > 1))
      shiny::validate(shiny::need(length(kk) > 0, "Number of samples filtered is too restrictive. Change the 'Filter samples'."))
      Y <- Y[, kk, drop = FALSE]
      dim(Y)

      ## ------------ Note: this takes a while. Maybe better precompute off-line...
      sdx <- apply(X, 1, sd, na.rm = TRUE)
      names(sdx) <- rownames(X)
      S <- matrix(NA, nrow = length(features), ncol = ncol(Y))
      rownames(S) <- names(features)
      colnames(S) <- colnames(Y)

      ## ------------ Create a Progress object
      if (!interactive()) {
        progress <- shiny::Progress$new()
        on.exit(progress$close())
        progress$set(message = "Calculating feature-set scores", value = 0)
      }

      gene.level <- TRUE
      gene.level <- (ft_level == "gene")

      for (i in 1:ncol(Y)) {
        if (!interactive()) progress$inc(1 / ncol(Y))

        grp <- Y[, i]
        grp <- as.character(grp)

        score <- rep(NA, length(features))
        names(score) <- names(features)
        for (j in 1:length(features)) {
          pp <- features[[j]]
          if (gene.level) {
            annot <- pgx$genes
            annot <- annot[match(unique(annot$symbol), annot$symbol), ]
            rownames(annot) <- annot$symbol
            pp <- playbase::filterProbes(annot, features[[j]])
            ## pp <- playbase::filterProbes(pgx$genes, features[[j]])
          }
          pp <- head(pp[order(-sdx[pp])], 1000) ## how many top SD??

          X1 <- playbase::rename_by(X, pgx$genes, "symbol")
          pp <- intersect(pp, rownames(X1))
          X1 <- X1[pp, , drop = FALSE]
          if (nrow(X1) == 0) {
            score[j] <- NA
            next
          }

          s1 <- s2 <- 1
          method <- input$clust_featureRank_method
          if (method %in% c("correlation", "meta")) {
            mx <- t(apply(X1, 1, function(x) tapply(x, grp, mean, na.rm = TRUE)))
            if (nrow(mx) == 0 || ncol(mx) == 0) next
            ## D <- 1 - cor(mx, use = "pairwise.complete.obs")
            D <- 1 - cor(mx, use = "pairwise")
            diag(D) <- NA
            s1 <- mean(D, na.rm = TRUE)
          }

          if (method %in% c("p-value", "meta")) {
            jj <- which(!is.na(grp))
            design <- model.matrix(~ grp[jj])
            fit <- tryCatch({
              suppressWarnings(limma::eBayes(limma::lmFit(X1[, jj, drop = FALSE], design)))
            }, error = function(w){NA})
            if (all(is.na(fit))) {
              score[j] <- NA
              next
            }
            suppressWarnings(suppressMessages(top <- limma::topTable(fit)))
            s2 <- mean(-log10(1e-99 + top$adj.P.Val), na.rm = TRUE)
          }

          f <- 1
          f <- (1 - exp(-(length(pp) / 20)**2)) ## penalize smaller sets
          score[j] <- f * (s1 * s2)**ifelse(method == "meta", 0.5, 1)
        }
        S[, i] <- score
      }

      S[is.na(S)] <- 0 ## missing values

      return(S)
    })

    render_featureRank <- function() {
      S <- calcFeatureRanking()
      if (is.null(S) || nrow(S) == 0 || ncol(S) == 0) {
        return(NULL)
      }

      ## top scoring
      S <- tail(S[order(rowSums(S, na.rm = TRUE)), , drop = FALSE], 25)
      rownames(S) <- paste(substring(rownames(S), 1, 50), "  ")

      playbase::pgx.stackedBarplot(
        x = t(S),
        showlegend = TRUE,
        xlab = "Discriminant score",
        ylab = "",
        horiz = TRUE
      )
    }

    clust_featureRank.RENDER <- function() {
      render_featureRank() %>%
        plotly_default() %>%
        plotly::layout(
          legend = list(orientation = "h")
        )
    }

    clust_featureRank.RENDER2 <- function() {
      render_featureRank() %>%
        plotly_modal_default()
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = clust_featureRank.RENDER,
      func2 = clust_featureRank.RENDER2,
      csvFunc = calcFeatureRanking, ##  *** downloadable data as CSV
      res = c(72, 90), ## resolution of plots
      pdf.width = 8, pdf.height = 10,
      add.watermark = watermark
    )
  })
}
