##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Drug Connectivity plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#'
#' @export
drugconnectivity_plot_cmap_dsea_ui <- function(
    id,
    title,
    info.text,
    caption,
    label = "",
    height) {
  ns <- shiny::NS(id)

  plot_opts <- shiny::tagList(
    withTooltip(
      shiny::radioButtons(
        ns("cmap_labeltype"), "label type:",
        c("drugs", "MOA", "target"),
        inline = TRUE
      ),
      "Label point with drugs, MOA terms or targets (if no drug selected).",
      placement = "left", options = list(container = "body")
    ),
    withTooltip(
      shiny::radioButtons(ns("cmap_nlabel"), "number of labels:", c(3, 10, 20, 100),
        selected = 10, inline = TRUE
      ),
      "Number of labels to show.",
      placement = "left", options = list(container = "body")
    ),
    withTooltip(
      shiny::checkboxGroupInput(
        ns("cmap_labeloptions"), "label options:",
        choices = c("show", "fixed"),
        selected = c("show"), inline = TRUE
      ), "Other labels options.",
      placement = "left", options = list(container = "body")
    )
  )

  PlotModuleUI(ns("plot"),
    title = title,
    caption = caption,
    label = label,
    plotlib = "base",
    info.text = info.text,
    options = plot_opts,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height,
    width = c("auto", "100%")
  )
}

#' Drug Connectivity plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
drugconnectivity_plot_cmap_dsea_server <- function(id,
                                                   pgx,
                                                   getActiveDSEA,
                                                   cmap_table,
                                                   getMOA.class,
                                                   getMOA.target,
                                                   dsea_method,
                                                   dsea_contrast,
                                                   watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      plotCMAP <- function(pgx, db, contr, moa.target, moa.class,
                           labtype = "drugs", showlabel = TRUE,
                           npoints = 100, nlabel = 10,
                           lab.wt = TRUE, lab.gamma = 1, lab.cex = 1,
                           opacity = 0.15, softmax = 1,
                           title = NULL, plotlib = "base") {
        if (!"drugs" %in% names(pgx)) {
          frame()
          text(0.5, 0.5, "Error: PGX object does not have CMAP results", col = "red3")
          return(NULL)
        }

        res <- pgx$drugs[[db]]
        if (!"clust" %in% names(res)) {
          frame()
          text(0.5, 0.5, "Error: PGX object does not have CMAP cluster positions", col = "red3")
          return(NULL)
        }

        pos <- res$clust
        var <- res$stats[, contr]
        ## highlight genes from table
        nes1 <- res$X[, contr] ## NES for drugs

        ## compute median position of drugs
        xdrugs <- gsub("[_@].*", "", rownames(pos))
        pos1 <- apply(pos, 2, function(x) tapply(x, xdrugs, median))
        pos1 <- pos1[names(nes1), ]
        var1 <- nes1

        ## compute median position of moa class
        xmoa <- res$annot[xdrugs, "moa"]
        xmoa <- strsplit(xmoa, split = "[\\|;,]")
        nmoa <- sapply(xmoa, length)
        ii <- unlist(mapply(rep, 1:nrow(pos), nmoa))
        pos2 <- apply(pos[ii, ], 2, function(x) tapply(x, unlist(xmoa), median))
        moa2 <- moa.class
        nes2 <- array(moa2$NES, dimnames = list(moa2$pathway))
        var2 <- nes2[rownames(pos2)]

        ## compute median position of targets
        xtarget <- res$annot[xdrugs, "target"]
        xtarget <- strsplit(xtarget, split = "[\\|;,]")
        ntarget <- sapply(xtarget, length)
        ii <- unlist(mapply(rep, 1:nrow(pos), ntarget))
        pos3 <- apply(pos[ii, ], 2, function(x) tapply(x, unlist(xtarget), median))
        #
        moa3 <- moa.target
        nes3 <- array(moa3$NES, dimnames = list(moa3$pathway))
        var3 <- nes3[rownames(pos3)]

        ## create extended positions and variable (drugs, moa, target)
        var <- var / max(abs(var), na.rm = TRUE) ## replicate level
        var1 <- var1 / max(abs(var1), na.rm = TRUE) ## drug level
        var2 <- var2 / max(abs(var2), na.rm = TRUE) ## MOA class level
        var3 <- var3 / max(abs(var3), na.rm = TRUE) ## target level

        rownames(pos) <- paste0("0:", rownames(pos))
        rownames(pos1) <- paste0("1:", rownames(pos1))
        rownames(pos2) <- paste0("2:", rownames(pos2))
        rownames(pos3) <- paste0("3:", rownames(pos3))
        names(var) <- paste0("0:", names(var))
        names(var1) <- paste0("1:", names(var1))
        names(var2) <- paste0("2:", names(var2))
        names(var3) <- paste0("3:", names(var3))

        xpos <- rbind(pos, pos1, pos2, pos3)
        xvar <- c(var, var1, var2, var3)
        sum(duplicated(names(xvar)))

        h1 <- h2 <- NULL
        ## limit number of labels/points
        if (labtype == "replicate") {
          xx <- res$stats[, contr]
          names(xx) <- paste0("0:", names(xx))
        } else if (labtype == "drugs") {
          xx <- res$X[, contr]
          names(xx) <- paste0("1:", names(xx))
        } else if (toupper(labtype) == "MOA") {
          xx <- moa.class$NES
          names(xx) <- moa.class$pathway
          names(xx) <- paste0("2:", names(xx))
        } else if (labtype == "target") {
          xx <- moa.target$NES
          names(xx) <- moa.target$pathway
          names(xx) <- paste0("3:", names(xx))
        } else {
          return(NULL)
        }
        xx <- xx[order(-abs(xx))]
        h1 <- head(names(xx), npoints)
        h2 <- head(h1, nlabel)
        if (!showlabel) {
          h2 <- NULL
        }

        if (is.null(title)) title <- contr

        plt <- NULL
        if (plotlib == "base") {
          wcex <- lab.cex
          if (lab.wt) {
            wcex <- 1.2 * (abs(xvar) / max(abs(xvar), na.rm = TRUE))**lab.gamma
            wcex <- lab.cex * wcex
            names(wcex) <- names(xvar)
          }
          wcex[is.na(wcex)] <- 1

          playbase::pgx.scatterPlotXY.BASE(
            xpos,
            var = xvar, title = title,
            xlab = "UMAP-x", ylab = "UMAP-y",
            labels = sub(".*:", "", rownames(xpos)),
            hilight = h1, hilight2 = h2, hilight.cex = 1.1,
            cex = 1, cex.lab = wcex, cex.title = 0.95,
            legend = TRUE, zsym = TRUE,
            rstep = 0.2, dlim = 0.01,
            softmax = softmax, opacity = opacity
          )
        } else {
          plt <- playbase::pgx.scatterPlotXY(
            xpos,
            var = xvar, plotlib = plotlib, title = title,
            xlab = "UMAP-x", ylab = "UMAP-y",
            hilight = h1, hilight2 = h2, hilight.cex = 1.1,
            cex = 1, cex.lab = cex.lab, cex.title = 1.0,
            legend = TRUE, zsym = TRUE,
            softmax = 1, opacity = opacity
          )
        }
        plt
      }

      plot_data <- shiny::reactive({
        res <- list(
          pgx = pgx,
          dsea = getActiveDSEA(),
          cmap_table = cmap_table,
          moa.class = getMOA.class(),
          moa.target = getMOA.target(),
          dsea_method = dsea_method(),
          dsea_contrast = dsea_contrast()
        )
        return(res)
      })

      plot.RENDER <- shiny::reactive({
        res <- plot_data()
        pgx <- res$pgx
        cmap_table <- res$cmap_table
        shiny::req(pgx$X)
        db <- "L1000/gene"
        contr <- "treatment:Gefitinib_vs_CT"

        db <- res$dsea_method
        contr <- res$dsea_contrast
        shiny::req(db)
        shiny::req(contr)

        dsea <- res$dsea

        ## get reactive values
        rows_selected <- 1
        rows_all <- 1:nrow(dsea$table)
        rows_selected <- cmap_table$rows_selected()
        rows_all <- cmap_table$rows_all()

        if (is.null(rows_all) || length(rows_all) == 0) {
          return(NULL)
        }

        drugs_all <- rownames(dsea$table)[rows_all]

        moa.class <- res$moa.class
        moa.target <- res$moa.target

        labtype <- input$cmap_labeltype
        nlabel <- as.integer(input$cmap_nlabel)
        showlab <- ("show" %in% input$cmap_labeloptions)
        lab.wt <- !("fixed" %in% input$cmap_labeloptions)

        ## ---------------  plot -------------------
        all.contr <- colnames(pgx$contrasts)
        contr <- all.contr[1]

        nr <- ceiling(sqrt(length(all.contr)))
        par(mfrow = c(nr, nr))
        for (contr in all.contr) {
          tt <- paste0(contr, " (", toupper(labtype), ")")
          #
          plotCMAP(pgx, db, contr,
            moa.target, moa.class,
            labtype = labtype, showlabel = showlab,
            lab.wt = lab.wt, lab.gamma = 1, lab.cex = 1.6,
            opacity = 0.15, softmax = 0,
            npoints = nlabel, nlabel = nlabel,
            title = tt, plotlib = "base"
          )
        }
      })

      PlotModuleServer(
        "plot",
        plotlib = "base",
        func = plot.RENDER,
        func2 = plot.RENDER,
        csvFunc = plot_data,
        res = c(80, 105),
        pdf.width = 8, pdf.height = 8,
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
