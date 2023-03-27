##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


clustering_plot_phenoplot_ui <- function(id,
                                         label = "",
                                         height) {
  ns <- shiny::NS(id)

  info <- tagsub("<strong>Phenotype distribution.</strong> This figure visualizes
  the distribution of the available phenotype data. You can choose to put the group labels in the
  figure or as separate legend in the {Label} setting, in the plot {{settings}}")

  phenoplot.opts <- shiny::tagList(
    shiny::radioButtons(ns("labelmode"), "Label", c("groups", "legend"), inline = TRUE)
  )

  PlotModuleUI(
    ns("pltmod"),
    title = "Phenotype distribution",
    label = label,
    plotlib = "base",
    info.text = info,
    options = phenoplot.opts,
    download.fmt = c("png", "pdf", "csv"),
    width = c("auto", "100%"),
    height = height
  )
}

clustering_plot_phenoplot_server <- function(id,
                                             pgx,
                                             hm_getClusterPositions,
                                             watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    plot_data <- reactive({
      pgx <- pgx
      shiny::req(pgx$Y)

      ## get t-SNE positions
      clust <- hm_getClusterPositions()
      ## pos = pgx$tsne2d
      pos <- clust$pos
      colnames(pos) <- c("x","y")
      
      Y <- pgx$Y[rownames(pos), , drop = FALSE]
      pheno <- colnames(Y)

      ## don't show these...
      pheno <- grep("batch|sample|donor|repl|surv", pheno,
        invert = TRUE, ignore.case = TRUE, value = TRUE
      )
      Y <- Y[,pheno] 
      
      ## complete dataframe for downloading
      df <- data.frame( pos, Y)
      
      return(
        list(
          df = df,
          pheno = pheno,
          labelmode = input$labelmode
        )
      )
    })

    renderPlots <- function( pd, mfrow, mar ) {

      dbg("[clustering_plot_phenoplot.R]")
      pheno <- pd[["pheno"]]
      Y <- pd[["df"]][,pheno]
      labelmode <- pd[["labelmode"]]
      pos <- pd[["df"]][,c("x","y")]

      ## set layout
      par(mfrow = mfrow, mar = mar)
      i <- 1

      cex1 <- 1.1 * c(1.8, 1.3, 0.8, 0.5)[cut(nrow(pos), breaks = c(-1, 40, 200, 1000, 1e10))]
      cex1 <- cex1 * ifelse(length(pheno) > 6, 0.8, 1)
      cex1 <- cex1 * ifelse(length(pheno) > 12, 0.8, 1)

      for (i in 1:min(20, length(pheno))) {
        ## ------- set colors
        colvar <- factor(Y[, 1])
        colvar <- factor(Y[, pheno[i]])
        colvar[which(colvar %in% c(NA, "", " ", "NA", "na"))] <- NA
        colvar <- factor(as.character(colvar))
        klrpal <- COLORS
        klr1 <- klrpal[colvar]
        klr1 <- paste0(gplots::col2hex(klr1), "99")
        jj <- which(is.na(klr1))
        if (length(jj)) klr1[jj] <- "#AAAAAA22"
        tt <- tolower(pheno[i])

        ## ------- start plot
        base::plot(pos[, ],
          pch = 19, cex = cex1, col = klr1,
          fg = gray(0.5), bty = "o", xaxt = "n", yaxt = "n",
          xlab = "tSNE1", ylab = "tSNE2"
        )
        title(tt, cex.main = 1.3, line = 0.5, col = "grey40")
        if (labelmode == "legend") {
          legend("bottomright",
            legend = levels(colvar), fill = klrpal,
            cex = 0.95, y.intersp = 0.85, bg = "white"
          )
        } else {
          grp.pos <- apply(pos, 2, function(x) tapply(x, colvar, mean, na.rm = TRUE))
          grp.pos <- apply(pos, 2, function(x) tapply(x, colvar, median, na.rm = TRUE))
          nvar <- length(setdiff(colvar, NA))
          if (nvar == 1) {
            grp.pos <- matrix(grp.pos, nrow = 1)
            rownames(grp.pos) <- setdiff(colvar, NA)[1]
          }
          labels <- rownames(grp.pos)
          boxes <- sapply(nchar(labels), function(n) paste(rep("\u2588", n), collapse = ""))
          cex2 <- 0.9 * cex1**0.33
          text(grp.pos, labels = boxes, cex = cex2 * 0.95, col = "#CCCCCC99")
          text(grp.pos, labels = labels, font = 2, cex = cex2)
        }
      }
    }

    plot.RENDER <- function() {

      pd <- plot_data()
      pheno <- pd[["pheno"]]
      
      ## layout
      mfrow = c(3,2)
      mar = c(0.3, 0.7, 2.8, 0.7)
      if (length(pheno) >= 6) {
        mfrow = c(4, 3)
        mar = c(0.3, 0.4, 2.8, 0.4) * 0.8
      }
      if (length(pheno) >= 12) {
        mfrow = c(5, 4)
        mar = c(0.2, 0.2, 2.5, 0.2) * 0.8
      }
      renderPlots( pd, mfrow, mar ) 
    }

    plot.RENDER2 <- function() {

      pd <- plot_data()
      pheno <- pd[["pheno"]]

      ## layout
      mfrow = c(2,3)
      mar = c(0.3, 0.7, 2.8, 0.7)
      if (length(pheno) >= 6) {
        mfrow = c(3, 4)
        mar = c(0.3, 0.4, 2.8, 0.4) * 0.8
      }
      if (length(pheno) >= 12) {
        mfrow = c(4, 5)
        mar = c(0.2, 0.2, 2.5, 0.2) * 0.8
      }
      renderPlots( pd, mfrow, mar ) 
    }
    
    
    PlotModuleServer(
      "pltmod",
      plotlib = "base",
      ## plotlib2 = "plotly",
      func = plot.RENDER,
      func2 = plot.RENDER2,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(85), ## resolution of plots
      pdf.width = 6, pdf.height = 9,
      add.watermark = watermark
    )
  })
}
