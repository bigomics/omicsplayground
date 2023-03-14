##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

#' Single cell plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#' @param width
#'
#' @export
singlecell_plot_phenoplot_ui <- function(id,
                                         label = "",
                                         height,
                                         width) {
  ns <- shiny::NS(id)

  phenoplot.opts <- shiny::tagList(
    withTooltip(
      shiny::radioButtons(ns("labelmode"), "Label:", c("groups", "legend"), inline = TRUE),
      "Select whether you want the group labels to be plotted inside the plots or in a seperate legend."
    )
  )

  phenoModule_info <- "<b>Phenotype plots.</b> The plots show the distribution of the phenotypes superposed on the t-SNE clustering. Often, we can expect the t-SNE distribution to be driven by the particular phenotype that is controlled by the experimental condition or unwanted batch effects."



  PlotModuleUI(ns("plot"),
    label = label,
    info.text = phenoModule_info,
    title = "Phenotypes",
    options = phenoplot.opts,
    download.fmt = c("png", "pdf", "csv"),
    height = height,
    width = width
  )
}

#' Single cell plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @export
singlecell_plot_phenoplot_server <- function(id,
                                             pgx,
                                             pfGetClusterPositions,
                                             watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    plot_data <- shiny::reactive({
      ## if(is.null(pgx)) return(NULL)
      shiny::req(pgx)
      clust.pos <- pfGetClusterPositions()
      if (is.null(clust.pos)) {
        return(NULL)
      }

      pos <- pgx$tsne2d
      pos <- clust.pos
      sel <- rownames(pos)
      pheno <- colnames(pgx$Y)
      return(list(
        pgx = pgx,
        pos = pos,
        sel = sel,
        pheno = pheno
      ))
    })

    plot.render <- function() {
      ## if(!input$tsne.all) return(NULL)

      pd <- plot_data()


      ## layout
      par(mfrow = c(2, 2), mar = c(0.3, 0.7, 2.8, 0.7))
      if (length(pd[["pheno"]]) > 4) par(mfrow = c(3, 2), mar = c(0.3, 0.7, 2.8, 0.7))
      if (length(pd[["pheno"]]) > 6) par(mfrow = c(4, 3), mar = c(0.3, 0.4, 2.8, 0.4) * 0.8)
      if (length(pd[["pheno"]]) > 12) par(mfrow = c(5, 4), mar = c(0.2, 0.2, 2.5, 0.2) * 0.8)

      cex1 <- 1.2 * c(1.8, 1.3, 0.8, 0.5)[cut(nrow(pd[["pos"]]), breaks = c(-1, 40, 200, 1000, 1e10))]
      cex1 <- cex1 * ifelse(length(pd[["pheno"]]) > 6, 0.8, 1)
      cex1 <- cex1 * ifelse(length(pd[["pheno"]]) > 12, 0.8, 1)

      ## is it a float/number???
      is.num <- function(y, fmin = 0.1) {
        suppressWarnings(numy <- as.numeric(as.character(y)))
        t1 <- !all(is.na(numy)) && is.numeric(numy)
        t2 <- (length(unique(y)) / length(y)) > fmin
        (t1 && t2)
      }

      i <- 6
      for (i in 1:min(20, length(pd[["pheno"]]))) {
        px <- 4
        px <- pd[["pheno"]][i]
        y <- pd[["pgx"]]$Y[pd[["sel"]], px]
        y[which(y %in% c(NA, "", " ", "NA", "na"))] <- NA
        if (sum(!is.na(y)) == 0) next

        if (is.num(y)) {
          klrpal <- colorRampPalette(c("grey90", "grey50", "red3"))(16)
          y <- rank(as.numeric(y))
          ny <- round(1 + 15 * (y - min(y)) / (max(y) - min(y)))
          klr0 <- klrpal[ny]
        } else {
          y <- factor(as.character(y))
          klrpal <- COLORS
          klrpal <- paste0(gplots::col2hex(klrpal), "99")
          klr0 <- klrpal[y]
        }

        jj <- which(is.na(klr0))
        if (length(jj)) klr0[jj] <- "#AAAAAA22"
        base::plot(pd[["pos"]],
          pch = 19, cex = cex1, col = klr0, fg = gray(0.5), bty = "o",
          xaxt = "n", yaxt = "n", xlab = "tSNE1", ylab = "tSNE2"
        )

        if (!is.num(y)) {
          if (input$labelmode == "legend") {
            legend("bottomright",
              legend = levels(y), fill = klrpal,
              cex = 0.9, y.intersp = 0.8, bg = "white"
            )
          } else {
            grp.pos <- apply(pd[["pos"]], 2, function(x) tapply(x, y, mean, na.rm = TRUE))
            grp.pos <- apply(pd[["pos"]], 2, function(x) tapply(x, y, median, na.rm = TRUE))
            nvar <- length(setdiff(y, NA))
            if (nvar == 1) {
              grp.pos <- matrix(grp.pos, nrow = 1)
              rownames(grp.pos) <- setdiff(y, NA)[1]
            }

            labels <- rownames(grp.pos)
            ## title("\u2591\u2592\u2593")
            boxes <- sapply(nchar(labels), function(n) paste(rep("\u2588", n), collapse = ""))
            boxes <- sapply(nchar(labels), function(n) paste(rep("█", n), collapse = ""))
            ## boxes = sapply(nchar(labels),function(n) paste(rep("#",n),collapse=""))
            cex2 <- c(1.3, 1.1, 0.9, 0.7)[cut(length(labels), breaks = c(-1, 5, 10, 20, 999))]
            text(grp.pos, labels = boxes, cex = cex2, col = "#CCCCCC99")
            text(grp.pos, labels = labels, font = 2, cex = 1.1 * cex2, col = "white")
            text(grp.pos, labels = labels, font = 2, cex = cex2)
            ## text( grp.pos[,], labels=rownames(grp.pos), font=2, cex=cex1**0.5)
          }
        }
        title(tolower(pd[["pheno"]][i]), cex.main = 1.3, line = 0.5, col = "grey40")
      }
    }

    PlotModuleServer(
      "plot",
      func = plot.render,
      res = c(85, 95),
      pdf.width = 12, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
