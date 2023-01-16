##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
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
drugconnectivity_plot_enplots_ui <- function(id,
                                            label = "",
                                            height = c(600, 800)) {
  ns <- shiny::NS(id)
  info_text <- strwrap("<strong>Drug connectivity</strong> correlates your
                       signature with known drug profiles from the L1000
                       database, and shows similar and opposite profiles by
                       running the GSEA algorithm on the drug profile
                       correlation space.")
  plot_opts <- shiny::tagList()

  PlotModuleUI(ns("plot"),
               title = "Drug connectivity",
               label = label,
               plotlib = "base",
               info.text = info_text,
               options = plot_opts,
               download.fmt = c("png", "pdf", "csv"),
               width = c("auto", "100%"),
               height = height
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
drugconnectivity_plot_enplots_server <- function(id,
                                                 pgx,
                                                 dsea_contrast,
                                                 dsea_method,
                                                 dsea_table,
                                                 getActiveDSEA,
                                                 watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {


      plot_data <- shiny::reactive({
        pgx <- pgx()
        dsea_contrast <- dsea_contrast()
        dsea_method <- dsea_method()
        shiny::req(pgx, dsea_contrast, dsea_method)
        shiny::validate(shiny::need("drugs" %in% names(pgx),
                                    "no 'drugs' in object."))
        if (is.null(pgx$drugs)) {
          return(NULL)
        }

        if (is.null(dsea_contrast)) {
          return(NULL)
        }

        res <- list(
          pgx = pgx,
          dsea_contrast = dsea_contrast,
          dsea_method = dsea_method
        )

        return(res)
      })

      plot.RENDER <- shiny::reactive({
        res <- plot_data()
        pgx <- res$pgx
        dsea_contrast <- res$dsea_contrast
        dsea_method <- res$dsea_method

        dsea <- getActiveDSEA()
        dt <- dsea$table

        ## filter with table selection/search
        ii <- dsea_table$rows_selected()
        jj <- dsea_table$rows_all()
        shiny::req(jj) ## must have non-empty table

        if (length(ii) > 0) {
          dt <- dt[ii, , drop = FALSE]
        }
        if (length(ii) == 0 && length(jj) > 0) {
          dt <- dt[jj, , drop = FALSE]
        }

        if (nrow(dt) == 0) {
          return(NULL)
        }

        ## rank vector for enrichment plots
        dmethod <- dsea_method
        rnk <- dsea$stats
        if (length(rnk) == 0) {
          return(NULL)
        }

        ## ENPLOT TYPE
        if (nrow(dt) == 1) {
          par(oma = c(1, 1, 1, 1))
          par(mfrow = c(1, 1), mar = c(4, 4, 1.1, 2), mgp = c(2.3, 0.9, 0))
          lab.cex <- 1
          xlab <- "Rank in ordered dataset"
          ylab <- "Rank metric"
          nc <- 1
        } else {
          dt <- head(dt, 16)
          lab.cex <- 0.75
          xlab <- ""
          ylab <- ""
          nc <- ceiling(sqrt(nrow(dt)))
          par(oma = c(0, 1.6, 0, 0))
          par(mfrow = c(nc, nc), mar = c(0.3, 1.0, 1.3, 0), mgp = c(1.9, 0.6, 0))
        }

        for (i in 1:nrow(dt)) {
          dx <- rownames(dt)[i]
          gmtdx <- grep(dx, names(rnk), fixed = TRUE, value = TRUE) ## L1000 naming
          dx1 <- substring(dx, 1, 26)
          par(cex.axis = 0.001)
          if (i %% nc == 1) par(cex.axis = 0.98)
          suppressWarnings(
            gsea.enplot(rnk, gmtdx,
                        main = dx1, cex.main = 1.2,
                        xlab = xlab, ylab = ylab
            )
          )
          nes <- round(dt$NES[i], 2)
          qv <- round(dt$padj[i], 3)
          tt <- c(paste("NES=", nes), paste("q=", qv))
          legend("topright", legend = tt, cex = 0.8, y.intersp = 0.85, bty = "n")
          if (i %% nc == 1 && nrow(dt) > 1) {
            mtext("rank metric", side = 2, line = 1.8, cex = lab.cex)
          }
        }
      })

      PlotModuleServer(
        "plot",
        plotlib = "base", # does not use plotly
        func = plot.RENDER,
        func2 = plot.RENDER,
        csvFunc = plot_data,
        res = c(78, 110),
        pdf.width = 6.5, pdf.height = 12.8,
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
