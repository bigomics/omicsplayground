##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Boxplots plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#'
#' @export
biomarker_plot_boxplots_ui <- function(
    id,
    title,
    info.text,
    info.methods,
    info.extra_link,
    caption,
    label = "",
    height,
    width) {
  ns <- shiny::NS(id)

  plot_options <- tagList(
    withTooltip(
      shiny::checkboxInput(
        ns("show_all"), "show all samples",
        value = FALSE
      ), "Show all samples or only selected.",
      placement = "right", options = list(container = "body")
    )
  )

  PlotModuleUI(ns("plot"),
    title = title,
    label = label,
    plotlib = "base",
    info.text = info.text,
    info.methods = info.methods,
    info.extra_link = info.extra_link,
    options = NULL,
    caption = caption,
    download.fmt = c("png", "pdf", "csv", "svg"),
    width = width,
    height = height
  )
}

#' Boxplots plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
biomarker_plot_boxplots_server <- function(id,
                                           pgx,
                                           calcVariableImportance,
                                           is_computed,
                                           watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      plot_data <- shiny::reactive({
        res <- calcVariableImportance()
        shiny::req(res)
        shiny::req(is_computed())
        
        ## get variables used in the tree solution
        leafs <- setdiff(res$rf$frame$var, "<leaf>")
        vars <- res$rf$orig.names[leafs]
        
        if (length(vars) == 0) {
          return(NULL)
        }
        vars <- intersect(vars, rownames(res$X))

        ## add some other variables
        if (length(vars) < 8) {
          jj <- order(-rowSums(res$R, na.rm = TRUE))
          other.vars <- rownames(res$R)[jj]
          vars <- head(unique(c(vars, other.vars)), 8)
        }

        y <- res$y
        is.surv <- grepl("Surv", res$rf$call)[2]
        is.surv
        if (is.surv) {
          y <- paste0("N", res$rf$where)
          names(y) <- colnames(res$X)
          table(y)
        }

        ## vars, X, y
        pdata <- list(
          X = res$X,
          vars = vars,
          y = y
        )
        return(pdata)
      })

      plot.RENDER <- function() {
        pdata <- plot_data()

        shiny::req(pdata)

        ## vars, X, y
        X <- pdata$X
        vars <- pdata$vars
        y <- pdata$y
        ny <- length(unique(y))

        par(
          mfrow = c(2, 4), mar = c(3.0, 3.0, 1.5, 0.5),
          mgp = c(1.6, 0.6, 0), oma = c(0.5, 0.5, 0.5, 0.5) * 0
        )
        if (length(vars) > 8) par(mfrow = c(3, 4), mar = c(2.8, 3, 2, 0.3) * 0.7)
        i <- 1
        for (i in 1:min(12, length(vars))) {
          g <- vars[i]
          gx <- X[g, ]
          g <- playbase::probe2symbol(g, pgx$genes, "gene_name", fill_na = TRUE)
          boxplot(
            gx ~ y,
            col = "grey85",
            ylim = range(gx),
            ylab = "expression",
            xlab = "",
            cex.axis = 0.001
          )
          axis(2, cex.axis = 0.9)

          cex1 <- ifelse(ny >= 8, 0.65, 0.8)
          title(g, cex.main = ifelse(nchar(g) > 20, 0.85, 1))
          nchar(y)
          too.big <- (max(nchar(y)) >= 8 && ny == 2) ||
            (max(nchar(y)) >= 5 && ny %in% c(3, 4) ||
              ny >= 5)
          if (too.big) {
            dy <- min(gx) - 0.12 * diff(range(gx))
            text(1:ny, dy, levels(factor(y)),
              xpd = NA,
              cex = cex1, srt = 30, adj = 1
            )
          } else {
            mtext(levels(factor(y)),
              side = 1, line = 0.7,
              cex = cex1 * 0.85, las = 1, at = 1:ny
            )
          }
        }
      }

      plot_data_csv <- function() {
        pdata <- plot_data()
        df <- rbind(pdata$y, pdata$X)
        return(df)
      }

      PlotModuleServer(
        "plot",
        plotlib = "base", # does not use plotly
        func = plot.RENDER,
        func2 = plot.RENDER, # no separate modal plot render
        csvFunc = plot_data_csv,
        res = c(90, 180),
        pdf.width = 10, pdf.height = 5.5,
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
