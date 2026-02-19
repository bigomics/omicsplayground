##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

preservationWGCNA_plot_moduletrait_ui <- function(
  id,
  title = "",
  info.text = "",
  caption = "",
  label = "",
  height = 400,
  width = 400
) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::checkboxInput(ns("weighted"), "Weight by preservation", FALSE),
    shiny::selectInput(
      inputId = ns("sortby"),
      label = "Sort modules by",
      choices = c("clust", "name", "zsummary"),
      selected = "clust"
    ),
    shiny::checkboxInput(ns("largemargin"), "Increase margin", FALSE),
  )

  PlotModuleUI(
    ns("heatmap"),
    title = title,
    label = label,
    info.text = info.text,
    options = options,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

preservationWGCNA_plot_moduletrait_barplot_ui <- function(
  id,
  title = "",
  info.text = "",
  caption = "",
  label = "",
  height = 400,
  width = 400
) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::checkboxInput(ns("colored"), "Colored", FALSE)
  )

  PlotModuleUI(
    ns("barplot"),
    title = title,
    label = label,
    info.text = info.text,
    options = options,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}


preservationWGCNA_plot_moduletrait_server <- function(id,
                                                      rwgcna,
                                                      rtrait) {
  moduleServer(id, function(input, output, session) {
    ## ----------------------------------------------------
    ## Heatmap
    ## ----------------------------------------------------
    heatmap.RENDER <- function() {
      res <- rwgcna()
      shiny::req(res)

      subplots <- c("zsummary", "consmt")
      if (input$weighted) subplots <- c("zsummary", "wt.consmt")

      par(mfrow = c(1, 2), mar = c(4, 9, 3, 2))
      if (input$largemargin) {
        par(mar = c(10, 9, 3, 2))
      }

      playbase::wgcna.plotPreservationModuleTraits(
        res,
        subplots = subplots,
        order.by = input$sortby,
        setpar = FALSE,
        rm.na = TRUE
      )
    }

    PlotModuleServer(
      "heatmap",
      func = heatmap.RENDER,
      pdf.width = 8,
      pdf.height = 12,
      res = c(72, 90),
      add.watermark = FALSE
    )

    ## ----------------------------------------------------
    ## Barplot
    ## ----------------------------------------------------
    barplot.RENDER <- function() {
      res <- rwgcna()
      trait <- rtrait()
      shiny::req(res)
      shiny::req(trait)

      ## Traits that are constant within a layer get NaN correlations in
      ## modTraits (cor() of constant vector). barplot() then receives all
      ## non-finite values and crashes with 'need finite ylim values'.
      ## Also guard against trait being absent from some layers entirely.
      m1 <- tryCatch(
        sapply(res$modTraits, function(x) {
          if (!trait %in% colnames(x)) return(rep(NA_real_, nrow(x)))
          x[, trait]
        }),
        error = function(e) NULL
      )
      shiny::validate(shiny::need(
        !is.null(m1) && any(is.finite(m1)),
        "Selected trait has no finite correlation values across condition groups. Please select a different trait."
      ))

      par(mfrow = c(1, 1), mar = c(8, 4, 2.5, 1))
      playbase::wgcna.plotTraitCorrelationBarPlots(
        res,
        trait = trait,
        multi = TRUE,
        colored = input$colored,
        beside = TRUE,
        setpar = FALSE
      )
    }

    PlotModuleServer(
      "barplot",
      func = barplot.RENDER,
      pdf.width = 10,
      pdf.height = 6,
      res = c(72, 100),
      add.watermark = FALSE
    )
  })
}
