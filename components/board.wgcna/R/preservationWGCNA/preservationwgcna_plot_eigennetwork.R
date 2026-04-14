##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

preservationWGCNA_plot_eigenNetwork_ui <- function(
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
    shiny::checkboxInput(
      inputId = ns("top20"),
      label = "Show top 20",
      value = FALSE
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    # options = options,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

preservationWGCNA_plot_eigenNetwork_server <- function(id,
                                                       rwgcna) {
  moduleServer(id, function(input, output, session) {
    plot.RENDER <- function() {
      res <- rwgcna()
      shiny::req(res)

      ## plotEigengeneNetworks calls cor(MEs) then hclust per set.
      ## It crashes with 'invalid dendrogram input' when any set has:
      ##   - fewer than 2 non-grey eigengenes  (hclust requires >= 2 obs)
      ##   - a constant or near-constant eigengene (cor → NaN → bad dist)
      ## Pre-check covers the most common cases; tryCatch handles anything
      ## else (WGCNA internals can fail for subtle reasons in small groups).
      mes_valid <- sapply(res$net$multiMEs, function(me) {
        non_grey <- me$data[, substring(colnames(me$data), 3) != "grey", drop = FALSE]
        if (ncol(non_grey) < 2) return(FALSE)
        all(apply(non_grey, 2, function(x) var(x, na.rm = TRUE) > 1e-10))
      })
      shiny::validate(shiny::need(
        all(mes_valid),
        paste0(
          "Cannot plot eigengene network: one or more condition groups has fewer than 2 ",
          "non-grey modules or has a constant eigengene. ",
          "Try reducing Min. module size or increasing Max. features."
        )
      ))

      par(mar = c(0, 3, 3, 1))
      tryCatch(
        WGCNA::plotEigengeneNetworks(
          res$net$multiMEs,
          setLabels = names(res$net$multiMEs),
          marHeatmap = c(1, 3, 3, 1)
        ),
        error = function(e) {
          shiny::validate(shiny::need(FALSE, paste0(
            "Cannot plot eigengene network: one or more condition groups has ",
            "near-constant eigengenes or an incompatible module structure. ",
            "Try reducing Min. module size or increasing Max. features."
          )))
        }
      )
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 8,
      pdf.height = 10,
      res = c(65, 90),
      add.watermark = FALSE
    )
  })
}
