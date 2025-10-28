##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


drugconnectivity_table_cmap_ui <- function(
  id,
  title,
  info.text,
  caption,
  width,
  height
) {
  ns <- shiny::NS(id)

  TableModuleUI(
    ns("datasets"),
    info.text = info.text,
    width = width,
    height = height,
    title = title,
    caption = caption,
    label = "b"
  )
}


drugconnectivity_table_cmap_server <- function(id,
                                               getActiveDSEA) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    table_data <- shiny::reactive({
      dsea <- getActiveDSEA()
      shiny::req(dsea)

      dt <- dsea$table
      return(dt)
    })

    table.RENDER <- function() {
      res <- table_data()
      res$moa <- playbase::shortstring(res$moa, 30)
      res$target <- playbase::shortstring(res$target, 80)
      res$drug <- playbase::shortstring(res$drug, 60)
      res$pval <- NULL
      res$padj <- NULL

      colnames(res) <- sub("moa", "MOA", colnames(res))
      DT::datatable(
        res,
        rownames = FALSE,
        class = "compact cell-border stripe hover",
        extensions = c("Scroller"),
        plugins = "scrollResize",
        selection = list(mode = "single", target = "row", selected = NULL),
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip",
          scrollX = TRUE,
          scrollY = 240, ## card is 380
          scrollResize = TRUE,
          scroller = TRUE,
          deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle("NES",
          background = color_from_middle(res[, "NES"], "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    }

    table.RENDER_modal <- shiny::reactive({
      dt <- table.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    cmap_table <- TableModuleServer(
      "datasets",
      func = table.RENDER,
      func2 = table.RENDER_modal,
      selector = "single"
    )

    return(cmap_table)
  }) ## end of moduleServer
} ## end of server
