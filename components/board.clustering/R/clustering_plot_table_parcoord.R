##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


clustering_table_parcoord_ui <- function(
  id,
  label = "",
  title,
  info.text,
  caption,
  height,
  width
) {
  ns <- shiny::NS(id)

  TableModuleUI(
    ns("datasets"),
    info.text = info.text,
    height = height,
    caption = caption,
    width = width,
    title = title,
    label = "b"
  )
}

clustering_plot_table_parcoord_server <- function(id,
                                                  getTopMatrix,
                                                  pgx,
                                                  watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    parcoord.matrix <- shiny::reactive({
      filt <- getTopMatrix()
      shiny::req(filt)
      zx <- filt$mat
      shiny::validate(shiny::need(
        ncol(zx) > 1,
        "Filter is too restrictive. Please change 'Filter samples:'."
      ))
      zx <- t(scale(t(zx)))
      zx <- round(zx, digits = 4)
      list(mat = zx, clust = filt$idx)
    })

    parcoord_table.RENDER <- function() {
      parcoord <- parcoord.matrix()

      mat <- parcoord$mat
      clust <- parcoord$clust
      df <- data.frame(gene.module = clust, mat, check.names = FALSE)
      row.names(df) <- playbase::probe2symbol(row.names(df), pgx$genes, "gene_name", fill_na = TRUE)
      numeric.cols <- 2:ncol(df)
      DT::datatable(
        df,
        rownames = TRUE,
        extensions = c("Buttons", "Scroller"),
        plugins = "scrollResize",
        selection = list(mode = "single", target = "row", selected = NULL),
        class = "compact hover",
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip",
          scrollX = TRUE,
          scrollY = "23vh",
          scrollResize = TRUE,
          scroller = TRUE,
          deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    }

    parcoord_table.RENDER_modal <- function() {
      dt <- parcoord_table.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    }

    TableModuleServer(
      "datasets",
      func = parcoord_table.RENDER,
      func2 = parcoord_table.RENDER_modal,
      selector = "none"
    )
  })
}
