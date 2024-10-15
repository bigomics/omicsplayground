##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_table_genes_ui <- function(
    id,
    label,
    title,
    info.text,
    caption,
    height,
    width) {
  ns <- shiny::NS(id)

  TableModuleUI(
    ns("datasets"),
    info.text = info.text,
    width = width,
    height = height,
    title = title,
    caption = caption,
    label = label
  )
}

wgcna_table_genes_server <- function(id,
                                     wgcna.compute,
                                     selected_module) {
  moduleServer(id, function(input, output, session) {
    geneTable.RENDER <- shiny::reactive({
      out <- wgcna.compute()

      k <- selected_module()
      genes <- out$me.genes[[k]]
      shiny::req(genes)
      tt <- playdata::GENE_TITLE[toupper(genes)]
      rho <- cor(out$datExpr[, genes], out$net$MEs[, k])[, 1]

      df <- data.frame(module = k, gene = genes, me.rho = rho, title = tt)
      numeric.cols <- grep("score|value|ratio|rho", colnames(df))

      DT::datatable(
        df,
        rownames = FALSE, #
        #
        extensions = c("Buttons", "Scroller"),
        selection = list(mode = "single", target = "row", selected = NULL),
        class = "compact cell-border stripe hover",
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip", #
          ## pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
          scrollX = TRUE, #
          #
          scrollY = "70vh",
          scroller = TRUE, deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    })

    geneTable.RENDER_modal <- shiny::reactive({
      dt <- geneTable.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    geneTable_module <- TableModuleServer(
      "datasets",
      func = geneTable.RENDER,
      func2 = geneTable.RENDER_modal,
      selector = "none"
    )

    return(geneTable_module)
  })
}
