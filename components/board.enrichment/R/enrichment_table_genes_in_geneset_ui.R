##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

enrichment_table_genes_in_geneset_ui <- function(id, width, height) {
  ns <- shiny::NS(id)

  info_text <- "By clicking on a gene set in the table <code>I</code>, it is possible to see the gene list of that gene set in this table. By clicking on a gene in this table, users can check the expression status of the gene for the selected contrast in the <code>Expression</code> barplot and its correlation to the gene set in the <code>Gene to gene set correlation</code> scatter plot under the <code>Plots</code> section."

  TableModuleUI(
    ns("datasets"),
    info.text = info_text,
    width = width,
    height = height,
    title = "Genes in gene set",
    label = "II"
  )

}

enrichment_table_genes_in_geneset_server <- function(id,
                                                     geneDetails) {
  moduleServer(id, function(input, output, session) {
    genetable.RENDER <- shiny::reactive({
      rpt <- geneDetails()
      if (is.null(rpt) || nrow(rpt) == 0) {
        shiny::validate(shiny::need(nrow(rpt) > 0, "warning. no genes."))
        return(NULL)
      }

      rpt$gene_title <- NULL
      if (!is.null(rpt) && nrow(rpt) > 0) {
        jj <- which(sapply(rpt, is.numeric))
        rpt[, jj] <- round(rpt[, jj], digits = 4)
        jj <- which(sapply(rpt, is.character) | sapply(rpt, is.factor))
        if (length(jj) > 0) rpt[, jj] <- apply(rpt[, jj, drop = FALSE], 2, shortstring, 60)
      } else {
        rpt <- data.frame("", 0, 0, 0)[0, ]
        colnames(rpt) <- c("gene_name", "fc", "p", "q")
      }

      colnames(rpt) <- sub("^GS$", "gene set", colnames(rpt))
      numeric.cols <- which(sapply(rpt, is.numeric))
      numeric.cols

      tbl <- DT::datatable(rpt,
        class = "compact cell-border stripe", rownames = FALSE,
        extensions = c("Scroller"),
        selection = list(mode = "single", target = "row", selected = 1),
        fillContainer = TRUE,
        options = list(
          dom = "frtip",
          #paging = TRUE,
          #pageLength = 15, ##  lengthMenu = c(20, 30, 40, 60, 100, 250),
          scrollX = TRUE,
          scrollY = 700,
          scroller = TRUE,
          deferRender = TRUE,
          search = list(
            regex = TRUE,
            caseInsensitive = TRUE
          )
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 4)

      if (nrow(rpt) > 0 && ("fc" %in% colnames(rpt))) {
        fx <- rpt[, "fc"]
        tbl <- tbl %>%
          DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
          DT::formatStyle("fc", background = color_from_middle(fx, "lightblue", "#f5aeae"))
      }
      tbl
    })

    genetable <- TableModuleServer(
      "datasets",
      func = genetable.RENDER,
      selector = "single"
    )

    return(genetable)
  })
}
