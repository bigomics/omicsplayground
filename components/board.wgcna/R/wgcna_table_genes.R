##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_table_genes_ui <- function(
    id,
    label = "",
    title = "",
    info.text = "",
    caption = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)

  options <- tagList(
    shiny::checkboxInput(ns("showpvalues"),"Show p-values", FALSE),
    shiny::checkboxInput(ns("showall"),"Show all modules", FALSE)
  )

  TableModuleUI(
    ns("datasets"),
    options = options,
    info.text = info.text,
    width = width,
    height = height,
    title = title,
    caption = caption,
    label = label
  )
}

wgcna_table_genes_server <- function(id,
                                     wgcna,
                                     pgx,
                                     selected_module,
                                     selected_trait) {
  moduleServer(id, function(input, output, session) {
    RENDER <- function(full = FALSE) {
      res <- wgcna()
      module <- selected_module()
      trait <- selected_trait()

      shiny::req(pgx$genes)
      shiny::req(res)
      shiny::req(module, trait)
      shiny::req(module != "" && trait != "")

      df <- playbase::wgcna.getGeneStats(
        res,
        module = module, trait = trait, plot = FALSE
      )
      symbol <- pgx$genes[rownames(df), "symbol"]
      feature <- rownames(df)
      feature1 <- sub(";.*", ";...", feature) ## take first
      df <- cbind(feature = feature1, symbol = symbol, df)
      if (all(symbol == feature)) {
        df$symbol <- NULL
      }

      if (input$showpvalues == FALSE) {
        df <- df[, grep("pvalue", colnames(df), invert = TRUE, ignore.case = TRUE), drop = FALSE]
      }

      ## only those in module
      if(!input$showall) {
        df <- df[ which(df$module == module), , drop=FALSE ]
      }

      numeric.cols <- grep("^module$|symbol|feature", colnames(df), invert=TRUE)
      colnames(df) <- sub("moduleMembership","MM",colnames(df))
      colnames(df) <- sub("traitSignificance","TS",colnames(df))
      colnames(df) <- sub("foldChange","logFC",colnames(df))

      # If trait is continuous, change logFC to rho
      # cont trait is in pgx$samples, binary ones are name-converted (not in pgx$samples)
      is.cont <- trait %in% colnames(pgx$samples)
      if (is.cont) {
        colnames(df) <- sub("logFC", "rho", colnames(df))
      }

      DT::datatable(
        df,
        rownames = FALSE,
        extensions = c("Buttons", "Scroller"),
        selection = list(mode = "single", target = "row", selected = NULL),
        plugins = "scrollResize",
        class = "compact cell-border stripe hover",
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip", #
          scrollX = TRUE, #
          scrollY = "70vh",
          scrollResize = TRUE,
          scroller = TRUE,
          deferRender = TRUE
        )
      ) %>%
        DT::formatSignif(numeric.cols, 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "10px", lineHeight = "70%") %>%
        DT::formatStyle(
          "score",
          background = color_from_middle(df$score, "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        ) 


    }

    RENDER_modal <- function() {
      dt <- RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    }

    geneTable_module <- TableModuleServer(
      "datasets",
      func = RENDER,
      func2 = RENDER_modal,
      selector = "none"
    )

    return(geneTable_module)
  })
}
