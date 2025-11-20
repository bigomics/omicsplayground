##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' @export
correlation_table_corr_ui <- function(
  id,
  title,
  info.text,
  caption,
  label = "",
  height,
  width
) {
  ns <- shiny::NS(id)

  TableModuleUI(
    ns("table"),
    info.text = info.text,
    height = height,
    width = width,
    title = title,
    caption = caption,
    label = label
  )
}


#' Expression plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
correlation_table_corr_server <- function(id,
                                          getPartialCorrelation,
                                          getGeneCorr,
                                          pgx,
                                          watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    # reactive function listeninng for changes in input
    plot_data <- shiny::reactive({
      R <- getGeneCorr()
      if (is.null(R)) {
        return(NULL)
      }

      P <- getPartialCorrelation()
      pcor <- P[match(rownames(R), rownames(P)), "pcor"]
      gene_table <- pgx$genes
      if (all(gene_table$human_ortholog == rownames(gene_table)) | all(is.na(gene_table$human_ortholog))) {
        gene_table_cols <- c("feature", "symbol", "gene_title")
      } else {
        gene_table_cols <- c("feature", "symbol", "human_ortholog", "gene_title")
      }

      tt <- gene_table[rownames(R), gene_table_cols]
      df <- data.frame(tt, cor = R[, "cor"], pcor = pcor)

      return(df)
    })

    ### TABLE
    cor_table.RENDER <- shiny::reactive({
      shiny::req(pgx$X)

      df <- plot_data()

      char_cols <- c("feature", "gene", "symbol", "human_ortholog", "gene_title")
      ## if (all(pgx$organism %in% c("Human", "human"))) {
      if (any(pgx$organism %in% c("Human", "human"))) {
        char_cols <- c("feature", "gene", "symbol", "human_ortholog", "gene_title")
      }
      if (sum(df$feature %in% df$symbol) > nrow(df) * .8) {
        df$feature <- NULL
      }

      if ("human_ortholog" %in% colnames(df)) df$human_ortholog <- NULL
      numeric.cols <- which(!colnames(df) %in% char_cols)

      DT::datatable(
        df,
        rownames = FALSE, #
        extensions = c("Buttons", "Scroller"),
        selection = list(mode = "single", target = "row", selected = 1),
        class = "compact cell-border stripe hover",
        plugins = "scrollResize",
        fillContainer = TRUE,
        options = list(
          dom = "lfrti",
          scrollX = TRUE,
          scrollY = 100,
          scrollResize = TRUE,
          scroller = TRUE,
          deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 4) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle("cor",
          background = color_from_middle(
            df[, "cor"], "lightblue", "#f5aeae"
          ),
          backgroundSize = "98% 88%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    })

    cor_table.RENDER_modal <- shiny::reactive({
      dt <- cor_table.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    TableModuleServer(
      "table",
      func = cor_table.RENDER,
      func2 = cor_table.RENDER_modal,
      selector = "none"
    )
  }) ## end of moduleServer
}
