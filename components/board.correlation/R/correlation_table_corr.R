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
  width) {
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
      if (is.null(R)) return(NULL)

      P <- getPartialCorrelation()
      pcor <- P[match(rownames(R), rownames(P)), "pcor"]

      title <- pgx$genes[rownames(R),"gene_title"]
      title <- substring(title, 1, 80)
      df <- data.frame(gene=rownames(R), title=title, cor=R[,"cor"], pcor=pcor)
      
      return(df)
    })

    ### TABLE
    cor_table.RENDER <- shiny::reactive({
      shiny::req(pgx)

      df <- plot_data()

      numeric.cols <- colnames(df)[3:ncol(df)]
      ## selectmode <- ifelse(input$corGSEAtable_multiselect,'multiple','single')

      DT::datatable(
        df,
        rownames = FALSE, ## escape = c(-1),
        extensions = c("Buttons", "Scroller"),
        selection = list(mode = "single", target = "row", selected = 1),
        class = "compact cell-border stripe hover",
        plugins = "scrollResize",
        fillContainer = TRUE,
        options = list(
          dom = "lfrti",
          ## pageLength = 20,
          ##  lengthMenu = c(20, 30, 40, 60, 100, 250),
          ##paging = FALSE,
          scrollX = FALSE,
          scrollY = 100,
          scrollResize = TRUE,
          scroller = TRUE,
          deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 4) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle("cor",
          background = playbase::color_from_middle(
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
