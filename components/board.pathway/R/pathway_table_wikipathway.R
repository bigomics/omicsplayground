##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


functional_table_wikipathway_ui <- function(
  id,
  title,
  info.text,
  caption,
  label,
  width,
  height
) {
  ns <- shiny::NS(id)

  TableModuleUI(
    ns("tablemodule"),
    info.text = info.text,
    width = width,
    caption = caption,
    height = height,
    title = title,
    label = label
  )
}

functional_table_wikipathway_server <- function(id,
                                                pgx,
                                                getFilteredWikipathwayTable,
                                                fa_contrast) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    table_data <- shiny::reactive({
      shiny::req(pgx$X)
      res <- list(
        pgx = pgx,
        df = getFilteredWikipathwayTable(),
        fa_contrast = fa_contrast()
      )
      return(res)
    })

    table_RENDER <- function() {
      res <- table_data()
      pgx <- res$pgx
      df <- res$df
      comparison <- res$fa_contrast

      if (is.null(comparison)) {
        return(NULL)
      }
      if (is.null(df)) {
        return(NULL)
      }
      if (nrow(df) == 0) {
        return(NULL)
      }

      ## add hyperlink
      url <- paste0("https://www.wikipathways.org/pathways/", df$pathway.id, ".html")
      pathway.id_link <- paste0(
        "<a href='", url, "' target='_blank'>",
        rep_len("<i class='fa-solid fa-arrow-up-right-from-square weblink'></i>", nrow(df)),
        "</a>"
      )
      df$pathway <- paste(df$pathway, pathway.id_link)
      df$pathway.id <- NULL

      numeric.cols <- colnames(df)[which(sapply(df, is.numeric))]

      DT::datatable(df,
        rownames = FALSE,
        escape = c(-1, -2),
        extensions = c("Scroller"),
        selection = list(
          mode = "single",
          target = "row",
          selected = 1
        ),
        fillContainer = TRUE,
        plugins = "scrollResize", ## resizes scrollable area
        options = list(
          dom = "lfrtip",
          scrollX = FALSE,
          scrollY = 800,
          scrollResize = TRUE,
          scroller = TRUE,
          deferRender = TRUE,
          autoWidth = TRUE,
          columnDefs = list(list(
            targets = 2, ## with no rownames column 1 is column 2
            render = DT::JS(
              "function(data, type, row, meta) {",
              "return type === 'display' && data.length > 50 ?",
              "'<span title=\"' + data + '\">' + data.substr(0, 50) + '...</span>' : data;",
              "}"
            )
          ))
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 4) %>%
        DT::formatStyle(
          0,
          target = "row",
          fontSize = "11px",
          lineHeight = "70%"
        ) %>%
        DT::formatStyle("logFC",
          background = color_from_middle(
            df[, "logFC"],
            "lightblue",
            "#f5aeae"
          ),
          backgroundSize = "98% 88%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    }

    table_RENDER_modal <- shiny::reactive({
      dt <- table_RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    my_table <- TableModuleServer(
      "tablemodule",
      func = table_RENDER,
      func2 = table_RENDER_modal,
      csvFunc = function() {
        table_data()$df
      },
      selector = "single"
    )

    return(my_table)
  }) ## end of moduleServer
} ## end of server
