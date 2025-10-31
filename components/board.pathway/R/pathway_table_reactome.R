##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


functional_table_reactome_ui <- function(
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

functional_table_reactome_server <- function(id,
                                             getFilteredReactomeTable,
                                             fa_contrast,
                                             scrollY) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    table_data <- shiny::reactive({
      df <- getFilteredReactomeTable()
      shiny::req(df, fa_contrast())
      res <- list(
        df = df,
        fa_contrast = fa_contrast
      )
      return(res)
    })

    table_RENDER <- function() {
      res <- table_data()
      df <- res$df
      if (nrow(df) == 0) {
        return(NULL)
      }
      comparison <- res$fa_contrast

      ## add hyperlink
      url <- paste0("https://reactome.org/content/detail/", df$reactome.id)
      reactome.id_link <- paste0(
        "<a href='", url, "' target='_blank'>",
        rep_len("<i class='fa-solid fa-arrow-up-right-from-square weblink'></i>", nrow(df)),
        "</a>"
      )
      df$pathway <- paste(df$pathway, reactome.id_link)
      df$reactome.id <- NULL

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
        plugins = "scrollResize",
        options = list(
          dom = "lfrtip",

          ## dom = "ft",
          scrollX = FALSE,
          scrollY = scrollY,
          scrollResize = TRUE,
          scroller = TRUE,
          deferRender = TRUE,
          autoWidth = TRUE,
          columnDefs = list(list(
            targets = 1, ## with no rownames column 1 is column 2
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

    table_RENDER_modal <- function() {
      dt <- table_RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    }

    my_table <- TableModuleServer(
      "tablemodule",
      func = table_RENDER,
      func2 = table_RENDER_modal,
      csvFunc = table_data,
      selector = "single"
    )

    return(my_table)
  }) ## end of moduleServer
} ## end of server
