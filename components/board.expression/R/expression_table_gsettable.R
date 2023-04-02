##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' UI code for table code: expression board
#'
#' @param id
#' @param label
#' @param height
#' @param width
#'
#' @export
expression_table_gsettable_ui <- function(id, width, height) {
  ns <- shiny::NS(id)

  gsettable_text <- "By clicking on a gene in the Table <code>I</code>, it is possible to see which genesets contain that gene in this table, and check the differential expression status in other comparisons from the <code>Gene in contrasts</code> plot under the <code>Plots</code> tab."

  TableModuleUI(
    ns("datasets"),
    info.text = gsettable_text,
    width = width,
    height = height,
    title = "Gene sets with gene",
    label = "II"
  )
}

#' Server side table code: expression board
#'
#' @param id
#' @param watermark
#'
#' @export
expression_table_gsettable_server <- function(id,
                                              gx_related_genesets,
                                              height,
                                              width,
                                              watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    gsettable.RENDER <- shiny::reactive({
      df <- gx_related_genesets()

      ##req(df)
      shiny::validate(shiny::need(!is.null(df),
        "Please select a gene in the table."))

      df$geneset <- wrapHyperLink(df$geneset, rownames(df))

      DT::datatable(df,
#        class = "compact",  ## not good!
        rownames = FALSE, escape = c(-1, -2),
        extensions = c("Scroller"),
        fillContainer = TRUE,
        options = list(
          dom = "frtip",
          scrollX = TRUE,
          scrollY = 240,
          scroller = TRUE,
          deferRender = TRUE,
          search = list(
            regex = TRUE,
            caseInsensitive = TRUE
          )
        ), ## end of options.list
        selection = list(mode = "single", target = "row", selected = NULL)
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle("fx", background = color_from_middle(df$fx, "lightblue", "#f5aeae"))
    })

    gsettable.RENDER_modal <- shiny::reactive({
      dt <- gsettable.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    gsettable <- TableModuleServer(
      "datasets",
      func = gsettable.RENDER,
      func2 = gsettable.RENDER_modal,
      selector = "single"
    )

    return(gsettable)
  }) # end module server
} # end server
