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
expression_table_gsettable_ui <- function(
    id,
    title,
    caption,
    info.text,
    width,
    height) {
  ns <- shiny::NS(id)

  TableModuleUI(
    ns("datasets"),
    info.text = info.text,
    width = width,
    height = height,
    caption = caption,
    title = title,
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
                                              scrollY,
                                              watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    gsettable.RENDER <- shiny::reactive({
      df <- gx_related_genesets()

      shiny::validate(shiny::need(
        !is.null(df),
        "Please select a gene in the table."
      ))

      external_links <- playbase::wrapHyperLink(
        rep_len("<i class='fa-solid fa-circle-info'></i>", nrow(df)),
        rownames(df)
      ) |> HandleNoLinkFound(
        NoLinkString = "<i class='fa-solid fa-circle-info'></i>",
        SubstituteString = "<i class='fa-solid fa-circle-info icon_container'></i><i class='fa fa-ban icon_nested'></i>"
      )
      DT::datatable(df,
        rownames = external_links,
        escape = c(-1, -2),
        extensions = c("Scroller"),
        plugins = "scrollResize",
        fillContainer = TRUE,
        options = list(
          dom = "frtip",
          scrollX = TRUE,
          scrollY = scrollY,
          scrollResize = TRUE,
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
      csvFunc = gx_related_genesets,
      selector = "single"
    )

    return(gsettable)
  }) # end module server
} # end server
