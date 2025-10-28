##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

connectivity_table_foldchange_ui <- function(
  id,
  title,
  info.text,
  caption,
  width,
  height,
  label = ""
) {
  ns <- shiny::NS(id)

  bslib::layout_columns(
    col_widths = 12,
    height = "35%",
    TableModuleUI(
      id = ns("table"),
      info.text = info.text,
      width = width,
      caption = caption,
      height = height,
      title = title,
      label = label
    )
  )
}

connectivity_table_foldchange_server <- function(id,
                                                 pgx,
                                                 getConnectivityScores,
                                                 columns,
                                                 getProfiles,
                                                 getConnectivityMatrix,
                                                 sigdb,
                                                 height) {
  moduleServer(id, function(input, output, session) {
    get_table <- reactive({
      F <- getProfiles()


      F <- F[order(rownames(F)), ]

      S <- getConnectivityScores()
      S1 <- S[match(colnames(F), S$pathway), c("score", "rho")]

      df <- data.frame(signature = colnames(F), S1, t(F))
      df
    })

    foldchangeTable.RENDER <- function() {
      df <- get_table()

      ## pathway is actually signature name
      df$signature <- playbase::shortstring(df$signature, 100)
      score.col <- which(colnames(df) == "score")
      numcols <- setdiff(colnames(df), grep("signature|dataset|pathway|dataset|contrast", colnames(df), value = TRUE))

      DT::datatable(df,
        rownames = FALSE,
        class = "compact cell-border stripe hover",
        extensions = c("Scroller"),
        selection = list(mode = "single", target = "row", selected = 1),
        plugins = "scrollResize",
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip",
          # dom = "lrtip",
          pageLength = 99999,
          scrollX = TRUE,
          scrollY = height,
          scrollResize = TRUE,
          scroller = TRUE,
          deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numcols, 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle("score",
          background = color_from_middle(df[, "score"], "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    }

    foldchangeTable.RENDER_modal <- function() {
      dt <- foldchangeTable.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    }

    foldchangeTable <- TableModuleServer(
      "table",
      func = foldchangeTable.RENDER,
      func2 = foldchangeTable.RENDER_modal,
      selector = "single"
    )

    return(foldchangeTable)
  })
}
