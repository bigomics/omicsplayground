##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

connectivity_table_similarity_scores_ui <- function(
  id,
  title,
  info.text,
  caption,
  width,
  height,
  label = ""
) {
  ns <- shiny::NS(id)

  TableModuleUI(
    id = ns("scores"),
    info.text = info.text,
    width = width,
    caption = caption,
    height = height,
    title = title,
    label = label
  )
}

connectivity_table_similarity_scores_server <- function(id,
                                                        getConnectivityScores,
                                                        columns,
                                                        height) {
  moduleServer(id, function(input, output, session) {
    get_datatable <- shiny::reactive({
      df <- getConnectivityScores()
      shiny::req(df)

      kk <- intersect(columns, colnames(df))
      df <- df[, kk]
      df <- df[abs(df$score) > 0, , drop = FALSE]

      df$pathway <- playbase::shortstring(df$pathway, 100)
      colnames(df) <- sub("pathway", "dataset/contrast", colnames(df))
      df
    })

    createDT <- function(df) {
      score.col <- which(colnames(df) == "score")
      numcols <- c("score", "pval", "padj", "NES.q", "odd.ratio", "ES", "NES", "rho", "R2", "tau")
      numcols <- intersect(numcols, colnames(df))

      feature_link <- playbase::wrapHyperLink(
        rep_len("<i class='fa-solid fa-arrow-up-right-from-square weblink'></i>", nrow(df)),
        df$`dataset/contrast`
      ) |> HandleNoLinkFound(
        NoLinkString = "<i class='fa-solid fa-arrow-up-right-from-square weblink'></i>",
        SubstituteString = "<i class='fa-solid fa-arrow-up-right-from-square blank_icon'></i>"
      )
      df[["dataset/contrast"]] <- paste(df[["dataset/contrast"]], "&nbsp;", feature_link)

      DT::datatable(df,
        rownames = FALSE,
        escape = c(-1),
        class = "compact cell-border stripe hover",
        extensions = c("Scroller"),
        selection = list(mode = "single", target = "row", selected = 1),
        plugins = "scrollResize",
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip",
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
          background = color_from_middle(
            df[, "score"], "lightblue", "#f5aeae"
          ),
          backgroundSize = "98% 88%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    }

    connectivityScoreTable.RENDER <- function() {
      df <- get_datatable()
      df <- df[, c("dataset/contrast", "score", "rho")]
      createDT(df)
    }

    connectivityScoreTable.RENDER_modal <- function() {
      df <- get_datatable()
      dt <- createDT(df)
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    }

    table_csv <- function() {
      df <- get_datatable()
      return(df)
    }

    connectivityScoreTable <- TableModuleServer(
      "scores",
      func = connectivityScoreTable.RENDER,
      func2 = connectivityScoreTable.RENDER_modal,
      csvFunc = table_csv,
      selector = "single"
    )

    return(connectivityScoreTable)
  })
}
