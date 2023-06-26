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
  label="") {
  ns <- shiny::NS(id)

  bslib::layout_column_wrap(
    width = 1,
    height = "35%",
    ##style = htmltools::css(grid_template_columns = "5fr 1fr"),
    TableModuleUI(
      id = ns("datasets"),
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
                                                  getTopProfiles,
                                                  getConnectivityMatrix,
                                                  sigdb,
                                                  height) {
  moduleServer(id, function(input, output, session) {

    get_table <- reactive({

      F <- getTopProfiles()
      F <- F[order(-rowMeans(F**2)), , drop = FALSE]
      F <- head(F,50)
      S <- getConnectivityScores()
      S1 <- S[match(colnames(F), S$pathway),c("score","rho")]

      df <- data.frame(signature=colnames(F), S1, t(F))
      df
    })
    
    connectivityScoreTable.RENDER <- function() {

      df <- get_table()
      
      ## pathway is actually signature name
      df$signature <- playbase::shortstring(df$signature, 100)
      score.col <- which(colnames(df) == "score")
      numcols <- setdiff(colnames(df),grep("signature|dataset|pathway|dataset|contrast",colnames(df),value=TRUE))

      DT::datatable(df,
        rownames = FALSE,
        class = "compact cell-border stripe hover",
        extensions = c("Scroller"),
        selection = list(mode = "single", target = "row", selected = 1),
        plugins = 'scrollResize',
        fillContainer = TRUE,
        options = list(
            #          dom = "lfrtip",
          dom = "lrtip",            
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
          background = playbase::color_from_middle(df[, "score"], "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    }

    connectivityScoreTable.RENDER_modal <- function() {
      dt <- connectivityScoreTable.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    }

    connectivityScoreTable <- TableModuleServer(
      "datasets",
      func = connectivityScoreTable.RENDER,
      func2 = connectivityScoreTable.RENDER_modal,
      selector = "single"
    )

    return(connectivityScoreTable)
  })
}
