##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

connectivity_table_similarity2_ui <- function(
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
    style = htmltools::css(grid_template_columns = "5fr 1fr"),
    TableModuleUI(
      id = ns("datasets"),
      info.text = info.text,
      width = width,
      caption = caption,
      height = height,
      title = title,
      label = label
    ),
    shiny::wellPanel(    
      shiny::radioButtons(
          inputId = ns("select_genes:"),
          label = "Select genes",
          choices = c("50","200","<custom>"),
          inline = TRUE,
      ),             
      shiny::conditionalPanel(
          "input.select_genes == '<custom>'",
          ns = ns,
          selectizeInput(
            inputId = ns("genes"),
            label = "Select genes:",
            choices = NULL,
            multiple = TRUE
          )
      ),
      selectizeInput(
        inputId = ns("datasets"),
        label = "Select datasets:",
        choices = NULL,
        multiple = TRUE
      )
    )
  )
  
}

connectivity_table_similarity2_server <- function(id,
                                                  pgx,
                                                  getConnectivityScores,
                                                  columns,
                                                  getConnectivityMatrix,
                                                  sigdb,
                                                  height) {
  moduleServer(id, function(input, output, session) {

    shiny::observe({
      shiny::req(pgx$X)
      genes <- rownames(pgx$X)
      updateSelectizeInput(session, "genes", choices=genes,
        selected=genes[1:3], server=TRUE)
    })

    shiny::observe({
      df <- getConnectivityScores()
      shiny::req(df)
      datasets <- sort(unique(gsub("^\\[|\\].*","",rownames(df))))
      datasets <- c("<all>",datasets)
      updateSelectizeInput(session, "datasets", choices=datasets,
        selected="<all>", server=TRUE)
    })

    get_table <- reactive({
      df <- getConnectivityScores()
      shiny::req(df)
      shiny::req(input$genes)
      shiny::req(input$datasets)      

      ##kk <- c("pathway", "score", "rho", "NES", "padj", "leadingEdge")
      kk <- intersect(columns, colnames(df))
      df <- df[, kk]
      df <- df[abs(df$score) > 0, , drop = FALSE]
      dbg("[connectivity_table_similarity2.R] dim(df) = ",dim(df))
      
      ## add genes
      genes <- input$genes
      dbg("[connectivity_table_similarity2.R] length(genes) = ",length(genes))
      datasets <- input$datasets
      dbg("[connectivity_table_similarity2.R] length(datasets) = ",length(datasets))
      
      if(length(genes)>0) {
        select <- df$pathway
        if("<all>" %in% datasets) {
          select <- df$pathway
        } else {
          pw.dataset <- gsub("^\\[|\\].*","",df$pathway)
          sel <- which(pw.dataset %in% datasets)
          select <- df$pathway[sel]
          df <- df[sel,]
        }
        M <- getConnectivityMatrix( sigdb(), select=select, genes=genes)
        dbg("[connectivity_table_similarity2.R] dim(M) = ",dim(M))
        df <- cbind(df, t(M))
      }
      
      df
    })
    
    connectivityScoreTable.RENDER <- function() {

      df <- get_table()
      
      ## pathway is actually signature name
      df$pathway <- playbase::shortstring(df$pathway, 100)

      colnames(df) <- sub("pathway", "dataset/contrast", colnames(df))
      score.col <- which(colnames(df) == "score")
      numcols <- c("score", "pval", "padj", "NES.q", "ES", "NES", "rho", "R2")
      numcols <- setdiff(colnames(df),grep("pathway|dataset|contrast",colnames(df),value=TRUE))

      DT::datatable(df,
        rownames = FALSE,
        class = "compact cell-border stripe hover",
        extensions = c("Scroller"),
        selection = list(mode = "single", target = "row", selected = 1),
        plugins = 'scrollResize',
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
          background = playbase::color_from_middle(
            df[, "score"], "lightblue", "#f5aeae"
          ),
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
