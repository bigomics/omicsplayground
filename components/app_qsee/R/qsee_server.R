##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## Was: board.upload/upload_module_batchcorrect.R
if(0) {
  pgx <- playdata::GEIGER_PGX
}

qsee_server <- function(id, pgx) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      
      message("[qsee_server:moduleServer] called ")
      message("[qsee_server:moduleServer] called ")
      pgx <- playdata::GEIGER_PGX
      X0 <- pgx$X
      Y0 <- pgx$samples[,-1]

      ## Holds user-uploaded matrices from the upload modal (NULL = use defaults)
      uploaded <- shiny::reactiveValues(X = NULL, Y = NULL)

      shiny::observeEvent(input$upload, {
        shiny::showModal(shiny::modalDialog(
          title = "Upload CSV files",
          shiny::fileInput(
            session$ns("counts_csv"),
            "Counts (counts.csv)",
            accept = c(".csv", "text/csv", "text/comma-separated-values")
          ),
          shiny::fileInput(
            session$ns("samples_csv"),
            "Samples (samples.csv)",
            accept = c(".csv", "text/csv", "text/comma-separated-values")
          ),
          easyClose = TRUE,
          footer = shiny::modalButton("Close")
        ))
      })

      shiny::observeEvent(input$counts_csv, {
        shiny::req(input$counts_csv)
        df <- playbase::read_counts(input$counts_csv$datapath)
        uploaded$X <- as.matrix(df)
      })

      shiny::observeEvent(input$samples_csv, {
        shiny::req(input$samples_csv)
        df <- playbase::read_samples(input$samples_csv$datapath)        
        uploaded$Y <- df
      })
      
      ## Clean expression matrix. Noise injection for the Normalization
      ## board is controlled there via the "Noise amount" slider.
      getX <- reactive({
        message("[qsee_server] getX ")
        X <- if (!is.null(uploaded$X)) uploaded$X else X0
        X[which(X == 0)] <- NA
        if(!playbase::is_logged(X)) X <- log2(X)
        X
      })

      getY <- reactive({
        if (!is.null(uploaded$Y)) uploaded$Y else Y0
      })
      
      qsee_normalization_server("normalize", rX=getX, rY=getY)
      qsee_imputation_server("impute", rX=getX, rY=getY)
      qsee_bsee_server("bsee", rX=getX, rY=getY)
      qsee_outlier_server("outlier", rX=getX, rY=getY)
      qsee_filtering_server("filtering", rX=getX, rY=getY)            
      
    } ## end-of-server
  )
}
