##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

qsee_server <- function(id, pgx=NULL) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      
      message("[qsee_server] called ")
      
      ## Holds user-uploaded matrices from the upload modal (NULL = use defaults)
      uploaded <- shiny::reactiveValues(
        X = NULL,
        Y = NULL,
        K = NULL
      )

      ## Track if we've already shown the "no dataset" popup (to avoid spamming)
      has_shown_no_data_popup <- shiny::reactiveVal(FALSE)

      shiny::observeEvent( list(pgx$X, pgx$samples), {
        uploaded$X <- pgx$X
        uploaded$Y <- pgx$samples
      })

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

      shiny::observeEvent(input$load_example, {
        example_pgx <- playdata::GEIGER_PGX
        uploaded$X <- example_pgx$X
        uploaded$Y <- example_pgx$samples[,-1]
        message("[qsee_server] loaded example data")
      })

      ## Handle popup buttons
      shiny::observeEvent(input$load_example_from_popup, {
        shiny::removeModal()
        example_pgx <- playdata::GEIGER_PGX
        uploaded$X <- example_pgx$X
        uploaded$Y <- example_pgx$samples[,-1]
        message("[qsee_server] loaded example data from popup")
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
        X <- uploaded$X 
        X[which(X == 0)] <- NA
        if(!playbase::is_logged(X)) X <- log2(X)
        X
      })

      getY <- reactive({
        uploaded$Y
      })
      
      qsee_normalization_server("normalize", rX=getX, rY=getY)
      qsee_imputation_server("impute", rX=getX, rY=getY)
      qsee_bsee_server("bsee", rX=getX, rY=getY)
      qsee_outlier_server("outlier", rX=getX, rY=getY)
      qsee_filtering_server("filtering", rX=getX, rY=getY)


      ## If no dataset is loaded (X is NULL/empty), show a helpful popup only
      ## after the Qsee module itself is visible.  In the main application Qsee
      ## is a hidden parent nav panel during startup.
      #shiny::observe({
      shiny::observeEvent( list(input$nav, input$is_visible), {      

        message("DBG [qsee_server] is_visible = ", input$is_visible)
        message("DBG [qsee_server] input$nav = ", input$nav)        
        
        shiny::req(isTRUE(input$is_visible))
        X <- uploaded$X
        noX <- is.null(X) || length(X) == 0 || nrow(X) == 0 || ncol(X) == 0
        message("DBG [qsee_server] is.null.X = ", is.null(X))
        message("DBG [qsee_server] noX = ", noX)        
        ##message("DBG [qsee_server] has_shown_no_data_popup = ", has_shown_no_data_popup())
        
        ##if (noX && !has_shown_no_data_popup()) {
        if (noX) {        
          #has_shown_no_data_popup(TRUE)
          shiny::showModal(
            shiny::modalDialog(
              title = "No dataset loaded",
              shiny::p(
                "No dataset has been loaded yet. Would you like to load ",
                "the example dataset to explore Qsee/Bsee?"
                ),
              footer = shiny::tagList(
                shiny::modalButton("Cancel"),
                shiny::actionButton(session$ns("load_example_from_popup"),
                    "Load example data", class = "btn-primary")
              ),
              size = "s",
              easyClose = FALSE
            )
          )
        }
      })
      

    } ## end-of-server
  )
}
