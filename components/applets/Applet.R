
runApplet <- function(applet, pgx, as="gadget", size='xl', ...) {

  ## random id
  id <- paste0("applet_",sample(1e8,1)) ## random

  app_inputs <- applet$inputs(id=id)
  if( inherits(app_inputs, "shiny.tag") &&
        inherits(app_inputs$children[[1]], "shiny.tag.list")) {
    ## strip away bigdash::tabSettings
    app_inputs <- app_inputs$children[[1]]
  }

  main_ui <- applet$ui(id=id)
  if(0) {
    ## get ui as list of modules
    main_ui <- applet$ui(id=id, as.taglist=TRUE)
    main_ui <- bslib::layout_columns(
      col_widths = 6,
      height = "calc(100vh - 200px)",
      heights_equal = "row",
      gap = "20px",
      !!!main_ui
    )
  }
  
  app_ui <- bslib::page_sidebar(
    title = NULL,
    sidebar = bslib::sidebar(app_inputs),
    main_ui
  )

  r_pgx <- pgx
  if(class(pgx) != "reactivevalues") {
    r_pgx <- do.call("reactiveValues",pgx)
  }
  
  if(as == "gadget") {    
    app_server <- function(input, output, session) {
      #applet$server(id=id, ...)
      #applet$server(id='app', pgx=pgx)
      applet$server(id=id, pgx = r_pgx)          
    }
    #shinyApp(app_ui, app_server)
    runGadget(app_ui, app_server)
  }

  if(as == "modal") {

    sizestyle <- switch( tolower(size),
      "m" = "height: calc(65vh - 180px); width: 50vw;",
      "l" = "height: calc(85vh - 180px); width: 75vw;",     
      "xl" = "height: calc(100vh - 180px); width: 95vw;"
    )

    ## run server
    applet$server(id=id, pgx=r_pgx, ... )
   
    ## Show applet as modal
    showModal(modalDialog2(
      title = applet$title,
#      style = sizestyle,
      size = size.
      easyClose = TRUE, 
      app_ui,
      tags$style(
        type = 'text/css',
        '.modal-dialog { width: fit-content !important; }',
        '.modal-content { width: 100vw; margin: -500px; }'
      )            
    ))
  }
  
}
