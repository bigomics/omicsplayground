
runApplet <- function(applet, pgx, as="gadget", layout=FALSE, size='m', ...) {

  ## random id
  id <- paste0("applet_",sample(1e8,1)) ## random

  app_inputs <- applet$inputs(id=id)
  ## strip away bigdash::tabSettings
  if( inherits(app_inputs, "shiny.tag") &&
        inherits(app_inputs$children[[1]], "shiny.tag.list")) {
    app_inputs <- app_inputs$children[[1]]
  }

  if(layout) {
    ## get ui as list of modules
    main_ui <- applet$ui(id=id, as="taglist")
    names(main_ui) <- NULL
    main_ui <- bslib::layout_columns(
      height = "calc(100vh - 200px)",
      heights_equal = "row",
      gap = "20px",
      col_widths = bslib::breakpoints(sm=12, md=12, lg=6, xl=6, xxl=4),
      !!!main_ui
    )
    
  } else {
    main_ui <- applet$ui(id=id)
  }

  theme <- bslib::bs_add_variables(bslib::bs_theme(),
    "grid-breakpoints" = # here e.g. with lg: 800px;
      "(xs: 0, sm: 576px, md: 768px, lg: 1200px, xl: 1600px, xxl: 2000px)",
    .where = "declarations"
  )
  
  app_ui <- bslib::page_sidebar(
    title = NULL,
    theme = theme,
    sidebar = bslib::sidebar(app_inputs),
    main_ui
  )
  
  r_pgx <- pgx
  if(class(pgx) != "reactivevalues") {
    r_pgx <- do.call("reactiveValues",pgx)
  }

  if(as == "gadget") {    
    app_server <- function(input, output, session) {
      applet$server(id=id, pgx = r_pgx)          
    }
    #shinyApp(app_ui, app_server)
    shiny::runGadget(app_ui, app_server)
  }

  if(as == "modal") {

    sizestyle <- switch( tolower(size),
      "m" = "height: calc(65vh - 180px); width: 50vw;",
      "l" = "height: calc(85vh - 180px); width: 75vw;",     
      "xl" = "height: calc(100vh - 180px); width: 95vw;",
      "fullscreen" = "height: 100vh; width: 100vw;"      
    )

    ## run server
    applet$server(id=id, pgx=r_pgx, ... )
   
    ## Show applet as modal
    showModal( modalDialog(
      title = applet$title,
      style = sizestyle,
      #size = size,
      easyClose = TRUE,
      app_ui,
      tags$style(
        type = 'text/css',
        '.modal-dialog { width: fit-content !important; }'
        ##'.modal-content { width: 100vw; }'
      )            
    ))
  }
  
}
