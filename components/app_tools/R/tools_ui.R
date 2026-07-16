##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

tools_ui <- function(id) {
  ns <- shiny::NS(id) ## namespace
  require(bslib)

  ui <- page_fluid(

    ## head
    bslib::layout_columns(
      style = "text-align: center; padding: 80px 0 50px 0;",
      col_widths = c(-4,4,-4),
      h1("Smart Tools"),
      p("Handy standalone utilities for your bioinformatics", style="margin-top: -20px;"),
      shiny::textInput(ns("tools_search"), NULL, placeholder="Search...")      
    ),

    ## cards
    bslib::layout_columns(
      col_widths = 4,
      style = "padding: 0 15%;",
      row_heights = "320px",
      bslib::card(
        class = "tools-card",        
        bslib::card_image(src=base64enc::dataURI(file = "www/applets/converter.png"),
          width=320, height=160, class="p-3" ),        
        bslib::card_header("ID Converter"),
        bslib::card_body(
          p("Convert and annotate your features using the latest databases."),
          actionButton(ns("runtool_idconvert"), "Run", class="action-pill")
        )
      ),
      bslib::card(
        class = "tools-card",
        bslib::card_image(src=base64enc::dataURI(file = "www/applets/qsee-bsee.png"),
          width=320, height=160, class="p-3" ),        
        bslib::card_header("Qsee/Bsee"),
        bslib::card_body(
          p("Visual QC analysis and check your data for batch effects."),
          actionButton(ns("runtool_qc"), "Run", class="action-pill")
        )
      ),
      bslib::card(
        class = "tools-card",
        bslib::card_image(src=base64enc::dataURI(file = "www/applets/smartprism.png"),
          width=320, height=160, class="p-3" ),        
        bslib::card_header("SmartPrism"),
        bslib::card_body(
          p("Create figures using AI by just saying what you want: 'Plot a fancy volcano'"),
          actionButton(ns("runtool_prism"), "Run", class="action-pill")
        )
      )
      ## --- end of cards ---
    )
  )
  
  return(ui)
}
