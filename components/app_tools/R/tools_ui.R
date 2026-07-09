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
        bslib::card_image(file="https://i0.wp.com/bigomics.ch/wp-content/uploads/2023/08/computer-with-op.webp?zoom=2&resize=280%2C160&ssl=1", width=280, height=160, class="p-3"),
        bslib::card_header("GraphSmart Prism"),
        bslib::card_body(
          p("Pimp up your plots using AI by just saying what you want: 'Plot a fancy volcano'"),
          actionButton(ns("runtool_prism"), "Run", class="action-pill")
        )
      ),
      bslib::card(
        class = "tools-card",
        #bslib::card_image(file="https://i0.wp.com/bigomics.ch/wp-content/uploads/2023/07/Omics-Playground-laptop-2023-e1689347474214.webp?zoom=2&resize=280%2C160&ssl=1", width=280, height=160, class="p-3" ),
        bslib::card_image(src=base64enc::dataURI(file = "www/qsee-bsee-logo.png"),
          width=280, height=160, class="p-3" ),
        bslib::card_header("Qsee/Bsee"),
        bslib::card_body(
          p("Visual quality control (QC) and batch effects analyzer for your raw data."),
          actionButton(ns("runtool_qc"), "Run", class="action-pill")
        )
      ),
      bslib::card(
        class = "tools-card",        
        bslib::card_image(file="https://i0.wp.com/bigomics.ch/wp-content/uploads/2023/07/Omics-Playground-laptop-2023-e1689347474214.webp?zoom=2&resize=280%2C160&ssl=1", width=280, height=160, class="p-3"),
        bslib::card_header("ID SmartConverter"),
        bslib::card_body(
          p("Automatically convert and annotate your transcriptomic, proteomics or lipidomics features using the latest databases."),
          actionButton(ns("runtool_idconvert"), "Run", class="action-pill")
        )
      )      
      ## --- end of cards ---
    )
  )
  
  return(ui)
}
