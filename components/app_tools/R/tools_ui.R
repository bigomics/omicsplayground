##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

tools_ui <- function(id) {
  ns <- shiny::NS(id) ## namespace
  require(bslib)

  ui <- page_fluid(

    ## head
    layout_columns(
      style = "text-align: center; padding: 80px 0 50px 0;",
      col_widths = c(-4,4,-4),
      h1("Smart Tools"),
      p("Handy standalone utilities for your bioinformatics", style="margin-top: -20px;"),
      shiny::textInput(ns("tools_search"), NULL, placeholder="Search...")      
    ),

    ## cards
    layout_columns(
      col_widths = 4,
      style = "padding: 0 15%;",
      row_heights = "320px",
      card(
        card_image(file="https://i0.wp.com/bigomics.ch/wp-content/uploads/2023/08/computer-with-op.webp?zoom=2&resize=280%2C160&ssl=1", width=280, height=160, class="p-3"),
        card_header("GraphSmart Prism"),
        card_body(
          p("Pimp up your plots using AI by just saying what you want: 'Plot a fancy volcano'"),
          actionButton(ns("runtool1"), "Run", class="action-pill")
        )
      ),
      card(
        card_image(file="https://i0.wp.com/bigomics.ch/wp-content/uploads/2023/07/Omics-Playground-laptop-2023-e1689347474214.webp?zoom=2&resize=280%2C160&ssl=1", width=280, height=160, class="p-3"),
        card_header("ID SmartConverter"),
        card_body(
          p("Automatically convert and annotate your transcriptomic, proteomics or lipidomics features using the latest databases."),
          actionButton(ns("runtool2"), "Run", class="action-pill")
        )
      ),
      card(
        card_image(file="https://i0.wp.com/bigomics.ch/wp-content/uploads/2023/07/Omics-Playground-laptop-2023-e1689347474214.webp?zoom=2&resize=280%2C160&ssl=1",
          width=280, height=160, class="p-3" ),
        card_header("BatchEffect SmartAnalyzer"),
        card_body(
          p("Analyze and correct your data for possible batch effects or unwanted covariates."),
          actionButton(ns("runtool3"), "Run", class="action-pill")
        )
      )
      ## --- end of cards ---
    )
  )
  
  return(ui)
}
