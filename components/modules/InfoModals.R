##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


show_upgrade_modal <- function(timeout.min = 40) {
  require(shiny)
  msg <- HTML(paste0("<center><h4><b>Ditch the ", timeout.min, "-minute limit</h4>
Upgrade today and get advanced analysis features of Playground Pro<br>without the time limit.</b></center>"))

  tbl <- fluidRow(
    column(1, ),
    column(
      5,
      div("\u2713 No 40-minute limit"),
      div("\u2713 Upload 20 datasets"),
      div("\u2713 No watermark")
    ),
    column(
      5,
      div("\u2713 Drug sensitivity analysis"),
      div("\u2713 WGCNA analysis"),
      div("\u2713 Experiment similarity analysis")
    ),
    column(1, )
  )

  shiny::showModal(shiny::modalDialog(
    msg,
    br(),
    tbl,
    footer = fillRow(
      height = 35, flex = c(NA, 1, NA),
      shiny::actionButton("upgrade", "Yes, upgrade me!"),
      br(),
      shiny::modalButton("Dismiss")
    ),
    size = "m",
    easyClose = TRUE
  ))
}
