##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

#' AI text card UI
#'
#' PlotModule wrapper for text-only AI output.
#'
#' @param id Shiny module namespace ID
#' @param title Card title
#' @param caption Caption text
#' @param info.text Help text
#' @param height Card height vector c(default, fullscreen)
#' @param width Card width vector c(default, fullscreen)
#'
#' @return Shiny tagList
AiTextCardUI <- function(id,
                         title = "AI Summary",
                         caption = "AI-generated summary.",
                         info.text = "",
                         height = c("100%", TABLE_HEIGHT_MODAL),
                         width = c("auto", "100%")) {
  ns <- shiny::NS(id)

  opts <- shiny::tagList(
    shiny::radioButtons(
      ns("style"),
      "AI summary:",
      choices = c("short", "long"),
      selected = "short",
      inline = TRUE
    ),
    shiny::checkboxInput(
      ns("show_prompt"),
      "Show prompt",
      FALSE
    ),
    shiny::actionButton(
      ns("generate"),
      "Generate!",
      icon = shiny::icon("refresh"),
      class = "btn-outline-primary"
    )
  )

  PlotModuleUI(
    ns("text"),
    outputFunc = shiny::htmlOutput,
    title = title,
    label = "",
    info.text = info.text,
    options = opts,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("pdf")
  )
}

#' AI diagram card UI
#'
#' PlotModule wrapper for visNetwork AI diagrams.
#'
#' @param id Shiny module namespace ID
#' @param title Card title
#' @param caption Caption text
#' @param info.text Help text
#' @param height Card height vector c(default, fullscreen)
#' @param width Card width vector c(default, fullscreen)
#'
#' @return Shiny tagList
AiDiagramCardUI <- function(id,
                            title = "AI Diagram",
                            caption = "AI-generated diagram.",
                            info.text = "",
                            height = c("100%", TABLE_HEIGHT_MODAL),
                            width = c("auto", "100%")) {
  ns <- shiny::NS(id)

  opts <- shiny::tagList(
    shiny::actionButton(
      ns("generate"),
      "Generate Diagram",
      icon = shiny::icon("refresh"),
      class = "btn-outline-primary"
    )
  )

  PlotModuleUI(
    ns("diagram"),
    plotlib = "visnetwork",
    title = title,
    label = "",
    info.text = info.text,
    options = opts,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "html")
  )
}

#' AI image card UI
#'
#' PlotModule wrapper for infographic image output.
#'
#' @param id Shiny module namespace ID
#' @param title Card title
#' @param caption Caption text
#' @param info.text Help text
#' @param height Card height vector c(default, fullscreen)
#' @param width Card width vector c(default, fullscreen)
#'
#' @return Shiny tagList
AiImageCardUI <- function(id,
                          title = "AI Image",
                          caption = "AI-generated infographic.",
                          info.text = "",
                          height = c("100%", TABLE_HEIGHT_MODAL),
                          width = c("auto", "100%")) {
  ns <- shiny::NS(id)

  opts <- shiny::tagList(
    shiny::actionButton(
      ns("generate"),
      "Generate Image",
      icon = shiny::icon("refresh"),
      class = "btn-outline-primary"
    )
  )

  PlotModuleUI(
    ns("image"),
    outputFunc = shiny::htmlOutput,
    title = title,
    label = "",
    info.text = info.text,
    options = opts,
    caption = caption,
    height = height,
    width = width
  )
}
