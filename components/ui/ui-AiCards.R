##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

#' AI text card UI
#'
#' PlotModule wrapper for on-demand AI text generation.
#'
#' @param id Shiny module namespace ID.
#' @param title Card title.
#' @param caption Caption text.
#' @param info.text Help text.
#' @param height Card height vector c(default, fullscreen).
#' @param width Card width vector c(default, fullscreen).
#'
#' @return Shiny UI tags.
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
      choices = c("Short" = "short_summary", "Long" = "long_summary"),
      selected = "short_summary",
      inline = TRUE
    ),
    shiny::checkboxInput(ns("show_prompt"), "Show prompt", FALSE),
    shiny::actionButton(
      ns("generate"),
      "Generate",
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
    no.download = TRUE
  )
}
