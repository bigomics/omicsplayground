##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

#' AI Summary Card UI
#'
#' Standard OmicsPlayground wrapper that bridges PlotModuleUI and
#' omicsai::omicsai_summary_card_ui. Provides consistent card styling
#' across all analysis boards.
#'
#' @param id Shiny module namespace ID
#' @param title Card title (default "AI Summary")
#' @param caption Board-specific caption text
#' @param info.text Help tooltip text (default empty)
#' @param height Card height, vector c(default, fullscreen)
#' @param width Card width, vector c(default, fullscreen)
#'
#' @return Shiny UI element
AISummaryCardUI <- function(id,
                            title = "AI Summary",
                            caption = "AI-generated summary.",
                            info.text = "",
                            height = c("100%", TABLE_HEIGHT_MODAL),
                            width = c("auto", "100%")) {
  ## --- Argument validation ---
  if (!is.character(id) || length(id) != 1 || !nzchar(id)) {
    stop("AISummaryCardUI: `id` must be a non-empty string.")
  }

  card_wrapper <- function(id, content, options, title, label, info.text,
                           caption, height, width, download.fmt, ...) {
    PlotModuleUI(
      id,
      outputFunc = shiny::htmlOutput,
      title = title,
      label = label,
      info.text = info.text,
      options = options,
      caption = caption,
      height = height,
      width = width,
      download.fmt = download.fmt
    )
  }

  omicsai::omicsai_summary_card_ui(
    id = id,
    card_wrapper = card_wrapper,
    title = title,
    label = "",
    info.text = info.text,
    caption = caption,
    height = height,
    width = width
  )
}
