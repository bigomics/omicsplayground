##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics Sagl. All rights reserved.
##


plotly_default <- function(e) {
  e %>%
    plotly::layout(
      font = list(family = "Lato")
    ) %>%
    ## plotly::config(displayModeBar = FALSE) %>%
    plotly::config(displaylogo = FALSE) %>%
    plotly::config(
      modeBarButtons = list(list("toImage", "resetScale2d")),
      toImageButtonOptions = list(format = "svg", height = 500, width = 900)
    )
}

plotly_modal_default <- function(e) {
  e %>%
    plotly::layout(
      font = list(family = "Lato", size = 18),
      legend = list(
        font = list(family = "Lato", size = 18)
      )
    ) %>%
    ## plotly::config(displayModeBar = FALSE) %>%
    plotly::config(displaylogo = FALSE) %>%
    plotly::config(
      modeBarButtons = list(list("toImage", "resetScale2d")),
      toImageButtonOptions = list(format = "svg", height = 500, width = 900)
    )
}
