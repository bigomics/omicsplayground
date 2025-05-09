##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## From https://github.com/plotly/plotly.js/blob/master/src/components/modebar/buttons.js
all.plotly.buttons <- c(
  "toImage",
  "senDataToCloud", "editInChartStudio", "zoom2d", "pan2d", "select2d",
  "lasso2d", "drawclosedpath", "drawopenpath", "drawline", "drawrect",
  "drawcircle", "eraseshape", "zoomIn2d", "zoomOut2d",
  "autoScale2d", "resetScale2d", "zoom3d", "pan3d",
  "orbitRotation", "tableRotation", "resetCameraDefault3d",
  "resetCameraLastSave3d", "hoverClosest3d", "zoomInGeo",
  "zoomOutGeo", "resetGeo", "hoverClosestGeo", "hoverClosestGl2d",
  "hoverClosestPie", "resetViewSankey", "toggleHover",
  "hoverClosestCartesian", "hoverCompareCartesian",
  "resetViews", "toggleSpikelines",
  "resetViewMapbox", "zoomInMapbox", "zoomOutMapbox"
)

DEFAULT_FONT <- "Cubano"
DEFAULT_FONT <- "Lato"

PLOTLY_COLORS <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
)

plotly_default <- function(e) {
  e %>%
    plotly::layout(
      font = list(family = DEFAULT_FONT),
      margin = list(l = 0, r = 0, t = 0, b = 0), ## remove margins
      legend = list(
        font = list(family = DEFAULT_FONT)
      ),
      hoverlabel = list(
        align = "left",
        font = list(family = DEFAULT_FONT)
      )
    ) %>%
    ## plotly::config(displayModeBar = FALSE) %>%
    plotly::config(
      displaylogo = FALSE,
      modeBarButtons = list(list("toImage", "resetScale2d")),
      toImageButtonOptions = list(format = "svg", height = 500, width = 900)
    )
}

plotly_modal_default <- function(e) {
  e %>%
    plotly::layout(
      font = list(family = DEFAULT_FONT, size = 18),
      margin = list(l = 0, r = 0, t = 0, b = 0), ## remove margins
      legend = list(
        font = list(family = DEFAULT_FONT, size = 18)
      ),
      hoverlabel = list(
        align = "left",
        font = list(family = DEFAULT_FONT)
      )
    ) %>%
    ## plotly::config(displayModeBar = FALSE) %>%
    plotly::config(
      displaylogo = FALSE,
      modeBarButtons = list(list("toImage", "resetScale2d")),
      toImageButtonOptions = list(format = "svg", height = 500, width = 900)
    )
}
