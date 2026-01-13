##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

withTooltip <- function(
  el,
  title,
  placement = "bottom",
  trigger = NULL,
  options = NULL
) {
  if (!is.null(trigger)) {}

  if (!is.null(options)) {}

  htmltools::tagAppendAttributes(
    el,
    title = title,
    `data-bs-placement` = placement,
    `data-bs-toggle` = "tooltip",
    `data-bs-trigger` = "hover"
  )
}
