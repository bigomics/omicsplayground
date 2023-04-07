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
  if(!is.null(trigger)) {
#    warning("`trigger` is ignored, used to be in shinyBS::tippify")
  }
  
  if(!is.null(options)) {
#    warning("`options` is ignored, used to be in shinyBS::tippify")
  }
  
  htmltools::tagAppendAttributes(
    el,
    title = title,
    `data-bs-placement` = placement,
    `data-bs-toggle` = "tooltip"
  )
}
