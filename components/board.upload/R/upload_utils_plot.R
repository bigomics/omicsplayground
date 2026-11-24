##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##


## =============================================================================
## ==================== PLOT OUTPUT ERROR HANDLING =============================
## =============================================================================


## Run a base graphics plot function and return an error plot when it fails.
## Shiny still serves the partially drawn device of a plot function that
## raised an error, so the user would see a stale or partial plot instead of
## the error. Drawing the error plot overrides that partial content; the plain
## plot function (without catch) is used by later renders to retry.
catch_plot_error <- function(func) {
  tryCatch(
    {
      func()
    },
    error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Error:", conditionMessage(e)))
      invisible(NULL)
    }
  )
}
