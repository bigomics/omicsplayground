##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##


record_plot_download <- function(plotName) {
  if (plotName %in% names(PLOT_DOWNLOAD_LOGGER$log)) {
    PLOT_DOWNLOAD_LOGGER$log[[plotName]] <- PLOT_DOWNLOAD_LOGGER$log[[plotName]] + 1
  } else {
    PLOT_DOWNLOAD_LOGGER$log[[plotName]] <- 1
  }
  # Convert the log to a string format
  PLOT_DOWNLOAD_LOGGERStr <- paste(names(PLOT_DOWNLOAD_LOGGER$log), PLOT_DOWNLOAD_LOGGER$log, sep = " = ", collapse = "; ")
  PLOT_DOWNLOAD_LOGGER$str <- PLOT_DOWNLOAD_LOGGERStr # Update the log with the string representation
}
