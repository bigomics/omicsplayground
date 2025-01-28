##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##


record_UPGRADE <- function(reportName = "UPGRADE_button") {
  if (reportName %in% names(UPGRADE_LOGGER$log)) {
    UPGRADE_LOGGER$log[[reportName]] <- UPGRADE_LOGGER$log[[reportName]] + 1
  } else {
    UPGRADE_LOGGER$log[[reportName]] <- 1
  }
  # Convert the log to a string format
  UPGRADE_LOGGERStr <- paste(names(UPGRADE_LOGGER$log), UPGRADE_LOGGER$log, sep = " = ", collapse = "; ")
  UPGRADE_LOGGER$str <- UPGRADE_LOGGERStr # Update the log with the string representation
}
