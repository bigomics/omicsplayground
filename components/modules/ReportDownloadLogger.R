##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##


record_report_download <- function(reportName) {
  if (reportName %in% names(REPORT_DOWNLOAD_LOGGER$log)) {
    REPORT_DOWNLOAD_LOGGER$log[[reportName]] <- REPORT_DOWNLOAD_LOGGER$log[[reportName]] + 1
  } else {
    REPORT_DOWNLOAD_LOGGER$log[[reportName]] <- 1
  }
  # Convert the log to a string format
  REPORT_DOWNLOAD_LOGGERStr <- paste(names(REPORT_DOWNLOAD_LOGGER$log), REPORT_DOWNLOAD_LOGGER$log, sep = " = ", collapse = "; ")
  REPORT_DOWNLOAD_LOGGER$str <- REPORT_DOWNLOAD_LOGGERStr # Update the log with the string representation
}
