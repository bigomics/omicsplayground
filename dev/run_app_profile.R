##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Run the app under R's sampling profiler (Rprof), started and stopped on
## demand, so a slow interaction can be attributed to actual R functions.
##
## Profiling is controlled by a sentinel file rather than a timer: get the app
## into the state you care about first (e.g. many boards open), then switch
## profiling on, reproduce the slowness, and switch it off. A timed window
## anchored to app start just profiles the boot.
##
##   Rscript dev/run_app_profile.R          # PROF_FLAG defaults to /tmp/opg-profile.on
##   touch  $PROF_FLAG                      # start profiling
##   rm     $PROF_FLAG                      # stop, flushing $PROF_FILE
##   Rscript -e 'summaryRprof("<PROF_FILE>")'
##

PROF_FILE <- Sys.getenv("PROF_FILE", "/tmp/opg-Rprof.out")
PROF_FLAG <- Sys.getenv("PROF_FLAG", "/tmp/opg-profile.on")

setwd("components/app/R")
source("global.R")
source("ui.R")
source("server.R")

local({
  running <- FALSE
  poll <- function() {
    on.exit(later::later(poll, delay = 1))
    wanted <- file.exists(PROF_FLAG)
    if (wanted && !running) {
      running <<- TRUE
      Rprof(PROF_FILE, interval = 0.02, line.profiling = TRUE)
      message("[PROFILE] >>> Rprof STARTED <<<")
    } else if (!wanted && running) {
      Rprof(NULL)
      running <<- FALSE
      message("[PROFILE] >>> Rprof STOPPED, wrote ", PROF_FILE, " <<<")
    }
  }
  later::later(poll, delay = 1)
})

shinyApp(
  ui = app_ui,
  server = app_server,
  uiPattern = ".*",
  options = list(
    launch.browser = FALSE,
    host = "0.0.0.0",
    port = 3838
  )
)
