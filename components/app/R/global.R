##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


message("\n\n\n")
message("  ___            _          ____  _                                              _ ")
message(" / _ \\ _ __ ___ (_) ___ ___|  _ \\| | __ _ _   _  __ _ _ __ ___  _   _ _ __    __| |")
message("| | | | '_ ` _ \\| |/ __/ __| |_) | |/ _` | | | |/ _` | '__/ _ \\| | | | '_ \\  / _` |")
message("| |_| | | | | | | | (__\\__ \\  __/| | (_| | |_| | (_| | | | (_) | |_| | | | || (_| |")
message(" \\___/|_| |_| |_|_|\\___|___/_|   |_|\\__,_|\\__, |\\__, |_|  \\___/ \\__,_|_| |_| \\__,_|")
message("                                          |___/ |___/                              ")
message("\n\n\n")


message("[GLOBAL] reading global.R ...")

if (Sys.info()["sysname"] != "Windows") {
  Sys.setlocale("LC_TIME", "en_US.UTF-8")
}

Sys.setenv("_R_CHECK_LENGTH_1_CONDITION_" = "true")


options(shiny.maxRequestSize = 999 * 1024^2) ## max 999Mb upload
options(shiny.fullstacktrace = TRUE)
# The following DT global options ensure
# 1. The header scrolls with the X scroll bar
# 2. The Y scroller works properly and no blank rows are displayed
options(DT.options = list(
  autoWidth = FALSE,
  scrollX = TRUE,
  fillContainer = FALSE
))

# Set global modal height values for tables.
# - The SCROLLY_MODAL defines the size of the scroll Y bar on the modals,
# this only defines the srollable part of the table, not the header height.
# - The TABLE_HEIGHT_MODAL defines the whole width of the table + header,
# this will define how close the caption is to the table.
SCROLLY_MODAL <<- "55vh"
TABLE_HEIGHT_MODAL <<- "75vh"


# Get the OPG root folder. Works only from inside the repo as it looks
# up to the closest parent folder matching 'omicsplayground'
get_opg_root <- function() {
  pwd <- getwd()
  dirs <- unlist(strsplit(pwd, "/"))
  root_dirs <- paste(dirs[1:max(grep("omicsplayground", dirs))], collapse = "/")
  root <- paste(root_dirs, collapse = "/")
  return(root)
}

## Set folders

OPG <- get_opg_root()
ETC <- file.path(OPG, "etc")
RDIR <- file.path(OPG, "components/base/R")
APPDIR <- file.path(OPG, "components/app/R")
ETC <- file.path(OPG, "etc") ## location of options, settings, DB files
FILES <- file.path(OPG, "lib")
FILESX <- file.path(OPG, "libx")
PGX.DIR <- file.path(OPG, "data")
SIGDB.DIR <- file.path(OPG, "libx/sigdb")

## like system.file()
pgx.system.file <- function(file = ".", package) {
  package <- sub("^board.", "", package)
  dir <- normalizePath(file.path(OPG, "components", paste0("board.", package)))
  file1 <- file.path(dir, "inst", file)
  if (file.exists(file1) && file != ".") {
    return(file1)
  }
  file.path(dir, file)
}

AUTHENTICATION <- "none"
WATERMARK <- FALSE
DEV <- FALSE
DEBUG <- TRUE
TIMEOUT <- 0

## Allow API like calls
ALLOW_URL_QUERYSTRING <- FALSE


## Determine if we are in ShinyProxy
SHINYPROXY <- (Sys.getenv("SHINYPROXY_USERNAME") != "" && "omicsplayground" %in% dir("/"))
USERNAME <- "anonymous"
if (SHINYPROXY) USERNAME <- Sys.getenv("SHINYPROXY_USERNAME")

main.start_time <- Sys.time()

if (DEV) {
  message("!!!!!!!!!!!!!!!!!!!! DEVELOPER MODE !!!!!!!!!!!!!!!!!!!!!!!!")
}

WORKDIR <- getwd()
message(">>>>> working directory = ", WORKDIR)
message(">>>>> LOADING INITIAL LIBS")

## some libraries that we often need and load fast
library(shiny)
library(shinyBS)
library(grid)
library(magrittr)


source(file.path(APPDIR, "utils/utils.R"), local = TRUE)

message("***********************************************")
message("***** RUNTIME ENVIRONMENT VARIABLES ***********")
message("***********************************************")

if (file.exists("Renviron.site")) {
  message("Loading local Renviron.site...")
  readRenviron("Renviron.site")
}

envcat("SHINYPROXY_USERNAME")
envcat("SHINYPROXY_USERGROUPS")
envcat("PLAYGROUND_AUTHENTICATION")
envcat("PLAYGROUND_USERNAME")

message("\n***********************************************")
message("*********** SETTING GLOBAL VARIABLES **********")
message("***********************************************")

message("OPG =", OPG)
message("RDIR =", RDIR)
message("FILES =", FILES)
message("FILESX =", FILESX)
message("PGX.DIR =", PGX.DIR)
message("APPDIR =", APPDIR)
message("SHINYPROXY = ", SHINYPROXY)

message("\n************************************************")
message("************* SOURCING FUNCTIONS ***************")
message("************************************************")

## MAIN SOURCING FUNCTION. SOURCES ALL R/SHINY CODE. ONLY SOURCE IF
## RUN IN SAME FOLDER.
if (file.exists("global.R")) {
  source(file.path(OPG, "components/00SourceAll.R"), chdir = TRUE)
}

message("\n************************************************")
message("************* PARSING OPTIONS ******************")
message("************************************************")

opt.file <- file.path(ETC, "OPTIONS")
if (!file.exists(opt.file)) stop("FATAL ERROR: cannot find OPTIONS file")
opt <- playbase::pgx.readOptions(file = opt.file) ## THIS IS GLOBAL!!!

## Check and set authentication method
if (Sys.getenv("PLAYGROUND_AUTHENTICATION") != "") {
  auth <- Sys.getenv("PLAYGROUND_AUTHENTICATION")
  message("[GLOBAL] overriding PLAYGROUND_AUTHENTICATION = ", auth)
  opt$AUTHENTICATION <- auth
}

## copy to global.R environment
WATERMARK <<- opt$WATERMARK
TIMEOUT <<- as.integer(opt$TIMEOUT) ## in seconds
PLOTLY_EDITOR <<- opt$PLOTLY_EDITOR

## show options
message("\n", paste(paste(names(opt), "\t= ", sapply(opt, paste, collapse = " ")), collapse = "\n"), "\n")


message("\n************************************************")
message("*********** READ MODULES/BOARDS ****************")
message("************************************************")


BOARDS <- c(
  "welcome", "load", "upload", "dataview", "clustersamples", "clusterfeatures",
  "diffexpr", "enrich", "isect", "pathway", "wordcloud", "drug", "sig", "cell",
  "corr", "bio", "cmap", "wgcna", "tcga", "comp", "user", "pcsf"
)
if (is.null(opt$BOARDS_ENABLED)) opt$BOARDS_ENABLED <- BOARDS
ENABLED <- array(rep(TRUE, length(BOARDS)), dimnames = list(BOARDS))
ENABLED <- array(BOARDS %in% opt$BOARDS_ENABLED, dimnames = list(BOARDS))

## --------------------------------------------------------------------
## --------------------- HANDLER MANAGER ------------------------------
## --------------------------------------------------------------------
## add handlerManager for log/crash reports

http.resp <- getFromNamespace("httpResponse", "shiny")

logHandler <- function(http.req) {
  if (!http.req$PATH_INFO == "/log") {
    return()
  }

  query <- shiny::parseQueryString(http.req$QUERY_STRING)

  if (is.null(query$msg)) {
    return(http.resp(400L, "application/json", jsonlite::toJSON(FALSE)))
  }

  if (query$msg == "") {
    return(http.resp(400L, "application/json", jsonlite::toJSON(FALSE)))
  }

  token <- Sys.getenv("HONCHO_TOKEN", "")
  if (token == "") {
    return(http.resp(403L, "application/json", jsonlite::toJSON(FALSE)))
  }

  uri <- sprintf("%s/log?token=%s", opt$HONCHO_URL, token)

  ## get the correct log file
  log.file <- NULL
  the.log <- "Could not find log file!"
  id <- query$msg ## use session id as message

  log.dirs <- "~/ShinyApps/log/*log /var/log/shiny-server/*log"
  suppressWarnings(log.file <- system(paste("grep -l -s", id, log.dirs), intern = TRUE))
  log.file <- tail(log.file, 1) ## take newest???

  if (length(log.file) == 0) {
    dbg("[GLOBAL.logHandler] could not resolve log file for session ID = ", id)
    return(http.resp(403L, "application/json", jsonlite::toJSON(FALSE)))
  }

  if (!is.null(log.file)) {
    ## truncate the log file
    the.log <- paste(system(paste("grep -B100 -A99999", id, log.file), intern = TRUE), collapse = "\n")
  }

  httr::POST(
    uri,
    body = list(
      msg = query$msg,
      log = "The log!",
      filename = "the_log.log"
    ),
    encode = "json"
  )

  http.resp(400L, "application/json", jsonlite::toJSON(TRUE))
}

## Are we ever going to use this??

message("\n\n")
message("=================================================================")
message("=================== GLOBAL INIT DONE ============================")
message("=================================================================")
message("\n\n")

## Calculate init time
main.init_time <- round(Sys.time() - main.start_time, digits = 4)
main.init_time
message("[GLOBAL] global init time = ", main.init_time, " ", attr(main.init_time, "units"))

shiny::addResourcePath("static", file.path(OPG, "components/app/R/www"))
