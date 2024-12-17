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

shiny::addResourcePath("custom", "www")



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
ETC <- file.path(OPG, "etc") ## location of options, settings, DB files
FILES <- file.path(OPG, "lib")
FILESX <- file.path(OPG, "libx")
APPDIR <- file.path(OPG, "components/app/R")
PGX.DIR <- file.path(OPG, "data")
SHARE.DIR <- file.path(OPG, "data_shared")
PUBLIC.DIR <- file.path(OPG, "data_public")
SIGDB.DIR <- file.path(OPG, "libx/sigdb")

## Set files
ACCESS_LOGFILE <- file.path(ETC, "access.log")

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
DEVMODE <- FALSE

## Allow API like calls
ALLOW_URL_QUERYSTRING <- FALSE

## Determine if we are in ShinyProxy
SHINYPROXY <- (Sys.getenv("SHINYPROXY_USERNAME") != "" && "omicsplayground" %in% dir("/"))
USERNAME <- "anonymous"
if (SHINYPROXY) USERNAME <- Sys.getenv("SHINYPROXY_USERNAME")

main.start_time <- Sys.time()

WORKDIR <- getwd()
message(">>>>> working directory = ", WORKDIR)
message(">>>>> LOADING INITIAL LIBS")

## some libraries that we often need and load fast
library(shiny)
library(shinyBS)
library(grid)
library(magrittr)
library(future)
library(promises)
future::plan(future::multisession)

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

opt.default <- list(
  TITLE = "Omics Playground",
  AUTHENTICATION = "none", ## none, password, login-code, login-code-redirect
  ALLOW_NEW_USERS = TRUE,
  ALLOW_PERSONAL_EMAIL = FALSE,
  USE_CREDENTIALS = FALSE,
  DOMAIN = NULL,
  BLOCKED_DOMAIN = "bigomics.com|massdynamics.com|pluto.bio|rosalind.bio",
  ## ENABLE_CHIRP         = TRUE,
  ENABLE_DELETE = TRUE,
  ENABLE_PGX_DOWNLOAD = TRUE,
  ENABLE_PUBLIC_SHARE = TRUE,
  ENABLE_UPLOAD = TRUE,
  ENABLE_USERDIR = TRUE,
  ENABLE_USER_SHARE = TRUE,
  ENABLE_USER_LOCK = TRUE,
  ENABLE_HEARTBEAT = TRUE,
  ENABLE_INACTIVITY = TRUE,
  ENABLE_ANNOT = FALSE,
  MAX_DATASETS = 25,
  MAX_SAMPLES = 1000,
  MAX_COMPARISONS = 20,
  MAX_GENES = 20000,
  MAX_GENESETS = 5000,
  MAX_SHARED_QUEUE = 3,
  MAX_SESSIONS = 2,
  TIMEOUT = 0,
  WATERMARK = TRUE,
  APACHE_COOKIE_PATH = OPG,
  DEVMODE = FALSE
)

opt.file <- file.path(ETC, "OPTIONS")
if (!file.exists(opt.file)) stop("FATAL ERROR: cannot find OPTIONS file")
opt <- playbase::pgx.readOptions(file = opt.file, default = opt.default) ## global!

message("\n************************************************")
message("************* SETTING DEFAULTS ***************")
message("************************************************")

defaults.file <- file.path(ETC, "DEFAULTS.yml")
if (!file.exists(defaults.file)) stop("FATAL ERROR: cannot find DEFAULTS.yml file")
DEFAULTS <<- yaml::read_yaml(defaults.file)

## Check and set authentication method
if (Sys.getenv("PLAYGROUND_AUTHENTICATION") != "") {
  auth <- Sys.getenv("PLAYGROUND_AUTHENTICATION")
  message("[GLOBAL] overriding PLAYGROUND_AUTHENTICATION = ", auth)
  opt$AUTHENTICATION <- auth
}
if (Sys.getenv("PLAYGROUND_APACHE_COOKIE_PATH") != "") {
  apache_cookie_path <- Sys.getenv("PLAYGROUND_APACHE_COOKIE_PATH")
  message("[GLOBAL] overriding PLAYGROUND_APACHE_COOKIE_PATH = ", apache_cookie_path)
  opt$APACHE_COOKIE_PATH <- apache_cookie_path
}

## copy to global.R environment
WATERMARK <<- opt$WATERMARK
## TIMEOUT <<- as.integer(opt$TIMEOUT) ## in seconds
PLOTLY_EDITOR <<- opt$PLOTLY_EDITOR
if (opt$DEVMODE) {
  message("!!!!!!!!!!!!!!!!!!!! DEVELOPER MODE !!!!!!!!!!!!!!!!!!!!!!!!")
}

## show options
message("\n", paste(paste(names(opt), "\t= ", sapply(opt, paste, collapse = " ")), collapse = "\n"), "\n")

## ------------------------------------------------
## Check HubSpot connection
## ------------------------------------------------
if (is.null(opt$HUBSPOT_CHECK)) opt$HUBSPOT_CHECK <- FALSE
if (opt$HUBSPOT_CHECK) {
  if (dir.exists(paste0(OPG, "/../omicsplayground-hubconnect"))) {
    dbg("[HubspotConnect]: Folder found, reading files")
    list_files <- list.files(paste0(OPG, "/../omicsplayground-hubconnect/R"), full.names = TRUE)
    sapply(list_files, source)
  } else {
    dbg("[HubspotConnect]: Folder NOT found, seting Hubspot check to FALSE")
    opt$HUBSPOT_CHECK <- FALSE
  }
}

## ------------------------------------------------
## ENABLE/DISABLE BOARDS
## ------------------------------------------------

BOARDS <- c(
  "welcome", "load", "upload", "dataview", "clustersamples", "clusterfeatures",
  "diffexpr", "enrich", "isect", "pathway", "wordcloud", "drug", "sig", "cell",
  "corr", "bio", "cmap", "wgcna", "tcga", "comp", "user", "pcsf",
  "multiomics"
)
##if (is.null(opt$BOARDS_ENABLED)) opt$BOARDS_ENABLED <- BOARDS
opt$BOARDS_ENABLED <- BOARDS
ENABLED <- array(BOARDS %in% opt$BOARDS_ENABLED, dimnames = list(BOARDS))

MODULES <- c('Welcome','Datasets','DataView','Clustering','Expression',
             'GeneSets','Compare','SystemsBio','MultiOmics')
if (is.null(opt$MODULES_ENABLED)) opt$MODULES_ENABLED <- MODULES
MODULES_ENABLED <- array(MODULES %in% opt$MODULES_ENABLED, dimnames = list(MODULES))


## ------------------------------------------------
## SESSION CONTROL
## ------------------------------------------------
if (is.null(opt$HOSTNAME) || opt$HOSTNAME == "") {
  opt$HOSTNAME <- toupper(system("hostname", intern = TRUE))
}
ACTIVE_SESSIONS <- c()
MAX_SESSIONS <- 3
if (!is.null(opt$MAX_SESSIONS)) MAX_SESSIONS <- opt$MAX_SESSIONS

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

## Initialize plot download logger
PLOT_DOWNLOAD_LOGGER <<- reactiveValues(log = list(), str = "")

## Initialize translator
library(shiny.i18n)
DICTIONARY <- file.path(FILES, "translation.json")
i18n <- shiny.i18n::Translator$new(translation_json_path = DICTIONARY)
i18n$set_translation_language("RNA-seq")

## Setup reticulate
##reticulate::use_virtualenv("reticulate")

