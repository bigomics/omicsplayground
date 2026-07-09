##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
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
## Make the PGX directory visible to omicsagentovi's disk-scanning tools (list_pgx, load_pgx)
options(omicspgxmcp.data_dir = PGX.DIR)
## Persistent copilot chat sessions + uploaded docs are anchored on
## auth$user_dir at session time (resolved inside CopilotBoardServer),
## matching the convention used by board.loading / compare / connectivity.
SHARE.DIR <- file.path(OPG, "data_shared")
PUBLIC.DIR <- file.path(OPG, "data_public")
SIGDB.DIR <- file.path(OPG, "libx/sigdb")

## Set files
ACCESS_LOGFILE <- file.path(ETC, "access.log")

# Fail fast at startup if the installed omicsai lacks the provider-catalog API
# the AI features depend on (model menus, live model discovery). A clear stop
# here beats a confusing "could not find function" surfacing deep inside a
# running Shiny session.
.opg_require_omicsai_catalog_api <- function() {
  required <- c(
    "ai_known_models",
    "ai_provider_catalog",
    "ai_select_model",
    "ai_validate_model",
    "ai.list_provider_models"
  )
  exports <- tryCatch(getNamespaceExports("omicsai"), error = function(e) character(0))
  missing <- setdiff(required, exports)
  if (!length(missing)) {
    return(invisible(TRUE))
  }

  stop(
    "Installed omicsai is missing required provider catalog API exports: ",
    paste(missing, collapse = ", "),
    ". Install omicsai >= 0.3.2 or update omicsai in renv.",
    call. = FALSE
  )
}

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

VERSION <- scan(file.path(OPG, "VERSION"), character())[1]
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
.opg_require_omicsai_catalog_api()
source(file.path(APPDIR, "ai_model_policy.R"), local = TRUE)

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
  ENABLE_DELETE = TRUE,
  ENABLE_PGX_DOWNLOAD = TRUE,
  ENABLE_PUBLIC_SHARE = TRUE,
  ENABLE_PUBLIC_LOAD = FALSE,
  ENABLE_PUBLIC_DELETE = FALSE,
  ENABLE_UPLOAD = TRUE,
  ENABLE_USERDIR = TRUE,
  ENABLE_USER_SHARE = TRUE,
  ENABLE_USER_LOCK = TRUE,
  ENABLE_HEARTBEAT = TRUE,
  ENABLE_INACTIVITY = TRUE,
  INACTIVITY_TIMEOUT = 1800,
  ENABLE_ANNOT = FALSE,
  ENABLE_METADATA = FALSE,
  ENABLE_UPGRADE = FALSE,
  ENCRYPTED_EMAIL = FALSE,
  MAX_DATASETS = 25,
  MAX_PUBLIC_DATASETS = NULL,
  MAX_SAMPLES = 1000,
  MAX_COMPARISONS = 20,
  MAX_GENES = 20000,
  MAX_GENESETS = 5000,
  MAX_METH_FEATURES = 450000,
  MAX_SHARED_QUEUE = 3,
  MAX_SESSIONS = 2,
  TIMEOUT = 0,
  WATERMARK = TRUE,
  APACHE_COOKIE_PATH = OPG,
  ALLOW_CUSTOM_FC = FALSE,
  DEVMODE = FALSE,
  USER_LEVEL = 'PRO',  
  ENABLE_MULTIOMICS = TRUE,
  ENABLE_ACROSS = FALSE,
  ENABLE_COOKIE_LOGIN = TRUE,
  PUBLIC_DATASETS_LABEL = "Public Datasets",
  LLM_MAXTURNS = 100,
  ENABLE_AI = FALSE,
  AI_PROVIDERS_ENABLED = c("bigomics", "openai", "anthropic", "google",
                           "github", "mistral", "custom"),
  AI_PROVIDER_LOCKED   = FALSE
)

opt.file <- file.path(ETC, "OPTIONS")
if (!file.exists(opt.file)) stop("FATAL ERROR: cannot find OPTIONS file")
opt <- playbase::pgx.readOptions(file = opt.file, default = opt.default) ## global!

message("\n************************************************")
message("************* SETTING DEFAULTS ***************")
message("************************************************")

defaults.file <- file.path(ETC, "DEFAULTS.yml")
if (file.exists(defaults.file)) {
  DEFAULTS <<- yaml::read_yaml(defaults.file)
} else {
  message("[GLOBAL] DEFAULTS.yml not found, using default configuration")
  DEFAULTS <<- list(
    computation_options = list(
      probe_filtering = list(
        default = c(
          "remove.notexpressed",
          "remove.unknown",
          "only.proteincoding"
        ),
        proteomics = list()
      )
    ),
    qc = list(
      impute = TRUE
    )
  )
}

## Load metadata options configuration
metadata.file <- file.path(ETC, "metadata_options.yml")
if (file.exists(metadata.file)) {
  METADATA_OPTIONS <<- yaml::read_yaml(metadata.file)
  message("[GLOBAL] Loaded metadata_options.yml")
} else {
  message("[GLOBAL] metadata_options.yml not found, metadata feature disabled")
  METADATA_OPTIONS <<- list(fields = list())
}

## Load OPG AI model selection policy. Provider/model facts are owned by
## omicsai::ai_provider_catalog()/ai_known_models(); this JSON file only
## overlays enablement, filtering, ordering, and defaults for OPG menus.
ai_model_policy.file <- file.path(ETC, "ai_model_policy.json")
AI_MODEL_POLICY <<- .opg_ai_read_policy(ai_model_policy.file)
message("[GLOBAL] Loaded AI model selection policy")

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
  "welcome", "summary", "load", "upload", "dataview", "clustersamples", "clusterfeatures",
  "diffexpr", "enrich", "isect", "pathway", "wordcloud", "drug", "sig", "cell",
  "corr", "bio", "cmap", "wgcna", "tcga", "comp", "user", "pcsf",
  "multiomics", "ideograms"
)
## if (is.null(opt$BOARDS_ENABLED)) opt$BOARDS_ENABLED <- BOARDS
opt$BOARDS_ENABLED <- BOARDS
ENABLED <- array(BOARDS %in% opt$BOARDS_ENABLED, dimnames = list(BOARDS))

MODULES <- c(
  "Welcome", "Summary", "Datasets", "DataView", "Clustering", "Expression",
  "GeneSets", "Compare", "SystemsBio", "MultiOmics", "WGCNA", "Epigenomics"
)
if (is.null(opt$MODULES_ENABLED)) opt$MODULES_ENABLED <- MODULES
if (is.null(opt$MODULES_MULTIOMICS)) opt$MODULES_MULTIOMICS <- MODULES
if (is.null(opt$MODULES_TRANSCRIPTOMICS)) opt$MODULES_TRANSCRIPTOMICS <- MODULES
if (is.null(opt$MODULES_METHYLOMICS)) opt$MODULES_METHYLOMICS <- setdiff(MODULES, "MultiOmics")
MODULES_ENABLED <- array(MODULES %in% opt$MODULES_ENABLED, dimnames = list(MODULES))
MODULES_MULTIOMICS <- array(MODULES %in% opt$MODULES_MULTIOMICS, dimnames = list(MODULES))
MODULES_TRANSCRIPTOMICS <- array(MODULES %in% opt$MODULES_TRANSCRIPTOMICS, dimnames = list(MODULES))
MODULES_METHYLOMICS <- array(MODULES %in% opt$MODULES_METHYLOMICS, dimnames = list(MODULES))
MODULES_LOADED <- array(rep(FALSE, length(MODULES)), dimnames = list(MODULES))

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

## Initialize report download logger
REPORT_DOWNLOAD_LOGGER <<- reactiveValues(log = list(), str = "")

## Initialize upgrade button logger
UPGRADE_LOGGER <<- reactiveValues(log = list(), str = "")

## Initialize translator
library(shiny.i18n)
DICTIONARY <- file.path(FILES, "translation.json")
i18n <- shiny.i18n::Translator$new(translation_json_path = DICTIONARY)
i18n$set_translation_language("RNA-seq")

## LLM model setup. Provider/model facts come from omicsai; OPG owns only the
## menu policy overlay in etc/ai_model_policy.json.
opt$AI_PROVIDERS <- unique(unlist(strsplit(as.character(opt$AI_PROVIDERS_ENABLED), ";")))
opt$AI_MODELS    <- .opg_ai_build_models(AI_MODEL_POLICY, opt$AI_PROVIDERS)
opt$AI_MENU_REPORTS          <- .opg_ai_menu_allowlist(opt$AI_MODELS, opt$AI_PROVIDERS, "reports")
opt$AI_MENU_IMAGES           <- .opg_ai_menu_allowlist(opt$AI_MODELS, opt$AI_PROVIDERS, "images")
opt$AI_MENU_COPILOT_DEEP     <- .opg_ai_menu_allowlist(opt$AI_MODELS, opt$AI_PROVIDERS, "copilot_deep")
opt$AI_MENU_COPILOT_BALANCED <- .opg_ai_menu_allowlist(opt$AI_MODELS, opt$AI_PROVIDERS, "copilot_balanced")
## LLM_MAXTURNS is read from etc/OPTIONS — single source of truth.

## Copilot tier selection — verify against the omicsagentovi registry
if (is.null(opt$COPILOT_MODEL)) {
  opt$COPILOT_MODEL <- "copilot-default"
}
if (requireNamespace("omicsagentovi", quietly = TRUE)) {
  .valid_tiers <- omicsagentovi::ovi_copilot_tiers()
  .bad_tiers <- setdiff(opt$COPILOT_MODEL, .valid_tiers)
  if (length(.bad_tiers) > 0L) {
    warning(sprintf(
      "[global] COPILOT_MODEL has unknown tiers: %s. Valid: %s. Dropping unknown entries.",
      paste(.bad_tiers, collapse = ", "), paste(.valid_tiers, collapse = ", ")
    ))
    opt$COPILOT_MODEL <- intersect(opt$COPILOT_MODEL, .valid_tiers)
  }
  if (length(opt$COPILOT_MODEL) == 0L) opt$COPILOT_MODEL <- "copilot-default"
  rm(.valid_tiers, .bad_tiers)
}

## Setup reticulate
## reticulate::use_virtualenv("reticulate")
