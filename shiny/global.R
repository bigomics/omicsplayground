##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

#########################################################################
##                                                                     ##
##              Main application for Omics Playground                  ##
##                                                                     ##
#########################################################################

message("\n\n")
message("###############################################################")
message("##################### OMICS PLAYGROUND ########################")
message("###############################################################")
message("\n")

## should we migrate all OPTIONS into this file??

if(Sys.info()["sysname"] != "Windows") {
    Sys.setlocale("LC_TIME","en_US.UTF-8")
}

##Sys.setlocale("LC_ALL", "C")  ## really??
Sys.setenv("_R_CHECK_LENGTH_1_CONDITION_" = "true")
##options(warnPartialMatchDollar = TRUE)
##options(warnPartialMatchArgs = TRUE)    

options(shiny.maxRequestSize = 999*1024^2)  ## max 999Mb upload
options(shiny.fullstacktrace = TRUE)
reticulate::use_miniconda('r-reticulate')

get_pkg_root <- function() {
    pwd <- strsplit(getwd(),split='/')[[1]]
    paste(pwd[1:max(grep("omicsplayground",pwd))],collapse='/')
}

message("[MAIN] reading global.R ...")

## Set folders
OPG       = ".."
OPG       = sub("/omicsplayground/.*","/omicsplayground",getwd())
OPG       = get_pkg_root()
RDIR      = file.path(OPG,"R")
APPDIR    = file.path(OPG,"shiny")
FILES     = file.path(OPG,"lib")
FILESX    = file.path(OPG,"libx")
PGX.DIR   = file.path(OPG,"data")
SIGDB.DIR = file.path(OPG,"libx/sigdb")

AUTHENTICATION = "none"
WATERMARK = FALSE
DEV       = FALSE
DEBUG     = TRUE
TIMEOUT   = 0

## Allow API like calls
ALLOW_URL_QUERYSTRING = FALSE
ALLOW_URL_QUERYSTRING = TRUE

## Determine if we are in ShinyProxy
SHINYPROXY = (Sys.getenv("SHINYPROXY_USERNAME")!="" && "omicsplayground" %in% dir("/"))
USERNAME   = "anonymous"
if(SHINYPROXY) USERNAME = Sys.getenv("SHINYPROXY_USERNAME")

main.start_time <- Sys.time()

if(DEV) {
    message('!!!!!!!!!!!!!!!!!!!! DEVELOPER MODE !!!!!!!!!!!!!!!!!!!!!!!!')
}

WORKDIR = getwd()
message(">>>>> working directory = ",WORKDIR)
message(">>>>> LOADING INITIAL LIBS")

## some libraries that we often need and load fast
library(shiny)
library(shinyBS)
library(grid)
library(ggplot2)
library(concaveman)

source(file.path(APPDIR,"utils/utils.R"), local = TRUE)

message("***********************************************")
message("***** RUNTIME ENVIRONMENT VARIABLES ***********")
message("***********************************************")

if(file.exists("Renviron.site")) {
    message("Loading local Renviron.site")    
    readRenviron("Renviron.site")
}

envcat("SHINYPROXY_USERNAME")
envcat("SHINYPROXY_USERGROUPS")
envcat("PLAYGROUND_AUTHENTICATION")
envcat("PLAYGROUND_USERID")
envcat("PLAYGROUND_EXPIRY")
envcat("PLAYGROUND_QUOTA")
envcat("PLAYGROUND_LEVEL")
envcat("PLAYGROUND_HELLO")

## --------------------------------------------------------------------
## -------------------------- INIT ------------------------------------
## --------------------------------------------------------------------

message("\n***********************************************")
message("*********** SETTING GLOBAL VARIABLES **********")
message("***********************************************")

message("OPG =",OPG)
message("RDIR =",RDIR)
message("FILES =",FILES)
message("FILESX =",FILESX)
message("PGX.DIR =",PGX.DIR)
message("SHINYPROXY = ",SHINYPROXY)

message("\n************************************************")
message("**************** READ FUNCTIONS ****************")
message("************************************************")

source(file.path(RDIR,"pgx-include.R"))    ## lots of libraries and source()
source(file.path(RDIR,"auth.R"))
source(file.path(RDIR,"ggplot-theme.R"))

message("\n************************************************")
message("**************** LOAD LIBRARIES ****************")
message("************************************************")

load(file.path(FILES,"sysdata.rda"),verbose=TRUE)


message("\n************************************************")
message("************* parsing OPTIONS file *************")
message("************************************************")

if(!file.exists("OPTIONS")) stop("FATAL ERROR: cannot find OPTIONS file")
opt <- pgx.readOptions(file="OPTIONS")

## Check and set authentication method
if(Sys.getenv("PLAYGROUND_AUTHENTICATION")!="") {
    auth <- Sys.getenv("PLAYGROUND_AUTHENTICATION")
    message("[ENV] overriding PLAYGROUND_AUTHENTICATION = ",auth)
    opt$AUTHENTICATION = auth
}
if(1 && opt$AUTHENTICATION=="shinyproxy" && !in.shinyproxy()) {
    Sys.setenv("SHINYPROXY_USERNAME"="Test Person")  ## only for testing!!
}
if(1 && opt$AUTHENTICATION=="firebase" && !file.exists("firebase.rds")) {
    message("[ENV] WARNING: Missing firebase.rds file!!! reverting authentication to 'none'")    
    opt$AUTHENTICATION = "none"
    ## opt$ENABLE_USERDIR = FALSE
    ## stop("[MAIN] FATAL Missing firebase.rds file")
}

## copy to global.R environment
WATERMARK <<- opt$WATERMARK
TIMEOUT   <<- as.integer(opt$TIMEOUT)  ## in seconds

## show options
message("\n",paste(paste(names(opt),"\t= ",sapply(opt,paste,collapse=" ")),collapse="\n"),"\n")


## --------------------------------------------------------------------
## --------------------- HANDLER MANAGER ------------------------------
## --------------------------------------------------------------------
## add handlerManager for log/crash reports


http.resp <- getFromNamespace("httpResponse", "shiny")

logHandler <- function(http.req){
  

    dbg("[MAIN.logHandler] >>>>> called! <<<<<")
    ##dbg("[MAIN.logHandler] names(http.req) = ",sort(names(http.req)))
    dbg("[MAIN.logHandler] http.req$PATH_INFO = ",http.req$PATH_INFO)
    
    if(!http.req$PATH_INFO == "/log") {
        return()
    }
    
    query <- shiny::parseQueryString(http.req$QUERY_STRING)
    dbg("[MAIN.logHandler] names(query) = ",names(query))
    dbg("[MAIN.logHandler] query$msg = ",query$msg)
    
    if(is.null(query$msg)) {
        dbg("[MAIN.logHandler] msg is NULL!")
        return(http.resp(400L, "application/json", jsonlite::toJSON(FALSE)))
    }

    if(query$msg == "") {
        dbg("[MAIN.logHandler] msg is empty!")        
        return(http.resp(400L, "application/json", jsonlite::toJSON(FALSE)))
    }
    
    token <- Sys.getenv("HONCHO_TOKEN", "")
    if(token == "") {
        dbg("[MAIN.logHandler] missing HONCHO_TOKEN!")        
        return(http.resp(403L, "application/json", jsonlite::toJSON(FALSE)))
    }
    
    uri <- sprintf("%s/log?token=%s", opt$HONCHO_URL, token)

    ## get the correct log file
    log.file = NULL
    the.log <- "Could not find log file!"
    id <- query$msg  ## use session id as message

    log.dirs <- "~/ShinyApps/log/*log /var/log/shiny-server/*log"
    suppressWarnings( log.file <- system(paste("grep -l -s",id,log.dirs),intern=TRUE) )
    log.file
    log.file <- tail(log.file,1)  ## take newest???
    log.file
    
    if(length(log.file)==0) {
        dbg("[MAIN.logHandler] could not resolve log file for session ID = ",id)
        return(http.resp(403L, "application/json", jsonlite::toJSON(FALSE)))
    }
    
    dbg("[logHandler] reading log.file = ",log.file)
    if(!is.null(log.file)) {
        ##the.log <- readr::read_file(log.file)
        ## truncate the log file
        the.log <- paste(system(paste("grep -B100 -A99999",id,log.file),intern=TRUE),collapse='\n')
    }

    dbg("[logHandler] sending log file... ")    
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

handlerManager <- getFromNamespace("handlerManager", "shiny")
handlerManager$removeHandler("/log")
handlerManager$addHandler(logHandler, "/log")

## --------------------------------------------------------------------
## ----------------- READ MODULES/BOARDS ------------------------------
## --------------------------------------------------------------------

BOARDS <- c("load","view","clust","expr","enrich","isect","func",
            "word","drug","sig","scell","cor","bio","cmap","ftmap",
            "wgcna", "tcga","multi","system","qa","corsa","comp","user")
if(is.null(opt$BOARDS_ENABLED))  opt$BOARDS_ENABLED = BOARDS
if(is.null(opt$BOARDS_DISABLED)) opt$BOARDS_DISABLED = NA

ENABLED  <- array(BOARDS %in% opt$BOARDS_ENABLED, dimnames=list(BOARDS))
DISABLED <- array(BOARDS %in% opt$BOARDS_DISABLED, dimnames=list(BOARDS))
ENABLED  <- ENABLED & !DISABLED
ENABLED

# load UI for each board
boards_ui <- dir(file.path(APPDIR,"boards/ui"), pattern="_ui.R$",full.names=TRUE)
for(b in boards_ui) {
    message("[MAIN] loading UI module ",basename(b))
    source(b, encoding = "UTF-8")
}

# load server for each board
boards_srv <- dir(file.path(APPDIR,"boards/server"), pattern="_server.R$",full.names=TRUE)
for(b in boards_srv) {
    message("[MAIN] loading server module ",basename(b))
    source(b, encoding = "UTF-8")
}

## loading modules
modules <- dir(file.path(APPDIR,"modules"), pattern="Module.R$", full.names=TRUE)
modules
for(m in modules) {
    message("[MAIN] loading module ",basename(m))
    source(m, encoding = "UTF-8")
}

## load plotModules
source(file.path(APPDIR,"modules/plotModules/PlotModule.R"), encoding = "UTF-8")
source(file.path(APPDIR,"modules/plotModules/dataviewTSNEPlotModule.R"), encoding = "UTF-8")

## disable connectivity map if we have no signature database folder
has.sigdb <- length(dir(SIGDB.DIR,pattern="sigdb.*h5"))>0
has.sigdb
if(has.sigdb==FALSE) ENABLED["cmap"] <- FALSE

## Main tab titles
MAINTABS = c("DataView","Clustering","Expression","Enrichment",
             "Signature","CellProfiling","DEV")

## Calculate init time
main.init_time <- round(Sys.time() - main.start_time,digits=4)
main.init_time
message("[MAIN] main init time = ",main.init_time," ",attr(main.init_time,"units"))


message("\n\n")
message("=================================================================")
message("=================== GLOBAL INIT DONE ============================")
message("=================================================================")
message("\n\n")
