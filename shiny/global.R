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

Sys.setlocale("LC_CTYPE","en_US.UTF-8") 
Sys.setlocale("LC_TIME","en_US.UTF-8")
##Sys.setlocale("LC_ALL", "C")  ## really??
Sys.setenv("_R_CHECK_LENGTH_1_CONDITION_" = "true")

options(shiny.maxRequestSize = 999*1024^2)  ## max 999Mb upload
reticulate::use_miniconda('r-reticulate')

message("[MAIN] reading global.R ...")
##OPG     = "~/Playground/omicsplayground"
OPG       = ".."
RDIR      = file.path(OPG,"R")
FILES     = file.path(OPG,"lib")
FILESX    = file.path(OPG,"libx")
PGX.DIR   = file.path(OPG,"data")
SIGDB.DIR = file.path(OPG,"libx/sigdb")

AUTHENTICATION = "none"
WATERMARK = FALSE
DEBUG     = FALSE
DEV       = dir.exists('/home/kwee')
##DEV     = FALSE
DEBUG     = TRUE
TIMEOUT   = 0

ALLOW_URL_QUERYSTRING = FALSE
ALLOW_URL_QUERYSTRING = TRUE

if(0 && DEV) {
  ## being pedantic... (https://adv-r.hadley.nz)
  options(warnPartialMatchDollar = TRUE)
  options(warnPartialMatchArgs = TRUE)    
  DEBUG  = TRUE
}

## Determine if we are in ShinyProxy
SHINYPROXY = (Sys.getenv("SHINYPROXY_USERNAME")!="" && "omicsplayground" %in% dir("/"))
USERNAME   = "anonymous"
if(SHINYPROXY) USERNAME = Sys.getenv("SHINYPROXY_USERNAME")

main.start_time <- Sys.time()

WORKDIR = getwd()
message(">>>>> working directory = ",WORKDIR)
message(">>>>> LOADING INITIAL LIBS")

## some libraries that we often need and load fast
library(shiny)
library(shinyBS)
library(grid)

source("utils/utils.R", local = TRUE)

message("***********************************************")
message("***** RUNTIME ENVIRONMENT VARIABLES ***********")
message("***********************************************")

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

DEV = (DEV && dir.exists("modulesx")) 
##DEV = FALSE
if(DEV) {
    message('!!!!!!!!!! DEVELOPER MODE !!!!!!!!!!!!!!')
}

message("\n************************************************")
message("**************** READ FUNCTIONS ****************")
message("************************************************")

source(file.path(RDIR,"pgx-include.R"))    ## lots of libraries and source()
source(file.path(RDIR,"pgx-functions.R")) ## functions...
source(file.path(RDIR,"pgx-files.R"))     ## file functions
source(file.path(RDIR,"pgx-init.R"))
source(file.path(RDIR,"auth.R"))


if(0) {    
    ## pgx.initDatasetFolder(PGX.DIR, force=TRUE, verbose=1)    
    load("../data/geiger2016-arginine.pgx")
    load("../data/GSE10846-dlbcl-nc.pgx")
    load("../data/bojkova2020-sarscov2-RC2.pgx")
    load("../data/gtex-aging-n40svaNnm.pgx")
    load("../data/axel-test3.pgx")        
    ngs = pgx.initialize(ngs)
}

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
## add handlerManager for log/crash reports
## --------------------------------------------------------------------


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

modules <- dir("modules", pattern="Module.R$")
modules
for(m in modules) {
    message("[MAIN] loading module ",m)
    source(paste0("modules/",m))
}

BOARDS <- c("load","view","clust","expr","enrich","isect","func",
            "word","drug","sig","scell","cor","bio","cmap","ftmap",
            "wgcna", "tcga","multi","system","qa","corsa","comp","user")
if(is.null(opt$BOARDS_ENABLED))  opt$BOARDS_ENABLED = BOARDS
if(is.null(opt$BOARDS_DISABLED)) opt$BOARDS_DISABLED = NA

ENABLED  <- array(BOARDS %in% opt$BOARDS_ENABLED, dimnames=list(BOARDS))
DISABLED <- array(BOARDS %in% opt$BOARDS_DISABLED, dimnames=list(BOARDS))
ENABLED  <- ENABLED & !DISABLED
ENABLED

# boards <- dir("boards", pattern="Board.R$")
# boards
# for(m in boards) {
#     message("[MAIN] loading board ",m)
#     source(paste0("boards/",m))
# }

# load ui for each board
source("./boards/ui/biomarker_ui.R", encoding = "UTF-8")
source("./boards/ui/clustering_ui.R", encoding = "UTF-8")
source("./boards/ui/compare_ui.R", encoding = "UTF-8")
source("./boards/ui/connectivity_ui.R", encoding = "UTF-8")
source("./boards/ui/correlation_ui.R", encoding = "UTF-8")
source("./boards/ui/dataview_ui.R", encoding = "UTF-8")
source("./boards/ui/drugconnectivity_ui.R", encoding = "UTF-8")
source("./boards/ui/enrichment_ui.R", encoding = "UTF-8")
source("./boards/ui/expression_ui.R", encoding = "UTF-8")
source("./boards/ui/featuremap_ui.R", encoding = "UTF-8")
source("./boards/ui/functional_ui.R", encoding = "UTF-8")
source("./boards/ui/intersection_ui.R", encoding = "UTF-8")
source("./boards/ui/loading_ui.R", encoding = "UTF-8")
source("./boards/ui/signature_ui.R", encoding = "UTF-8")
source("./boards/ui/singlecell_ui.R", encoding = "UTF-8")
source("./boards/ui/tcga_ui.R", encoding = "UTF-8")
source("./boards/ui/user_ui.R", encoding = "UTF-8")
source("./boards/ui/wgcna_ui.R", encoding = "UTF-8")
source("./boards/ui/wordcloud_ui.R", encoding = "UTF-8")


# load server for each board
source("./boards/server/biomarker_server.R", encoding = "UTF-8")
source("./boards/server/clustering_server.R", encoding = "UTF-8")
source("./boards/server/compare_server.R", encoding = "UTF-8")
source("./boards/server/connectivity_server.R", encoding = "UTF-8")
source("./boards/server/correlation_server.R", encoding = "UTF-8")
source("./boards/server/dataview_server.R", encoding = "UTF-8")
source("./boards/server/drugconnectivity_server.R", encoding = "UTF-8")
source("./boards/server/enrichment_server.R", encoding = "UTF-8")
source("./boards/server/expression_server.R", encoding = "UTF-8")
source("./boards/server/featuremap_server.R", encoding = "UTF-8")
source("./boards/server/functional_server.R", encoding = "UTF-8")
source("./boards/server/intersection_server.R", encoding = "UTF-8")
source("./boards/server/loading_server.R", encoding = "UTF-8")
source("./boards/server/signature_server.R", encoding = "UTF-8")
source("./boards/server/singlecell_server.R", encoding = "UTF-8")
source("./boards/server/tcga_server.R", encoding = "UTF-8")
source("./boards/server/user_server.R", encoding = "UTF-8")
source("./boards/server/wgcna_server.R", encoding = "UTF-8")
source("./boards/server/wordcloud_server.R", encoding = "UTF-8")


##ENABLED[c("wgcna","system","multi")] <- FALSE
ENABLED[c("system","multi","corsa")] <- FALSE 
if(0 && DEV && dir.exists("modulesx")) {
    ## Very early development modules/boards (ALWAYS SHOW FOR DEV)
    ##
    xboards <- dir("modulesx", pattern="Board.R$")
    xboards
    m=xboards[1]
    for(m in xboards) {
        message("[MAIN] loading DEVELOPMENT modules ",m)
        source(paste0("modulesx/",m))
    }
    ENABLED[] <- TRUE  ## enable all modules
    boards <- unique(c(boards, xboards))
}
ENABLED

## disable connectivity map if we have no signature database folder
has.sigdb <- length(dir(SIGDB.DIR,pattern="sigdb.*h5"))>0
has.sigdb
if(has.sigdb==FALSE) ENABLED["cmap"] <- FALSE

MAINTABS = c("DataView","Clustering","Expression","Enrichment",
             "Signature","CellProfiling","DEV")

main.init_time <- round(Sys.time() - main.start_time,digits=4)
main.init_time
message("[MAIN] main init time = ",main.init_time," ",attr(main.init_time,"units"))