##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2021 BigOmics Analytics Sagl. All rights reserved.
##

#########################################################################
##                                                                     ##
##              Main application for Omics Playground                  ##
##                                                                     ##
#########################################################################

library(shiny)
library(shinyjs)
library(shinyWidgets)
library(plotly)
library(shinybusy)

message("\n\n")
message("###############################################################")
message("##################### OMICS PLAYGROUND ########################")
message("###############################################################")
message("\n")
message("************************************************")
message("********* RUNTIME ENVIRONMENT VARIABLES ********")
message("************************************************")

##Sys.setenv("S HINYPROXY_USERNAME"="Test Person")
main.start_time <- Sys.time()

envcat <- function(var) message(var," = ",Sys.getenv(var))
envcat("SHINYPROXY_USERNAME")
envcat("SHINYPROXY_USERGROUPS")
envcat("PLAYGROUND_AUTHENTICATION")
envcat("PLAYGROUND_USERID")
envcat("PLAYGROUND_EXPIRY")
envcat("PLAYGROUND_LEVEL")
envcat("PLAYGROUND_HELLO")

## --------------------------------------------------------------------
## -------------------------- INIT ------------------------------------
## --------------------------------------------------------------------

message("\n")
message("***********************************************")
message("*********** SETTING GLOBAL VARIABLES **********")
message("***********************************************")

Sys.setlocale("LC_CTYPE","en_US.UTF-8") 
Sys.setlocale("LC_TIME","en_US.UTF-8")
##Sys.setlocale("LC_ALL", "C")  ## really??
Sys.setenv("_R_CHECK_LENGTH_1_CONDITION_" = "true")

## being pedantic... (https://adv-r.hadley.nz)
options(warnPartialMatchDollar = TRUE)
options(warnPartialMatchArgs = TRUE)    
options(shiny.maxRequestSize = 999*1024^2)  ## max 999Mb upload

message("[MAIN] reading global.R ...")

##OPG     = "~/Playground/omicsplayground"
OPG     = ".."
RDIR    = file.path(OPG,"R")
FILES   = file.path(OPG,"lib")
FILESX  = file.path(OPG,"libx")
PGX.DIR = file.path(OPG,"data")

WATERMARK = FALSE
USER_MODE = "pro"
DEV     = FALSE
##DEV     = TRUE
DEBUG   = TRUE

## Determine if we are in ShinyProxy
SHINYPROXY = (Sys.getenv("SHINYPROXY_USERNAME")!="" && "omicsplayground" %in% dir("/"))
USERNAME = "anonymous"
if(SHINYPROXY) USERNAME = Sys.getenv("SHINYPROXY_USERNAME")

## dbg <- function(msg) if(DEBUG) message(cat(msg))
dbg <- function(...) {
    if(DEBUG) {
        msg <- list(...)
        msg <- paste(sapply(msg, function(s) paste(s,collapse=" ")),collapse=" ")
        message(msg)
    }
}

message("OPG =",OPG)
message("RDIR =",RDIR)
message("FILES =",FILES)
message("FILESX =",FILESX)
message("PGX.DIR =",PGX.DIR)
message("SHINYPROXY = ",SHINYPROXY)

src.local=TRUE  ## local or not-local, that's the question...
src.local=FALSE ## local or not-local, that's the question...
source(file.path(RDIR,"pgx-include.R"),local=src.local)    ## lots of libraries and source()
source(file.path(RDIR,"pgx-functions.R"), local=src.local) ## functions...
source(file.path(RDIR,"pgx-files.R"), local=src.local)     ## file functions
source(file.path(RDIR,"pgx-init.R"),local=src.local)       ## global variables

if(0) {
    save.image(file="../cache/image.RData")
    system.time( load(file="../cache/image.RData") )
}

message("\n")
message("************************************************")
message("************* parsing OPTIONS file *************")
message("************************************************")

if(!file.exists("OPTIONS")) stop("FATAL ERROR: cannot find OPTIONS file")
opt <- pgx.readOptions(file="OPTIONS")

## over-ride options (for DEBUGGING)
## opt$AUTHENTICATION = "none"
## opt$AUTHENTICATION = "password"
## opt$AUTHENTICATION = "register"
## opt$AUTHENTICATION = "firebase"

if(Sys.getenv("PLAYGROUND_AUTHENTICATION")!="") {
    auth <- Sys.getenv("PLAYGROUND_AUTHENTICATION")
    message("[ENV] overriding PLAYGROUND_AUTHENTICATION = ",auth)
    opt$AUTHENTICATION = auth
}

## copy to global environment
SHOW_QUESTIONS = FALSE
AUTHENTICATION = opt$AUTHENTICATION

DEV = (DEV && dir.exists("modulesx")) 
##DEV = FALSE
if(DEV) {
    message('****************** DEVELOPER MODE ********************')
}

## show options
message("\n",paste(paste(names(opt),"\t= ",sapply(opt,paste,collapse=" ")),collapse="\n"),"\n")

## --------------------------------------------------------------------
## ------------------------ READ FUNCTIONS ----------------------------
## --------------------------------------------------------------------


source("app-init.R", local=FALSE)
##pgx.initDatasetFolder(PGX.DIR, force=TRUE, verbose=1)
pgx.initDatasetFolder(PGX.DIR, force=FALSE, verbose=1)

if(0) {    
    ##PGX.DIR="../test/"
    ##pgx.initDatasetFolder(PGX.DIR, force=TRUE, verbose=1)    
    load("../data/geiger2016-arginine-test.pgx")
    load("../data/GSE10846-dlbcl-nc.pgx")
    load("../data/GSE22886-immune.pgx")
    load("../data/gtex-aging-n40svaNnm.pgx")
    ngs = pgx.initialize(ngs)
}

## --------------------------------------------------------------------
## ----------------- READ MODULES/BOARDS ------------------------------
## --------------------------------------------------------------------

source("modules/AuthenticationModule.R",local=src.local)
source("modules/ComputePgxModule.R",local=src.local)
source("modules/MakeContrastModule.R",local=src.local)
source("modules/NormalizeCountsModule.R",local=src.local)
source("modules/BatchCorrectModule.R",local=src.local)
source("modules/UploadModule.R",local=src.local)
##source("modules/UsersMapModule.R_")

BOARDS <- c("load","view","clust","expr","enrich","isect","func",
            "word","drug","sig","scell","cor","bio","cmap",
            "wgcna", "tcga","multi","system","qa","corsa","comp")
if(is.null(opt$BOARDS_ENABLED)) opt$BOARDS_ENABLED = BOARDS
if(is.null(opt$BOARDS_DISABLED)) opt$BOARDS_DISABLED = NA

ENABLED  <- array(BOARDS %in% opt$BOARDS_ENABLED, dimnames=list(BOARDS))
DISABLED <- array(BOARDS %in% opt$BOARDS_DISABLED, dimnames=list(BOARDS))
ENABLED  <- ENABLED & !DISABLED
ENABLED

boards <- dir("boards", pattern="Board.R$")
boards
for(m in boards) {
    message("[MAIN] loading board ",m)
    source(paste0("boards/",m), local=src.local)
    ##source(paste0("boards/",m), local=FALSE)
}

##ENABLED[c("wgcna","system","multi")] <- FALSE
ENABLED[c("system","multi","corsa")] <- FALSE
if(1 && DEV && dir.exists("modulesx")) {
    ## Very early development modules/boards
    ##
    xboards <- dir("modulesx", pattern="Board.R$")
    xboards
    m=xboards[1]
    for(m in xboards) {
        message("[MAIN] loading DEVELOPMENT board ",m)
        source(paste0("modulesx/",m), local=src.local)
    }
    ENABLED[] <- TRUE  ## enable all modules
    boards <- unique(c(boards, xboards))
}
ENABLED

## disable connectivity map if we have no signature database folder
has.sigdb <- length(dir(FILESX,pattern="sigdb.*h5")>0)
has.sigdb
if(has.sigdb==FALSE) ENABLED["cmap"] <- FALSE

MAINTABS = c("DataView","Clustering","Expression","Enrichment",
             "Signature","CellProfiling","DEV")

main.init_time <- round(Sys.time() - main.start_time,digits=4)
main.init_time
message("[MAIN] total init time = ",main.init_time," ",attr(main.init_time,"units"))