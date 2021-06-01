##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

## we may eventually migrate all OPTIONS into this file

message("[MAIN] reading global.R ...")

##OPG     = "~/Playground/omicsplayground"
OPG     = ".."
RDIR    = file.path(OPG,"R")
FILES   = file.path(OPG,"lib")
FILESX  = file.path(OPG,"libx")
PGX.DIR = file.path(OPG,"data")

USER_MODE = "pro"
DEV       = FALSE
WATERMARK = FALSE
DEBUG     = FALSE
DEBUG     = TRUE

## Determine if we are in ShinyProxy
SHINYPROXY = (Sys.getenv("SHINYPROXY_USERNAME")!="" && "omicsplayground" %in% dir("/"))
USERNAME = "anonymous"
if(SHINYPROXY) USERNAME = Sys.getenv("SHINYPROXY_USERNAME")

if(0) {
    TITLE           = "Omics Playground"
    AUTHENTICATION  = "none"
    ##AUTHENTICATION = password
    ENABLE_UPLOAD   = TRUE
    ENABLE_SAVE     = FALSE
    ENABLE_DELETE   = FALSE
    MAX_SAMPLES     = 20
    MAX_COMPARISONS = 5
    MAX_GENES       = 19999
    WATERMARK       = TRUE
    ##BOARDS_ENABLED = load,view,clust,expr,enrich,isect,func,word,drug,sig,scell,cor,bio,cmap
    ##BOARDS_ENABLED = load,view,clust,expr,cor,enrich,func
    BOARDS_DISABLED = c("system","multi")
    ##BOARDS_DISABLED = tcga,multi,sig,isect,bio
}

## dbg <- function(msg) if(DEBUG) message(cat(msg))
dbg <- function(...) {
    if(DEBUG) {
        msg <- list(...)
        msg <- paste(sapply(msg, function(s) paste(s,collapse=" ")),collapse=" ")
        message(msg)
    }
}

