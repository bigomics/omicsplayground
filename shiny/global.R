##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

## should we migrate all OPTIONS into this file??

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

