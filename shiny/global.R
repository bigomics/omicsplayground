##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

## we may eventually migrate all OPTIONS into this file

message("[] reading global.R ...")

##OPG     = "~/Playground/omicsplayground"
OPG     = ".."
RDIR    = file.path(OPG,"R")
FILES   = file.path(OPG,"lib")
FILESX  = file.path(OPG,"libx")
PGX.DIR = file.path(OPG,"data")

USER_MODE = "pro"
DEV       = FALSE
WATERMARK = FALSE
DEBUG     = TRUE
