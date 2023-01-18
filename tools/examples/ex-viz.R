##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

library(shiny)

## RUN FROM root folder!
##setwd(pkgload::pkg_path())
setwd("~/Playground/omicsplayground")
source("components/00SourceAll.R",chdir=TRUE)  ## global variable

load("data/example-data.pgx",verbose=1)
pgx=ngs;remove(ngs)
viz.PhenoMaps(pgx)
