##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


# This file is supposed to run from the root Playground folder
if (basename(getwd()) != "omicsplayground") {
    stop("Please run from the OmicsPlayground root folder")
}

if(!require("renv")) install.packages("renv")
if(!require("remotes")) install.packages("remotes")
if(!require("devtools")) install.packages("devtools")

source("dev/functions.R")
pkg <- scan_packages(path='components') 

pkg.remotes <- pkg$remotes
pkg.imports <- pkg$imports

## determine packages in playbase
playbase.imports <- strsplit(packageDescription("playbase")$Imports,split=',')[[1]]
playbase.imports <- trimws(gsub("[\n]","",playbase.imports))
playbase.remotes <- strsplit(packageDescription("playbase")$Remotes,split=',')[[1]]
playbase.remotes <- gsub("[\n]","",playbase.remotes)
playbase.remotes <- gsub(".*github.com/|/archive.*","",playbase.remotes)
playbase.remotes <- gsub(".*contrib/|_[0-9].*$","",playbase.remotes)
playbase.remotes <- gsub(".*/|@.*","",playbase.remotes)
playbase.remotes <- trimws(playbase.remotes)
playbase.pkg <- sort(unique(c(playbase.imports, playbase.remotes)))

bigomics.pkg <- c("playbase","playdata","bigdash","bigLoaders","wizardR")

## we only report packages that are not in playbase
pkg.imports <- setdiff(pkg.imports, playbase.pkg)

sel1 <- !names(pkg.remotes) %in% playbase.pkg
sel2 <- names(pkg.remotes) %in% bigomics.pkg
pkg.remotes <- sort(pkg.remotes[sel1 | sel2])

pkg.imports
names(pkg.remotes)

if(file.exists("DESCRIPTION")) {
  hash <- substring(tempfile(),nchar(tempfile())-3, nchar(tempfile()))
  file.copy("DESCRIPTION",paste0("DESCRIPTION.save.",hash))
  message("WARNING: overwriting existing DESCRIPTION file.")
}

desc.header <- readLines("dev/description.app")
desc.file <- "DESCRIPTION"
write(desc.header, file=desc.file)
write("Imports:", file=desc.file, append=TRUE)
write(paste0("    ",pkg.imports,","), file=desc.file, append=TRUE)
write("Remotes:", file=desc.file, append=TRUE)
write(paste0("    ",pkg.remotes,","), file=desc.file, append=TRUE)



