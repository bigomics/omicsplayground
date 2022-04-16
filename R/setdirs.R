get_pkg_root <- function() {
    pwd <- strsplit(getwd(),split='/')[[1]]
    paste(pwd[1:max(grep("omicsplayground",pwd))],collapse='/')
}

## Set folders
OPG       <<- get_pkg_root()
RDIR      <<- file.path(OPG,"R")
APPDIR    <<- file.path(OPG,"shiny")
FILES     <<- file.path(OPG,"lib")
FILESX    <<- file.path(OPG,"libx")
PGX.DIR   <<- file.path(OPG,"data")
SIGDB.DIR <<- file.path(OPG,"libx/sigdb")
