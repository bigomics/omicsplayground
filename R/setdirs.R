
## Set folders in global environment
OPG       <<- pkgload::pkg_path()  ## search DESCRIPTION file
RDIR      <<- file.path(OPG,"R")
APPDIR    <<- file.path(OPG,"shiny")
FILES     <<- file.path(OPG,"lib")
FILESX    <<- file.path(OPG,"libx")
PGX.DIR   <<- file.path(OPG,"data")
SIGDB.DIR <<- file.path(OPG,"libx/sigdb")
