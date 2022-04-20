# Building a Prod-Ready, Robust Shiny Application.
#
# README: each step of the dev files is optional, and you don't have to
# fill every dev scripts before getting started.
# 01_start.R should be filled at start.
# 02_dev.R should be used to keep track of your development during the project.
# 03_deploy.R should be used once you need to deploy your app.
#
#
######################################
#### CURRENT FILE: DEPLOY SCRIPT #####
######################################


# Detach all loaded components and clean your environment
golem::detach_all_attached()
rm(list=ls(all.names = TRUE))

## set to package root!!
source("dev-utils.R")
appdir <- get_appdir() 
setwd(appdir)
appdir
wd="components/base"

remove.packages("omics.base")
remove.packages("omics.app")
remove.packages("omics.board.biomarker")

## Create ---------------------------------------------------------------

system("make clean.force")
pkg.build("components/base", appdir)
pkg.build("components/app", appdir)
pkg.build("components/board.biomarker", appdir)

setwd(appdir)
install.packages("components/omics.base_0.0.0.9000.tar.gz")
install.packages("components/omics.app_0.0.0.9000.tar.gz")
install.packages("components/omics.board.biomarker_0.0.0.9000.tar.gz")
