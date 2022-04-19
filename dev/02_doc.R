# Building a Prod-Ready, Robust Shiny Application.
#
# README: each step of the dev files is optional, and you don't have to
# fill every dev scripts before getting started.
# 01_start.R should be filled at start.
# 02_dev.R should be used to keep track of your development during the project.
# 03_deploy.R should be used once you need to deploy your app.
#
#
###################################
#### CURRENT FILE: DEV SCRIPT #####
###################################

# Detach all loaded components and clean your environment
golem::detach_all_attached()
rm(list=ls(all.names = TRUE))

## Create documentaion ----------------------------------
source("dev-utils.R")
appdir <- get_appdir() 
setwd(appdir)
appdir
wd="components/base"

system("make clean.force")
pkg.doc("components/base", appdir)
pkg.doc("components/app", appdir)
pkg.doc("components/board.biomarker", appdir)

