# Building a Prod-Ready, Robust Shiny Application.
#
# README: each step of the dev files is optional, and you don't have to
# fill every dev scripts before getting started.
# 01_start.R should be filled at start.
# 02_dev.R should be used to keep track of your development during the project.
# 03_deploy.R should be used once you need to deploy your app.
#
#
########################################
#### CURRENT FILE: ON START SCRIPT #####
########################################


## Copy DESCRIPTION ----
## Add meta data about your application
##
## /!\ Note: if you want to change the name of your app during development,
## either re-run this function, call golem::set_golem_name(), or don't forget
## to change the name in the app_sys() function in app_config.R /!\
##
golem::detach_all_attached()
rm(list=ls(all.names = TRUE))

source("dev-utils.R")
appdir <- get_app_root() 
setwd(appdir)
appdir

dir.pkg <- dir("components", full.names=TRUE)
dir.pkg

desc <- readLines("dev/description.component")

for(d in dir.pkg) {
  setwd(file.path(appdir,d))
  desc[1] <- paste("Package:", paste0("omics.",basename(d)))
  desc[2] <- paste("Title: OmicsPlayground",toupper(basename(d)),"component package")  
  write( desc, file="DESCRIPTION")  
}


## Create global "source_all" file
setwd(appdir)
create_SourceAll('components',add.comments=FALSE)
##create_headers('R', add.source=FALSE)
##create_allcode('R') 
source("components/00SourceAll.R")
source("components/00SourceAll.R",chdir=TRUE)
