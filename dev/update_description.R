## Copy DESCRIPTION ----
## Add meta data about your application
##
## /!\ Note: if you want to change the name of your app during development,
## either re-run this function, call golem::set_golem_name(), or don't forget
## to change the name in the app_sys() function in app_config.R /!\
##

source("dev-utils.R")
appdir <- get_appdir() 
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


