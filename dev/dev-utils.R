##appdir="~/R/golem/omics.app";wd="packages/board1"

get_appdir <- function() {
    wd <- strsplit(getwd(),split='/')[[1]]
    root.files <- c(".git",".dockerignore")
    chk <- rep(NA,length(wd))
    for(i in 1:length(wd)) {
        dd <- paste(wd[1:i],collapse='/')
        chk[i] <- all(root.files %in% dir(dd,all.files=TRUE))
    }
    if(!any(chk)) {
        return(NULL)
    } 
    paste(head(wd,tail(which(chk==TRUE),1)),collapse='/')
}

pkg.doc <- function(wd, appdir) {

  olddir <- getwd()
  wd <- file.path(appdir,wd)
  setwd(wd)
  getwd()
  ## usethis::proj_set(wd,force=TRUE);usethis::proj_get()  ## does not work!!! goes to root
  pkgload::pkg_path()
  
  ##usethis::create_package(path='.')
  getwd()
  ## usethis::use_package("shiny")  ## not correct, goes to app root!!

  ## Documentation --------------------------------------------------------
  ## ERROR!!! seems running out of disk space...
  ## devtools::build_vignettes()  ## builds vignettes, copy to doc
  devtools::document(wd)       ## builds Rd in man folder
  devtools::build_manual(pkg=wd, path=file.path(appdir,'doc'))   ## builds PDF reference manual
  roxygen2::roxygenise(wd)
  setwd(olddir)
}


pkg.build <- function(wd, appdir) {

  olddir <- getwd()
  wd <- file.path(appdir,wd)
  setwd(wd)
  getwd()
  ## usethis::proj_set(wd,force=TRUE);usethis::proj_get()  ## does not work!!! goes to root
  pkgload::pkg_path()
  
  ##usethis::create_package(path='.')
  getwd()
  ## usethis::use_package("shiny")  ## not correct, goes to app root!!

  ## Documentation --------------------------------------------------------
  devtools::document(wd)  ## builds NAMESPACE and Rd in man folder
  ##roxygen2::roxygenise()
  
  ## Deploy ---------------------------------------------------------------
  file.copy(file.path(appdir,'.Rbuildignore'),file.path(wd,".Rbuildignore"),overwrite=TRUE)
  devtools::build(pkg=wd, path=file.path(appdir,"packages"))
  ## devtools::install(pkg=wd, quick=FALSE)
  setwd(olddir)
}

