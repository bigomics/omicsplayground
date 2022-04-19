##appdir="~/R/golem/omics.app";wd="packages/board1"

get_app_root <- function() {
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


##path='R';hdr.file="R/00Headers.R";excl.files=NULL
create_SourceAll <- function(path='R',
                             add.comments=TRUE,
                             add.source=TRUE,
                             excl.files=NULL)
{
  
  out.file = file.path(path[1],"00SourceAll.R")
  message("Creating SourceAll file: ",out.file)
  rfiles <- list.files(path, pattern='.[rR]$', recursive=TRUE, full.names=TRUE)

  ## exclude some files
  rfiles <- setdiff(rfiles,out.file)    
  excl.files <- c(excl.files,"global.R","app.R")   ## skip 'executable' files
  rfiles <- rfiles[!basename(rfiles) %in% excl.files]

  rfiles0 <- list.files('R', pattern='.[rR]$', recursive=FALSE, full.names=FALSE)
  
  write('## Generated automatically: do not edit by hand\n', out.file, append=FALSE)    

  write("message('source all called from wd = ',getwd())", out.file, append=TRUE)            
  write(paste("if(!file.exists('00SourceAll.R')) {"), out.file, append=TRUE)        
  write("  message('WARNING: not in source folder. skipping.')", out.file, append=TRUE)
  write("} else {", out.file, append=TRUE)        
  write("  message('Note: sourcing all code...')", out.file, append=TRUE)
  i=1
  for(i in 1:length(rfiles)) {

    ## if the file is in the R folder, we don't need roxygen comments
    if(add.comments && !(basename(rfiles[i]) %in% rfiles0)) {
      code <- readLines(rfiles[i])
      code <- code[grep("^$|^[ ]*$",code,invert=TRUE)]  ## drop empty lines
      code <- code[grep("^#[a-z,A_Z #]",code,invert=TRUE)]  ## drop comments
      nn <- grep('^#\'',code)
      if(length(nn)) {
        nn <- unique(c(nn,nn+1))
        nn <- nn[which(nn >=1 & nn <= length(code))]
        hdr <- code[sort(nn)]
        hdr <- gsub("function\\(.*","function(){}",hdr)
        hdr <- paste0("  ",hdr)
        hdr <- c("",hdr)
        write(hdr, out.file, append=TRUE)
      }
    }
    if(add.source) {
      rm.path <- paste0("^",path,"|^",path,"/")
      fn  <- paste0("\'",sub(rm.path,"",rfiles[i]),"\'")      
      src <- paste0("  source(",fn,",encoding='UTF-8')")
      write(src, out.file, append=TRUE)
    }
  }
  write("  message('done! (sourcing all code)')", out.file, append=TRUE)    
  write("}", out.file, append=TRUE)      
}
