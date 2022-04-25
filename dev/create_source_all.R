
## Create global "source_all" file

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
  rfiles <- setdiff(rfiles, out.file)    
  excl.files <- c(excl.files,"global.R","app.R","server.R","ui.R")
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


if(1) {

    source("dev/dev-utils.R")
    appdir <- get_appdir() 
    setwd(appdir)
    appdir
    
    setwd(appdir)
    create_SourceAll('components',add.comments=FALSE)
    
}
