
##path='R';hdr.file="R/00Headers.R";excl.files=NULL
create_headers <- function(path='R', add.header=TRUE, add.source=TRUE, excl.files=NULL) {
    
    out.file = file.path(path[1],"00Headers.R")
    message("Creating header file: ",out.file)
    rfiles <- list.files(path, pattern='.[rR]$', recursive=TRUE, full.names=TRUE)

    ## exclude some files
    rfiles <- setdiff(rfiles,out.file)    
    excl.files <- c(excl.files,"global.R","app.R")   ## skip 'executable' files
    rfiles <- rfiles[!basename(rfiles) %in% excl.files]

    rfiles0 <- list.files('R', pattern='.[rR]$', recursive=FALSE, full.names=FALSE)
    
    write('## Generated automatically: do not edit by hand\n', out.file, append=FALSE)    
    write("## find package root folder", out.file, append=TRUE)
    write("PKG <<- pkgload::pkg_path()\n", out.file, append=TRUE)
    write("message('package root = ',PKG)\n", out.file, append=TRUE)            
    
    i=1
    for(i in 1:length(rfiles)) {

        ## if the file is in the R folder, we don't need roxygen comments
        if(add.header && !(basename(rfiles[i]) %in% rfiles0)) {
            code <- readLines(rfiles[i])
            code <- code[grep("^$|^[ ]*$",code,invert=TRUE)]  ## drop empty lines
            code <- code[grep("^#[a-z,A_Z #]",code,invert=TRUE)]  ## drop comments
            nn <- grep('^#\'',code)
            if(length(nn)) {
                nn <- unique(c(nn,nn+1))
                nn <- nn[which(nn >=1 & nn <= length(code))]
                hdr <- code[sort(nn)]
                hdr <- gsub("function\\(.*","function(){}",hdr)
                hdr <- c("",hdr)
                write(hdr, out.file, append=TRUE)
            }
        }
        if(add.source) {
            fn  <- paste0("file.path(PKG,\'",rfiles[i],"\')")
            src <- paste0("source(",fn,",local=TRUE)")
            write(src, out.file, append=TRUE)
        }
    }
}

create_allcode <- function(path='R') {

    out.file = file.path(path[1],"00AllCode.R")
    message("Creating all code file: ",out.file)
    rfiles <- list.files(path,pattern='.[rR]$', recursive=TRUE, full.names=TRUE)

    ## exclude some files
    rfiles <- setdiff(rfiles,out.file)    
    excl.files <- c("global.R","app.R")   ## skip 'executable' files
    rfiles <- rfiles[!basename(rfiles) %in% excl.files]
    
    write('# Generated automatically: do not edit by hand', out.file, append=FALSE)    
    i=1
    for(i in 1:length(rfiles)) {
        code <- readLines(rfiles[i])
        code <- code[grep("^$|^[ ]*$",code,invert=TRUE)]  ## drop empty lines
        code <- code[grep("^#[a-z,A_Z #]",code,invert=TRUE)]  ## drop comments        
        write(paste("\n## --------------------- ",rfiles[i]," -----------------------"),
              out.file, append=TRUE)
        write(code, out.file, append=TRUE)        
    }
}
