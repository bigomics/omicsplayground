
path='R';hdr.file="R/00Headers.R"
create_headers <- function(path='R', hdr.file="R/00Headers.R") {

    message("Creating header file: ",hdr.file)
    rfiles <- list.files(path,pattern='.[rR]$', recursive=TRUE, full.names=TRUE)
    ##hdr.file <- file.path(path,'R','00Headers.R')

    ## exclude some files
    rfiles <- setdiff(rfiles,hdr.file)
    excl.files <- c("global.R","00Headers.R","requirements.R",
                    "app.R","server.R","ui.R",
                    "pgx-init.R","pgx-include.R")
    rfiles <- grep(paste(excl.files,collapse="|"),rfiles,invert=TRUE,value=TRUE)
    
    write('# Generated automatically: do not edit by hand!\n', hdr.file, append=FALSE)    
    i=1
    for(i in 1:length(rfiles)) {
        code <- readLines(rfiles[i])
        nn <- grep('^#\'',code)
        nn <- unique(c(nn,nn+1))
        nn <- nn[which(nn >=1 & nn <= length(code))]
        hdr <- code[sort(nn)]
        hdr <- gsub("function\\(.*","function(){}",hdr)
        write(hdr, hdr.file, append=TRUE)
        src <- paste0("source(\'",rfiles[i],"\',local=TRUE)\n")
        write(src, hdr.file, append=TRUE)
    }
}

create_headers.SAVE <- function(path='.') {
    message("Creating header file R/00Headers.R")
    rfiles <- list.files(file.path(path,'R'),pattern='.[rR]$', recursive=TRUE, full.names=TRUE)
    hdr.file <- file.path(path,'R','00Headers.R')
    rfiles <- setdiff(rfiles,hdr.file)
    write('# Generated automatically: do not edit by hand', hdr.file, append=FALSE)    
    i=1
    for(i in 1:length(rfiles)) {
        code <- readLines(rfiles[i])
        nn <- grep('^[#\']',code)
        nn <- unique(c(nn,nn+1))
        nn <- nn[which(nn >=1 & nn <= length(code))]
        hdr <- code[sort(nn)]
        hdr <- gsub("function\\(.*","function(){}",hdr)
        write(hdr, hdr.file, append=TRUE)
        src <- paste0("source(\'",rfiles[i],"\',local=TRUE)\n")
        write(src, hdr.file, append=TRUE)
    }
}


