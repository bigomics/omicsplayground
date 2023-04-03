mem.proc <- function(digits=0) {
  mem <- "[? MB]"
  if(Sys.info()["sysname"] %in% c("Linux")) {
    file <- paste("/proc", Sys.getpid(), "stat", sep = "/")
    what <- vector("character", 52)
    ## In your logging routine
    vsz <- as.numeric(scan(file, what = what, quiet = TRUE)[23])
    vsz <- vsz / (1024**2) ## MB
    ##cat("Virtual size: ", vsz, " MB\n", sep = "")
    mem <- paste0(round(vsz,digits),"MB")
  }
  mem
}

info <- function(..., type="INFO") {
  dd <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg = "some message"
  msg = sapply( list(...),paste,collapse=" ")
  dd <- paste0("[",dd,"]")
  mm <- paste0("[",mem.proc(),"]")
  type <- paste0("[",type,"]")
  message(paste0(type,dd,mm," --- ",sub("\n$","",paste(msg,collapse=" "))))
}

dbg <- function(...) info(..., type="DBUG")
