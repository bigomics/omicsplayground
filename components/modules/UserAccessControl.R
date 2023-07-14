##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

ACCESS_LOGFILE = file.path(ETC,"access.log")
## unlink(ACCESS_LOGFILE)

record_access <- function(user, action, host='', client='', session_id='',
                          access.file=ACCESS_LOGFILE) {
  time.now <- as.POSIXct(Sys.time())
  client.ip <- system("curl -s http://ipinfo.io",intern=TRUE)
  client.ip <- stringr::str_extract_all( client.ip[2], '[0-9][0-9.]*[0-9]')[[1]]
  hostname <- system("cat /etc/hostname",intern=TRUE)  
  login_data <- data.frame(
    user=user, session=substring(session_id,1,8),
    action=action, time=time.now,
    host=hostname, client=client.ip)
  do.append <- file.exists(access.file)
  data.table::fwrite(login_data, file=access.file, quote=TRUE, append=do.append)
}

path="~/Playground/omicsplayground/data/ivo.kwee@bigomics.ch"

read_lock <- function(path) {  
  
  lock_file <- dir(path,"^LOCK__.*",full.name=FALSE)
  lock_file
  dbg("[read_lock] lock_file = ", lock_file)
  
  if(length(lock_file)==0) {
    message("UNLOCKED: no lock file")
    return(NULL)
  }

  if(length(lock_file)>1) {
    message("WARNING: multiple lock files in folder!")
    mtimes <- sapply(lock_file, function(f) file.mtime(file.path(path,f)))
    most_recent <- which.max(mtimes)
    ## remove older lock files
    for(i in setdiff(1:length(lock_file),most_recent)) {
      file.remove(file.path(path,lock_file[i]))
    }
    lock_file <- lock_file[most_recent]
  }

  ## ok LOCK file exists
  lock_user <- strsplit(lock_file, split='LOCK__')[[1]][2]
  lock_time <- file.mtime(file.path(path,lock_file))
  
  ## if the lock is older than say 5minutes, release it
  delta <- (Sys.time() - as.POSIXct(lock_time))
  delta <- round(as.numeric(delta, units='secs'),digits=2)

  info <- list(
    user = lock_user,
    time = lock_time,
    delta = delta,
    file = lock_file,
    path = path
  )
  info
  return(info)
}

write_lock <- function(user, path, force=FALSE, update=TRUE, max_idle=30) {

  dbg("[write_lock] user = ", user)
  dbg("[write_lock] path = ", path)
  
  lock <- read_lock(path=path)
  lock
  
  has.lockfile <- !is.null(lock)
  mylock_file <- paste0("LOCK__",user)
  mylock_file

  ## determine is the lock is stale
  is_stale <- FALSE
  if(!is.null(lock)) {
    is_stale <- lock$delta > max_idle
  }
  
  has.lockfile  
  lock$file
  mylock_file
  
  if(!has.lockfile || is_stale || force) {
    ## remove any previous locks
    if(is_stale) {
      dbg("[write_lock] STALE LOCK by ",lock$user)
      user <- strsplit(lock$user,split='__')[[1]]
      record_access(user=user[1], action='logout.stale', session_id=user[2])      
      file.remove(file.path(path,lock$file))
    }
    if(!has.lockfile) {
      dbg("[write_lock] no lock file")
    }
    write(NULL, file.path(path, mylock_file))
    lock$status <- TRUE  ## successful own lock
    return(lock)
  }
  
  ## if there is a lock file but it is from the current user, then
  ## update the time stamp of file (touch or rewrite)
  is_mylock <- (has.lockfile && lock$file == mylock_file)
  is_mylock
  if(has.lockfile && is_mylock) {
    dbg("[write_lock] UNLOCKED: because you are lock owner: ", lock$user)
    if(update) {
      write(NULL, file.path(path, mylock_file))
      ## system(paste("touch", mylock_file))
    }
    lock$status <- TRUE  ## successful own lock    
    return(lock)
  }
  
  ## if session id and user are different then return as locked
  if(has.lockfile && !is_mylock) {
    dbg("[write_lock] LOCKED by ",lock$user)
    lock$status <- FALSE  ## failed to get lock
    return(lock)
  } 

  ## should not come here
  dbg("[write_lock] WARNING: locked for unknown reason!!")
  lock$status <- FALSE  ## failed to get lock
  return(lock)
}

## write_lock <- function(user, path) {
##   mylock_file <- file.path(path,paste0("LOCK__",user))
##   if(!file.exists(mylock_file)) {
##     write_lock(user=user, path=path, force=TRUE) 
##   } else {
##     system(paste("touch", mylock_file))
##   }
## }

remove_lock <- function(user, path) {
  mylock_file <- file.path(path,paste0("LOCK__",user))
  file.remove(mylock_file)
}
