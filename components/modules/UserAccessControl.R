##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

ACCESS_LOGFILE = file.path(ETC,"access.log")
## unlink(ACCESS_LOGFILE)

pgx.record_access <- function(user, action, session=session,
                              access.file=ACCESS_LOGFILE) {

  if(is.null(user) || is.null(action)) return(NULL)
  if(length(user)==0 || length(action)==0) return(NULL)
  if(is.na(user) || is.na(action)) return(NULL)
  if(user=='' || action=='') return(NULL)    

  time.now <- as.POSIXct(Sys.time())
  public.ip <- system("curl -s http://api.ipify.org",intern=TRUE)  
#  public.ip <- system("curl -s http://ipinfo.io",intern=TRUE)    
#  public.ip <- stringr::str_extract_all( public.ip[2], '[0-9][0-9.]*[0-9]')[[1]]
  hostname <- system("hostname",intern=TRUE)  
  session_id <- substring(session$token,1,8)
  remote.addr <- session$request$REMOTE_ADDR
  
  dbg("[pgx.record_access] session_id = ",session_id)
  dbg("[pgx.record_access] hostname = ",hostname)
  dbg("[pgx.record_access] public.ip = ",public.ip)
  dbg("[pgx.record_access] remote.addr = ",remote.addr)    
  dbg("[pgx.record_access] time = ",time.now)        

  login_data <- data.frame(
    user = user,
    session = session_id,
    action = action,
    time = time.now,
    host = hostname,
    public = public.ip,
    client = remote.addr
  )
  do.append <- file.exists(access.file)
  data.table::fwrite(login_data, file=access.file, quote=TRUE, append=do.append)
}

## pgx.record_access(user=user[1], action='logout.stale', session_id=user[2])      

FolderLock <- R6::R6Class("FolderLock",
  private = list(
    heartbeat_secs = 15*1000
  ),
  active = list(
    heartbeat = function(secs) {
      if (missing(secs)) {
        self$heartbeat_secs
      } else {
        self$heartbeat_secs <- secs
      }
    }
  ),
  public = list(
    is_locked = FALSE,
    path = NULL,
    user = NULL,
    max_idle = 60,
    #' Initialize the R6 Object
    #'
    #' @param rds_path The path to the rds file.
    #'
    initialize = function(max_idle=60 ) {
      self$max_idle <- max_idle
      invisible(self)
    },
    set_user = function(user, path='.') {
      self$path <- path
      self$user <- user
    },
    remove_all_locks = function() {
      other_lock_files <- dir(self$path,"^LOCK__.*",full.name=TRUE)
      if(length(other_lock_files)>0) {
        lapply( other_lock_files, file.remove)
      }
    },
    reset = function() {      
      self$path <- NULL
      self$user <- NULL
      is_locked = FALSE
    },
    lockfile = function(full.name=FALSE) {
      if(is.null(self$path)) return(NULL)
      f <- paste0("LOCK__",self$user)
      if(full.name) f <- file.path( self$path, f)
      return(f)
    },
    read_lock = function() {
      if(is.null(self$path)) return(NULL)      
      lock_file <- dir(self$path,"^LOCK__.*",full.name=FALSE)
      lock_file
      message("[pgx.read_lock] lock_file = ", lock_file)
      
      if(length(lock_file)==0) {
        message("UNLOCKED: no lock file")
        return(NULL)
      }
      if(length(lock_file)>1) {
        message("WARNING: multiple lock files in folder!")
        mtimes <- sapply(lock_file, function(f) file.mtime(file.path(self$path,f)))
        most_recent <- which.max(mtimes)
        ## remove older lock files
        for(i in setdiff(1:length(lock_file),most_recent)) {
          file.remove(file.path(self$path,lock_file[i]))
        }
        lock_file <- lock_file[most_recent]
      }

      ## ok LOCK file exists
      lock_user <- strsplit(lock_file, split='LOCK__')[[1]][2]
      lock_time <- file.mtime(file.path(self$path,lock_file))
      delta_secs <- (Sys.time() - as.POSIXct(lock_time))
      delta_secs <- round(as.numeric(delta_secs, units='secs'),digits=2)
      
      ## compute status here?
      is_locked <- delta_secs < self$max_idle  
      
      info <- list(
        is_locked = is_locked,
        user = lock_user,
        time = lock_time,
        delta_secs = delta_secs,
        file = lock_file,
        path = self$path
      )
      return(info)
    },
    remove_lock = function() {
      if(is.null(self$path)) return(NULL)
      file.remove(self$lockfile())
      self$is_locked = FALSE
    },
    write_lock = function(force=FALSE) {
      if(is.null(self$path)) return(NULL)
      other_lock_files <- dir(self$path,"^LOCK__.*",full.name=TRUE)
      if(length(other_lock_files)>0) {
        lapply( other_lock_files, file.remove)
      }
      mylock_file <- self$lockfile(full.name=TRUE)
      write(NULL, mylock_file)
      self$is_locked = TRUE
    },
    #' Show shinyalert popup message that the lock has been succesful
    #' 
    shinyalert_unlocked = function(lock) {
      id <- strsplit(self$user,split="__")[[1]]
      shinyalert::shinyalert(
        title = "SUCCESS!",
        text = paste(
          "successfully locked by you",
          "<br><br>name =",id[1],                            
          "<br>session =",id[2],
          "<br>lock_time =",lock$time
        ),
        ##closeOnEsc = FALSE, showConfirmButton = FALSE,
        animation = FALSE,
        html = TRUE, immediate=TRUE
      )
    },
    #' Show shinyalert popup message that the folder is locked by
    #' someone else
    shinyalert_locked = function(lock, session) {
      
      id <- strsplit(lock$user,split="__")[[1]]
      my_id <- strsplit(self$user,split="__")[[1]]

      msg.text <- paste(
        "This account is locked by someone else. Please try again later.",
        "<br><br>name =",id[1],                            
        "<br>session =",id[2],
        "<br>lock_time =",lock$time,
#        "<br>delta =",paste0(lock$delta_secs,"sec"),
#        "<br><br>your name =",my_id[1],                            
#        "<br>your session =",my_id[2],
        NULL
      )

      shinyalert::shinyalert(
        title = "LOCKED!",
        text = msg.text,
        confirmButtonText = "Close window",
        callbackR = function(x) { if(x) session$close() },
        #callbackJS = "function(x) { if(x) {setTimeout(function(){window.close();},500);} }",
        callbackJS = "function(x) { if(x) { window.close() }}",              
        ##closeOnEsc = FALSE, 
        animation = FALSE,
        html = TRUE, immediate=TRUE            
      )
    },
    start_shiny_observer = function(auth, access_id, session, poll_secs=15) {
      
      observe({
        if(auth$logged) {
          shiny::req(auth$user_dir)

          self$set_user(user=access_id(), path=auth$user_dir)          
          cur <- self$read_lock()

          if(is.null(cur) || !cur$is_locked) {
            pgx.record_access(self$user, "login", session=session)
            cur <- list( user= self$user, is_locked = TRUE)
          }

          is_mylock <- cur$user == access_id()           
          if( is_mylock )  {
            self$write_lock()
            if(is.null(cur$time)) cur <- self$read_lock()
            self$shinyalert_unlocked(lock=cur)
            invalidateLater(poll_secs*1000)
          } else {
            self$shinyalert_locked( lock = cur, session )
            invalidateLater(Inf)
          }
        } else {
          if(is.null(self$user)) return(NULL)
          user <- sub("__.*","",self$user)  ## strip postfix
          pgx.record_access(user, "logout", session=session)
          self$remove_lock()
          self$reset()
        }
      })
    } 
  )  
) ## end of R6 class
    


