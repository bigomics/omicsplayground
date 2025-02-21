##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## ACCESS_LOGFILE = file.path(ETC,"access.log")
## unlink(ACCESS_LOGFILE)

pgx.record_access <- function(user,
                              action,
                              comment = "",
                              comment2 = "",
                              comment3 = "",
                              comment4 = "",
                              session = session,
                              num_datasets = "",
                              time = Sys.time(),
                              access.file = ACCESS_LOGFILE,
                              ip = "") {
  if (is.null(user) || is.null(action)) {
    return(NULL)
  }
  if (length(user) == 0 || length(action) == 0) {
    return(NULL)
  }
  if (is.na(user) || is.na(action)) {
    return(NULL)
  }
  if (user == "" || action == "") {
    return(NULL)
  }
  if (is.null(ip)) {
    ip <- ""
  }

  user <- sub("__.*", "", user) ## strip postfix
  time <- as.POSIXct(time)
  session_id <- substring(session$token, 1, 16)
  hostname <- opt$HOSTNAME
  remote_addr <- session$request$REMOTE_ADDR

  if (1) {
    dbg("[pgx.record_access] action = ", action)
    dbg("[pgx.record_access] user = ", user)
    dbg("[pgx.record_access] session_id = ", session_id)
    dbg("[pgx.record_access] hostname = ", hostname)
    dbg("[pgx.record_access] remote_addr = ", remote_addr)
    dbg("[pgx.record_access] time = ", time)
  }

  login_data <- data.frame(
    user = user,
    action = action,
    time = time,
    session = session_id,
    hostname = hostname,
    client.ip = remote_addr,
    comment = comment,
    comment2 = comment2,
    comment3 = comment3,
    comment4 = comment4,
    num_datasets = num_datasets,
    ip = ip
  )
  do.append <- file.exists(access.file)
  data.table::fwrite(login_data, file = access.file, quote = TRUE, append = do.append)
}

## pgx.record_access(user=user[1], action='logout.stale', session_id=user[2])

FolderLock <- R6::R6Class("FolderLock",
  private = list(
    poll_secs = 15,
    max_idle = 60,
    show_success = FALSE,
    show_details = FALSE
  ),
  active = list(),
  public = list(
    path = NULL,
    user = NULL,
    user_id = NULL,
    #' Initialize the R6 Object
    #'
    #' @param rds_path The path to the rds file.
    #'
    initialize = function(poll_secs = 15,
                          max_idle = 60,
                          show_success = FALSE,
                          show_details = FALSE) {
      private$max_idle <- max_idle
      private$poll_secs <- poll_secs
      private$show_success <- show_success
      private$show_details <- show_details
      invisible(self)
    },
    set_user = function(user, session, path) {
      self$path <- path
      self$user <- paste0(user, "__", substring(session$token, 1, 16))
      self$user_id <- user
    },
    remove_all_locks = function() {
      other_lock_files <- dir(self$path, "^LOCK__.*", full.names = TRUE)
      if (length(other_lock_files) > 0) {
        lapply(other_lock_files, file.remove)
      }
    },
    reset = function() {
      self$path <- NULL
      self$user <- NULL
    },
    lockfile = function(full.path = FALSE) {
      if (is.null(self$path)) {
        return(NULL)
      }
      f <- paste0("LOCK__", self$user)
      if (full.path) f <- file.path(self$path, f)
      return(f)
    },
    read_lock = function() {
      if (is.null(self$path)) {
        return(NULL)
      }
      lock_file <- dir(self$path, paste0("^LOCK__", self$user_id, ".*"), full.name = FALSE)
      if (length(lock_file) == 0) {
        message("UNLOCKED: no lock file")
        return(NULL)
      }
      if (length(lock_file) > 1) {
        message("WARNING: multiple lock files in folder!")
        mtimes <- sapply(lock_file, function(f) file.mtime(file.path(self$path, f)))
        most_recent <- which.max(mtimes)
        ## remove older lock files
        for (i in setdiff(1:length(lock_file), most_recent)) {
          file.remove(file.path(self$path, lock_file[i]))
        }
        lock_file <- lock_file[most_recent]
      }

      ## ok LOCK file exists
      lock_user <- strsplit(lock_file, split = "LOCK__")[[1]][2]
      lock_time <- file.mtime(file.path(self$path, lock_file))
      delta_secs <- (Sys.time() - as.POSIXct(lock_time))
      delta_secs <- round(as.numeric(delta_secs, units = "secs"), digits = 2)

      ## compute status here?
      is_locked <- delta_secs < private$max_idle

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
      if (is.null(self$path)) {
        return(NULL)
      }
      f <- self$lockfile(full.path = TRUE)
      if (file.exists(f)) {
        dbg("[FileLock] removing lockfile", self$lockfile())
        file.remove(f)
      }
    },
    write_lock = function(force = FALSE) {
      if (is.null(self$path)) {
        return(NULL)
      }
      other_lock_files <- dir(self$path, paste0("^LOCK__", self$user_id, ".*"), full.names = TRUE)
      if (length(other_lock_files) > 0) {
        lapply(other_lock_files, file.remove)
      }
      mylock_file <- self$lockfile(full.path = TRUE)
      write(NULL, mylock_file)
    },
    #' Show shinyalert popup message that the lock has been succesful
    #'
    shinyalert_success = function(lock) {
      id <- strsplit(self$user, split = "__")[[1]]

      msg.text <- paste(
        "successfully locked by you",
        "<br><br>name =", id[1],
        "<br>session =", id[2],
        "<br>lock_time =", lock$time,
        "<br>delta =", paste0(lock$delta_secs, "sec")
      )

      shinyalert::shinyalert(
        title = "SUCCESS!",
        text = HTML(msg.text),
        ## closeOnEsc = FALSE, showConfirmButton = FALSE,
        animation = FALSE,
        html = TRUE,
        immediate = TRUE
      )
    },
    #' Show shinyalert popup message that the folder is locked by
    #' someone else
    shinyalert_locked = function(lock, session) {
      id <- strsplit(lock$user, split = "__")[[1]]
      my_id <- strsplit(self$user, split = "__")[[1]]

      msg.text <- paste(
        "This account is locked by someone else. <br>Please try again later.",
        "<br><br>name =", id[1],
        "<br>session =", id[2]
      )

      if (private$show_details) {
        msg.text <- paste(
          msg.text,
          "<br>lock_time =", lock$time,
          "<br>delta =", paste0(lock$delta_secs, "sec"),
          "<br><br>your name =", my_id[1],
          "<br>your session =", my_id[2]
        )
      }

      shinyalert::shinyalert(
        title = "LOCKED!",
        text = HTML(msg.text),
        animation = FALSE,
        html = TRUE,
        immediate = TRUE,
        confirmButtonText = "Close window",
        callbackR = function(x) {
          ## if (x) session$close()
          shinyalert::shinyalert(
            title = "Close this session?",
            text = "Please log out from any other browser or browser tab. Please try again later.",
            confirmButtonText = "Close session",
            callbackR = function(x) {
              dbg("[FolderLock] call closing session. x= ", x)
              if (x) {
                dbg("[FolderLock] closing session")
                session$close()
              }
            },
            animation = FALSE,
            html = TRUE,
            immediate = TRUE
          )
        }
      )
    },
    start_shiny_observer = function(auth, session) {
      observe({
        if (auth$logged) {
          shiny::req(auth$user_dir)

          no_email <- is.null(auth$email) || is.na(auth$email) || auth$email == ""
          user_id <- ifelse(no_email, auth$username, auth$email)
          self$set_user(user = user_id, session, path = auth$user_dir)
          cur <- self$read_lock()

          if (is.null(cur) || !cur$is_locked) {
            ## this is probably a clean login, either with no lockfile
            ## present or the lockfile is stale. Note: normal login is
            ## recorded in server.R
            if (!is.null(cur) && !cur$is_locked) {
              ## lockfile is stale. NEED RETHINK!! should record the
              ## correct host/client IP and time.
              pgx.record_access(cur$user, action = "stale.logout", time = cur$time, session = session)
            }
            cur <- list(user = self$user, is_locked = TRUE)
          }

          is_mylock <- (cur$user == self$user)
          if (is_mylock) {
            self$write_lock()
            if (private$show_success) {
              if (is.null(cur$time)) cur <- self$read_lock()
              self$shinyalert_success(lock = cur)
            }
            invalidateLater(private$poll_secs * 1000)
          } else {
            pgx.record_access(
              user = self$user,
              action = "locked.login",
              comment = "login attempt while locked by other user",
              session = session
            )
            self$shinyalert_locked(lock = cur, session)
            invalidateLater(Inf)
          }
        } else {
          ## !auth$logged
          if (is.null(self$user)) {
            invalidateLater(Inf) ## really needed?
            return(NULL)
          }
          ## at logout we record the logout action and remove the
          ## lockfile, then reset the user/path
          self$remove_lock()
          self$reset()
          invalidateLater(Inf) ## really needed?
        }
      })
    }
  )
) ## end of R6 class


pgx.start_heartbeat <- function(auth, session, online_dir, delta = 60) {
  reactive({
    ## shiny::req(auth$email)
    user <- auth$email

    session_id <- substring(session$token, 1, 16)
    hostname <- opt$HOSTNAME
    if (is.null(hostname) || hostname == "") {
      hostname <- system("hostname", intern = TRUE)
    }
    hostname2 <- isolate(session$clientData$url_hostname)
    ip <- session$request$REMOTE_ADDR
    online_id <- paste0("ONLINE__", user, "__", session_id, "__", hostname, ":", ip)
    online_file <- file.path(online_dir, online_id)

    if (auth$logged) {
      if (!dir.exists(online_dir)) dir.create(online_dir)
      write(NULL, file = online_file)
    } else {
      if (file.exists(online_file)) {
        file.remove(online_file)
      }
    }

    ## remove old ONLINE files
    online_tags <- dir(online_dir, pattern = "^ONLINE", full.name = FALSE)
    files <- file.path(online_dir, online_tags)
    if (length(files) > 0) {
      mtimes <- sapply(files, function(f) as.POSIXct(file.mtime(f)))
      lapsed <- (Sys.time() - mtimes)
      lapsed <- round(as.numeric(lapsed, units = "secs"), digits = 2)
      is_stale <- which(lapsed > 3 * delta)
      if (length(is_stale) > 0) {
        file.remove(files[is_stale])
      }
    }

    ## invalidate-rinse-repeat
    ##    dbg("[pgx.start_heartbeat] ticking...")
    invalidateLater(delta * 1000)

    ## return filename
    invisible(online_file)
  })
}
