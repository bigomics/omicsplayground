##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#########################################################################
##                                                                     ##
##              Utility Functions for Omics Playground                 ##
##                                                                     ##
#########################################################################

is.symlink <- function(paths) isTRUE(nzchar(Sys.readlink(paths), keepNA = TRUE))

getAppVersion <- function(add.auth = FALSE) {
  version <- scan(file.path(OPG, "VERSION"), character())[1]
  if (add.auth) {
    auth.method <- paste0(opt$AUTHENTICATION, ifelse(opt$USE_CREDENTIALS, "+cred", ""))
    version <- paste0(version, "/", auth.method)
  }
  version
}

getFirstName <- function(name) {
  first.name <- strsplit(name, split = "[@ .]")[[1]][1]
  first.name <- paste0(
    toupper(substring(first.name, 1, 1)),
    substring(first.name, 2, nchar(first.name))
  )
  first.name
}

req2 <- function(x) {
  if ("reactivevalues" %in% class(x)) {
    if (length(names(x)) == 0) {
      return(req(FALSE))
    }
    return(req(all(!sapply(x, is.null))))
  }
  if ("reactive" %in% class(x)) {
    return(req(x, x()))
  }
  req(x)
}

envcat <- function(var) {
  message(var, " = ", Sys.getenv(var))
}

mem.vmrss <- function(digits = 0) {
  mem <- "[? MB]"
  if (Sys.info()["sysname"] %in% c("Linux")) {
    proc <- paste("/proc", Sys.getpid(), "status", sep = "/")
    rss <- system(paste("grep -i vmrss", proc), intern = TRUE)
    rss <- gsub("VmRSS:[\t ]+| kB", "", rss)
    rss <- as.numeric(rss) / (1024) ## MB
    mem <- paste0(round(rss, digits), "MB")
  }
  mem
}

mem.proc <- function(digits = 0) {
  mem <- "[? MB]"
  if (Sys.info()["sysname"] %in% c("Linux")) {
    file <- paste("/proc", Sys.getpid(), "stat", sep = "/")
    what <- vector("character", 52)
    ## In your logging routine
    vsz <- as.numeric(scan(file, what = what, quiet = TRUE)[23])
    vsz <- vsz / (1024**2) ## MB
    mem <- paste0(round(vsz, digits), "MB")
  }
  mem
}

info <- function(..., type = "INFO") {
  dd <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- "some message"
  msg <- sapply(list(...), paste, collapse = " ")
  dd <- paste0("[", dd, "]")
  mm <- paste0("[", mem.proc(), "/", mem.vmrss(), "]")
  type <- paste0("[", type, "]")
  message(paste0(type, dd, mm, " --- ", sub("\n$", "", paste(msg, collapse = " "))))
}

dbg <- function(...) info(..., type = "DBUG")

## Parse access logs
ACCESS.LOG <- NULL

tipify2 <- function(...) {
  withTooltip(..., placement = "top", options = list(container = "body"))
}
tipifyL <- function(...) {
  withTooltip(..., placement = "left", options = list(container = "body"))
}
tipifyR <- function(...) {
  withTooltip(..., placement = "right", options = list(container = "body"))
}
tipifyT <- function(...) {
  withTooltip(..., placement = "top", options = list(container = "body"))
}
tipifyB <- function(...) {
  withTooltip(..., placement = "bottom", options = list(container = "body"))
}

# TODO: this function can be a variable
## in.shinyproxy <- function() {
##   return(Sys.getenv("SHINYPROXY_USERNAME") != "")
## }

tabRequire <- function(pgx, session, tabname, slot, enable = TRUE) {
  has.slot <- (slot %in% names(pgx)) && !is.null(pgx[[slot]])
  if (has.slot && enable) {
    bigdash.showTab(session, tabname)
  } else {
    bigdash.hideTab(session, tabname)
  }
}

fileRequire <- function(file, tabname, subtab) {
  file1 <- playbase::search_path(c(FILES, FILESX), file)
  has.file <- !is.null(file1) && file.exists(file1)
  if (!has.file) {
    message(paste("[MAIN] file ", file, " not found. Hiding", subtab, "\n"))
    hideTab(tabname, subtab)
  } else {
    message(paste("[MAIN] file ", file, " available. Showing", subtab, "\n"))
    showTab(tabname, subtab)
  }
}

# theming for the loading spinner
waiter::waiter_set_theme(html = waiter::spin_3(), color = waiter::transparent(.5))

tabView <- function(title, tab.inputs, tab.ui, id = title) {
  shiny::tabPanel(title,
    id = id,
    # enabkle loading spinner for plots etc.
    waiter::autoWaiter(color = waiter::transparent()),
    shiny::sidebarLayout(
      shiny::sidebarPanel(width = 2, tab.inputs, id = "sidebar"),
      shiny::mainPanel(width = 10, tab.ui)
    )
  )
}

toggleTab <- function(inputId, target, do.show, req.file = NULL, session = session) {
  if (!is.null(req.file)) {
    file1 <- playbase::search_path(c(FILES, FILESX), req.file)
    has.file <- !is.null(file1[1])
    do.show <- do.show && has.file
  }
  if (do.show) {
    shiny::showTab(inputId, target)
  } else {
    shiny::hideTab(inputId, target)
  }
}


## ======================================================================
## ======================== SEVER =======================================
## ======================================================================

sever_disconnected <- function() {
  sever_crash(error = NULL)
}

sendErrorLogToCustomerSuport <- function(user_email, pgx_name, raw_dir, error, path_to_creds = "hubspot_creds") {
  if (!file.exists(path_to_creds)) {
    message("[sendErrorMessageToCustomerSuport] WARNING : ticket not opened. cannot get credential =", path_to_creds)
    return(NULL)
  }

  user_email <- trimws(user_email)

  if (user_email == "") {
    user_email <- "No email (dev?)"
  }

  message <- glue::glue(
    "
        The user name is : {user_email}

        The ds name is: {pgx_name}

        Upload folder is: {raw_dir}

          The error is:

          {error}"
  )

  # Define the payload
  payload <- list(
    fields = list(
      list(
        objectTypeId = "0-1",
        name = "email",
        value = user_email
      ),
      list(
        objectTypeId = "0-5",
        name = "subject",
        value = "Blue Screen Crash Ticket"
      ),
      list(
        objectTypeId = "0-5",
        name = "content",
        value = message
      )
    )
  )

  # Convert the payload to JSON
  json_payload <- jsonlite::toJSON(payload, auto_unbox = TRUE)

  bearer_token <- readLines(path_to_creds)

  message("[sendErrorLogToCustomerSuport] Sending error log to customer support")


  # Send the POST request to HubSpot
  response <- httr::POST(
    url = "https://api.hsforms.com/submissions/v3/integration/secure/submit/24974201/9485b387-8cbf-4e1b-b93c-21b0d04956f3",
    httr::add_headers(
      `Content-Type` = "application/json",
      `Authorization` = paste("Bearer", bearer_token)
    ),
    body = json_payload,
    encode = "json"
  )
}


sever_crash <- function(error = NULL) {
  err_message <- NULL
  err_traceback <- NULL
  if (!is.null(error)) {
    err_traceback <- capture.output(
      printStackTrace(
        error,
        full = get_devmode_option(
          "shiny.fullstacktrace",
          FALSE
        ),
        offset = getOption("shiny.stacktraceoffset", TRUE)
      ),
      type = "message"
    )
    err_message <- error$message
  }

  shiny::tagList(
    shiny::div(
      style = "
          width: 100vw;
          height: 100vh;
        ",
      shiny::div(
        style = "
            transform: translateY(50%);
            background-color: #004c7d;
          ",
        shiny::tags$h1(
          "Woops!",
          style = "color:white;font-family:lato;"
        ),
        shiny::p("You have been disconnected", style = "font-size:15px;"),
        shiny::br(),
        shiny::div(shiny::img(
          src = base64enc::dataURI(file = "www/disconnected.png"),
          width = 540, height = 300
        )),
        shiny::br(),
        sever::reload_button("Relaunch", class = "info"),
        if (!is.null(error)) {
          tags$button("Show error",
            id = "showModalBtn",
            onClick = "document.getElementById('crashModal').style.display = 'block';",
            class = "btn btn-danger"
          )
        }
      ),
      tags$div(
        id = "crashModal",
        class = "modal",
        style = "
                      display: none;
                      position: fixed;
                      z-index: 1;
                      left: 0;
                      top: 0;
                      width: 100%;
                      height: 100%;
                      overflow: auto;
                      background-color: rgba(0,0,0,0.4);
                     ",
        tags$div(
          class = "modal-content",
          style = "
              background-color: #fefefe;
              margin: 15% auto;
              padding: 20px;
              border: 1px solid #888;
              width: 45%;
              color: black;
              text-align: left;
              height: 70vh;
              overflow-y: auto;
            ",
          tags$button(
            class = "btn btn-info", HTML("&times;"),
            onClick = "document.getElementById('crashModal').style.display = 'none';",
            style = "
                                  position: absolute;
                                  top: 5px;
                                  right: 5px;
                                 "
          ),
          tags$h3("Error Message"),
          tags$p(err_message),
          tags$h3("Error traceback"),
          tags$p(HTML(paste(err_traceback, collapse = "<br>")))
        )
      )
    )
  )
}

sever_ciao <- function(title = "We hope you enjoyed your stay!",
                       msg = "user session ended due to inactivity") {
  shiny::tagList(
    shiny::div(
      style = "
          width: 100vw;
          height: 100vh;
        ",
      shiny::div(
        style = "
            transform: translateY(30%);
            background-color: #004c7d;
          ",
        shiny::tags$h1(
          title,
          style = "color:white;font-family:lato;"
        ),
        shiny::p(msg, style = "font-size:18px;"),
        shiny::br(),
        shiny::div(
          shiny::img(
            src = base64enc::dataURI(file = "www/monster-ciao.png"),
            ## width = 362,
            height = 350
          ),
          style = "transform: translateX(50px);"
        ),
        shiny::br(),
        shiny::br(),
        shiny::br(),
        sever::reload_button("Relaunch", class = "info"),
      )
    )
  )
}

sever_serverfull <- function(srv) {
  shiny::tagList(
    shiny::tags$h1("Sorry, the Playground is full!", style = "color:white;font-family:lato;"),
    shiny::p("Our servers are at capacity. Please try again later.", style = "font-size:15px;"),
    shiny::br(),
    shiny::div(paste("server =", srv), style = "font-size:11px;text-align:center;"),
    shiny::br(), shiny::br(),
    sever::reload_button("Relaunch", class = "info")
  )
}

## From https://github.com/plotly/plotly.js/blob/master/src/components/modebar/buttons.js
all.plotly.buttons <- c(
  "toImage",
  "sendDataToCloud", "editInChartStudio", "zoom2d", "pan2d", "select2d",
  "lasso2d", "drawclosedpath", "drawopenpath", "drawline", "drawrect",
  "drawcircle", "eraseshape", "zoomIn2d", "zoomOut2d",
  "autoScale2d", "resetScale2d", "zoom3d", "pan3d",
  "orbitRotation", "tableRotation", "resetCameraDefault3d",
  "resetCameraLastSave3d", "hoverClosest3d", "zoomInGeo",
  "zoomOutGeo", "resetGeo", "hoverClosestGeo", "hoverClosestGl2d",
  "hoverClosestPie", "resetViewSankey", "toggleHover",
  "hoverClosestCartesian", "hoverCompareCartesian",
  "resetViews", "toggleSpikelines",
  "resetViewMapbox", "zoomInMapbox", "zoomOutMapbox"
)
