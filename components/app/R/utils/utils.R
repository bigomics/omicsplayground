##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#########################################################################
##                                                                     ##
##              Utility Functions for Omics Playground                 ##
##                                                                     ##
#########################################################################


req2 <- function(x) {
    if("reactivevalues" %in% class(x)) {
        if(length(names(x))==0) return(req(FALSE))
        return(req(all(!sapply(x,is.null))))
    }
    if("reactive" %in% class(x)) {
        return(req(x,x()))
    }
    req(x)
}

envcat <- function(var) {
	message(var," = ",Sys.getenv(var))
}

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

## Parse access logs
ACCESS.LOG <- NULL
if(0) {
    access.dirs = c("/var/www/html/logs", "/var/log/apache2","/var/log/apache",
                    "../logs","/var/log/httpd","/var/log/nginx")
    access.dirs <- access.dirs[dir.exists(access.dirs)]
    access.dirs
    ##ACCESS.LOG <- playbase::pgx.parseAccessLogs(access.dirs[], filter.get=NULL)
    ACCESS.LOG <- playbase::pgx.parseAccessLogs(access.dirs[], filter.get="playground")
    names(ACCESS.LOG)
    sum(ACCESS.LOG$visitors$count)
}


tipify2 <- function(...) {
	withTooltip(..., placement="top", options = list(container = "body"))
}
tipifyL <- function(...) {
	withTooltip(..., placement="left", options = list(container = "body"))
}
tipifyR <- function(...) {
	withTooltip(..., placement="right", options = list(container = "body"))
}
tipifyT <- function(...) {
	withTooltip(..., placement="top", options = list(container = "body"))
}
tipifyB <- function(...) {
	withTooltip(..., placement="bottom", options = list(container = "body"))
}

# TODO: this function isn't being used
premium.feature <- function(...) {
	message("[premium.feature] USER_MODE = ",USER_MODE)
	message("[premium.feature] DEV = ",DEV)
	el <- list(...)
	if(USER_MODE %in% c("pro","premium","dev")) return(el)
	tipify( disabled(...),
		 "This is a Premium feature. Upgrade to enable this feature."
	)
}

# TODO: this function can be a variable
in.shinyproxy <- function() {
        return(Sys.getenv("SHINYPROXY_USERNAME") != "")
}

tabRequire <- function(pgx, session, tabname, slot) {
        has.slot <- (slot %in% names(pgx))
        if(!has.slot) {
          cat(paste("[MAIN] object has no ",slot," results. hiding tab.\n"))
          ##hideTab(tabname, subtab)
          bigdash.hideTab(session, tabname)
	} else {
          ##showTab(tabname, subtab)
          bigdash.showTab(session, tabname)
	}
}

fileRequire <- function(file, tabname, subtab) {
	file1 <- playbase::search_path(c(FILES,FILESX),file)
	has.file <- !is.null(file1) && file.exists(file1)
	if(!has.file) {
		message(paste("[MAIN] file ",file," not found. Hiding",subtab,"\n"))
		hideTab(tabname, subtab)
	} else {
		message(paste("[MAIN] file ",file," available. Showing",subtab,"\n"))
		showTab(tabname, subtab)
	}
}

# theming for the loading spinner
waiter::waiter_set_theme(html = waiter::spin_3(), color = waiter::transparent(.5))

tabView <- function(title, tab.inputs, tab.ui, id=title) {
    shiny::tabPanel(title, id=id,
        # enabkle loading spinner for plots etc.
        waiter::autoWaiter(color = waiter::transparent()),
             shiny::sidebarLayout(
                 shiny::sidebarPanel( width=2, tab.inputs, id="sidebar"),
                 shiny::mainPanel( width=10, tab.ui)
             ))
}

toggleTab <- function(inputId, target, do.show, req.file=NULL, session=session ) {
    if(!is.null(req.file)) {
        file1 <- playbase::search_path(c(FILES,FILESX),req.file)
        has.file <- !is.null(file1[1])
        do.show <- do.show && has.file
    }
    if(do.show) {
      shiny::showTab(inputId, target)
    } else {
      shiny::hideTab(inputId, target)
    }
}


##======================================================================
##======================== SEVER =======================================
##======================================================================

sever_screen <- function() {
    shiny::tagList(
               shiny::tags$h1(
                               "Woops!", style = "color:white;font-family:lato;"
                           ),
               shiny::p("You have been disconnected", style="font-size:15px;"),
               shiny::br(),
               shiny::div(shiny::img(src=base64enc::dataURI(file="www/disconnected.png"),
                                     width=450,height=250)),
               shiny::div(
                          id="logSub",
                          ##        shiny::textAreaInput(
                          ##               inputId = "logMsg",
                          ##               label = "",
                          ##               width = "100%", height="80px",
                          ##               value = "If this was a crash, please help and describe here the last thing you did."
                          ##        ),
                          shiny::br(),
                          shiny::tags$a(
                                          onClick = "sendLog()",
                                          class = "btn btn-sm btn-warning",
                                          "Send error to developers"
                                      )
                      ),
               shiny::div(
                          id="logSubbed",
                          style="display:none;",
                          shiny::p("Mission Control has been notified. Thank you!", style="font-size:15px;")
                      ),
               shiny::br(),
               shiny::div(
                          id="sever-reload-btn",
                          sever::reload_button("Relaunch", class = "info"),
                          style="display:none;"
                      )
           )
}

sever_screen0 <- function(error = NULL) {
  err_message <- NULL
  err_traceback <- NULL
  if(!is.null(error)){
    err_traceback <- capture.output(
      printStackTrace(
        error,
        full = get_devmode_option("shiny.fullstacktrace",
                                  FALSE),
        offset = getOption("shiny.stacktraceoffset", TRUE)),
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
          ",
          shiny::tags$h1(
            "Woops!", style="color:white;font-family:lato;"
          ),
          shiny::p("You have been disconnected", style="font-size:15px;"),
          shiny::br(),
          shiny::div(shiny::img(src=base64enc::dataURI(file="www/disconnected.png"),
                                width=540,height=300)),
          shiny::br(),
          sever::reload_button("Relaunch", class = "info"),
          if(!is.null(error)){
            tags$button("Show error",
                        id = "showModalBtn",
                        onClick = "document.getElementById('crashModal').style.display = 'block';",
                        class = "btn btn-danger")
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
            tags$button(class = "btn btn-info", HTML("&times;"),
                        onClick = "document.getElementById('crashModal').style.display = 'none';",
                        style = "
                                  position: absolute;
                                  top: 5px;
                                  right: 5px;
                                 "),
            tags$h3("Error Message"),
            tags$p(err_message),
            tags$h3("Error traceback"),
            tags$p(HTML(paste(err_traceback, collapse = "<br>")))
          )
        )
      )
    )
}


sever_screen2 <- function(session_id) {
  shiny::tagList(
    shiny::tags$h1(
      "Woops!", style = "color:white;font-family:lato;"
    ),
    shiny::p("You have been disconnected", style="font-size:15px;"),
    shiny::br(),
    shiny::div(shiny::img(src=base64enc::dataURI(file="www/disconnected.png"),
                          width=450,height=250)),
    shiny::div(
      id="logSub",
      ##        shiny::textAreaInput(
      ##               inputId = "logMsg",
      ##               label = "",
      ##               width = "100%", height="80px",
      ##               value = "If this was a crash, please help and describe here the last thing you did."
      ##        ),
      shiny::br(),
      shiny::tags$a(
        onClick = HTML(paste0("sendLog2('",session_id,"')")),
        class = "btn btn-sm btn-warning",
        "Send error to developers"
      )
    ),
    shiny::div(
      id="logSubbed",
      style="display:none;",
      shiny::p("Mission Control has been notified. Thank you!", style="font-size:15px;")
    ),
    shiny::br(),
    shiny::div(
      id="sever-reload-btn",
      sever::reload_button("Relaunch", class = "info"),
      style="display:none;"
    )
  )
}


## From https://github.com/plotly/plotly.js/blob/master/src/components/modebar/buttons.js
all.plotly.buttons = c(
	"toImage",
	"senDataToCloud","editInChartStudio","zoom2d","pan2d","select2d",
	"lasso2d","drawclosedpath","drawopenpath","drawline","drawrect",
	"drawcircle","eraseshape","zoomIn2d","zoomOut2d",
	"autoScale2d","resetScale2d","zoom3d","pan3d",
	"orbitRotation","tableRotation","resetCameraDefault3d",
	"resetCameraLastSave3d","hoverClosest3d","zoomInGeo",
	"zoomOutGeo","resetGeo","hoverClosestGeo","hoverClosestGl2d",
	"hoverClosestPie","resetViewSankey","toggleHover",
	"hoverClosestCartesian","hoverCompareCartesian",
	"resetViews","toggleSpikelines",
	"resetViewMapbox","zoomInMapbox","zoomOutMapbox"
)
