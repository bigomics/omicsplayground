##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
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

# TODO: this is a version that respects the DEBUG flag, but it's not being used
dbg <- function(...) {
	if(DEBUG) {
		msg <- list(...)
		msg <- paste(sapply(msg, function(s) paste(s,collapse=" ")),collapse=" ")
		message(msg)
	}
}

dbg <- function(...) {
	##msg = paste0(ifelse(is.null(module),"",paste0("<",module,"> ")),msg)
	msg = sapply( list(...),paste,collapse=" ")
	message(paste0("DBG ",sub("\n$","",paste(msg,collapse=" "))))
}

## Parse access logs
ACCESS.LOG <- NULL
if(0) {
    access.dirs = c("/var/www/html/logs", "/var/log/apache2","/var/log/apache",
                    "../logs","/var/log/httpd","/var/log/nginx")
    access.dirs <- access.dirs[dir.exists(access.dirs)]
    access.dirs
    ##ACCESS.LOG <- pgx.parseAccessLogs(access.dirs[], filter.get=NULL)
    ACCESS.LOG <- pgx.parseAccessLogs(access.dirs[], filter.get="playground")
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
	tipify(disabled(...),
				 "This is a Premium feature. Upgrade to enable this feature."
	)    
	
}

# TODO: this function can be a variable
in.shinyproxy <- function() {
	## Determine if we are in ShinyProxy
	##
	vars <- c("SHINYPROXY_USERNAME","SHINYPROXY_USERGROUPS",
						"PLAYGROUND_USERID","PLAYGROUND_LEVEL")
	vars <- c("SHINYPROXY_USERNAME","SHINYPROXY_USERGROUPS")
	vals <- sapply(vars,Sys.getenv)
	all(vals!="")
}

tabRequire <- function(pgx, slot, tabname, subtab) {
	if(!slot %in% names(pgx)) {
		cat(paste("[MAIN] object has no ",slot," results. hiding tab.\n"))
		hideTab(tabname, subtab)
	} else {
		showTab(tabname, subtab)
	}
}

fileRequire <- function(file, tabname, subtab) {
	file1 <- search_path(c(FILES,FILESX),file)
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

toggleTab <- function(inputId, target, do.show, req.file=NULL ) {
    if(!is.null(req.file)) {
        file1 <- search_path(c(FILES,FILESX),req.file)
        has.file <- !is.null(file1[1])
        do.show <- do.show && has.file
    }
    if(do.show) {
        shiny::showTab(inputId, target)
    }
    if(!do.show) {
        shiny::hideTab(inputId, target)
    }
}


##======================================================================
##======================== SEVER =======================================
##======================================================================

sever_screen <- function() {
    shiny::tagList(
               shiny::tags$h1(
                               "Houston, we have a problem", style = "color:white;font-family:lato;"
                           ),
               shiny::p("You have been disconnected!", style="font-size:15px;"),
               shiny::br(),
               shiny::div(shiny::img(src=base64enc::dataURI(file="www/lost-in-space.gif"),
                                     width=500,height=250)),
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

sever_screen0 <- function() {
    shiny::tagList(
               shiny::tags$h1(
                               "Houston, we have a problem", style="color:white;font-family:lato;"
                           ),
               shiny::p("You have been disconnected!", style="font-size:15px;"),
               shiny::br(),
               shiny::div(shiny::img(src=base64enc::dataURI(file="www/lost-in-space.gif"),
                          width=500,height=250)),
               shiny::br(),
               sever::reload_button("Relaunch", class = "info")
           )
}


sever_screen2 <- function(session_id) {
  shiny::tagList(
    shiny::tags$h1(
      "Houston, we have a problem", style = "color:white;font-family:lato;"
    ),
    shiny::p("You have been disconnected!", style="font-size:15px;"),
    shiny::br(),
    shiny::div(shiny::img(src=base64enc::dataURI(file="www/lost-in-space.gif"),
                          width=500,height=250)),
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
