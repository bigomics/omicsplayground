##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##


message("===============================================================")
message("======================= init.R ================================")
message("===============================================================")

## Parse access logs
access.dirs = c("/var/www/html/logs", "/var/log/apache2","/var/log/apache",
                "/var/log/httpd","/var/log/nginx","../logs")
access.dirs <- access.dirs[dir.exists(access.dirs)]
access.dirs
##ACCESS.LOG <- pgx.parseAccessLogs(access.dirs[], filter.get=NULL)
ACCESS.LOG <- pgx.parseAccessLogs(access.dirs[], filter.get="playground")
names(ACCESS.LOG)
sum(ACCESS.LOG$visitors$count)

##-----------------------------------------------------
## Orca server
##-----------------------------------------------------

message("*****************************************")
message("***** starting local ORCA server ********")
message("*****************************************\n")

## see: pgx-module.R
ORCA <- initOrca(launch=TRUE) 
class(ORCA)
if(is.null(ORCA)) {
    warning("##### FATAL:: Could not connect to ORCA server. Please start ORCA. #####")
    stop()
}

##======================================================================
##==================== FUNCTIONS =======================================
##======================================================================

tipify2 <- function(...) {
    shinyBS::tipify(..., placement="top", options = list(container = "body"))
}

premium.feature <- function(...) {
    message("[premium.feature] USER_MODE = ",USER_MODE)
    message("[premium.feature] DEV = ",DEV)        
    el <- list(...)
    if(USER_MODE %in% c("pro","premium","dev")) return(el)
    tipify(disabled(...),
           "This is a Premium feature. Upgrade to enable this feature."
           )    

}

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

tabView <- function(title, tab.inputs, tab.ui) {
    tabPanel(title, ## id=title,
             sidebarLayout(
                 sidebarPanel( width=2, tab.inputs, id="sidebar"),
                 mainPanel( width=10, tab.ui)
             ))
}

## dev.tabView <- function(title, tab.inputs, tab.ui) {
##     if(!DEV.MODE) return(NULL)
##     tabView(title, tab.inputs, tab.ui)
## }
## dev.tabPanel <- function(id, ui) {
##     if(!DEV.MODE) return(NULL)
##     tabPanel(id, ui)
## }

social_buttons <- function() {
    div(
        id="social-buttons",
        tagList(
            tipify( tags$a( href="https://omicsplayground.readthedocs.io", icon("book"), target="_blank"),
                   "Read our online documentation at Read-the-docs", placement="top"),
            tipify( tags$a( href="https://www.youtube.com/watch?v=_Q2LJmb2ihU&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-",
                           icon("youtube"), target="_blank"),
                   "Watch our tutorials on YouTube", placement="top"),
            tipify( tags$a( href="https://github.com/bigomics/omicsplayground",
                           icon("github"), target="_blank"),
                   "Get the source code or report a bug at GitHub", placement="top"),
            tipify( tags$a( href="https://hub.docker.com/r/bigomics/omicsplayground",
                           icon("docker"), target="_blank"),
                   "Pull our docker from Docker", placement="top"),
            tipify( tags$a( href="https://groups.google.com/d/forum/omicsplayground",
                           icon("users"), target="_blank"),
                   "Get help at our user forum", placement="top")            
        )
    )
}

TAGS.JSSCRIPT =
    ## https://stackoverflow.com/questions/36995142/get-the-size-of-the-window-in-shiny
    tags$head(tags$script('
    var dimension = [0, 0];
    $(document).on("shiny:connected", function(e) {
        dimension[0] = window.innerWidth;
        dimension[1] = window.innerHeight;
        Shiny.onInputChange("dimension", dimension);
    });
    $(window).resize(function(e) {
        dimension[0] = window.innerWidth;
        dimension[1] = window.innerHeight;
        Shiny.onInputChange("dimension", dimension);
    });
')) 


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
    "resetViewMapbox","zoomInMapbox","zoomOutMapbox")




##======================================================================
##==================== END-OF-FILE =====================================
##======================================================================
