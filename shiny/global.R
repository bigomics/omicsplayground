
tabView <- function(title, tab.inputs, tab.ui) {
    tabPanel(title, ## id=title,
             sidebarLayout(
                 sidebarPanel( width=2, tab.inputs, id="sidebar"),
                 mainPanel( width=10, tab.ui)
             ))
}

## dev.tabView <- function(title, tab.inputs, tab.ui) {
##     if(!DEV.VERSION) return(NULL)
##     tabView(title, tab.inputs, tab.ui)
## }
## dev.tabPanel <- function(id, ui) {
##     if(!DEV.VERSION) return(NULL)
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
