
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

